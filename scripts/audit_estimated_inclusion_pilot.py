"""Independent consistency audit for estimated-inclusion pilot output files.

This auditor reads only persisted CSV/JSON artifacts and independently
reconstructs arm summaries and paired absolute-error contrasts. It does not
import the simulation or aggregation implementation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import minimize
from scipy.special import expit, log_expit
from scipy.stats import norm


def _wilson(successes: int, trials: int) -> tuple[float, float]:
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("invalid binomial counts")
    z = 1.959963984540054
    p = successes / trials
    denominator = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / denominator
    half = z * np.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / denominator
    return max(0.0, float(center - half)), min(1.0, float(center + half))


def _close(actual: object, expected: object, label: str, *, atol: float = 1e-9) -> None:
    if expected is None:
        if actual is not None:
            raise AssertionError(f"{label}: expected null, got {actual!r}")
        return
    if actual is None or not np.isfinite(float(actual)) or not np.isclose(
        float(actual), float(expected), rtol=1e-9, atol=atol
    ):
        raise AssertionError(f"{label}: stored={actual!r}, reconstructed={expected!r}")



def _strict_bool(value: object, label: str) -> bool:
    """Parse persisted Boolean values without Python's truthiness trap."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1"}:
            return True
        if normalized in {"false", "0"}:
            return False
    raise AssertionError(f"{label}: expected an explicit Boolean/0/1 value, got {value!r}")


def _required_finite_number(value: object, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AssertionError(f"{label}: expected a finite number") from exc
    if not np.isfinite(number):
        raise AssertionError(f"{label}: expected a finite number")
    return number

def _integer_values(values: pd.Series, label: str) -> np.ndarray:
    """Parse integer-valued CSV columns without silently truncating decimals."""
    numeric = pd.to_numeric(values, errors="raise").to_numpy(dtype=float)
    if not np.all(np.isfinite(numeric)) or not np.all(numeric == np.floor(numeric)):
        raise AssertionError(f"{label}: expected finite integer-valued entries")
    return numeric.astype(np.int64)



def _audit_inference_payload(row: pd.Series, h0_true: float) -> None:
    """Reconstruct the toy selected-event posterior from persisted numerical inputs.

    This intentionally does not import the production likelihood or posterior
    helpers. It uses direct event-by-event quadrature as an independent
    numerical reference and also reconstructs the reported density ISE.
    """
    schema = row.get("inference_audit_schema")
    if schema != "toy_selection_posterior_reconstruction_v1":
        raise AssertionError("successful row lacks the supported inference-audit schema")
    raw = row.get("inference_audit_payload_json")
    if not isinstance(raw, str) or not raw:
        raise AssertionError("successful row lacks inference_audit_payload_json")
    try:
        payload = json.loads(raw)
        scoring_raw = row.get("scoring_audit_payload_json")
        scoring = json.loads(scoring_raw) if isinstance(scoring_raw, str) and scoring_raw else {}
        z = np.asarray(payload["redshift_grid"], dtype=float)
        density = np.asarray(payload["host_density"], dtype=float)
        target = np.asarray(scoring["target_density_for_ise"], dtype=float)
        h0 = np.asarray(payload["h0_grid"], dtype=float)
        distances = np.asarray(payload["observed_distances_mpc"], dtype=float)
        sigma = float(payload["sigma_mpc"])
        threshold = float(payload["detection_threshold_mpc"])
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise AssertionError("inference audit payload is malformed") from exc

    if (z.ndim != 1 or z.size < 21 or density.shape != z.shape
            or target.shape != z.shape or h0.ndim != 1 or h0.size < 31
            or distances.ndim != 1 or distances.size < 1):
        raise AssertionError("inference audit payload arrays have invalid shapes")
    if (not all(np.all(np.isfinite(x)) for x in (z, density, target, h0, distances))
            or np.any(np.diff(z) <= 0) or np.any(np.diff(h0) <= 0)
            or np.any(density < 0) or np.any(target <= 0)
            or np.any(distances < threshold) or sigma <= 0 or threshold <= 0):
        raise AssertionError("inference audit payload contains invalid numerical inputs")
    if not np.isclose(trapezoid(density, z), 1.0, rtol=1e-7, atol=1e-9):
        raise AssertionError("persisted host density does not integrate to one")

    # Refit the estimated logistic inclusion model from the saved independent
    # validation records. This is a separate optimization implementation and
    # never reads hidden parent truth or the scoring-only probability MAE.
    arm_name = str(row.get("arm", ""))
    if arm_name.startswith("estimated_"):
        try:
            validation_z = np.asarray(payload["validation_redshift"], dtype=float)
            validation_y_raw = np.asarray(payload["validation_included"])
            stored_model = payload["fitted_inclusion_model"]
            stored_probability = np.asarray(payload["density_inclusion_probability"], dtype=float)
        except (KeyError, TypeError, ValueError) as exc:
            raise AssertionError("estimated-arm validation audit inputs are malformed") from exc
        if (validation_z.ndim != 1 or validation_z.size < 10
                or validation_y_raw.shape != validation_z.shape
                or not isinstance(stored_model, dict)
                or stored_probability.ndim != 1):
            raise AssertionError("estimated-arm validation audit arrays have invalid shapes")
        if (not np.all(np.isfinite(validation_z))
                or not np.all(np.isin(validation_y_raw, [0, 1, False, True]))
                or np.unique(validation_y_raw).size != 2):
            raise AssertionError("estimated-arm validation records are invalid")
        validation_y = validation_y_raw.astype(float)
        z_min, z_max = float(z[0]), float(z[-1])
        if np.any((validation_z < z_min) | (validation_z > z_max)):
            raise AssertionError("validation redshifts fall outside the declared fitting domain")
        if validation_z.size != int(row.get("validation_size")):
            raise AssertionError("validation record count disagrees with arm validation_size")
        x_validation = (validation_z - (z_min + z_max)/2.0) / ((z_max-z_min)/2.0)

        degree = 2 if arm_name.startswith("estimated_quadratic_") else 1
        quadratic_basis = x_validation**2 - 1.0 / 3.0

        def independent_eta(theta):
            eta = theta[0] + theta[1] * x_validation
            if degree == 2:
                eta = eta + theta[2] * quadratic_basis
            return eta

        def independent_nll(theta):
            eta = independent_eta(theta)
            return float(-np.sum(
                validation_y*log_expit(eta) + (1.0-validation_y)*log_expit(-eta)
            ))

        def independent_gradient(theta):
            residual = expit(independent_eta(theta)) - validation_y
            gradient = [np.sum(residual), np.dot(residual, x_validation)]
            if degree == 2:
                gradient.append(np.dot(residual, quadratic_basis))
            return np.asarray(gradient, dtype=float)

        bounds = ((-25.0, 25.0), (-50.0, 50.0)) + (
            ((-50.0, 50.0),) if degree == 2 else ()
        )
        refit = minimize(
            independent_nll, np.zeros(degree + 1), jac=independent_gradient,
            method="L-BFGS-B", bounds=bounds,
            options={"maxiter": 10000, "maxfun": 50000, "maxls": 100,
                     "ftol": 1e-15, "gtol": 1e-8},
        )
        # L-BFGS-B can report an abnormal line-search termination for a
        # nearly separated sample even when it has reached the constrained
        # optimum. Independently verify the box-constrained KKT condition
        # rather than accepting a failed status or blindly relaxing tolerances.
        if np.all(np.isfinite(refit.x)):
            gradient = independent_gradient(refit.x)
            projected_gradient = gradient.copy()
            for j, (lower, upper) in enumerate(bounds):
                if refit.x[j] <= lower + 1e-7 and gradient[j] > 0:
                    projected_gradient[j] = 0.0
                elif refit.x[j] >= upper - 1e-7 and gradient[j] < 0:
                    projected_gradient[j] = 0.0
            kkt_residual = float(np.max(np.abs(projected_gradient)))
        else:
            kkt_residual = float("inf")
        if (not np.all(np.isfinite(refit.x)) or not np.isfinite(refit.fun)
                or kkt_residual > 1e-5):
            raise AssertionError(
                "independent logistic inclusion refit did not converge "
                f"(status={refit.message}; projected_gradient={kkt_residual:.3g})"
            )
        stored_intercept = _required_finite_number(stored_model.get("intercept"), "fitted intercept")
        stored_slope = _required_finite_number(stored_model.get("slope"), "fitted slope")
        stored_quadratic = _required_finite_number(
            stored_model.get("quadratic", 0.0), "fitted quadratic coefficient"
        )
        expected_coefficients = [stored_intercept, stored_slope]
        if degree == 2:
            expected_coefficients.append(stored_quadratic)
        if not np.allclose(refit.x, expected_coefficients, rtol=2e-5, atol=2e-5):
            raise AssertionError("stored fitted inclusion coefficients disagree with independent refit")
        x_catalogue = (np.asarray(payload["catalogue_redshift"], dtype=float)
                       - (z_min + z_max)/2.0) / ((z_max-z_min)/2.0)
        eta_catalogue = refit.x[0] + refit.x[1] * x_catalogue
        if degree == 2:
            eta_catalogue += refit.x[2] * (x_catalogue**2 - 1.0 / 3.0)
        refit_probability = expit(eta_catalogue)
        if (stored_probability.shape != refit_probability.shape
                or not np.allclose(stored_probability, refit_probability, rtol=5e-5, atol=5e-5)):
            raise AssertionError("stored estimated inclusion probabilities disagree with independent refit")
        expected_model_name = (
            "quadratic_logistic_inclusion_v1" if degree == 2 else "logistic_inclusion_v1"
        )
        if (stored_model.get("model") != expected_model_name
                or int(stored_model.get("degree", 1)) != degree
                or int(stored_model.get("n_validation", -1)) != validation_z.size
                or int(stored_model.get("n_included", -1)) != int(validation_y.sum())
                or int(stored_model.get("n_not_included", -1)) != int(validation_y.size-validation_y.sum())):
            raise AssertionError("stored inclusion-model metadata disagrees with validation records")
        if row.get("inclusion_model_degree") != degree:
            raise AssertionError("row inclusion_model_degree disagrees with arm family")
    elif payload.get("fitted_inclusion_model") is not None:
        raise AssertionError("non-estimated arm unexpectedly contains a fitted inclusion model")

    # Enforce a strict boundary between inference inputs and scoring truth.
    forbidden_truth_keys = {
        "true_inclusion_probability", "parent_redshift", "hidden_parent_redshift",
        "true_probability", "parent_catalogue", "target_density_for_ise",
    }
    if forbidden_truth_keys.intersection(payload):
        raise AssertionError("inference audit payload contains scoring-only or hidden-parent truth")
    if scoring.get("schema") != "inclusion_probability_scoring_only_v1":
        raise AssertionError("successful row lacks supported scoring-only payload schema")
    if target.shape != z.shape or not np.all(np.isfinite(target)) or np.any(target <= 0):
        raise AssertionError("scoring-only target density is malformed")
    if arm_name.startswith("estimated_"):
        try:
            parent_z = np.asarray(scoring["parent_redshift"], dtype=float)
            true_parent_p = np.asarray(scoring["true_inclusion_probability"], dtype=float)
            estimated_parent_p = np.asarray(scoring["estimated_inclusion_probability"], dtype=float)
        except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
            raise AssertionError("estimated-arm scoring-only payload is missing or malformed") from exc
        if (scoring.get("schema") != "inclusion_probability_scoring_only_v1"
                or parent_z.ndim != 1 or parent_z.size < 1
                or true_parent_p.shape != parent_z.shape
                or estimated_parent_p.shape != parent_z.shape):
            raise AssertionError("scoring-only payload arrays or schema are invalid")
        if (not all(np.all(np.isfinite(x)) for x in (parent_z, true_parent_p, estimated_parent_p))
                or np.any((true_parent_p <= 0) | (true_parent_p > 1))
                or np.any((estimated_parent_p <= 0) | (estimated_parent_p > 1))):
            raise AssertionError("scoring-only inclusion probabilities are outside (0, 1]")
        if np.any((parent_z < z[0]) | (parent_z > z[-1])):
            raise AssertionError("scoring parent redshifts fall outside the simulation domain")
        x_parent = (parent_z - (z_min + z_max)/2.0) / ((z_max-z_min)/2.0)
        eta_parent = refit.x[0] + refit.x[1] * x_parent
        if degree == 2:
            eta_parent += refit.x[2] * (x_parent**2 - 1.0 / 3.0)
        refit_parent_probability = expit(eta_parent)
        if not np.allclose(
            estimated_parent_p, refit_parent_probability, rtol=5e-5, atol=5e-5
        ):
            raise AssertionError(
                "scoring-only estimated probabilities disagree with independent logistic refit"
            )
        reconstructed_mae = float(np.mean(np.abs(estimated_parent_p - true_parent_p)))
        _close(row.get("inclusion_probability_mae"), reconstructed_mae,
               "independent_scoring.inclusion_probability_mae", atol=2e-10)
    elif any(key in scoring for key in (
        "parent_redshift", "true_inclusion_probability", "estimated_inclusion_probability"
    )):
        raise AssertionError("non-estimated arm unexpectedly contains parent-probability scoring truth")

    # Reconstruct the reflected, inverse-inclusion-weighted KDE from its
    # primitive inputs. This deliberately does not call estimate_host_density.
    try:
        catalogue_z = np.asarray(payload["catalogue_redshift"], dtype=float)
        inclusion_p = np.asarray(payload["density_inclusion_probability"], dtype=float)
        merger_rate = np.asarray(payload["merger_weight_grid"], dtype=float)
        raw_bandwidth_factor = payload["bandwidth_factor"]
    except (KeyError, TypeError, ValueError) as exc:
        raise AssertionError("KDE reconstruction inputs are missing or malformed") from exc
    if (catalogue_z.ndim != 1 or catalogue_z.size < 5
            or inclusion_p.shape != catalogue_z.shape
            or merger_rate.shape != z.shape):
        raise AssertionError("KDE reconstruction arrays have invalid shapes")
    if (not all(np.all(np.isfinite(x)) for x in (catalogue_z, inclusion_p, merger_rate))
            or np.any((catalogue_z < z[0]) | (catalogue_z > z[-1]))
            or np.unique(catalogue_z).size < 3
            or np.any((inclusion_p <= 0) | (inclusion_p > 1))
            or np.any(merger_rate <= 0)):
        raise AssertionError("KDE reconstruction inputs contain invalid values")
    if raw_bandwidth_factor is None:
        bandwidth_factor = None
    else:
        bandwidth_factor = _required_finite_number(raw_bandwidth_factor, "bandwidth_factor")
        if bandwidth_factor <= 0:
            raise AssertionError("KDE bandwidth factor must be positive")
    raw_weights = 1.0 / inclusion_p
    scaled_weights = raw_weights / np.max(raw_weights)
    weight_sum = float(np.sum(scaled_weights))
    effective_n = weight_sum**2 / float(np.sum(scaled_weights**2))
    weighted_mean = float(np.average(catalogue_z, weights=scaled_weights))
    weighted_variance = float(np.average(
        (catalogue_z - weighted_mean)**2, weights=scaled_weights
    ))
    if not np.isfinite(weighted_variance) or weighted_variance <= 0:
        raise AssertionError("KDE weighted sample variance is invalid")
    factor = effective_n**(-1.0 / 5.0) if bandwidth_factor is None else bandwidth_factor
    bandwidth = factor * np.sqrt(weighted_variance)
    reflected = np.concatenate((
        catalogue_z, 2.0*z[0] - catalogue_z, 2.0*z[-1] - catalogue_z
    ))
    reflected_weights = np.tile((scaled_weights / weight_sum) / 3.0, 3)
    kernel_norm = 1.0 / (np.sqrt(2.0*np.pi) * bandwidth)
    reconstructed_galaxy_density = np.empty_like(z)
    for start in range(0, z.size, 256):
        stop = min(start + 256, z.size)
        standardized = (z[start:stop, None] - reflected[None, :]) / bandwidth
        reconstructed_galaxy_density[start:stop] = (
            np.exp(-0.5*standardized**2) @ reflected_weights
        ) * kernel_norm
    unnormalized_density = reconstructed_galaxy_density * merger_rate
    kde_area = float(trapezoid(unnormalized_density, z))
    if not np.isfinite(kde_area) or kde_area <= 0:
        raise AssertionError("independent KDE normalization failed")
    reconstructed_density = unnormalized_density / kde_area
    if not np.allclose(density, reconstructed_density, rtol=2e-8, atol=2e-10):
        max_difference = float(np.max(np.abs(density - reconstructed_density)))
        raise AssertionError(
            "persisted host density disagrees with independent reflected-KDE reconstruction "
            f"(max_abs_difference={max_difference:.3g})"
        )

    # Independent direct quadrature: integrate each event's Gaussian mixture
    # separately, then apply the fixed-N detection-conditioning correction.
    event_log_terms = []
    beta_values = []
    for candidate in h0:
        means = 299792.458 * z / candidate
        event_integrals = [
            trapezoid(density * norm.pdf(distance, loc=means, scale=sigma), z)
            for distance in distances
        ]
        if np.any(np.asarray(event_integrals) <= 0):
            raise AssertionError("independent event likelihood integral is nonpositive")
        event_log_terms.append(float(np.sum(np.log(event_integrals))))
        beta = trapezoid(
            density * norm.cdf((means - threshold) / sigma), z
        )
        if not np.isfinite(beta) or beta <= 0 or beta > 1.0 + 1e-12:
            raise AssertionError("independent detection efficiency is invalid")
        beta_values.append(float(beta))

    selected_loglike = np.asarray(event_log_terms) - distances.size * np.log(beta_values)
    posterior = np.exp(selected_loglike - np.max(selected_loglike))
    posterior_area = trapezoid(posterior, h0)
    if not np.isfinite(posterior_area) or posterior_area <= 0:
        raise AssertionError("independent posterior normalization failed")
    posterior /= posterior_area
    cdf = np.r_[0.0, cumulative_trapezoid(posterior, h0)]
    cdf /= cdf[-1]
    q16, q50, q84 = (
        float(np.interp(q, cdf, h0)) for q in (0.16, 0.50, 0.84)
    )
    reconstructed = {
        "posterior_mean": float(trapezoid(h0 * posterior, h0)),
        "posterior_median": q50,
        "ci68_low": q16,
        "ci68_high": q84,
        "ci68_width": q84 - q16,
        "density_ise": float(trapezoid((density - target) ** 2, z)),
        "covered_68": bool(q16 <= h0_true <= q84),
        "absolute_error": abs(float(trapezoid(h0 * posterior, h0)) - h0_true),
    }
    for key, expected in reconstructed.items():
        actual = row.get(key)
        if isinstance(expected, bool):
            if _strict_bool(actual, key) is not expected:
                raise AssertionError(f"{key} disagrees with independent inference reconstruction")
        else:
            _close(actual, expected, f"independent_inference.{key}", atol=2e-8)

def audit_pilot_artifacts(output_dir: str | Path) -> dict:
    root = Path(output_dir)
    csv_path, json_path = root / "replicates.csv", root / "summary.json"
    if not csv_path.is_file() or not json_path.is_file():
        raise FileNotFoundError("pilot output must contain replicates.csv and summary.json")
    frame = pd.read_csv(csv_path)
    summary = json.loads(json_path.read_text(encoding="utf-8"))
    required = {
        "replicate", "replicate_seed", "catalogue_seed", "audit_seed", "arm",
        "validation_size", "status", "error", "fit_status", "fit_boundary_limited",
        "n_catalogue", "n_detected",
    }
    missing = required - set(frame.columns)
    if missing:
        raise AssertionError(f"replicate table missing columns: {sorted(missing)}")
    if frame.empty or frame.duplicated(["replicate", "arm"]).any():
        raise AssertionError("replicate table is empty or has duplicate replicate/arm rows")
    if not set(frame["status"]).issubset({"success", "failed"}):
        raise AssertionError("unknown fit status in replicate table")

    # Reconstruct the per-replicate positivity/overlap diagnostics independently.
    positivity_columns = [
        "true_probability_min", "true_probability_max", "true_probability_p01",
        "true_probability_p05", "true_probability_below_005",
        "true_probability_below_010", "true_max_inverse_probability",
        "true_retention_rate", "audit_included_count",
    ]
    missing_positivity = set(positivity_columns) - set(frame.columns)
    if missing_positivity:
        raise AssertionError(f"replicate table missing positivity diagnostics: {sorted(missing_positivity)}")
    pmin = pd.to_numeric(frame["true_probability_min"], errors="raise").to_numpy(dtype=float)
    pmax = pd.to_numeric(frame["true_probability_max"], errors="raise").to_numpy(dtype=float)
    p01 = pd.to_numeric(frame["true_probability_p01"], errors="raise").to_numpy(dtype=float)
    p05 = pd.to_numeric(frame["true_probability_p05"], errors="raise").to_numpy(dtype=float)
    below05 = pd.to_numeric(frame["true_probability_below_005"], errors="raise").to_numpy(dtype=float)
    below10 = pd.to_numeric(frame["true_probability_below_010"], errors="raise").to_numpy(dtype=float)
    retention = pd.to_numeric(frame["true_retention_rate"], errors="raise").to_numpy(dtype=float)
    max_inverse = pd.to_numeric(frame["true_max_inverse_probability"], errors="raise").to_numpy(dtype=float)
    for label, values in (("true_probability_min", pmin), ("true_probability_max", pmax),
                          ("true_probability_p01", p01), ("true_probability_p05", p05),
                          ("true_probability_below_005", below05),
                          ("true_probability_below_010", below10),
                          ("true_retention_rate", retention),
                          ("true_max_inverse_probability", max_inverse)):
        if not np.all(np.isfinite(values)):
            raise AssertionError(f"{label} must be finite on every replicate row")
    if np.any((pmin <= 0) | (pmax >= 1) | (pmin > p01) | (p01 > p05) | (p05 > pmax)):
        raise AssertionError("true inclusion-probability extrema/quantiles are invalid")
    if np.any((below05 < 0) | (below05 > 1) | (below10 < 0) | (below10 > 1)
              | (retention < 0) | (retention > 1) | (below05 > below10)):
        raise AssertionError("positivity fractions or retention rates are outside valid bounds")
    if not np.allclose(max_inverse, 1.0 / pmin, rtol=1e-9, atol=1e-9):
        raise AssertionError("true_max_inverse_probability disagrees with true_probability_min")
    for replicate, group in frame.groupby("replicate"):
        for column in positivity_columns[:-1]:
            values = pd.to_numeric(group[column], errors="raise").to_numpy(dtype=float)
            if not np.allclose(values, values[0], rtol=1e-9, atol=1e-9):
                raise AssertionError(f"{column} varies across arms in replicate {replicate}")
        for _, row in group.iterrows():
            arm = str(row["arm"])
            count = row["audit_included_count"]
            if arm.startswith("estimated_"):
                m = int(arm.rsplit("_", 1)[1])
                if pd.isna(count) or float(count) != int(float(count)) or not 0 <= int(float(count)) <= m:
                    raise AssertionError(f"{arm}: invalid audit_included_count")
                if int(float(row["validation_size"])) != m:
                    raise AssertionError(f"{arm}: audit sample size differs from arm definition")
            elif not pd.isna(count):
                raise AssertionError(f"{arm}: non-estimated arm must not store an audit_included_count")
    success = frame["status"].eq("success")
    failed_errors = frame.loc[~success, "error"].fillna("").astype(str)
    if not failed_errors.str.len().gt(0).all():
        raise AssertionError("failed rows must preserve a nonempty error reason")
    numeric = [
        "posterior_mean", "posterior_median", "ci68_low", "ci68_high", "ci68_width",
        "absolute_error", "density_ise", "weight_ess", "weight_ess_fraction",
        "max_inverse_weight",
    ]
    if success.any():
        for column in [*numeric, "covered_68", "inclusion_probability_mae"]:
            if column not in frame.columns:
                raise AssertionError(f"replicate table missing column {column}")
    for column in numeric:
        if column in frame.columns:
            values = frame.loc[success, column].to_numpy(dtype=float)
            if values.size and not np.all(np.isfinite(values)):
                raise AssertionError(f"successful rows contain non-finite values in {column}")
    if success.any():
        if (frame.loc[success, "ci68_low"] > frame.loc[success, "ci68_high"]).any():
            raise AssertionError("credible interval endpoints are reversed")
        if (frame.loc[success, "absolute_error"] < 0).any():
            raise AssertionError("absolute posterior-mean errors cannot be negative")
    boundary_flags = np.array([_strict_bool(v, "fit_boundary_limited") for v in frame.loc[success, "fit_boundary_limited"]], dtype=bool)
    boundary_status = frame.loc[success, "fit_status"].eq("boundary_limited").to_numpy(dtype=bool)
    if not np.array_equal(boundary_flags, boundary_status):
        raise AssertionError("fit boundary flags disagree with fit_status")

    h0_true = _required_finite_number(summary.get("h0_true"), "summary.h0_true")
    requested_replicates = summary.get("n_replicates_requested")
    if isinstance(requested_replicates, bool) or not isinstance(requested_replicates, int) or requested_replicates < 1:
        raise AssertionError("summary.n_replicates_requested must be a positive integer")
    validation_sizes = summary.get("validation_sizes")
    if not isinstance(validation_sizes, list) or not validation_sizes:
        raise AssertionError("summary.validation_sizes must be a nonempty list")
    if any(isinstance(v, bool) or not isinstance(v, int) or v < 1 for v in validation_sizes):
        raise AssertionError("summary.validation_sizes must contain positive integers")

    # Cross-check saved population counts against the declared simulation design.
    n_parent = summary.get("n_parent")
    n_detected_declared = summary.get("n_detected")
    if isinstance(n_parent, bool) or not isinstance(n_parent, int) or n_parent < 100:
        raise AssertionError("summary.n_parent must be an integer >= 100")
    if isinstance(n_detected_declared, bool) or not isinstance(n_detected_declared, int) or n_detected_declared < 1:
        raise AssertionError("summary.n_detected must be a positive integer")
    catalogue_counts = _integer_values(frame["n_catalogue"], "n_catalogue")
    detected_counts = _integer_values(frame["n_detected"], "n_detected")
    if np.any((catalogue_counts < 0) | (catalogue_counts > n_parent)):
        raise AssertionError("n_catalogue must lie within the declared parent population")
    if not np.all(detected_counts == n_detected_declared):
        raise AssertionError("n_detected disagrees with the declared simulation design")
    if not np.allclose(retention, catalogue_counts / n_parent, rtol=1e-9, atol=1e-9):
        raise AssertionError("true_retention_rate disagrees with n_catalogue / n_parent")
    validation_counts = pd.to_numeric(frame["audit_included_count"], errors="coerce")
    if validation_counts.loc[frame["arm"].str.startswith("estimated_")].isna().any():
        raise AssertionError("estimated arms must preserve their validation inclusion count")
    estimated = frame.loc[frame["arm"].str.startswith("estimated_")].copy()
    if not estimated.empty:
        parsed_counts = _integer_values(
            estimated["audit_included_count"], "audit_included_count"
        )
        if np.any((parsed_counts < 0) | (parsed_counts > estimated["validation_size"].to_numpy(dtype=int))):
            raise AssertionError("audit_included_count must lie within its validation sample")
        estimated["_audit_count"] = parsed_counts
        for replicate, group in estimated.groupby("replicate"):
            ordered = group.sort_values("validation_size")
            if np.any(np.diff(ordered["_audit_count"].to_numpy(dtype=int)) < 0):
                raise AssertionError(
                    f"nested validation inclusion counts decrease with sample size in replicate {replicate}"
                )
    # Independently validate observable-only pre-posterior feature records.
    estimated_success = success & frame["arm"].str.startswith("estimated_")
    if estimated_success.any():
        required_features = ["diagnostic_schema", "n_validation", "n_validation_included", "validation_inclusion_fraction", "catalogue_outside_validation_range_fraction", "catalogue_near_validation_boundary_fraction", "fitted_probability_min", "fitted_probability_p01", "fitted_probability_p05", "fitted_probability_median", "fitted_probability_max", "inverse_weight_ess", "inverse_weight_sum", "inverse_weight_sum_squares", "inverse_weight_ess_fraction", "inverse_weight_entropy_ess_fraction", "inverse_weight_max_share", "inverse_weight_top_1pct_mass", "inverse_weight_top_5pct_mass", "inverse_weight_max", "inverse_weight_p95", "inverse_weight_p99"]
        missing_features = [name for name in required_features if name not in frame.columns]
        if missing_features:
            raise AssertionError(f"replicate table missing pre-posterior features: {missing_features}")
    if "diagnostic_schema" in frame.columns:
        if frame.loc[~frame["arm"].str.startswith("estimated_"), "diagnostic_schema"].notna().any():
            raise AssertionError("non-estimated arms must not contain pre-posterior diagnostic records")
        diagnostic_rows = frame.loc[frame["diagnostic_schema"].notna()].copy()
        if not diagnostic_rows.empty:
            if not diagnostic_rows["diagnostic_schema"].eq("preposterior_weight_stability_v1").all():
                raise AssertionError("unknown pre-posterior diagnostic schema")
            if diagnostic_rows["fit_status"].isna().any():
                raise AssertionError("pre-posterior fit status contains missing values")
            if not diagnostic_rows["fit_status"].eq(frame.loc[diagnostic_rows.index, "fit_status"]).all():
                raise AssertionError("pre-posterior fit status disagrees with replicate record")
            if not np.array_equal(diagnostic_rows["fit_boundary_limited"].astype(bool).to_numpy(), frame.loc[diagnostic_rows.index, "fit_boundary_limited"].astype(bool).to_numpy()):
                raise AssertionError("pre-posterior boundary flag disagrees with replicate record")
            if not np.array_equal(diagnostic_rows["n_catalogue"].to_numpy(dtype=int), frame.loc[diagnostic_rows.index, "n_catalogue"].to_numpy(dtype=int)):
                raise AssertionError("pre-posterior catalogue count disagrees with replicate record")
            numeric_features = ["validation_inclusion_fraction", "catalogue_outside_validation_range_fraction", "catalogue_near_validation_boundary_fraction", "fitted_probability_min", "fitted_probability_p01", "fitted_probability_p05", "fitted_probability_median", "fitted_probability_max", "inverse_weight_ess", "inverse_weight_sum", "inverse_weight_sum_squares", "inverse_weight_ess_fraction", "inverse_weight_entropy_ess_fraction", "inverse_weight_max_share", "inverse_weight_top_1pct_mass", "inverse_weight_top_5pct_mass", "inverse_weight_max", "inverse_weight_p95", "inverse_weight_p99"]
            for column in numeric_features:
                values = pd.to_numeric(diagnostic_rows[column], errors="coerce").to_numpy(dtype=float)
                if not np.all(np.isfinite(values)):
                    raise AssertionError(f"pre-posterior feature {column} contains non-finite values")
            n_validation_feature = _integer_values(diagnostic_rows["n_validation"], "n_validation")
            n_included_feature = _integer_values(diagnostic_rows["n_validation_included"], "n_validation_included")
            if not np.array_equal(n_validation_feature, diagnostic_rows["validation_size"].to_numpy(dtype=int)):
                raise AssertionError("pre-posterior validation size disagrees with arm design")
            if not np.array_equal(n_included_feature, diagnostic_rows["audit_included_count"].to_numpy(dtype=int)):
                raise AssertionError("pre-posterior validation inclusion count disagrees with audit count")
            expected_fraction = n_included_feature / n_validation_feature
            if not np.allclose(diagnostic_rows["validation_inclusion_fraction"], expected_fraction, rtol=1e-9, atol=1e-9):
                raise AssertionError("pre-posterior validation inclusion fraction is inconsistent")
            probability_columns = ["fitted_probability_min", "fitted_probability_p01", "fitted_probability_p05", "fitted_probability_median", "fitted_probability_max"]
            probability_matrix = diagnostic_rows[probability_columns].to_numpy(dtype=float)
            if np.any((probability_matrix <= 0.0) | (probability_matrix >= 1.0)) or np.any(np.diff(probability_matrix, axis=1) < -1e-12):
                raise AssertionError("pre-posterior fitted probability quantiles are invalid or unordered")
            support_gap = diagnostic_rows["catalogue_outside_validation_range_fraction"].to_numpy(dtype=float)
            if np.any((support_gap < 0.0) | (support_gap > 1.0)):
                raise AssertionError("pre-posterior catalogue support-gap fraction is outside [0, 1]")
            support_edge = diagnostic_rows["catalogue_near_validation_boundary_fraction"].to_numpy(dtype=float)
            if np.any((support_edge < 0.0) | (support_edge > 1.0)):
                raise AssertionError("pre-posterior validation-boundary fraction is outside [0, 1]")
            concentration = diagnostic_rows[["inverse_weight_entropy_ess_fraction", "inverse_weight_max_share", "inverse_weight_top_1pct_mass", "inverse_weight_top_5pct_mass"]].to_numpy(dtype=float)
            if np.any((concentration < 0.0) | (concentration > 1.0 + 1e-10)):
                raise AssertionError("pre-posterior weight concentration metrics are outside [0, 1]")
            if np.any(concentration[:, 1] > concentration[:, 2] + 1e-10) or np.any(concentration[:, 2] > concentration[:, 3] + 1e-10):
                raise AssertionError("pre-posterior top-tail weight mass is not nested")
            if np.any(concentration[:, 0] + 1e-10 < diagnostic_rows["inverse_weight_ess_fraction"].to_numpy(dtype=float)):
                raise AssertionError("entropy ESS fraction is below quadratic ESS fraction")
            ess = diagnostic_rows["inverse_weight_ess"].to_numpy(dtype=float)
            ess_fraction = diagnostic_rows["inverse_weight_ess_fraction"].to_numpy(dtype=float)
            catalogue_size = diagnostic_rows["n_catalogue"].to_numpy(dtype=float)
            if np.any((ess <= 0.0) | (ess > catalogue_size + 1e-9)):
                raise AssertionError("pre-posterior inverse-weight ESS is outside valid bounds")
            if not np.allclose(ess_fraction, ess / catalogue_size, rtol=1e-8, atol=1e-10):
                raise AssertionError("pre-posterior ESS fraction is inconsistent")
            if "weight_ess" in diagnostic_rows.columns:
                if not np.allclose(ess, diagnostic_rows["weight_ess"], rtol=1e-8, atol=1e-8):
                    raise AssertionError("pre-posterior ESS disagrees with saved estimator diagnostic")
            weight_sum = diagnostic_rows["inverse_weight_sum"].to_numpy(dtype=float)
            weight_sum_squares = diagnostic_rows["inverse_weight_sum_squares"].to_numpy(dtype=float)
            if np.any((weight_sum <= 0.0) | (weight_sum_squares <= 0.0)):
                raise AssertionError("pre-posterior inverse-weight aggregates must be positive")
            reconstructed_ess = (weight_sum / np.sqrt(weight_sum_squares)) ** 2
            if not np.allclose(ess, reconstructed_ess, rtol=1e-8, atol=1e-8):
                raise AssertionError("pre-posterior ESS does not reconstruct from weight aggregates")
            expected_max_share = diagnostic_rows["inverse_weight_max"].to_numpy(dtype=float) / weight_sum
            if not np.allclose(diagnostic_rows["inverse_weight_max_share"], expected_max_share, rtol=1e-8, atol=1e-10):
                raise AssertionError("pre-posterior maximum normalized weight share is inconsistent")
            if "max_inverse_weight" in diagnostic_rows.columns:
                if not np.allclose(diagnostic_rows["inverse_weight_max"], diagnostic_rows["max_inverse_weight"], rtol=1e-8, atol=1e-8):
                    raise AssertionError("pre-posterior maximum weight disagrees with saved estimator diagnostic")
            weight_quantiles = diagnostic_rows[["inverse_weight_max", "inverse_weight_p95", "inverse_weight_p99"]].to_numpy(dtype=float)
            if np.any(np.diff(weight_quantiles[:, 1:], axis=1) < -1e-12) or np.any(weight_quantiles[:, 1:] > weight_quantiles[:, :1] + 1e-9):
                raise AssertionError("pre-posterior inverse-weight quantiles are inconsistent")
    expected_arms = {"catalogue_only", "oracle"} | {
        f"estimated_{family}_{v}"
        for v in validation_sizes for family in ("linear", "quadratic")
    }
    arms_in_rows = set(frame["arm"])
    if arms_in_rows != expected_arms:
        raise AssertionError(f"replicate arm set differs from declared design: expected={sorted(expected_arms)}, observed={sorted(arms_in_rows)}")
    expected_replicates = set(range(requested_replicates))
    replicate_values = _integer_values(frame["replicate"], "replicate")
    observed_replicates = set(np.unique(replicate_values).tolist())
    if observed_replicates != expected_replicates:
        raise AssertionError(f"replicate IDs differ from declared design: expected={sorted(expected_replicates)}, observed={sorted(observed_replicates)}")
    expected_pairs = {(replicate, arm) for replicate in expected_replicates for arm in expected_arms}
    observed_pairs = set(zip(replicate_values.tolist(), frame["arm"]))
    if observed_pairs != expected_pairs:
        missing_pairs = sorted(expected_pairs - observed_pairs)[:5]
        extra_pairs = sorted(observed_pairs - expected_pairs)[:5]
        raise AssertionError(f"incomplete replicate-arm matrix; missing={missing_pairs}, extra={extra_pairs}")
    expected_validation_size = {"catalogue_only": 0, "oracle": 0}
    expected_validation_size.update({
        f"estimated_{family}_{v}": v
        for v in validation_sizes for family in ("linear", "quadratic")
    })
    for arm, expected_size in expected_validation_size.items():
        sizes = set(_integer_values(frame.loc[frame["arm"] == arm, "validation_size"], f"{arm}.validation_size").tolist())
        if sizes != {expected_size}:
            raise AssertionError(f"{arm}: validation_size does not match declared design")
    if not isinstance(summary.get("arm_summaries"), list):
        raise AssertionError("summary.arm_summaries must be a list")
    if any(not isinstance(item, dict) or "arm" not in item for item in summary["arm_summaries"]):
        raise AssertionError("each summary arm entry must be an object with an arm name")
    summary_arm_names = [item["arm"] for item in summary["arm_summaries"]]
    if len(summary_arm_names) != len(set(summary_arm_names)):
        raise AssertionError("summary contains duplicate arm summaries")
    arms_in_summary = set(summary_arm_names)
    if arms_in_rows != arms_in_summary:
        raise AssertionError("summary arm set does not match replicate table")
    for replicate, group in frame.groupby("replicate"):
        for seed_column in ("replicate_seed", "catalogue_seed", "audit_seed"):
            if group[seed_column].nunique(dropna=False) != 1:
                raise AssertionError(f"{seed_column} varies across arms in replicate {replicate}")
        if group["n_catalogue"].nunique(dropna=False) != 1 or group["n_detected"].nunique(dropna=False) != 1:
            raise AssertionError(f"catalogue/event counts differ across arms in replicate {replicate}")

    summary_by_arm = {item["arm"]: item for item in summary["arm_summaries"]}
    for arm in arms_in_rows:
        all_rows = frame.loc[frame["arm"] == arm]
        group = all_rows.loc[all_rows["status"] == "success"]
        stored = summary_by_arm[arm]
        _close(stored["n_attempted"], len(all_rows), f"{arm}.n_attempted")
        _close(stored["n_success"], len(group), f"{arm}.n_success")
        _close(stored["n_failed"], len(all_rows) - len(group), f"{arm}.n_failed")
        if group.empty:
            if stored.get("status") != "no_successful_fits":
                raise AssertionError(f"{arm}: zero-success arm must be marked no_successful_fits")
            continue
        for _, inference_row in group.iterrows():
            _audit_inference_payload(inference_row, h0_true)
        if stored.get("status") == "no_successful_fits":
            raise AssertionError(f"{arm}: summary marks an arm non-evaluable despite successful rows")
        posterior_mean = group["posterior_mean"].to_numpy(dtype=float)
        ci_low = group["ci68_low"].to_numpy(dtype=float)
        ci_high = group["ci68_high"].to_numpy(dtype=float)
        ci_width = group["ci68_width"].to_numpy(dtype=float)
        row_abs_error = group["absolute_error"].to_numpy(dtype=float)
        error = posterior_mean - h0_true
        if not np.allclose(ci_width, ci_high - ci_low, rtol=1e-9, atol=1e-9):
            raise AssertionError(f"{arm}: ci68_width disagrees with interval endpoints")
        if not np.allclose(row_abs_error, np.abs(error), rtol=1e-9, atol=1e-9):
            raise AssertionError(f"{arm}: absolute_error disagrees with posterior_mean and h0_true")
        if (group["density_ise"].to_numpy(dtype=float) < 0).any():
            raise AssertionError(f"{arm}: density_ise cannot be negative")
        ess = group["weight_ess"].to_numpy(dtype=float)
        ess_fraction = group["weight_ess_fraction"].to_numpy(dtype=float)
        n_catalogue = group["n_catalogue"].to_numpy(dtype=float)
        if (ess <= 0).any() or (ess > n_catalogue + 1e-9).any():
            raise AssertionError(f"{arm}: weight_ess must be positive and no greater than n_catalogue")
        if not np.allclose(ess_fraction, ess / n_catalogue, rtol=1e-9, atol=1e-9):
            raise AssertionError(f"{arm}: weight_ess_fraction disagrees with weight_ess / n_catalogue")
        if (group["max_inverse_weight"].to_numpy(dtype=float) < 1.0 - 1e-9).any():
            raise AssertionError(f"{arm}: max_inverse_weight cannot be below one")
        mae = group["inclusion_probability_mae"].to_numpy(dtype=float)
        if arm != "catalogue_only" and not np.all(np.isfinite(mae)):
            raise AssertionError(f"{arm}: inclusion_probability_mae must be finite for successful rows")
        if np.any(np.isfinite(mae) & (mae < 0)):
            raise AssertionError(f"{arm}: inclusion_probability_mae cannot be negative")
        covered = np.array([_strict_bool(v, f"{arm}.covered_68") for v in group["covered_68"]], dtype=bool)
        reconstructed_coverage = ((ci_low <= h0_true) & (h0_true <= ci_high))
        if not np.array_equal(covered, reconstructed_coverage):
            raise AssertionError(f"{arm}: covered_68 flags disagree with interval endpoints and h0_true")
        successes, trials = int(np.sum(covered)), int(len(covered))
        low, high = _wilson(successes, trials)
        reconstructed = {
            "mean_bias": float(np.mean(error)),
            "rmse": float(np.sqrt(np.mean(error**2))),
            "coverage_68": float(np.mean(covered)),
            "coverage_numerator": successes,
            "coverage_denominator": trials,
            "coverage_wilson_95_low": low,
            "coverage_wilson_95_high": high,
            "mean_interval_width_68": float(group["ci68_width"].mean()),
            "mean_density_ise": float(group["density_ise"].mean()),
            "mean_weight_ess_fraction": float(group["weight_ess_fraction"].mean()),
            "mean_max_inverse_weight": float(group["max_inverse_weight"].mean()),
        }
        for key, value in reconstructed.items():
            _close(stored.get(key), value, f"{arm}.{key}")
        mae_values = group["inclusion_probability_mae"].dropna().to_numpy(dtype=float)
        reconstructed_mae = None if mae_values.size == 0 else float(np.mean(mae_values))
        _close(stored.get("mean_inclusion_probability_mae"), reconstructed_mae,
               f"{arm}.mean_inclusion_probability_mae")

    successful = frame.loc[success]
    raw_contrasts = summary.get("contrasts", [])
    if not isinstance(raw_contrasts, list):
        raise AssertionError("summary.contrasts must be a list")
    if any(not isinstance(item, dict) or "comparison_arm" not in item for item in raw_contrasts):
        raise AssertionError("each contrast entry must be an object with comparison_arm")
    contrast_names = [item["comparison_arm"] for item in raw_contrasts]
    if len(contrast_names) != len(set(contrast_names)):
        raise AssertionError("summary contains duplicate comparison contrasts")
    stored_contrasts = {item["comparison_arm"]: item for item in raw_contrasts}
    expected_comparisons = [arm for arm in arms_in_rows if arm != "catalogue_only"]
    unexpected_contrasts = set(stored_contrasts) - set(expected_comparisons)
    if unexpected_contrasts:
        raise AssertionError(f"unexpected stored contrasts: {sorted(unexpected_contrasts)}")
    if successful.empty and stored_contrasts:
        raise AssertionError("all-failed artifact must not contain paired contrasts")
    baseline = (
        successful.loc[successful["arm"] == "catalogue_only"].set_index("replicate")
        if not successful.empty else None
    )
    for arm in expected_comparisons:
        if successful.empty:
            continue
        comparison = successful.loc[successful["arm"] == arm].set_index("replicate")
        paired = baseline[["absolute_error"]].merge(
            comparison[["absolute_error"]], left_index=True, right_index=True,
            how="inner", suffixes=("_baseline", "_comparison")
        )
        if paired.empty:
            if arm in stored_contrasts:
                raise AssertionError(f"{arm}: contrast exists but has no paired successful replicates")
            continue
        values = (
            paired["absolute_error_baseline"] - paired["absolute_error_comparison"]
        ).to_numpy(dtype=float)
        stored = stored_contrasts.get(arm)
        if stored is None:
            raise AssertionError(f"{arm}: missing paired contrast")
        if stored.get("reference_arm") != "catalogue_only":
            raise AssertionError(f"{arm}: contrast reference_arm must be catalogue_only")
        _close(stored.get("n_paired_replicates"), values.size, f"{arm}.contrast.n")
        mean = float(np.mean(values))
        _close(stored.get("mean_absolute_error_reduction"), mean, f"{arm}.contrast.mean")
        if values.size < 2:
            se = low = high = None
        else:
            se = float(np.std(values, ddof=1) / np.sqrt(values.size))
            low, high = mean - 1.96 * se, mean + 1.96 * se
        _close(stored.get("standard_error"), se, f"{arm}.contrast.se")
        _close(stored.get("normal_95_low"), low, f"{arm}.contrast.low")
        _close(stored.get("normal_95_high"), high, f"{arm}.contrast.high")

    # Independently reconstruct the primary paired linear-versus-quadratic
    # contrast at each validation size. Positive values favor the quadratic arm.
    raw_model_contrasts = summary.get("paired_model_contrasts")
    if not isinstance(raw_model_contrasts, list):
        raise AssertionError("summary.paired_model_contrasts must be a list")
    stored_by_size = {}
    for item in raw_model_contrasts:
        if not isinstance(item, dict) or "validation_size" not in item:
            raise AssertionError("each paired model contrast must identify validation_size")
        size = int(item["validation_size"])
        if size in stored_by_size:
            raise AssertionError("duplicate paired model contrast validation_size")
        stored_by_size[size] = item
    expected_model_contrast_sizes = set()
    for size in validation_sizes:
        linear_name = f"estimated_linear_{size}"
        quadratic_name = f"estimated_quadratic_{size}"
        linear = successful.loc[successful["arm"] == linear_name].set_index("replicate")
        quadratic = successful.loc[successful["arm"] == quadratic_name].set_index("replicate")
        paired_models = linear[["absolute_error"]].merge(
            quadratic[["absolute_error"]], left_index=True, right_index=True,
            how="inner", suffixes=("_linear", "_quadratic")
        )
        if paired_models.empty:
            continue
        expected_model_contrast_sizes.add(int(size))
        values = (
            paired_models["absolute_error_linear"]
            - paired_models["absolute_error_quadratic"]
        ).to_numpy(dtype=float)
        stored = stored_by_size.get(int(size))
        if stored is None:
            raise AssertionError(f"{size}: missing paired linear-versus-quadratic contrast")
        if (stored.get("reference_arm") != linear_name
                or stored.get("comparison_arm") != quadratic_name
                or stored.get("positive_favors_quadratic") is not True):
            raise AssertionError(f"{size}: paired model contrast labels are inconsistent")
        _close(stored.get("n_paired_replicates"), values.size, f"{size}.model_contrast.n")
        mean = float(np.mean(values))
        _close(stored.get("mean_absolute_error_reduction"), mean, f"{size}.model_contrast.mean")
        if values.size < 2:
            se = low = high = None
        else:
            se = float(np.std(values, ddof=1) / np.sqrt(values.size))
            low, high = mean - 1.96 * se, mean + 1.96 * se
        _close(stored.get("standard_error"), se, f"{size}.model_contrast.se")
        _close(stored.get("normal_95_low"), low, f"{size}.model_contrast.low")
        _close(stored.get("normal_95_high"), high, f"{size}.model_contrast.high")
    if set(stored_by_size) != expected_model_contrast_sizes:
        raise AssertionError("paired model contrast sizes do not match evaluable paired arms")

    # Independently reconstruct all paired secondary outcomes. The first
    # three compare each estimated arm with catalogue-only; inclusion-probability
    # MAE compares the two fitted model families because catalogue-only has no
    # estimated inclusion probabilities.
    raw_secondary = summary.get("paired_secondary_contrasts")
    if not isinstance(raw_secondary, list):
        raise AssertionError("summary.paired_secondary_contrasts must be a list")
    stored_secondary = {}
    for item in raw_secondary:
        if not isinstance(item, dict):
            raise AssertionError("each paired secondary contrast must be an object")
        key = (
            item.get("metric"), item.get("reference_arm"), item.get("comparison_arm")
        )
        if key in stored_secondary:
            raise AssertionError("duplicate paired secondary contrast")
        stored_secondary[key] = item

    expected_secondary = {}
    def _independent_secondary(metric, reference_arm, comparison_arm, column):
        reference = successful.loc[successful["arm"] == reference_arm].set_index("replicate")
        comparison = successful.loc[successful["arm"] == comparison_arm].set_index("replicate")
        paired_rows = reference[[column]].merge(
            comparison[[column]], left_index=True, right_index=True,
            how="inner", suffixes=("_reference", "_comparison"),
        )
        if paired_rows.empty:
            return
        values = (
            paired_rows[f"{column}_reference"] - paired_rows[f"{column}_comparison"]
        ).to_numpy(dtype=float)
        if not np.all(np.isfinite(values)):
            raise AssertionError(f"{metric}: paired values must be finite")
        expected_secondary[(metric, reference_arm, comparison_arm)] = values

    if not successful.empty:
        for size in validation_sizes:
            for family in ("linear", "quadratic"):
                arm = f"estimated_{family}_{size}"
                for metric, column in (
                    ("coverage_68", "covered_68"),
                    ("interval_width_68", "ci68_width"),
                    ("host_density_ise", "density_ise"),
                ):
                    _independent_secondary(metric, "catalogue_only", arm, column)
            _independent_secondary(
                "inclusion_probability_mae",
                f"estimated_linear_{size}", f"estimated_quadratic_{size}",
                "inclusion_probability_mae",
            )

    if set(stored_secondary) != set(expected_secondary):
        missing = set(expected_secondary) - set(stored_secondary)
        extra = set(stored_secondary) - set(expected_secondary)
        raise AssertionError(
            f"paired secondary contrast keys disagree; missing={sorted(missing)}, extra={sorted(extra)}"
        )
    for key, values in expected_secondary.items():
        metric, reference_arm, comparison_arm = key
        stored = stored_secondary[key]
        label = f"{metric}.{reference_arm}->{comparison_arm}"
        if (stored.get("metric") != metric
                or stored.get("reference_arm") != reference_arm
                or stored.get("comparison_arm") != comparison_arm
                or stored.get("difference_direction") != "reference_minus_comparison"):
            raise AssertionError(f"{label}: paired secondary contrast labels are inconsistent")
        _close(stored.get("n_paired_replicates"), values.size, f"{label}.n")
        mean = float(np.mean(values))
        _close(
            stored.get("mean_difference_reference_minus_comparison"),
            mean, f"{label}.mean",
        )
        if values.size < 2:
            se = low = high = None
        else:
            se = float(np.std(values, ddof=1) / np.sqrt(values.size))
            low, high = mean - 1.96 * se, mean + 1.96 * se
        _close(stored.get("standard_error"), se, f"{label}.se")
        _close(stored.get("normal_95_low"), low, f"{label}.low")
        _close(stored.get("normal_95_high"), high, f"{label}.high")

    boundary_success = np.array([
        _strict_bool(v, "fit_boundary_limited")
        for v in frame.loc[success, "fit_boundary_limited"]
    ], dtype=bool)
    per_arm_counts = {}
    for arm in sorted(arms_in_rows):
        arm_rows = frame.loc[frame["arm"] == arm]
        arm_success = arm_rows["status"].eq("success")
        per_arm_counts[arm] = {
            "n_attempted": int(len(arm_rows)),
            "n_success": int(arm_success.sum()),
            "n_failed": int((~arm_success).sum()),
            "failure_fraction": float((~arm_success).mean()),
            "status": "evaluable" if arm_success.any() else "non_evaluable_no_successful_fits",
        }
    n_boundary = int(boundary_success.sum())
    comparison_diagnostics = {}
    for arm in sorted(arms_in_rows - {"catalogue_only"}):
        baseline_all = frame.loc[frame["arm"] == "catalogue_only"].set_index("replicate")
        comparison_all = frame.loc[frame["arm"] == arm].set_index("replicate")
        paired_status = baseline_all[["status"]].merge(
            comparison_all[["status"]], left_index=True, right_index=True,
            how="inner", suffixes=("_baseline", "_comparison")
        )
        baseline_ok = paired_status["status_baseline"].eq("success")
        comparison_ok = paired_status["status_comparison"].eq("success")
        comparison_diagnostics[arm] = {
            "n_attempted_common_replicates": int(len(paired_status)),
            "n_paired_successful_replicates": int((baseline_ok & comparison_ok).sum()),
            "n_baseline_success_comparison_failure": int((baseline_ok & ~comparison_ok).sum()),
            "n_baseline_failure_comparison_success": int((~baseline_ok & comparison_ok).sum()),
            "n_both_failed": int((~baseline_ok & ~comparison_ok).sum()),
            "successful_pair_contrast_conditions_on_fit_success": True,
        }
    return {
        "status": "PASS",
        "n_rows": int(len(frame)),
        "n_success": int(success.sum()),
        "n_failed": int((~success).sum()),
        "n_arms": int(len(arms_in_rows)),
        "n_replicates_expected": int(requested_replicates),
        "n_replicate_arm_pairs_expected": int(len(expected_pairs)),
        "n_replicate_arm_pairs_observed": int(len(observed_pairs)),
        "failure_fraction": float((~success).mean()),
        "n_boundary_limited_successful_fits": n_boundary,
        "boundary_limited_fraction_of_successful_fits": (
            None if not success.any() else float(n_boundary / int(success.sum()))
        ),
        "per_arm_counts": per_arm_counts,
        "comparison_diagnostics": comparison_diagnostics,
        "audit_scope": "persisted-table integrity, truth-separated inclusion scoring, independent validation-logistic refit, reflected weighted-KDE reconstruction, direct-quadrature posterior reconstruction, density-ISE reconstruction, and summary reconstruction",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Independently audit estimated-inclusion pilot artifacts.")
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    report = audit_pilot_artifacts(args.output_dir)
    (args.output_dir / "independent_audit.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
