"""Independently audit persisted validation-bootstrap posterior propagation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import minimize
from scipy.special import expit, log_expit, logsumexp
from scipy.stats import norm

C_KM_S = 299_792.458


def _quad_weights(grid: np.ndarray) -> np.ndarray:
    delta = np.diff(grid)
    weights = np.empty_like(grid)
    weights[0] = delta[0] / 2.0
    weights[-1] = delta[-1] / 2.0
    weights[1:-1] = (delta[:-1] + delta[1:]) / 2.0
    return weights


def _posterior_from_density(
    z: np.ndarray, density: np.ndarray, h0: np.ndarray,
    distances: np.ndarray, sigma: float, threshold: float,
) -> tuple[float, float, float, float]:
    q = _quad_weights(z)
    valid = (density > 0) & (q > 0)
    log_mass = np.log(density[valid]) + np.log(q[valid])
    means = C_KM_S * z[valid][None, :] / h0[:, None]
    event_terms = norm.logpdf(
        distances[None, :, None],
        loc=means[:, None, :],
        scale=sigma,
    ) + log_mass[None, None, :]
    joint = np.sum(logsumexp(event_terms, axis=2), axis=1)
    detection = norm.cdf(
        (C_KM_S * z[None, :] / h0[:, None] - threshold) / sigma
    )
    beta = np.sum(density[None, :] * q[None, :] * detection, axis=1)
    if np.any(~np.isfinite(beta)) or np.any(beta <= 0) or np.any(beta > 1 + 1e-12):
        raise AssertionError("independent reconstruction found invalid detection efficiency")
    selected_loglike = joint - distances.size * np.log(beta)
    posterior = np.exp(selected_loglike - np.max(selected_loglike))
    posterior /= trapezoid(posterior, h0)
    cdf = np.r_[0.0, cumulative_trapezoid(posterior, h0)]
    cdf /= cdf[-1]
    low, median, high = (float(np.interp(p, cdf, h0)) for p in (0.16, 0.50, 0.84))
    mean = float(trapezoid(h0 * posterior, h0))
    return mean, median, low, high


def _independent_refit_density(
    grid: np.ndarray, observed: np.ndarray, validation_z: np.ndarray,
    validation_y: np.ndarray, merger_rate: np.ndarray,
    bandwidth_factor: float | None,
) -> tuple[np.ndarray, float, float, float, float, float, bool]:
    """Reconstruct the bounded-logit fit and reflected weighted KDE independently."""
    midpoint = 0.5 * (float(grid[0]) + float(grid[-1]))
    half_range = 0.5 * (float(grid[-1]) - float(grid[0]))
    x = (validation_z - midpoint) / half_range
    y = validation_y.astype(float)

    def objective(theta: np.ndarray) -> float:
        eta = theta[0] + theta[1] * x
        return float(-np.sum(y * log_expit(eta) + (1.0 - y) * log_expit(-eta)))

    fit = minimize(
        objective, np.zeros(2), method="L-BFGS-B",
        bounds=((-25.0, 25.0), (-50.0, 50.0)),
        options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-8},
    )
    if not fit.success or not np.all(np.isfinite(fit.x)):
        raise RuntimeError("independent logistic fit failed")
    intercept, slope = (float(v) for v in fit.x)
    probability = expit(intercept + slope * ((observed - midpoint) / half_range))
    inverse = 1.0 / probability
    scaled = inverse / np.max(inverse)
    ess = float(np.sum(scaled) ** 2 / np.sum(scaled ** 2))
    maximum = float(np.max(inverse))
    p99 = float(np.quantile(inverse, 0.99))

    weighted_mean = float(np.average(observed, weights=scaled))
    weighted_variance = float(np.average((observed - weighted_mean) ** 2, weights=scaled))
    if not np.isfinite(weighted_variance) or weighted_variance <= 0:
        raise ArithmeticError("independent weighted variance is invalid")
    factor = ess ** (-0.2) if bandwidth_factor is None else float(bandwidth_factor)
    bandwidth = factor * np.sqrt(weighted_variance)
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        raise ArithmeticError("independent KDE bandwidth is invalid")
    reflected = np.concatenate((observed, 2.0 * grid[0] - observed, 2.0 * grid[-1] - observed))
    weights = np.tile((scaled / np.sum(scaled)) / 3.0, 3)
    kernel = np.exp(
        -0.5 * ((grid[:, None] - reflected[None, :]) / bandwidth) ** 2
    ) @ weights
    galaxy_density = kernel / (np.sqrt(2.0 * np.pi) * bandwidth)
    raw = galaxy_density * merger_rate
    density = raw / trapezoid(raw, grid)
    boundary = bool(
        np.isclose(abs(intercept), 25.0, rtol=0.0, atol=1e-5)
        or np.isclose(abs(slope), 50.0, rtol=0.0, atol=1e-5)
    )
    return density, intercept, slope, ess, maximum, p99, boundary


def _assert_close(actual: np.ndarray, expected: np.ndarray, label: str) -> None:
    if actual.shape != expected.shape or not np.allclose(
        actual, expected, rtol=2e-7, atol=2e-8, equal_nan=False
    ):
        raise AssertionError(f"{label} does not match independent reconstruction")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--pilot-csv", type=Path, required=True)
    args = parser.parse_args()

    summary = json.loads((args.artifact_dir / "summary.json").read_text(encoding="utf-8"))
    inputs = json.loads((args.artifact_dir / "inputs.json").read_text(encoding="utf-8"))
    if summary.get("schema") != "validation_bootstrap_artifact_v1":
        raise AssertionError("unsupported bootstrap artifact schema")
    if inputs.get("schema") != "validation_bootstrap_inputs_v1":
        raise AssertionError("unsupported bootstrap input schema")
    if summary.get("finite_population_correction_applied") is not False:
        raise AssertionError("artifact must explicitly disclose absent finite-population correction")
    payload = inputs.get("inference_payload")
    if not isinstance(payload, dict):
        raise AssertionError("missing inference-only input payload")
    forbidden = {"parent_redshift", "true_inclusion_probability", "estimated_inclusion_probability"}
    if forbidden.intersection(payload):
        raise AssertionError("hidden-parent/scoring truth leaked into bootstrap inference inputs")

    pilot = pd.read_csv(args.pilot_csv)
    rows = pilot.loc[
        (pilot["replicate"] == int(inputs["replicate"]))
        & (pilot["arm"] == str(inputs["arm"]))
        & (pilot["status"] == "success")
    ]
    if len(rows) != 1:
        raise AssertionError("pilot CSV does not contain exactly one matching successful row")
    pilot_payload = json.loads(rows.iloc[0]["inference_audit_payload_json"])
    if pilot_payload != payload:
        raise AssertionError("persisted bootstrap inputs differ from the selected pilot inference payload")

    with np.load(args.artifact_dir / "bootstrap_draws.npz", allow_pickle=False) as data:
        arrays = {name: np.asarray(data[name]) for name in data.files}
    required = {
        "density_draws", "resample_indices", "success_mask",
        "intercept_draws", "slope_draws",
        "coefficient_boundary_limited_draws", "inverse_weight_ess_draws",
        "inverse_weight_maximum_draws", "inverse_weight_p99_draws",
        "posterior_mean_draws", "posterior_median_draws",
        "posterior_ci68_low_draws", "posterior_ci68_high_draws",
    }
    if set(arrays) != required:
        raise AssertionError("bootstrap NPZ has missing or unexpected arrays")

    n_attempted = int(summary["n_attempted"])
    n_successful = int(summary["n_successful"])
    n_failed = int(summary["n_failed"])
    if n_attempted < 2 or n_successful + n_failed != n_attempted:
        raise AssertionError("attempted/successful/failed bootstrap counts do not reconcile")
    if len(summary.get("failure_reasons", [])) != n_failed:
        raise AssertionError("failure reason count does not match failed resamples")
    if summary.get("seed") is None or not isinstance(summary["seed"], int):
        raise AssertionError("bootstrap seed is missing")
    z = np.asarray(payload["redshift_grid"], dtype=float)
    h0 = np.asarray(payload["h0_grid"], dtype=float)
    density_draws = arrays["density_draws"]
    resample_indices = arrays["resample_indices"]
    success_mask = arrays["success_mask"]
    validation_z = np.asarray(payload["validation_redshift"], dtype=float)
    validation_y = np.asarray(payload["validation_included"], dtype=int)
    merger_rate = np.asarray(payload["merger_weight_grid"], dtype=float)
    if resample_indices.shape != (n_attempted, validation_z.size):
        raise AssertionError("resample-index matrix has incorrect dimensions")
    if success_mask.shape != (n_attempted,) or int(np.sum(success_mask)) != n_successful:
        raise AssertionError("bootstrap success mask does not reconcile with counts")
    if np.any((resample_indices < 0) | (resample_indices >= validation_z.size)):
        raise AssertionError("bootstrap resample index lies outside validation records")
    replay_rng = np.random.default_rng(int(summary["seed"]))
    expected_indices = np.stack([
        replay_rng.integers(0, validation_z.size, size=validation_z.size)
        for _ in range(n_attempted)
    ], axis=0)
    if not np.array_equal(resample_indices, expected_indices):
        raise AssertionError("stored resample indices do not reproduce from the declared seed")
    if density_draws.shape != (n_successful, z.size):
        raise AssertionError("density draw dimensions disagree with declared success count")
    # Refit every attempted resample with a separate implementation. Failed
    # resamples must align exactly with the persisted success mask.
    reconstructed_density = []
    reconstructed_scalars = []
    bandwidth_factor = payload.get("bandwidth_factor")
    for attempt, indices in enumerate(resample_indices):
        try:
            rebuilt = _independent_refit_density(
                z, np.asarray(payload["catalogue_redshift"], dtype=float),
                validation_z[indices], validation_y[indices], merger_rate,
                bandwidth_factor,
            )
            if not success_mask[attempt]:
                raise AssertionError("stored failure mask marks a reproducibly successful refit as failed")
            reconstructed_density.append(rebuilt[0])
            reconstructed_scalars.append(rebuilt[1:])
        except (ValueError, RuntimeError, ArithmeticError, FloatingPointError):
            if success_mask[attempt]:
                raise AssertionError("stored success mask marks an independently failed refit as successful")
    if len(reconstructed_density) != n_successful:
        raise AssertionError("independent refit count disagrees with successful draw count")
    if n_successful:
        _assert_close(density_draws, np.stack(reconstructed_density), "bootstrap host densities")
        scalars = np.asarray(reconstructed_scalars, dtype=float)
        for name, column in (
            ("intercept_draws", 0), ("slope_draws", 1),
            ("inverse_weight_ess_draws", 2), ("inverse_weight_maximum_draws", 3),
            ("inverse_weight_p99_draws", 4),
        ):
            _assert_close(arrays[name], scalars[:, column], name)
        if not np.array_equal(
            arrays["coefficient_boundary_limited_draws"].astype(bool),
            scalars[:, 5].astype(bool),
        ):
            raise AssertionError("coefficient boundary flags do not match independent refits")
    if not np.all(np.isfinite(density_draws)) or np.any(density_draws < 0):
        raise AssertionError("density draws contain invalid values")
    if n_successful and not np.allclose(
        trapezoid(density_draws, z, axis=1), 1.0, rtol=1e-7, atol=1e-8
    ):
        raise AssertionError("one or more bootstrap densities fail unit normalization")

    scalar_draw_names = (
        "intercept_draws", "slope_draws", "coefficient_boundary_limited_draws",
        "inverse_weight_ess_draws", "inverse_weight_maximum_draws",
        "inverse_weight_p99_draws", "posterior_mean_draws",
        "posterior_median_draws", "posterior_ci68_low_draws",
        "posterior_ci68_high_draws",
    )
    for name in scalar_draw_names:
        if arrays[name].shape != (n_successful,):
            raise AssertionError(f"{name} length disagrees with successful draw count")
    for name in scalar_draw_names:
        if name != "coefficient_boundary_limited_draws" and not np.all(np.isfinite(arrays[name])):
            raise AssertionError(f"{name} contains non-finite values")
    if n_successful and (
        np.any(arrays["inverse_weight_ess_draws"] <= 0)
        or np.any(arrays["inverse_weight_maximum_draws"] <= 0)
        or np.any(arrays["inverse_weight_p99_draws"] <= 0)
        or np.any(arrays["posterior_ci68_low_draws"] >= arrays["posterior_ci68_high_draws"])
    ):
        raise AssertionError("positive-weight or interval-order invariant failed")

    reconstructed = [
        _posterior_from_density(
            z, density, h0,
            np.asarray(payload["observed_distances_mpc"], dtype=float),
            float(payload["sigma_mpc"]),
            float(payload["detection_threshold_mpc"]),
        )
        for density in density_draws
    ]
    if n_successful:
        values = np.asarray(reconstructed, dtype=float)
        _assert_close(arrays["posterior_mean_draws"], values[:, 0], "posterior means")
        _assert_close(arrays["posterior_median_draws"], values[:, 1], "posterior medians")
        _assert_close(arrays["posterior_ci68_low_draws"], values[:, 2], "posterior lower limits")
        _assert_close(arrays["posterior_ci68_high_draws"], values[:, 3], "posterior upper limits")
    expected_percentiles = (
        np.quantile(arrays["posterior_mean_draws"], [0.025, 0.975]).tolist()
        if n_successful >= 2 else None
    )
    stored_percentiles = summary.get("posterior_mean_bootstrap_percentile_95")
    if expected_percentiles is None:
        if stored_percentiles is not None:
            raise AssertionError("percentile interval should be null with fewer than two successful draws")
    elif not np.allclose(stored_percentiles, expected_percentiles, rtol=1e-10, atol=1e-10):
        raise AssertionError("stored posterior-mean bootstrap percentiles do not match draws")

    print(json.dumps({
        "status": "PASS",
        "replicate": int(inputs["replicate"]),
        "arm": inputs["arm"],
        "n_attempted": n_attempted,
        "n_successful": n_successful,
        "n_failed": n_failed,
        "posterior_draws_independently_reconstructed": n_successful,
        "scope": "independent validation-resample refit, weighted-KDE, and posterior-propagation reconstruction",
    }, indent=2))


if __name__ == "__main__":
    main()
