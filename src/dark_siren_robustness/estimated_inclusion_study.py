"""Non-confirmatory pilot for inclusion probabilities estimated from an external audit.

This is a controlled synthetic study only. It does not use real galaxy surveys
or gravitational-wave observations and is not an H0 measurement.
"""
from __future__ import annotations

from collections.abc import Sequence
import json

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.special import expit

from dark_siren_robustness.catalogue_density import estimate_host_density
from dark_siren_robustness.estimated_inclusion import (
    bootstrap_validation_inclusion_density,
    estimate_host_density_from_validation,
    predict_inclusion_probability,
)
from dark_siren_robustness.generative_selection import (
    generate_synthetic_population,
    latent_galaxy_intensity,
)
from dark_siren_robustness.host_density_benchmark import posterior_summary
from dark_siren_robustness.preposterior_diagnostics import build_preposterior_diagnostics
from dark_siren_robustness.selection_toy import selection_aware_toy_likelihood


def _wilson(successes: int, trials: int) -> tuple[float, float]:
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("invalid binomial counts")
    z = 1.959963984540054
    p = successes / trials
    denom = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / denom
    half = z * np.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / denom
    return max(0.0, float(center - half)), min(1.0, float(center + half))


def _contrast_summary(values: np.ndarray) -> dict[str, float | int | None]:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("contrast values must be a nonempty finite vector")
    mean = float(np.mean(values))
    if values.size < 2:
        se = low = high = None
    else:
        se = float(np.std(values, ddof=1) / np.sqrt(values.size))
        low, high = mean - 1.96 * se, mean + 1.96 * se
    return {
        "n_paired_replicates": int(values.size),
        "mean_absolute_error_reduction": mean,
        "standard_error": se,
        "normal_95_low": low,
        "normal_95_high": high,
        "positive_favors_comparison_arm": True,
    }


def run_estimated_inclusion_pilot(
    *,
    n_replicates: int = 10,
    n_parent: int = 4000,
    n_detected: int = 30,
    validation_sizes: Sequence[int] = (100, 300, 1000),
    seed: int = 20261030,
    h0_true: float = 70.0,
    sigma_mpc: float = 140.0,
    detection_threshold_mpc: float = 250.0,
    h0_min: float = 50.0,
    h0_max: float = 90.0,
    h0_points: int = 81,
    z_points: int = 401,
    true_intercept: float = 0.25,
    true_slope: float = -1.2,
    true_quadratic: float = 0.0,
    bandwidth_factor: float | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Compare catalogue-only, oracle, and audit-estimated inclusion weighting.

    The validation sample is a simple random sample of the parent frame,
    independent of catalogue retention. Validation samples are nested by size
    within replicate. Non-oracle fitting sees only validation redshifts and
    retention indicators; true probabilities are used only for scoring.

    The pilot records failed arms instead of silently dropping them. Its small
    default replicate count is for runtime/variance planning, not confirmatory
    scientific claims.
    """
    raw_sizes = tuple(validation_sizes)
    if any(
        isinstance(m, (bool, np.bool_)) or not isinstance(m, (int, np.integer))
        for m in raw_sizes
    ):
        raise ValueError("validation sizes must be integer counts, not booleans or floats")
    sizes = tuple(sorted(set(int(m) for m in raw_sizes)))
    if n_replicates < 1 or n_parent < 100 or n_detected < 1:
        raise ValueError("n_replicates/n_detected must be positive and n_parent >= 100")
    if not sizes or any(m < 10 for m in sizes) or max(sizes) > n_parent:
        raise ValueError("validation_sizes must contain integers in [10, n_parent]")
    if h0_points < 31 or z_points < 101 or not (0 < h0_min < h0_max):
        raise ValueError("invalid H0 or redshift grid settings")
    if not np.isfinite(h0_true) or not h0_min < h0_true < h0_max:
        raise ValueError("h0_true must be finite and inside the H0 grid")
    if not all(np.isfinite(value) for value in (true_intercept, true_slope, true_quadratic)):
        raise ValueError("true inclusion-model coefficients must be finite")
    if bandwidth_factor is not None and (
        not np.isfinite(bandwidth_factor) or bandwidth_factor <= 0
    ):
        raise ValueError("bandwidth_factor must be finite and positive")

    z_grid = np.linspace(0.01, 0.40, z_points)
    h0_grid = np.linspace(h0_min, h0_max, h0_points)
    z_min, z_max = float(z_grid[0]), float(z_grid[-1])
    z_mid, z_half = (z_min + z_max) / 2.0, (z_max - z_min) / 2.0
    standardized_grid = (z_grid - z_mid) / z_half
    merger_rate = 1.0 + 1.5 * (z_grid - z_min) / (z_max - z_min)
    target_density = latent_galaxy_intensity(z_grid) * merger_rate
    target_density = target_density / trapezoid(target_density, z_grid)

    rows: list[dict] = []
    for replicate in range(n_replicates):
        seed_sequence = np.random.SeedSequence([int(seed), int(replicate)])
        sim_ss, catalogue_ss, audit_ss = seed_sequence.spawn(3)
        sim_seed = int(sim_ss.generate_state(1, dtype=np.uint32)[0])
        catalogue_seed = int(catalogue_ss.generate_state(1, dtype=np.uint32)[0])
        audit_seed = int(audit_ss.generate_state(1, dtype=np.uint32)[0])

        sim = generate_synthetic_population(
            seed=sim_seed,
            n_parent=n_parent,
            n_detected=n_detected,
            h0_true=h0_true,
            sigma_mpc=sigma_mpc,
            detection_threshold_mpc=detection_threshold_mpc,
            z_min=z_min,
            z_max=z_max,
        )
        parent_z = sim["parent_truth"]["redshift"].to_numpy(dtype=float)
        distances = sim["gw_events_observed"]["distance_mpc"].to_numpy(dtype=float)

        # Use a prespecified logistic law so the fitted family is correctly
        # specified in the primary arms. The generator's default catalogue
        # table is deliberately ignored; this experiment uses its own thinning.
        standardized_parent = (parent_z - z_mid) / z_half
        # A nonzero quadratic term creates a deliberately misspecified stress
        # condition for the linear-logit validation estimator. Subtracting 1/3
        # centers x^2 under a *uniform reference measure* on [-1, 1] in logit
        # space. It does not make the intercept the midpoint log-odds (at x=0,
        # eta = intercept - quadratic/3), nor does it generally center under
        # the non-uniform simulated parent-redshift distribution.
        true_probability = expit(
            true_intercept + true_slope * standardized_parent
            + true_quadratic * (standardized_parent**2 - 1.0 / 3.0)
        )
        if (
            not np.all(np.isfinite(true_probability))
            or np.any((true_probability <= 0.0) | (true_probability >= 1.0))
        ):
            raise ValueError(
                "true inclusion probabilities numerically saturated at 0 or 1; "
                "use less extreme true selection coefficients"
            )
        catalogue_rng = np.random.default_rng(catalogue_seed)
        retained = catalogue_rng.random(n_parent) < true_probability
        observed_z = parent_z[retained]
        # Record overlap/positivity diagnostics for every replicate, including
        # those that fail before any posterior can be constructed.
        true_probability_min = float(np.min(true_probability))
        true_probability_max = float(np.max(true_probability))
        true_probability_p01 = float(np.quantile(true_probability, 0.01))
        true_probability_p05 = float(np.quantile(true_probability, 0.05))
        true_probability_below_005 = float(np.mean(true_probability < 0.05))
        true_probability_below_010 = float(np.mean(true_probability < 0.10))
        true_max_inverse_probability = float(1.0 / true_probability_min)
        true_retention_rate = float(np.mean(retained))
        audit_rng = np.random.default_rng(audit_seed)
        audit_order = audit_rng.permutation(n_parent)
        audit_records = {
            "redshift": parent_z[audit_order],
            "included": retained[audit_order].astype(int),
        }
        if observed_z.size < 5 or np.unique(observed_z).size < 3:
            # Preserve the replicate failure explicitly for all arms.
            arms = ["catalogue_only", "oracle"] + [f"estimated_{m}" for m in sizes]
            for arm in arms:
                rows.append({
                    "replicate": replicate, "replicate_seed": sim_seed,
                    "catalogue_seed": catalogue_seed, "audit_seed": audit_seed,
                    "arm": arm, "validation_size": 0 if arm in ("catalogue_only", "oracle") else int(arm.split("_")[1]),
                    "status": "failed", "error": "too few retained catalogue objects",
                    # No inclusion model or pre-posterior features were fitted:
                    # the catalogue failed the minimum-size gate before arm setup.
                    "fit_status": "not_attempted",
                    "fit_boundary_limited": False,
                    "n_catalogue": int(observed_z.size), "n_detected": int(distances.size),
                    "true_probability_min": true_probability_min,
                    "true_probability_max": true_probability_max,
                    "true_probability_p01": true_probability_p01,
                    "true_probability_p05": true_probability_p05,
                    "true_probability_below_005": true_probability_below_005,
                    "true_probability_below_010": true_probability_below_010,
                    "true_max_inverse_probability": true_max_inverse_probability,
                    "true_retention_rate": true_retention_rate,
                    "audit_included_count": (
                        int(np.sum(audit_records["included"][:int(arm.split("_")[1])]))
                        if arm.startswith("estimated_") else None
                    ),
                })
            continue

        arm_specs: list[tuple[str, int, object]] = [
            ("catalogue_only", 0, None),
            ("oracle", 0, "oracle"),
        ]
        arm_specs.extend((f"estimated_{m}", m, m) for m in sizes)
        for arm, audit_size, spec in arm_specs:
            try:
                fit_status = "not_applicable"
                fit_boundary_limited = False
                preposterior_features: dict[str, object] = {}
                density_inclusion_probability = np.ones(observed_z.size, dtype=float)
                audit_validation_redshift = np.array([], dtype=float)
                audit_validation_included = np.array([], dtype=int)
                fitted_inclusion_model = None
                scoring_audit_payload: dict[str, object] = {
                    "schema": "inclusion_probability_scoring_only_v1",
                    "target_density_for_ise": target_density.tolist(),
                }
                weight_ess = float("nan")
                max_inverse_weight = float("nan")
                if arm == "catalogue_only":
                    density = estimate_host_density(
                        z_grid, observed_z, merger_weight=merger_rate,
                        bandwidth_factor=bandwidth_factor,
                    )
                    probability_mae = None
                    weight_ess = float(observed_z.size)
                    max_inverse_weight = 1.0
                elif arm == "oracle":
                    fit_status = "oracle"
                    observed_p = true_probability[retained]
                    density_inclusion_probability = observed_p
                    density = estimate_host_density(
                        z_grid, observed_z, inclusion_probability=observed_p,
                        merger_weight=merger_rate, bandwidth_factor=bandwidth_factor,
                    )
                    scaled_w = (1.0 / observed_p) / np.max(1.0 / observed_p)
                    weight_ess = float(np.sum(scaled_w) ** 2 / np.sum(scaled_w**2))
                    max_inverse_weight = float(np.max(1.0 / observed_p))
                    probability_mae = 0.0
                else:
                    m = int(spec)
                    validation_z = audit_records["redshift"][:m]
                    validation_included = audit_records["included"][:m]
                    audit_validation_redshift = np.asarray(validation_z, dtype=float)
                    audit_validation_included = np.asarray(validation_included, dtype=int)
                    result = estimate_host_density_from_validation(
                        z_grid, observed_z, validation_z, validation_included,
                        merger_weight=merger_rate, bandwidth_factor=bandwidth_factor,
                    )
                    density = result["density"]
                    density_inclusion_probability = np.asarray(
                        result["observed_inclusion_probability"], dtype=float
                    )
                    fitted_inclusion_model = dict(result["model"])
                    fit_status = str(result["model"]["fit_status"])
                    fit_boundary_limited = bool(result["model"]["coefficient_bound_hit"])
                    weight_ess = float(result["inverse_weight_effective_sample_size"])
                    max_inverse_weight = float(result["inverse_weight_maximum"])
                    # Score probability estimation over the full synthetic parent
                    # population, not only retained objects. Parent truth is used
                    # strictly for simulation evaluation, never for fitting/weighting.
                    estimated_parent_probability = predict_inclusion_probability(
                        parent_z, result["model"]
                    )
                    probability_mae = float(np.mean(np.abs(
                        estimated_parent_probability - true_probability
                    )))
                    # Keep scoring truth in a distinct, explicitly scoring-only
                    # payload. It is never passed to the estimator or posterior.
                    scoring_audit_payload.update({
                        "parent_redshift": parent_z.tolist(),
                        "true_inclusion_probability": true_probability.tolist(),
                        "estimated_inclusion_probability": estimated_parent_probability.tolist(),
                    })
                    # Freeze observable-only warning features before constructing
                    # the GW likelihood or inspecting any posterior outcome.
                    preposterior_features = build_preposterior_diagnostics(
                        result["model"], validation_z, validation_included,
                        observed_z, result["observed_inclusion_probability"],
                    )

                likelihood = selection_aware_toy_likelihood(
                    z_grid, density, h0_grid, distances, sigma=sigma_mpc,
                    detection_threshold=detection_threshold_mpc,
                )
                post = posterior_summary(h0_grid, likelihood["selected_log_likelihood"])
                rows.append({
                    "replicate": replicate, "replicate_seed": sim_seed,
                    "catalogue_seed": catalogue_seed, "audit_seed": audit_seed,
                    "arm": arm, "validation_size": int(audit_size),
                    "status": "success", "error": "",
                    "fit_status": fit_status, "fit_boundary_limited": fit_boundary_limited,
                    "n_catalogue": int(observed_z.size), "n_detected": int(distances.size),
                    "true_probability_min": true_probability_min,
                    "true_probability_max": true_probability_max,
                    "true_probability_p01": true_probability_p01,
                    "true_probability_p05": true_probability_p05,
                    "true_probability_below_005": true_probability_below_005,
                    "true_probability_below_010": true_probability_below_010,
                    "true_max_inverse_probability": true_max_inverse_probability,
                    "true_retention_rate": true_retention_rate,
                    "audit_included_count": (
                        int(np.sum(audit_records["included"][:audit_size]))
                        if audit_size else None
                    ),
                    "posterior_mean": post["posterior_mean"],
                    "posterior_median": post["posterior_median"],
                    "ci68_low": post["ci68_low"], "ci68_high": post["ci68_high"],
                    "ci68_width": post["ci68_width"],
                    "covered_68": bool(post["ci68_low"] <= h0_true <= post["ci68_high"]),
                    "absolute_error": float(abs(post["posterior_mean"] - h0_true)),
                    "density_ise": float(trapezoid((density - target_density) ** 2, z_grid)),
                    # Persist the numerical inputs needed for a genuinely
                    # independent posterior/selection reconstruction. This
                    # payload contains no hidden parent truth or inclusion
                    # probabilities; target_density is included only to audit
                    # the reported density ISE.
                    "inference_audit_schema": "toy_selection_posterior_reconstruction_v1",
                    "inference_audit_payload_json": json.dumps({
                        "redshift_grid": z_grid.tolist(),
                        "host_density": np.asarray(density, dtype=float).tolist(),
                        "catalogue_redshift": observed_z.tolist(),
                        "density_inclusion_probability": density_inclusion_probability.tolist(),
                        "merger_weight_grid": merger_rate.tolist(),
                        "bandwidth_factor": (
                            None if bandwidth_factor is None else float(bandwidth_factor)
                        ),
                        "validation_redshift": audit_validation_redshift.tolist(),
                        "validation_included": audit_validation_included.tolist(),
                        "fitted_inclusion_model": fitted_inclusion_model,
                        "h0_grid": h0_grid.tolist(),
                        "observed_distances_mpc": distances.tolist(),
                        "sigma_mpc": float(sigma_mpc),
                        "detection_threshold_mpc": float(detection_threshold_mpc),
                    }, separators=(",", ":")),
                    "scoring_audit_payload_json": json.dumps(
                        scoring_audit_payload, separators=(",", ":")
                    ),
                    "weight_ess": weight_ess,
                    "weight_ess_fraction": weight_ess / observed_z.size,
                    "max_inverse_weight": max_inverse_weight,
                    "inclusion_probability_mae": probability_mae,
                    **preposterior_features,
                })
            except Exception as exc:  # failed fits are retained as explicit records
                rows.append({
                    "replicate": replicate, "replicate_seed": sim_seed,
                    "catalogue_seed": catalogue_seed, "audit_seed": audit_seed,
                    "arm": arm, "validation_size": int(audit_size),
                    "status": "failed", "error": f"{type(exc).__name__}: {exc}",
                    "fit_status": "failed", "fit_boundary_limited": False,
                    "n_catalogue": int(observed_z.size), "n_detected": int(distances.size),
                    "true_probability_min": true_probability_min,
                    "true_probability_max": true_probability_max,
                    "true_probability_p01": true_probability_p01,
                    "true_probability_p05": true_probability_p05,
                    "true_probability_below_005": true_probability_below_005,
                    "true_probability_below_010": true_probability_below_010,
                    "true_max_inverse_probability": true_max_inverse_probability,
                    "true_retention_rate": true_retention_rate,
                    "audit_included_count": (
                        int(np.sum(audit_records["included"][:audit_size]))
                        if audit_size else None
                    ),
                    "weight_ess": weight_ess,
                    "weight_ess_fraction": (
                        weight_ess / observed_z.size if np.isfinite(weight_ess) else np.nan
                    ),
                    "max_inverse_weight": max_inverse_weight,
                    **preposterior_features,
                })

    frame = pd.DataFrame(rows)
    summaries: list[dict] = []
    successful = frame.loc[frame["status"] == "success"].copy()
    for arm in frame["arm"].drop_duplicates():
        group = successful.loc[successful["arm"] == arm]
        n_success = int(len(group))
        if not n_success:
            summaries.append({
                "arm": arm,
                "n_attempted": int((frame["arm"] == arm).sum()),
                "n_success": 0,
                "n_failed": int((frame["arm"] == arm).sum()),
                "status": "no_successful_fits",
            })
            continue
        covered = group["covered_68"].astype(bool).to_numpy()
        n_covered = int(np.sum(covered))
        w_low, w_high = _wilson(n_covered, n_success)
        if n_success:
            summaries.append({
                "arm": arm,
                "n_attempted": int((frame["arm"] == arm).sum()),
                "n_success": n_success,
                "n_failed": int(((frame["arm"] == arm) & (frame["status"] != "success")).sum()),
                "mean_bias": float(np.mean(group["posterior_mean"] - h0_true)),
                "rmse": float(np.sqrt(np.mean((group["posterior_mean"] - h0_true) ** 2))),
                "coverage_68": float(np.mean(covered)),
                "coverage_numerator": n_covered,
                "coverage_denominator": n_success,
                "coverage_wilson_95_low": w_low,
                "coverage_wilson_95_high": w_high,
                "mean_interval_width_68": float(group["ci68_width"].mean()),
                "mean_density_ise": float(group["density_ise"].mean()),
                "mean_weight_ess_fraction": float(group["weight_ess_fraction"].mean()),
                "mean_max_inverse_weight": float(group["max_inverse_weight"].mean()),
                "mean_inclusion_probability_mae": (
                    None if group["inclusion_probability_mae"].isna().all()
                    else float(group["inclusion_probability_mae"].dropna().mean())
                ),
            })
        else:
            summaries.append({
                "arm": arm, "n_attempted": int((frame["arm"] == arm).sum()),
                "n_success": 0, "n_failed": int((frame["arm"] == arm).sum()),
                "status": "no_successful_fits",
            })

    contrasts: list[dict] = []
    if not successful.empty:
        baseline = successful.loc[successful["arm"] == "catalogue_only"].set_index("replicate")
        for arm in ["oracle"] + [f"estimated_{m}" for m in sizes]:
            comparison = successful.loc[successful["arm"] == arm].set_index("replicate")
            paired = baseline[["absolute_error"]].merge(
                comparison[["absolute_error"]], left_index=True, right_index=True,
                how="inner", suffixes=("_baseline", "_comparison")
            )
            values = (
                paired["absolute_error_baseline"] - paired["absolute_error_comparison"]
            ).to_numpy(dtype=float)
            if values.size:
                item = _contrast_summary(values)
                item.update({"comparison_arm": arm, "reference_arm": "catalogue_only"})
                contrasts.append(item)

    summary = {
        "study": "estimated_catalogue_inclusion_probability_pilot",
        "status": "NON_CONFIRMATORY_SYNTHETIC_PILOT",
        "seed": int(seed), "n_replicates_requested": int(n_replicates),
        "n_parent": int(n_parent), "n_detected": int(n_detected),
        "h0_true": float(h0_true), "sigma_mpc": float(sigma_mpc),
        "detection_threshold_mpc": float(detection_threshold_mpc),
        "h0_grid": [float(h0_min), float(h0_max), int(h0_points)],
        "redshift_grid": [float(z_min), float(z_max), int(z_points)],
        "bandwidth_factor": None if bandwidth_factor is None else float(bandwidth_factor),
        "validation_sizes": list(sizes),
        "true_inclusion_model": {
            "family": "quadratic_logistic" if true_quadratic != 0.0 else "linear_logistic",
            "logit_probability": "intercept + slope*x + quadratic*(x^2 - 1/3)",
            "intercept": float(true_intercept), "slope": float(true_slope),
            "quadratic": float(true_quadratic),
            "standardization": "(z - domain_midpoint) / domain_half_range",
            "quadratic_centering": (
                "x^2 - 1/3; centered in logit space under uniform x on [-1, 1] "
                "only, not necessarily under the simulated parent distribution"
            ),
        },
        "primary_contrast": "abs_error_catalogue_only - abs_error_comparison_arm",
        "contrasts": contrasts, "arm_summaries": summaries,
        "interpretation_limit": (
            "Controlled synthetic toy only; not a real galaxy-survey or LVK validation. "
            "Pilot uncertainty is not confirmatory evidence."
        ),
    }
    return frame, summary


def run_validation_bootstrap_posterior_sensitivity(
    *,
    redshift_grid: np.ndarray,
    observed_redshift: np.ndarray,
    validation_redshift: np.ndarray,
    validation_included: np.ndarray,
    h0_grid: np.ndarray,
    observed_gw_distances_mpc: np.ndarray,
    sigma_mpc: float,
    detection_threshold_mpc: float,
    n_bootstrap: int = 200,
    seed: int = 20261030,
    merger_weight: np.ndarray | None = None,
    bandwidth_factor: float | None = None,
) -> dict:
    """Propagate validation-record bootstrap draws through the fixed toy GW likelihood.

    The observed catalogue and GW-event distances remain fixed. Each successful
    validation bootstrap refits inclusion probabilities, rebuilds the host
    density, and recomputes the selection-conditioned H0 posterior. Percentile
    ranges across posterior summaries quantify validation-sample sensitivity;
    they are not posterior credible intervals.
    """
    z_grid = np.asarray(redshift_grid, dtype=float)
    h0_values = np.asarray(h0_grid, dtype=float)
    distances = np.asarray(observed_gw_distances_mpc, dtype=float)
    if (
        h0_values.ndim != 1 or h0_values.size < 31
        or not np.all(np.isfinite(h0_values)) or np.any(np.diff(h0_values) <= 0)
        or np.any(h0_values <= 0)
    ):
        raise ValueError("h0_grid must be finite, positive, strictly increasing, and have >=31 points")
    if (
        distances.ndim != 1 or distances.size < 1
        or not np.all(np.isfinite(distances)) or np.any(distances <= 0)
    ):
        raise ValueError("observed_gw_distances_mpc must be finite and positive")
    if not np.isfinite(sigma_mpc) or sigma_mpc <= 0:
        raise ValueError("sigma_mpc must be finite and positive")
    if not np.isfinite(detection_threshold_mpc) or detection_threshold_mpc <= 0:
        raise ValueError("detection_threshold_mpc must be finite and positive")
    if np.any(distances < detection_threshold_mpc):
        raise ValueError("all supplied GW distances must satisfy the detection threshold")

    bootstrap = bootstrap_validation_inclusion_density(
        z_grid, observed_redshift, validation_redshift, validation_included,
        n_bootstrap=n_bootstrap, seed=seed,
        merger_weight=merger_weight, bandwidth_factor=bandwidth_factor,
    )
    posterior_means: list[float] = []
    posterior_medians: list[float] = []
    credible_lows: list[float] = []
    credible_highs: list[float] = []
    for density in bootstrap["density_draws"]:
        likelihood = selection_aware_toy_likelihood(
            z_grid, density, h0_values, distances,
            sigma=float(sigma_mpc),
            detection_threshold=float(detection_threshold_mpc),
        )
        posterior = posterior_summary(h0_values, likelihood["selected_log_likelihood"])
        posterior_means.append(float(posterior["posterior_mean"]))
        posterior_medians.append(float(posterior["posterior_median"]))
        credible_lows.append(float(posterior["ci68_low"]))
        credible_highs.append(float(posterior["ci68_high"]))

    means = np.asarray(posterior_means, dtype=float)
    if means.size >= 2:
        mean_interval = np.quantile(means, [0.025, 0.975]).tolist()
    else:
        mean_interval = None
    return {
        "schema": "validation_bootstrap_posterior_sensitivity_v1",
        "seed": int(seed),
        "n_attempted": int(bootstrap["n_attempted"]),
        "n_successful": int(bootstrap["n_successful"]),
        "n_failed": int(bootstrap["n_failed"]),
        "posterior_mean_draws": means,
        "posterior_median_draws": np.asarray(posterior_medians, dtype=float),
        "posterior_ci68_low_draws": np.asarray(credible_lows, dtype=float),
        "posterior_ci68_high_draws": np.asarray(credible_highs, dtype=float),
        "posterior_mean_bootstrap_percentile_95": mean_interval,
        "density_draws": bootstrap["density_draws"],
        "intercept_draws": bootstrap["intercept_draws"],
        "slope_draws": bootstrap["slope_draws"],
        "coefficient_boundary_limited_draws": bootstrap[
            "coefficient_boundary_limited_draws"
        ],
        "inverse_weight_ess_draws": bootstrap["inverse_weight_ess_draws"],
        "inverse_weight_maximum_draws": bootstrap["inverse_weight_maximum_draws"],
        "inverse_weight_p99_draws": bootstrap["inverse_weight_p99_draws"],
        "failure_reasons": bootstrap["failure_reasons"],
        "interpretation": (
            "validation-bootstrap sensitivity of posterior summaries conditional on "
            "the fixed catalogue and GW distances; percentile range is not an H0 "
            "posterior credible interval or a full uncertainty interval"
        ),
    }
