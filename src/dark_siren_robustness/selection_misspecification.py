"""Selection-function misspecification experiment for the synthetic host-density toy.

All results are synthetic. The true catalogue selection law generates the
observed catalogue; analysts receive only the assumed law assigned to their arm.
"""
from __future__ import annotations

import json
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid

from dark_siren_robustness.catalogue_density import estimate_host_density
from dark_siren_robustness.generative_selection import (
    catalogue_inclusion_probability,
    generate_synthetic_population,
    latent_galaxy_intensity,
)
from dark_siren_robustness.host_density_benchmark import posterior_summary
from dark_siren_robustness.selection_toy import selection_aware_toy_likelihood


def _normalized(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    area = float(trapezoid(values, grid))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("Target density normalization failed")
    return values / area


def _wilson_interval(successes: int, trials: int, z_value: float = 1.959963984540054) -> tuple[float, float]:
    """95% Wilson score interval for a binomial proportion."""
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("invalid binomial count")
    p = successes / trials
    z2 = z_value * z_value
    denominator = 1.0 + z2 / trials
    center = (p + z2 / (2.0 * trials)) / denominator
    half = z_value * np.sqrt(p * (1.0 - p) / trials + z2 / (4.0 * trials * trials)) / denominator
    return max(0.0, float(center - half)), min(1.0, float(center + half))


def _paired_mean_contrast(values: np.ndarray) -> dict:
    mean = float(np.mean(values))
    if len(values) < 2:
        return {"mean": mean, "standard_error": None, "normal_95_low": None, "normal_95_high": None}
    standard_error = float(np.std(values, ddof=1) / np.sqrt(len(values)))
    return {"mean": mean, "standard_error": standard_error,
            "normal_95_low": mean - 1.96 * standard_error,
            "normal_95_high": mean + 1.96 * standard_error}


def run_selection_misspecification(
    *,
    n_replicates: int = 30,
    n_parent: int = 6000,
    n_detected: int = 40,
    seed: int = 20261029,
    h0_true: float = 70.0,
    sigma_mpc: float = 140.0,
    detection_threshold_mpc: float = 250.0,
    h0_min: float = 50.0,
    h0_max: float = 90.0,
    h0_points: int = 161,
    z_points: int = 801,
    bandwidth_factor: float = 0.8,
) -> tuple[pd.DataFrame, dict]:
    """Run paired arms under one fixed true catalogue-selection law.

    Arms are oracle, catalogue-only, correct inverse-inclusion weighting,
    shifted-transition weighting, and shallow-slope weighting. The assumed
    selection model is only used to calculate weights; it never changes the
    generated catalogue or GW-like event distances.
    """
    if n_replicates < 1 or n_parent < 100 or n_detected < 1:
        raise ValueError("n_replicates/n_detected must be positive and n_parent >= 100")
    if h0_points < 31 or z_points < 101 or not (0 < h0_min < h0_max):
        raise ValueError("invalid H0 or redshift grid settings")
    if not np.isfinite(bandwidth_factor) or bandwidth_factor <= 0:
        raise ValueError("bandwidth_factor must be finite and positive")

    z = np.linspace(0.01, 0.40, z_points)
    h0 = np.linspace(h0_min, h0_max, h0_points)
    merger = 1.0 + 1.5 * (z - z[0]) / (z[-1] - z[0])
    target_density = _normalized(latent_galaxy_intensity(z) * merger, z)
    # Fixed truth law used by the generator: floor=.10, midpoint=.22, width=.05.
    true_selection = {"floor": 0.10, "midpoint": 0.22, "width": 0.05}
    analyst_models: dict[str, dict | None] = {
        "oracle": None,
        "catalogue_only": None,
        "correct_weighting": true_selection,
        "shifted_transition": {"floor": 0.10, "midpoint": 0.27, "width": 0.05},
        "shallow_slope": {"floor": 0.10, "midpoint": 0.22, "width": 0.075},
    }
    records: list[dict] = []
    for rep in range(n_replicates):
        rep_seed = seed + rep
        sim = generate_synthetic_population(
            seed=rep_seed, n_parent=n_parent, n_detected=n_detected,
            h0_true=h0_true, sigma_mpc=sigma_mpc,
            detection_threshold_mpc=detection_threshold_mpc,
            z_min=float(z[0]), z_max=float(z[-1]),
        )
        parent_z = sim["parent_truth"]["redshift"].to_numpy(dtype=float)
        observed_z = sim["catalogue_observed"]["redshift"].to_numpy(dtype=float)
        distances = sim["gw_events_observed"]["distance_mpc"].to_numpy(dtype=float)
        for arm, assumed in analyst_models.items():
            if arm == "oracle":
                sample = parent_z
                probability = None
            else:
                sample = observed_z
                probability = None if assumed is None else catalogue_inclusion_probability(
                    observed_z, **assumed
                )
            density = estimate_host_density(
                z, sample, inclusion_probability=probability,
                merger_weight=merger, bandwidth_factor=bandwidth_factor,
            )
            likelihood = selection_aware_toy_likelihood(
                z, density, h0, distances, sigma=sigma_mpc,
                detection_threshold=detection_threshold_mpc,
            )
            post = posterior_summary(h0, likelihood["selected_log_likelihood"])
            weights = np.ones(len(sample), dtype=float) if probability is None else 1.0 / probability
            normalized_weights = weights / weights.sum()
            ess = float(weights.sum() ** 2 / np.sum(weights ** 2))
            records.append({
                "replicate": rep, "seed": rep_seed, "arm": arm,
                "h0_true": h0_true, "n_parent": len(parent_z),
                "n_catalogue_observed": len(observed_z), "n_detected": len(distances),
                "assumed_selection_model": json.dumps(assumed, sort_keys=True),
                "density_ise": float(trapezoid((density - target_density) ** 2, z)),
                "posterior_mean": post["posterior_mean"],
                "posterior_median": post["posterior_median"],
                "ci68_low": post["ci68_low"], "ci68_high": post["ci68_high"],
                "ci68_width": post["ci68_width"],
                "posterior_integral": post["posterior_integral"],
                "posterior_density_grid_json": post["posterior_density_grid_json"],
                "covered_68": bool(post["ci68_low"] <= h0_true <= post["ci68_high"]),
                "weight_ess": ess,
                "weight_ess_fraction": ess / len(sample),
                "max_normalized_weight": float(np.max(normalized_weights)),
                "q95_normalized_weight": float(np.quantile(normalized_weights, 0.95)),
                "min_assumed_inclusion_probability": (
                    1.0 if probability is None else float(np.min(probability))
                ),
                "min_detection_efficiency": float(np.min(likelihood["detection_efficiency"])),
            })
    frame = pd.DataFrame(records)
    aggregate = []
    for arm, group in frame.groupby("arm", sort=False):
        errors = group["posterior_mean"].to_numpy(dtype=float) - h0_true
        covered = group["covered_68"].astype(bool).to_numpy()
        coverage_low, coverage_high = _wilson_interval(int(np.sum(covered)), len(covered))
        aggregate.append({
            "arm": arm, "n_replicates": len(group),
            "mean_bias": float(np.mean(errors)),
            "rmse": float(np.sqrt(np.mean(errors ** 2))),
            "empirical_68_coverage": float(np.mean(covered)),
            "coverage_numerator": int(np.sum(covered)),
            "coverage_denominator": int(len(covered)),
            "coverage_wilson_95_low": coverage_low,
            "coverage_wilson_95_high": coverage_high,
            "mean_68_interval_width": float(group["ci68_width"].mean()),
            "mean_density_ise": float(group["density_ise"].mean()),
            "mean_weight_ess_fraction": float(group["weight_ess_fraction"].mean()),
            "mean_max_normalized_weight": float(group["max_normalized_weight"].mean()),
        })
    baseline = frame[frame["arm"] == "catalogue_only"].set_index("replicate")
    paired_contrasts = []
    for arm in analyst_models:
        if arm == "catalogue_only":
            continue
        comparison = frame[frame["arm"] == arm].set_index("replicate")
        squared_error_difference = (
            (comparison["posterior_mean"] - h0_true) ** 2
            - (baseline["posterior_mean"] - h0_true) ** 2
        ).to_numpy(dtype=float)
        density_ise_difference = (
            comparison["density_ise"] - baseline["density_ise"]
        ).to_numpy(dtype=float)
        paired_contrasts.append({
            "arm": arm,
            "reference_arm": "catalogue_only",
            "negative_values_favor_comparison_arm": True,
            "squared_error_difference": _paired_mean_contrast(squared_error_difference),
            "density_ise_difference": _paired_mean_contrast(density_ise_difference),
        })
    summary = {
        "status": "synthetic_selection_misspecification_pilot_not_real_H0_measurement",
        "experiment": "selection_function_misspecification_v1",
        "seed_start": seed, "n_replicates": n_replicates,
        "n_parent": n_parent, "n_detected": n_detected,
        "h0_true_km_s_mpc": h0_true, "sigma_mpc": sigma_mpc,
        "detection_threshold_mpc": detection_threshold_mpc,
        "h0_grid": [h0_min, h0_max, h0_points],
        "redshift_grid": [float(z[0]), float(z[-1]), z_points],
        "bandwidth_factor": bandwidth_factor,
        "true_catalogue_selection": true_selection,
        "assumed_models": {key: value for key, value in analyst_models.items()},
        "density_target": "normalized latent_galaxy_intensity(z) times declared linear merger weight",
        "arms": list(analyst_models),
        "aggregate_metrics": aggregate,
        "paired_contrasts_vs_catalogue_only": paired_contrasts,
        "paired_contrast_interval_note": "Normal-approximation intervals use the replicate-level paired-difference standard error; for one replicate only, uncertainty bounds are null. Use more replicates for research claims.",
        "limitations": [
            "Synthetic low-redshift Euclidean toy, not real GW or LVK inference.",
            "Oracle density uses parent redshifts and is not observable in practice.",
            "Correct weighting assumes the true selection law is known.",
            "The target density is the declared analytic toy intensity times a merger weight.",
            "Coverage is imprecise at pilot replicate counts; no physical H0 uncertainty is represented.",
            "No event-rate marginalization, realistic detector response, peculiar velocities, or redshift uncertainty.",
        ],
    }
    return frame, summary
