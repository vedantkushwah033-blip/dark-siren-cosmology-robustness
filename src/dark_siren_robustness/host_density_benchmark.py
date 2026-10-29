"""Exploratory synthetic comparison of host-density estimates in toy H0 inference."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from dark_siren_robustness.catalogue_density import estimate_host_density
from dark_siren_robustness.generative_selection import catalogue_inclusion_probability, generate_synthetic_population, latent_galaxy_intensity
from dark_siren_robustness.selection_toy import selection_aware_toy_likelihood

def posterior_summary(grid, loglike):
    """Normalize a log-likelihood on an explicit, finite H0 grid.

    The function uses an implicit flat prior over the supplied grid. It rejects
    malformed inputs rather than allowing NaNs, duplicate grid points, or
    mismatched arrays to contaminate posterior summaries.
    """
    grid = np.asarray(grid, dtype=float)
    loglike = np.asarray(loglike, dtype=float)
    if grid.ndim != 1 or grid.size < 2 or not np.all(np.isfinite(grid)):
        raise ValueError("H0 grid must be a finite one-dimensional array with at least two points")
    if np.any(np.diff(grid) <= 0):
        raise ValueError("H0 grid must be strictly increasing")
    if loglike.ndim != 1 or loglike.shape != grid.shape or not np.all(np.isfinite(loglike)):
        raise ValueError("loglike must be finite, one-dimensional, and match the H0 grid")
    shifted = loglike - np.max(loglike)
    density = np.exp(shifted)
    area = float(trapezoid(density, grid))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("H0 posterior normalization failed")
    density /= area
    cdf = np.r_[0.0, cumulative_trapezoid(density, grid)]
    cdf /= cdf[-1]
    q16, q50, q84 = (float(np.interp(q, cdf, grid)) for q in (0.16, 0.50, 0.84))
    return {"posterior_mean": float(trapezoid(grid*density, grid)),
            "posterior_median": q50, "ci68_low": q16, "ci68_high": q84,
            "ci68_width": q84-q16, "posterior_integral": float(trapezoid(density, grid)),
            # Persist the actual normalized grid so the auditor can reconstruct
            # summaries rather than trusting a scalar integral claim.
            "posterior_density_grid_json": json.dumps(density.tolist(), separators=(",", ":"))}

def run_benchmark(*, n_replicates=10, n_parent=6000, n_detected=40, seed=20261028,
                  h0_true=70.0, sigma_mpc=140.0, detection_threshold_mpc=250.0,
                  h0_min=50.0, h0_max=90.0, h0_points=161, z_points=801,
                  bandwidth_factor=None, supplied_midpoint=0.22, supplied_width=0.05):
    if n_replicates < 1 or n_parent < 100 or n_detected < 1:
        raise ValueError("n_replicates/n_detected must be positive and n_parent >= 100")
    if h0_points < 31 or z_points < 101 or not (0 < h0_min < h0_max):
        raise ValueError("invalid H0 or redshift grid settings")
    if bandwidth_factor is not None:
        bandwidth_factor = float(bandwidth_factor)
        if not np.isfinite(bandwidth_factor) or bandwidth_factor <= 0:
            raise ValueError("bandwidth_factor must be finite and positive when supplied")
    if not np.isfinite(supplied_midpoint) or supplied_midpoint <= 0:
        raise ValueError("supplied_midpoint must be finite and positive")
    if not np.isfinite(supplied_width) or supplied_width <= 0:
        raise ValueError("supplied_width must be finite and positive")
    z = np.linspace(0.01, 0.40, z_points)
    h0 = np.linspace(h0_min, h0_max, h0_points)
    merger = 1.0 + 1.5*(z-z[0])/(z[-1]-z[0])
    true_density = latent_galaxy_intensity(z) * merger
    true_density /= trapezoid(true_density, z)
    rows = []
    for rep in range(n_replicates):
        sim = generate_synthetic_population(seed=seed+rep, n_parent=n_parent, n_detected=n_detected,
            h0_true=h0_true, sigma_mpc=sigma_mpc, detection_threshold_mpc=detection_threshold_mpc,
            z_min=float(z[0]), z_max=float(z[-1]))
        parent_z = sim["parent_truth"]["redshift"].to_numpy(dtype=float)
        observed_z = sim["catalogue_observed"]["redshift"].to_numpy(dtype=float)
        # Known selection is evaluated from the declared law using observable redshifts.
        p_include = catalogue_inclusion_probability(
            observed_z, midpoint=supplied_midpoint, width=supplied_width
        )
        distances = sim["gw_events_observed"]["distance_mpc"].to_numpy(dtype=float)
        arms = {"oracle": (parent_z, None), "catalogue_only": (observed_z, None),
                "known_completeness": (observed_z, p_include)}
        for arm, (sample, probabilities) in arms.items():
            host = estimate_host_density(z, sample, inclusion_probability=probabilities,
                                         merger_weight=merger, bandwidth_factor=bandwidth_factor)
            like = selection_aware_toy_likelihood(z, host, h0, distances, sigma=sigma_mpc,
                                                  detection_threshold=detection_threshold_mpc)
            s = posterior_summary(h0, like["selected_log_likelihood"])
            estimation_weights = (np.ones(sample.size, dtype=float) if probabilities is None
                                  else 1.0 / np.asarray(probabilities, dtype=float))
            ess = float(estimation_weights.sum()**2 / np.sum(estimation_weights**2))
            rows.append({"replicate": rep, "seed": seed+rep, "arm": arm, "h0_true": h0_true,
                "n_parent": len(parent_z), "n_catalogue_observed": len(observed_z),
                "n_detected": len(distances), "min_detection_efficiency": float(np.min(like["detection_efficiency"])),
                "l1_host_density_error": float(trapezoid(np.abs(host-true_density), z)),
                "max_weight_ratio": float(np.max(estimation_weights) / np.mean(estimation_weights)),
                "ess_fraction": float(ess / estimation_weights.size),
                **s, "covered_68": bool(s["ci68_low"] <= h0_true <= s["ci68_high"])})
    frame = pd.DataFrame(rows)
    aggregate = []
    for arm, group in frame.groupby("arm", sort=False):
        errors = group.posterior_mean.to_numpy() - h0_true
        aggregate.append({"arm": arm, "n_replicates": len(group), "mean_bias": float(np.mean(errors)),
            "rmse": float(np.sqrt(np.mean(errors**2))), "empirical_68_coverage": float(group.covered_68.mean()),
            "mean_68_interval_width": float(group.ci68_width.mean()),
            "mean_posterior_mean": float(group.posterior_mean.mean())})
    summary = {"status": "synthetic_exploratory_benchmark_not_real_H0_measurement",
        "method": "paired seeds; fixed H0; Gaussian distance noise; fixed-threshold detection",
        "n_replicates": n_replicates, "n_parent": n_parent, "n_detected": n_detected,
        "seed_start": seed, "h0_true_km_s_mpc": h0_true, "sigma_mpc": sigma_mpc,
        "detection_threshold_mpc": detection_threshold_mpc,
        "h0_grid": [h0_min, h0_max, h0_points], "redshift_grid": [float(z[0]), float(z[-1]), z_points],
        "bandwidth_factor": None if bandwidth_factor is None else float(bandwidth_factor),
        "bandwidth_rule": ("weighted Scott rule (h = weighted_SD * n_eff**(-1/5))"
                           if bandwidth_factor is None else
                           "explicit physical bandwidth h = factor * weighted_SD of original sample"),
        "supplied_inclusion_midpoint": float(supplied_midpoint),
        "supplied_inclusion_width": float(supplied_width),
        "arms": {
            "oracle": "parent truth redshifts; benchmark reference only",
            "catalogue_only": "observed catalogue without completeness correction",
            "known_completeness": "observed catalogue with known declared inclusion probabilities"},
        "limitations": ["No real detector response or event-rate inference.",
            "Known-completeness assumes the selection law is exactly known.",
            "Bandwidth and grid sensitivity are evaluated in separate paired campaigns.",
            "The supplied inclusion law may deliberately differ from the true generator law.",
            "Density L1 error and weight diagnostics are scoring/diagnostic metrics, not real-survey guarantees.",
            "Coverage is noisy for small replicate counts.",
            "Truth is used for scoring and oracle only, never non-oracle estimation."],
        "aggregate_metrics": aggregate}
    return frame, summary

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replicates", type=int, default=10); p.add_argument("--n-parent", type=int, default=6000)
    p.add_argument("--n-detected", type=int, default=40); p.add_argument("--seed", type=int, default=20261028)
    p.add_argument("--h0-points", type=int, default=161); p.add_argument("--z-points", type=int, default=801)
    p.add_argument("--output-dir", type=Path, default=Path("artifacts/host_density_benchmark"))
    a = p.parse_args()
    rows, summary = run_benchmark(n_replicates=a.replicates, n_parent=a.n_parent,
        n_detected=a.n_detected, seed=a.seed, h0_points=a.h0_points, z_points=a.z_points)
    a.output_dir.mkdir(parents=True, exist_ok=True)
    rows.to_csv(a.output_dir/"replicate_metrics.csv", index=False)
    (a.output_dir/"summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
if __name__ == "__main__": main()
