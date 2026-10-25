"""Known-host calibration surface over true H0 and distance-noise levels.

This isolates the one-host likelihood from candidate-mixture uncertainty. The
simulation uses the same low-redshift d_L = c z / H0 model in generation and
inference. Truths at 50 and 90 deliberately sit on the prior boundaries.
This remains a synthetic diagnostic, not an observational analysis.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

H0_TRUTHS = (50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0)
FRACTIONAL_ERRORS = (0.05, 0.10, 0.15, 0.20)
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0


def wilson_interval(successes: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion."""
    if n <= 0 or successes < 0 or successes > n:
        raise ValueError("require n > 0 and 0 <= successes <= n")
    p = successes / n
    z2 = z * z
    denom = 1.0 + z2 / n
    center = (p + z2 / (2.0 * n)) / denom
    half = z * math.sqrt(p * (1.0 - p) / n + z2 / (4.0 * n * n)) / denom
    # Avoid a tiny positive lower bound from floating-point cancellation at p=0,
    # and make the complementary all-success boundary exact as well.
    lower = 0.0 if successes == 0 else center - half
    upper = 1.0 if successes == n else center + half
    return max(0.0, lower), min(1.0, upper)


def summarize_cell(rows: list[dict], true_h0: float, fractional_error: float) -> dict:
    """Recompute coverage, error, width, pull and boundary diagnostics for one cell."""
    if not rows:
        raise ValueError("rows must not be empty")
    medians = [float(row["median"]) for row in rows]
    means = [float(row["mean"]) for row in rows]
    lo68 = [float(row["lower_68"]) for row in rows]
    hi68 = [float(row["upper_68"]) for row in rows]
    lo95 = [float(row["lower_95"]) for row in rows]
    hi95 = [float(row["upper_95"]) for row in rows]
    sd = [float(row["posterior_sd"]) for row in rows]
    if not all(math.isfinite(x) for xs in (medians, means, lo68, hi68, lo95, hi95, sd) for x in xs):
        raise ValueError("all posterior summaries must be finite")
    if any(a > b for a, b in zip(lo68, hi68)) or any(a > b for a, b in zip(lo95, hi95)):
        raise ValueError("credible interval bounds are reversed")
    if any(lo95_i > lo68_i or hi68_i > hi95_i for lo95_i, lo68_i, hi68_i, hi95_i in zip(lo95, lo68, hi68, hi95)):
        raise ValueError("68% intervals must be nested inside 95% intervals")
    if any(
        x < PRIOR_MIN or x > PRIOR_MAX
        for values in (medians, means, lo68, hi68, lo95, hi95)
        for x in values
    ):
        raise ValueError("posterior summaries and interval bounds must lie within prior support")
    if any(x <= 0.0 for x in sd):
        raise ValueError("posterior standard deviations must be positive")
    if (
        not math.isfinite(true_h0)
        or true_h0 < PRIOR_MIN
        or true_h0 > PRIOR_MAX
        or not math.isfinite(fractional_error)
        or fractional_error <= 0
    ):
        raise ValueError("truth must be finite and within prior support; fractional distance error must be positive and finite")
    n = len(rows)
    c68n = sum(a <= true_h0 <= b for a, b in zip(lo68, hi68))
    c95n = sum(a <= true_h0 <= b for a, b in zip(lo95, hi95))
    median_errors = [x - true_h0 for x in medians]
    mean_errors = [x - true_h0 for x in means]
    pulls = [error / max(s, 1e-12) for error, s in zip(median_errors, sd)]
    w68 = [b - a for a, b in zip(lo68, hi68)]
    w95 = [b - a for a, b in zip(lo95, hi95)]
    ci68 = wilson_interval(c68n, n)
    ci95 = wilson_interval(c95n, n)
    return {
        "true_h0": true_h0,
        "fractional_distance_error": fractional_error,
        "n": n,
        "coverage_68": c68n / n,
        "coverage_68_successes": c68n,
        "coverage_68_wilson_low": ci68[0],
        "coverage_68_wilson_high": ci68[1],
        "coverage_68_mc_se": math.sqrt((c68n / n) * (1.0 - c68n / n) / n),
        "coverage_95": c95n / n,
        "coverage_95_successes": c95n,
        "coverage_95_wilson_low": ci95[0],
        "coverage_95_wilson_high": ci95[1],
        "coverage_95_mc_se": math.sqrt((c95n / n) * (1.0 - c95n / n) / n),
        "median_bias": float(sorted(median_errors)[n // 2] if n % 2 else (sorted(median_errors)[n // 2 - 1] + sorted(median_errors)[n // 2]) / 2.0),
        "mean_median_error": sum(median_errors) / n,
        "mean_posterior_mean_error": sum(mean_errors) / n,
        "rmse_median": math.sqrt(sum(x * x for x in median_errors) / n),
        "rmse_mean": math.sqrt(sum(x * x for x in mean_errors) / n),
        "mean_width_68": sum(w68) / n,
        "mean_width_95": sum(w95) / n,
        "mean_posterior_sd": sum(sd) / n,
        "mean_median_pull": sum(pulls) / n,
        "sd_median_pull": math.sqrt(sum((x - sum(pulls) / n) ** 2 for x in pulls) / (n - 1)) if n > 1 else 0.0,
        "median_near_prior_edge_fraction": sum(
            (x - PRIOR_MIN <= 0.5) or (PRIOR_MAX - x <= 0.5) for x in medians
        ) / n,
        "truth_on_prior_boundary": true_h0 in (PRIOR_MIN, PRIOR_MAX),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", default="original-project")
    parser.add_argument("--replicates", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261025)
    parser.add_argument("--output-dir", default="results/known-host-calibration-surface")
    args = parser.parse_args()
    if not 100 <= args.replicates <= 2000:
        parser.error("--replicates must be between 100 and 2000")
    if args.seed < 0:
        parser.error("--seed must be non-negative")

    import numpy as np

    original = Path(args.original_dir).resolve()
    if not (original / "src" / "experiments.py").is_file():
        raise SystemExit(f"Could not find original simulator under {original}")
    sys.path.insert(0, str(original))
    from src.simulations_v2 import simulate_galaxy_catalogue, choose_host, simulate_gw_event
    from src.bayesian_inference import posterior_h0, posterior_summary

    grid = np.linspace(PRIOR_MIN, PRIOR_MAX, 1601)
    rows: list[dict] = []
    for truth_index, true_h0 in enumerate(H0_TRUTHS):
        for replicate in range(args.replicates):
            # Keep noise-level pairing within each truth/replicate, but use a
            # disjoint seed block across truths to avoid accidental cross-truth
            # reuse of identical catalogues and host draws.
            seed = args.seed + truth_index * args.replicates + replicate
            # Reuse the same redshifts, luminosities, host draw and standardized
            # noise across the four error levels to improve paired comparisons.
            cat = simulate_galaxy_catalogue(n_galaxies=500, H0_true=true_h0, seed=seed)
            host = choose_host(cat, weighting="luminosity", seed=seed + 100000)
            for frac_error in FRACTIONAL_ERRORS:
                gw = simulate_gw_event(
                    cat["true_distance_mpc"][host],
                    fractional_distance_error=frac_error,
                    seed=seed + 200000,
                )
                posterior = posterior_h0(
                    grid,
                    np.asarray([cat["redshift"][host]], dtype=float),
                    gw["observed_distance_mpc"],
                    gw["sigma_mpc"],
                    weights=np.ones(1),
                    prior_min=PRIOR_MIN,
                    prior_max=PRIOR_MAX,
                )
                summary = posterior_summary(grid, posterior)
                rows.append({
                    "run": replicate,
                    "seed": seed,
                    "true_h0": true_h0,
                    "fractional_distance_error": frac_error,
                    "host_redshift": float(cat["redshift"][host]),
                    "true_distance_mpc": float(cat["true_distance_mpc"][host]),
                    "observed_distance_mpc": float(gw["observed_distance_mpc"]),
                    "distance_sigma_mpc": float(gw["sigma_mpc"]),
                    "median": float(summary["median"]),
                    "mean": float(summary["mean"]),
                    "posterior_sd": float(summary["posterior_sd"]),
                    "lower_68": float(summary["lower_68"]),
                    "upper_68": float(summary["upper_68"]),
                    "lower_95": float(summary["lower_95"]),
                    "upper_95": float(summary["upper_95"]),
                    "width_68": float(summary["upper_68"] - summary["lower_68"]),
                    "width_95": float(summary["upper_95"] - summary["lower_95"]),
                    "pull_median": float((summary["median"] - true_h0) / max(summary["posterior_sd"], 1e-12)),
                })

    summaries = [
        summarize_cell(
            [row for row in rows if row["true_h0"] == truth and row["fractional_distance_error"] == frac],
            truth, frac
        )
        for truth in H0_TRUTHS for frac in FRACTIONAL_ERRORS
    ]
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for filename, data in (("per_run_results.csv", rows), ("cell_summary.csv", summaries)):
        with (out / filename).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    metadata = {
        "title": "Known-host calibration surface",
        "original_repository": "vedantkushwah033-blip/bayesian-dark-siren-cosmology",
        "original_commit": "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f",
        "true_h0_values": list(H0_TRUTHS),
        "fractional_distance_errors": list(FRACTIONAL_ERRORS),
        "replicates_per_cell": args.replicates,
        "cell_count": len(summaries),
        "raw_row_count": len(rows),
        "prior": [PRIOR_MIN, PRIOR_MAX],
        "posterior_grid_points": len(grid),
        "seed": args.seed,
        "paired_design": "Within each truth/replicate, reuse the catalogue seed, host selection, and standardized distance-noise draw across four error levels; use disjoint seed blocks across truths: base_seed + truth_index * replicates + replicate_index.",
        "interpretation_warning": "Synthetic one-host diagnostic under the matched low-redshift distance model. Boundary truths are intentionally included to expose prior-edge effects; not a real GW measurement.",
        "cells": summaries,
    }
    (out / "cell_summary.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps({k: metadata[k] for k in ("title", "cell_count", "raw_row_count", "replicates_per_cell", "prior")}, indent=2))


if __name__ == "__main__":
    main()
