"""Compare fixed-H0 repeated-sampling coverage under two prior ranges."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np

LEVELS = (0.50, 0.68, 0.80, 0.90, 0.95)
DEFAULT_H0_VALUES = (55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0)
PRIOR_RANGES = ((50.0, 90.0), (30.0, 120.0))


def interval(grid: np.ndarray, posterior: np.ndarray, level: float) -> tuple[float, float]:
    cdf = np.concatenate(
        [[0.0], np.cumsum((posterior[1:] + posterior[:-1]) * np.diff(grid) / 2)]
    )
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ValueError("Posterior CDF normalization failed")
    cdf /= cdf[-1]
    alpha = (1.0 - level) / 2.0
    return float(np.interp(alpha, cdf, grid)), float(np.interp(1.0 - alpha, cdf, grid))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", default="original-project")
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--seed", type=int, default=73021)
    parser.add_argument("--output", default="results/fixed-h0-coverage-grid.csv")
    parser.add_argument("--h0-values", type=float, nargs="+", default=list(DEFAULT_H0_VALUES))
    args = parser.parse_args()
    if args.replicates < 2:
        parser.error("--replicates must be at least 2")
    if any(not 30.0 <= h <= 120.0 for h in args.h0_values):
        parser.error("Every H0 value must be inside both comparison priors [30, 120]")

    original = Path(args.original_dir).resolve()
    sys.path.insert(0, str(original))
    from src.simulations_v2 import simulate_galaxy_catalogue, choose_host, simulate_gw_event
    from src.bayesian_inference import posterior_h0, posterior_summary

    rows: list[dict] = []
    for prior_min, prior_max in PRIOR_RANGES:
        grid = np.linspace(
            prior_min, prior_max, int(round((prior_max - prior_min) / 0.025)) + 1
        )
        for h0_true in args.h0_values:
            for scenario in ("known_host", "unknown_host"):
                interval_data = {level: [] for level in LEVELS}
                medians, means, lower_95, upper_95 = [], [], [], []
                for r in range(args.replicates):
                    seed = args.seed + int(round(h0_true * 1000)) * 100000 + r * 10
                    cat = simulate_galaxy_catalogue(
                        n_galaxies=500, H0_true=h0_true, seed=seed
                    )
                    host = choose_host(cat, weighting="uniform", seed=seed + 1)
                    gw = simulate_gw_event(
                        cat["true_distance_mpc"][host], seed=seed + 2
                    )
                    if scenario == "known_host":
                        redshifts = np.array([cat["redshift"][host]])
                        weights = np.array([1.0])
                    else:
                        redshifts = cat["redshift"]
                        weights = np.ones(len(redshifts))
                    posterior = posterior_h0(
                        grid,
                        redshifts,
                        gw["observed_distance_mpc"],
                        gw["sigma_mpc"],
                        weights=weights,
                        prior_min=prior_min,
                        prior_max=prior_max,
                    )
                    summary = posterior_summary(grid, posterior)
                    medians.append(float(summary["median"]))
                    means.append(float(summary["mean"]))
                    lower_95.append(float(summary["lower_95"]))
                    upper_95.append(float(summary["upper_95"]))
                    for level in LEVELS:
                        low, high = interval(grid, posterior, level)
                        interval_data[level].append((low, high))

                row = {
                    "prior_min": prior_min,
                    "prior_max": prior_max,
                    "scenario": scenario,
                    "H0_true": h0_true,
                    "replicates": args.replicates,
                    "mean_median_error": float(np.mean(np.asarray(medians) - h0_true)),
                    "mean_posterior_mean_error": float(np.mean(np.asarray(means) - h0_true)),
                    "mean_width_68": float(
                        np.mean([hi - lo for lo, hi in interval_data[0.68]])
                    ),
                    "mean_width_95": float(
                        np.mean([hi - lo for lo, hi in interval_data[0.95]])
                    ),
                    "fraction_95_interval_touches_lower_prior_bound": float(
                        np.mean(np.asarray(lower_95) <= prior_min + 0.05)
                    ),
                    "fraction_95_interval_touches_upper_prior_bound": float(
                        np.mean(np.asarray(upper_95) >= prior_max - 0.05)
                    ),
                }
                for level in LEVELS:
                    pairs = interval_data[level]
                    coverage = float(
                        np.mean([low <= h0_true <= high for low, high in pairs])
                    )
                    tag = int(round(level * 100))
                    row[f"coverage_{tag}"] = coverage
                    row[f"coverage_{tag}_mc_se"] = float(
                        np.sqrt(coverage * (1.0 - coverage) / args.replicates)
                    )
                rows.append(row)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} H0/prior/scenario rows to {output}")


if __name__ == "__main__":
    main()
