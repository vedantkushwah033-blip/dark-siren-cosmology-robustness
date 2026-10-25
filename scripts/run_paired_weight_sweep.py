"""Paired host-weight stress test across host separation, noise, and prior support.

Each synthetic observation is evaluated with correct and reversed host weights
under both a narrow and a broad H0 prior. The broad-prior arm checks whether an
apparent lack of weight sensitivity is caused by excluding alternate-host H0 modes.
This is a synthetic low-redshift diagnostic, not real-data validation.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary

H0_TRUE = 70.0
Z1 = 0.02
Z2_VALUES = (0.025, 0.04, 0.06)
NOISE_FRACTIONS = (0.03, 0.10, 0.20)
TRUE_WEIGHTS = np.array([0.7, 0.3])
WRONG_WEIGHTS = np.array([0.3, 0.7])
PRIOR_SUPPORTS = {
    "narrow_50_90": (50.0, 90.0),
    "broad_20_220": (20.0, 220.0),
}
# 0.025 H0 spacing is retained while widening support enough to contain the
# alternate-host H0 modes implied by z=0.02 and z=0.06 at injected H0=70.
GRID = np.linspace(10.0, 250.0, 9601)


def _percentile_ci(values, rng, resamples, statistic="mean"):
    """Percentile bootstrap CI, resampling replicate indices as paired units."""
    values = np.asarray(values, dtype=float)
    draws = np.empty(resamples, dtype=float)
    for start in range(0, resamples, 250):
        stop = min(start + 250, resamples)
        idx = rng.integers(0, len(values), size=(stop - start, len(values)))
        if statistic == "mean":
            draws[start:stop] = values[idx].mean(axis=1)
        elif statistic == "rmse_delta":
            draws[start:stop] = (
                np.sqrt(np.mean(values[idx, 1] ** 2, axis=1))
                - np.sqrt(np.mean(values[idx, 0] ** 2, axis=1))
            )
        else:
            raise ValueError(f"Unknown bootstrap statistic: {statistic}")
    low, high = np.quantile(draws, [0.025, 0.975])
    return float(low), float(high)


def _infer(distance, sigma, redshifts, weights, prior_min, prior_max):
    density = posterior_h0(
        GRID, redshifts, distance, sigma, weights=weights,
        prior_min=prior_min, prior_max=prior_max,
    )
    return posterior_summary(GRID, density)


def run_paired_sweep(replicates=1000, seed=20261028, bootstrap_resamples=4000):
    if bootstrap_resamples < 1000:
        raise ValueError("Use at least 1000 bootstrap resamples.")
    if replicates < 200:
        raise ValueError("Use at least 200 replicates per design cell.")

    rng = np.random.default_rng(seed)
    rows, summaries = [], []
    cell_index = 0
    for z2 in Z2_VALUES:
        redshifts = np.array([Z1, z2])
        for noise_fraction in NOISE_FRACTIONS:
            sigma = noise_fraction * C_KM_S * 0.03 / H0_TRUE
            cell_rows = []
            for replicate in range(replicates):
                host_index = int(rng.choice(2, p=TRUE_WEIGHTS))
                distance_true = C_KM_S * redshifts[host_index] / H0_TRUE
                distance_observed = float(rng.normal(distance_true, sigma))
                while distance_observed <= 0:
                    distance_observed = float(rng.normal(distance_true, sigma))
                record = {
                    "z1": Z1, "z2": z2, "redshift_separation": z2 - Z1,
                    "noise_fraction": noise_fraction, "replicate": replicate,
                    "true_host_index": host_index,
                    "true_host_redshift": float(redshifts[host_index]),
                    "true_H0": H0_TRUE,
                    "observed_distance_mpc": distance_observed, "sigma_mpc": sigma,
                }
                for support_index, (support_name, (prior_min, prior_max)) in enumerate(PRIOR_SUPPORTS.items()):
                    posteriors = {}
                    for label, weights in (("correct", TRUE_WEIGHTS), ("reversed", WRONG_WEIGHTS)):
                        post = _infer(distance_observed, sigma, redshifts, weights, prior_min, prior_max)
                        posteriors[label] = post
                        prefix = f"{support_name}_{label}"
                        record[f"{prefix}_median"] = post["median"]
                        record[f"{prefix}_median_error"] = post["median"] - H0_TRUE
                        record[f"{prefix}_lower_68"] = post["lower_68"]
                        record[f"{prefix}_upper_68"] = post["upper_68"]
                        record[f"{prefix}_covered_68"] = int(post["lower_68"] <= H0_TRUE <= post["upper_68"])
                        record[f"{prefix}_width_68"] = post["upper_68"] - post["lower_68"]
                    prefix = support_name
                    record[f"paired_abs_error_difference_{prefix}"] = (
                        abs(posteriors["reversed"]["median"] - H0_TRUE)
                        - abs(posteriors["correct"]["median"] - H0_TRUE)
                    )
                rows.append(record)
                cell_rows.append(record)

            for support_index, (support_name, (prior_min, prior_max)) in enumerate(PRIOR_SUPPORTS.items()):
                correct_error = np.array([r[f"{support_name}_correct_median_error"] for r in cell_rows])
                reversed_error = np.array([r[f"{support_name}_reversed_median_error"] for r in cell_rows])
                abs_error_diffs = np.array([r[f"paired_abs_error_difference_{support_name}"] for r in cell_rows])
                coverage_diffs = np.array([
                    r[f"{support_name}_reversed_covered_68"] - r[f"{support_name}_correct_covered_68"]
                    for r in cell_rows
                ], dtype=float)
                bootstrap_seed = seed + 104729 * cell_index + 1009 * support_index
                abs_ci = _percentile_ci(abs_error_diffs, np.random.default_rng(bootstrap_seed), bootstrap_resamples)
                coverage_ci = _percentile_ci(coverage_diffs, np.random.default_rng(bootstrap_seed + 1), bootstrap_resamples)
                rmse_ci = _percentile_ci(
                    np.column_stack([correct_error, reversed_error]),
                    np.random.default_rng(bootstrap_seed + 2), bootstrap_resamples,
                    statistic="rmse_delta",
                )
                correct_coverage = np.array([r[f"{support_name}_correct_covered_68"] for r in cell_rows])
                reversed_coverage = np.array([r[f"{support_name}_reversed_covered_68"] for r in cell_rows])
                correct_width = np.array([r[f"{support_name}_correct_width_68"] for r in cell_rows])
                reversed_width = np.array([r[f"{support_name}_reversed_width_68"] for r in cell_rows])
                summaries.append({
                    "prior_support": support_name, "prior_min": prior_min, "prior_max": prior_max,
                    "z1": Z1, "z2": z2, "redshift_separation": z2 - Z1,
                    "noise_fraction": noise_fraction, "replicates": replicates,
                    "correct_mean_median_error": float(correct_error.mean()),
                    "reversed_mean_median_error": float(reversed_error.mean()),
                    "correct_rmse": float(np.sqrt(np.mean(correct_error ** 2))),
                    "reversed_rmse": float(np.sqrt(np.mean(reversed_error ** 2))),
                    "rmse_diff_reversed_minus_correct": float(np.sqrt(np.mean(reversed_error ** 2)) - np.sqrt(np.mean(correct_error ** 2))),
                    "rmse_diff_ci95_low": rmse_ci[0], "rmse_diff_ci95_high": rmse_ci[1],
                    "correct_coverage_68": float(correct_coverage.mean()),
                    "reversed_coverage_68": float(reversed_coverage.mean()),
                    "coverage_68_diff_reversed_minus_correct": float(coverage_diffs.mean()),
                    "coverage_68_diff_ci95_low": coverage_ci[0], "coverage_68_diff_ci95_high": coverage_ci[1],
                    "correct_mean_width_68": float(correct_width.mean()),
                    "reversed_mean_width_68": float(reversed_width.mean()),
                    "mean_paired_abs_error_diff_reversed_minus_correct": float(abs_error_diffs.mean()),
                    "median_paired_abs_error_diff_reversed_minus_correct": float(np.median(abs_error_diffs)),
                    "fraction_reversed_has_larger_abs_error": float(np.mean(abs_error_diffs > 0)),
                    "paired_abs_error_diff_ci95_low": abs_ci[0], "paired_abs_error_diff_ci95_high": abs_ci[1],
                    "bootstrap_resamples": bootstrap_resamples,
                    "bootstrap_seed_abs_error": bootstrap_seed,
                    "bootstrap_seed_coverage": bootstrap_seed + 1,
                    "bootstrap_seed_rmse": bootstrap_seed + 2,
                })
            cell_index += 1
    return rows, summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261028)
    parser.add_argument("--bootstrap-resamples", type=int, default=4000)
    parser.add_argument("--output-dir", default="results/paired-weight-sweep")
    args = parser.parse_args()
    rows, summaries = run_paired_sweep(args.replicates, args.seed, args.bootstrap_resamples)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, data in (("per_run.csv", rows), ("summary.csv", summaries)):
        with (out / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    print(
        f"Completed {len(rows)} synthetic observations with 4 posterior fits each "
        f"({4 * len(rows)} fits), across {len(summaries)} prior-support × design cells."
    )
    for summary in summaries:
        print(summary)


if __name__ == "__main__":
    main()
