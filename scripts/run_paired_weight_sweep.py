"""Paired stress test of host-weight misspecification across separation and noise.

Each replicate's observation is evaluated under both correct and reversed
weights, ensuring the contrast is paired. Synthetic low-z model only.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary

H0_TRUE = 70.0
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
Z1 = 0.02
Z2_VALUES = (0.025, 0.04, 0.06)
NOISE_FRACTIONS = (0.03, 0.10, 0.20)
TRUE_WEIGHTS = np.array([0.7, 0.3])
WRONG_WEIGHTS = np.array([0.3, 0.7])
GRID = np.linspace(25.0, 125.0, 4001)


def run_paired_sweep(replicates=1000, seed=20261028):
    if replicates < 200:
        raise ValueError("Use at least 200 replicates per design cell.")
    rng = np.random.default_rng(seed)
    summaries, rows = [], []
    for z2 in Z2_VALUES:
        redshifts = np.array([Z1, z2])
        for noise_fraction in NOISE_FRACTIONS:
            sigma = noise_fraction * C_KM_S * 0.03 / H0_TRUE
            cell = []
            for i in range(replicates):
                host = int(rng.choice(2, p=TRUE_WEIGHTS))
                d_true = C_KM_S * redshifts[host] / H0_TRUE
                d_obs = float(rng.normal(d_true, sigma))
                while d_obs <= 0:
                    d_obs = float(rng.normal(d_true, sigma))
                record = {
                    "z1": Z1, "z2": z2, "noise_fraction": noise_fraction,
                    "replicate": i, "true_host_index": host,
                    "observed_distance_mpc": d_obs, "sigma_mpc": sigma,
                }
                estimates = {}
                for label, weights in (("correct", TRUE_WEIGHTS), ("reversed", WRONG_WEIGHTS)):
                    density = posterior_h0(
                        GRID, redshifts, d_obs, sigma, weights=weights,
                        prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
                    )
                    post = posterior_summary(GRID, density)
                    estimates[label] = post
                    record[f"{label}_median"] = post["median"]
                    record[f"{label}_median_error"] = post["median"] - H0_TRUE
                    record[f"{label}_lower_68"] = post["lower_68"]
                    record[f"{label}_upper_68"] = post["upper_68"]
                    record[f"{label}_covered_68"] = int(post["lower_68"] <= H0_TRUE <= post["upper_68"])
                    record[f"{label}_width_68"] = post["upper_68"] - post["lower_68"]
                record["paired_abs_error_difference_reversed_minus_correct"] = (
                    abs(estimates["reversed"]["median"] - H0_TRUE)
                    - abs(estimates["correct"]["median"] - H0_TRUE)
                )
                rows.append(record)
                cell.append(record)
            diffs = np.array([r["paired_abs_error_difference_reversed_minus_correct"] for r in cell])
            err_correct = np.array([r["correct_median_error"] for r in cell])
            err_wrong = np.array([r["reversed_median_error"] for r in cell])
            summaries.append({
                "z1": Z1, "z2": z2, "redshift_separation": z2 - Z1,
                "noise_fraction": noise_fraction, "replicates": replicates,
                "correct_mean_median_error": float(err_correct.mean()),
                "reversed_mean_median_error": float(err_wrong.mean()),
                "correct_rmse": float(np.sqrt(np.mean(err_correct**2))),
                "reversed_rmse": float(np.sqrt(np.mean(err_wrong**2))),
                "correct_coverage_68": float(np.mean([r["correct_covered_68"] for r in cell])),
                "reversed_coverage_68": float(np.mean([r["reversed_covered_68"] for r in cell])),
                "mean_paired_abs_error_diff_reversed_minus_correct": float(diffs.mean()),
                "median_paired_abs_error_diff_reversed_minus_correct": float(np.median(diffs)),
                "fraction_reversed_has_larger_abs_error": float(np.mean(diffs > 0)),
                "correct_mean_width_68": float(np.mean([r["correct_width_68"] for r in cell])),
                "reversed_mean_width_68": float(np.mean([r["reversed_width_68"] for r in cell])),
            })
    return rows, summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261028)
    parser.add_argument("--output-dir", default="results/paired-weight-sweep")
    args = parser.parse_args()
    rows, summaries = run_paired_sweep(args.replicates, args.seed)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, data in (("per_run.csv", rows), ("summary.csv", summaries)):
        with (out / name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    print(f"Completed {len(rows)} paired observations × 2 inference weightings = {2*len(rows)} posterior fits across {len(summaries)} design cells.")
    for row in summaries:
        print(row)


if __name__ == "__main__":
    main()
