"""Paired fixed-truth comparison of toy selection-corrected and uncorrected inference.

The same detected mock observation is analysed twice, once with the ordinary
host-mixture likelihood and once with the likelihood divided by beta(H0).
This is a controlled toy-model diagnostic, not a realistic GW selection study.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary
from dark_siren_robustness.selection import posterior_h0_selected

PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
HOST_REDSHIFTS = np.array([0.02, 0.04])
HOST_WEIGHTS = np.array([0.70, 0.30])
SIGMA_MPC = 20.0
DETECTION_DISTANCE_MAX = 180.0
TRUTHS = (55.0, 70.0, 85.0)


def draw_detected_observation(h0_true, rng, max_attempts=100000):
    """Draw host and noisy distance jointly, conditional on 0 < distance < threshold."""
    for _ in range(max_attempts):
        host_index = int(rng.choice(len(HOST_REDSHIFTS), p=HOST_WEIGHTS))
        mean = C_KM_S * HOST_REDSHIFTS[host_index] / h0_true
        distance = float(rng.normal(mean, SIGMA_MPC))
        if 0.0 < distance < DETECTION_DISTANCE_MAX:
            return host_index, distance
    raise RuntimeError("Exceeded attempts while sampling a detected observation.")


def _metrics(rows, key):
    truth = np.asarray([float(row["true_H0"]) for row in rows])
    estimates = np.asarray([float(row[f"{key}_median"]) for row in rows])
    lower68 = np.asarray([float(row[f"{key}_lower_68"]) for row in rows])
    upper68 = np.asarray([float(row[f"{key}_upper_68"]) for row in rows])
    lower95 = np.asarray([float(row[f"{key}_lower_95"]) for row in rows])
    upper95 = np.asarray([float(row[f"{key}_upper_95"]) for row in rows])
    error = estimates - truth
    return {
        "n": len(rows),
        "bias": float(error.mean()),
        "rmse": float(np.sqrt(np.mean(error**2))),
        "coverage_68": float(np.mean((truth >= lower68) & (truth <= upper68))),
        "coverage_95": float(np.mean((truth >= lower95) & (truth <= upper95))),
        "mean_width_68": float(np.mean(upper68 - lower68)),
        "mean_width_95": float(np.mean(upper95 - lower95)),
        "mean_absolute_error": float(np.mean(np.abs(error))),
    }


def run_comparison(replicates_per_truth=600, seed=20261028, grid_size=2001):
    if replicates_per_truth < 100:
        raise ValueError("Use at least 100 replicates per fixed truth.")
    if grid_size < 1001:
        raise ValueError("grid_size must be at least 1001.")
    rng = np.random.default_rng(seed)
    grid = np.linspace(25.0, 125.0, grid_size)
    rows = []
    for h0_true in TRUTHS:
        for replicate in range(replicates_per_truth):
            host_index, distance = draw_detected_observation(h0_true, rng)
            base = posterior_h0(
                grid, HOST_REDSHIFTS, distance, SIGMA_MPC, weights=HOST_WEIGHTS,
                prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
            )
            corrected = posterior_h0_selected(
                grid, HOST_REDSHIFTS, distance, SIGMA_MPC, DETECTION_DISTANCE_MAX,
                weights=HOST_WEIGHTS, prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
            )
            base_summary = posterior_summary(grid, base)
            corrected_summary = posterior_summary(grid, corrected)
            row = {
                "true_H0": h0_true,
                "replicate": replicate,
                "host_index": host_index,
                "observed_distance_mpc": distance,
            }
            for label, summary in (("uncorrected", base_summary), ("corrected", corrected_summary)):
                for metric in ("median", "lower_68", "upper_68", "lower_95", "upper_95"):
                    row[f"{label}_{metric}"] = summary[metric]
            row["uncorrected_error"] = base_summary["median"] - h0_true
            row["corrected_error"] = corrected_summary["median"] - h0_true
            row["paired_abs_error_change_corrected_minus_uncorrected"] = (
                abs(row["corrected_error"]) - abs(row["uncorrected_error"])
            )
            rows.append(row)

    summary_rows = []
    for truth in TRUTHS:
        cell = [row for row in rows if float(row["true_H0"]) == truth]
        for method in ("uncorrected", "corrected"):
            summary_rows.append({
                "truth_H0": truth,
                "method": method,
                **_metrics(cell, method),
                "seed": seed,
                "replicates_per_truth": replicates_per_truth,
                "grid_size": grid_size,
                "sigma_mpc": SIGMA_MPC,
                "detection_threshold_mpc": DETECTION_DISTANCE_MAX,
            })
    for method in ("uncorrected", "corrected"):
        summary_rows.append({
            "truth_H0": "ALL_TRUTHS_DESCRIPTIVE_ONLY",
            "method": method,
            **_metrics(rows, method),
            "seed": seed,
            "replicates_per_truth": replicates_per_truth,
            "grid_size": grid_size,
            "sigma_mpc": SIGMA_MPC,
            "detection_threshold_mpc": DETECTION_DISTANCE_MAX,
        })
    return rows, summary_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates-per-truth", type=int, default=600)
    parser.add_argument("--seed", type=int, default=20261028)
    parser.add_argument("--grid-size", type=int, default=2001)
    parser.add_argument("--output-dir", default="results/selection-paired-comparison")
    args = parser.parse_args()
    rows, summaries = run_comparison(args.replicates_per_truth, args.seed, args.grid_size)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, data in (("per_run.csv", rows), ("summary.csv", summaries)):
        with (out / name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    print(f"PASS: {len(rows)} paired fixed-truth mock observations analysed twice.")
    for row in summaries:
        print(row)


if __name__ == "__main__":
    main()
