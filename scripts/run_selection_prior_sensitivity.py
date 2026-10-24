"""Prior-sensitivity audit for paired toy selection-normalized inference."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary
from dark_siren_robustness.selection import posterior_h0_selected

HOST_REDSHIFTS = np.array([0.02, 0.04])
HOST_WEIGHTS = np.array([0.70, 0.30])
SIGMA_MPC = 20.0
DETECTION_DISTANCE_MAX = 180.0
TRUTHS = (55.0, 70.0, 85.0)
PRIORS = ((50.0, 90.0), (40.0, 100.0), (30.0, 120.0))


def draw_detected_observation(h0_true, rng, max_attempts=100000):
    for _ in range(max_attempts):
        host = int(rng.choice(len(HOST_REDSHIFTS), p=HOST_WEIGHTS))
        mean = C_KM_S * HOST_REDSHIFTS[host] / h0_true
        distance = float(rng.normal(mean, SIGMA_MPC))
        if 0.0 < distance < DETECTION_DISTANCE_MAX:
            return host, distance
    raise RuntimeError("Exceeded attempts while sampling detected observation.")


def summarize(rows, method):
    truth = np.asarray([float(r["true_H0"]) for r in rows])
    med = np.asarray([float(r[f"{method}_median"]) for r in rows])
    lo68 = np.asarray([float(r[f"{method}_lower_68"]) for r in rows])
    hi68 = np.asarray([float(r[f"{method}_upper_68"]) for r in rows])
    lo95 = np.asarray([float(r[f"{method}_lower_95"]) for r in rows])
    hi95 = np.asarray([float(r[f"{method}_upper_95"]) for r in rows])
    err = med - truth
    return {
        "n": len(rows),
        "bias": float(np.mean(err)),
        "rmse": float(np.sqrt(np.mean(err**2))),
        "mae": float(np.mean(np.abs(err))),
        "coverage_68": float(np.mean((truth >= lo68) & (truth <= hi68))),
        "coverage_95": float(np.mean((truth >= lo95) & (truth <= hi95))),
        "mean_width_68": float(np.mean(hi68 - lo68)),
        "mean_width_95": float(np.mean(hi95 - lo95)),
    }


def run_audit(replicates_per_truth=600, seed=20261029, grid_size=4001):
    if replicates_per_truth < 100:
        raise ValueError("Use at least 100 replicates per truth.")
    if grid_size < 1001:
        raise ValueError("grid_size must be at least 1001.")
    rng = np.random.default_rng(seed)
    grid = np.linspace(25.0, 125.0, grid_size)
    rows = []
    for truth in TRUTHS:
        for replicate in range(replicates_per_truth):
            host, distance = draw_detected_observation(truth, rng)
            for prior_min, prior_max in PRIORS:
                base = posterior_h0(grid, HOST_REDSHIFTS, distance, SIGMA_MPC,
                                    weights=HOST_WEIGHTS, prior_min=prior_min, prior_max=prior_max)
                corrected = posterior_h0_selected(
                    grid, HOST_REDSHIFTS, distance, SIGMA_MPC, DETECTION_DISTANCE_MAX,
                    weights=HOST_WEIGHTS, prior_min=prior_min, prior_max=prior_max)
                row = {
                    "true_H0": truth, "replicate": replicate, "host_index": host,
                    "observed_distance_mpc": distance, "prior_min": prior_min,
                    "prior_max": prior_max,
                }
                for name, posterior in (("uncorrected", base), ("corrected", corrected)):
                    sm = posterior_summary(grid, posterior)
                    for key in ("median", "lower_68", "upper_68", "lower_95", "upper_95"):
                        row[f"{name}_{key}"] = sm[key]
                row["paired_abs_error_change"] = abs(row["corrected_median"] - truth) - abs(row["uncorrected_median"] - truth)
                rows.append(row)
    summaries = []
    for truth in TRUTHS:
        for prior_min, prior_max in PRIORS:
            cell = [r for r in rows if r["true_H0"] == truth and r["prior_min"] == prior_min]
            for method in ("uncorrected", "corrected"):
                summaries.append({
                    "true_H0": truth, "prior_min": prior_min, "prior_max": prior_max,
                    "method": method, **summarize(cell, method), "seed": seed,
                    "replicates_per_truth": replicates_per_truth, "grid_size": grid_size,
                    "sigma_mpc": SIGMA_MPC, "detection_distance_max_mpc": DETECTION_DISTANCE_MAX,
                })
            changes = np.asarray([r["paired_abs_error_change"] for r in cell])
            summaries.append({
                "true_H0": truth, "prior_min": prior_min, "prior_max": prior_max,
                "method": "paired_corrected_minus_uncorrected_absolute_error",
                "n": len(cell), "mean_paired_change": float(changes.mean()),
                "median_paired_change": float(np.median(changes)),
                "paired_change_sd": float(changes.std(ddof=1)),
                "paired_change_mcse": float(changes.std(ddof=1) / np.sqrt(len(changes))),
                "seed": seed, "replicates_per_truth": replicates_per_truth,
                "grid_size": grid_size, "sigma_mpc": SIGMA_MPC,
                "detection_distance_max_mpc": DETECTION_DISTANCE_MAX,
            })
    return rows, summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates-per-truth", type=int, default=600)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--grid-size", type=int, default=4001)
    parser.add_argument("--output-dir", default="results/selection-prior-sensitivity")
    args = parser.parse_args()
    rows, summaries = run_audit(args.replicates_per_truth, args.seed, args.grid_size)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, data in (("per_run.csv", rows), ("summary.csv", summaries)):
        with (out / name).open("w", newline="", encoding="utf-8") as f:
            fieldnames = list(dict.fromkeys(key for row in data for key in row))
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="raise")
            writer.writeheader()
            writer.writerows(data)
    print(f"PASS: {len(rows)} prior-matched paired rows; {len(summaries)} summaries.")
    for row in summaries:
        print(row)


if __name__ == "__main__":
    main()
