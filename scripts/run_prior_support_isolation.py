"""Isolate prior-boundary effects in a matched single-host synthetic model.

The data-generating process and likelihood both use d ~ Normal(c z / H0, sigma).
This is a controlled statistical diagnostic, not a realistic GW event simulation.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary

H0_VALUES = (55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0)
PRIORS = {"narrow_50_90": (50.0, 90.0), "broad_30_120": (30.0, 120.0)}


def run_audit(replicates: int = 500, seed: int = 20261024, grid_size: int = 5001):
    if replicates < 100:
        raise ValueError("Use at least 100 replicates per truth/prior cell.")
    if grid_size < 1001:
        raise ValueError("grid_size must be at least 1001.")
    grid = np.linspace(25.0, 125.0, grid_size)
    rng = np.random.default_rng(seed)
    rows = []
    z = 0.03
    for h0_true in H0_VALUES:
        true_distance = C_KM_S * z / h0_true
        sigma = 0.10 * true_distance
        for replicate in range(replicates):
            observed_distance = float(rng.normal(true_distance, sigma))
            while observed_distance <= 0:
                observed_distance = float(rng.normal(true_distance, sigma))
            for prior_name, (prior_min, prior_max) in PRIORS.items():
                posterior = posterior_h0(
                    grid, [z], observed_distance, sigma,
                    prior_min=prior_min, prior_max=prior_max,
                )
                summary = posterior_summary(grid, posterior)
                lo68, hi68 = summary["lower_68"], summary["upper_68"]
                lo95, hi95 = summary["lower_95"], summary["upper_95"]
                rows.append({
                    "true_H0": h0_true,
                    "replicate": replicate,
                    "prior": prior_name,
                    "observed_distance_mpc": observed_distance,
                    "posterior_median": summary["median"],
                    "median_error": summary["median"] - h0_true,
                    "lower_68": lo68,
                    "upper_68": hi68,
                    "width_68": hi68 - lo68,
                    "covers_68": int(lo68 <= h0_true <= hi68),
                    "lower_95": lo95,
                    "upper_95": hi95,
                    "width_95": hi95 - lo95,
                    "covers_95": int(lo95 <= h0_true <= hi95),
                    "touches_prior_boundary_68": int(
                        lo68 <= prior_min + 0.15 or hi68 >= prior_max - 0.15
                    ),
                    "touches_prior_boundary_95": int(
                        lo95 <= prior_min + 0.15 or hi95 >= prior_max - 0.15
                    ),
                })
    summary_rows = []
    for h0_true in H0_VALUES:
        for prior_name, (prior_min, prior_max) in PRIORS.items():
            cell = [r for r in rows if r["true_H0"] == h0_true and r["prior"] == prior_name]
            n = len(cell)
            c68 = float(np.mean([r["covers_68"] for r in cell]))
            c95 = float(np.mean([r["covers_95"] for r in cell]))
            errors = np.asarray([r["median_error"] for r in cell])
            summary_rows.append({
                "true_H0": h0_true,
                "prior": prior_name,
                "prior_min": prior_min,
                "prior_max": prior_max,
                "replicates": n,
                "mean_signed_error": float(errors.mean()),
                "rmse_median": float(np.sqrt(np.mean(errors ** 2))),
                "coverage_68": c68,
                "coverage_68_mc_se": float(np.sqrt(c68 * (1 - c68) / n)),
                "coverage_95": c95,
                "coverage_95_mc_se": float(np.sqrt(c95 * (1 - c95) / n)),
                "mean_width_68": float(np.mean([r["width_68"] for r in cell])),
                "mean_width_95": float(np.mean([r["width_95"] for r in cell])),
                "boundary_touch_rate_68": float(np.mean([r["touches_prior_boundary_68"] for r in cell])),
                "boundary_touch_rate_95": float(np.mean([r["touches_prior_boundary_95"] for r in cell])),
            })
    return rows, summary_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261024)
    parser.add_argument("--grid-size", type=int, default=5001)
    parser.add_argument("--output-dir", default="results/prior-support-isolation")
    args = parser.parse_args()
    rows, summaries = run_audit(args.replicates, args.seed, args.grid_size)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, data in (("per_run.csv", rows), ("summary.csv", summaries)):
        with (out / name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    print(f"Validated design: {len(H0_VALUES)} truths × {len(PRIORS)} priors × {args.replicates} replicates = {len(rows)} fits.")
    for row in summaries:
        print(row)


if __name__ == "__main__":
    main()
