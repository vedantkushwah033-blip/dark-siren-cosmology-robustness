#!/usr/bin/env python3
"""Aggregate paired A-minus-B contrasts over interior truths with a stratified bootstrap.

Reads the raw CSV produced by run_distance_uncertainty_sensitivity.py. For each
fractional-error setting, the estimand is the equally weighted average over
interior truth cells H0={55,60,65,70,75,80,85}; replicate rows are resampled
within each truth cell, preserving the fixed simulation design.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

A = "A_fixed_realized_sigma"
B = "B_fractional_sigma_H0"
INTERIOR_TRUTHS = (55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0)
METRICS = ("coverage68_diff", "absolute_median_error_diff", "squared_median_error_diff",
           "width68_diff", "width95_diff")


def summarize(path: Path, seed: int = 4182026, bootstrap_resamples: int = 10000) -> dict:
    if seed < 0 or bootstrap_resamples < 100:
        raise ValueError("seed must be non-negative and bootstrap_resamples >= 100")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    indexed = {}
    for row in rows:
        key = (float(row["true_h0"]), float(row["fractional_error"]),
               int(row["replicate"]), row["arm"])
        if key in indexed:
            raise ValueError(f"duplicate raw record: {key}")
        indexed[key] = row
    truths = sorted({key[0] for key in indexed})
    fractions = sorted({key[1] for key in indexed})
    if not set(INTERIOR_TRUTHS).issubset(truths):
        raise ValueError("raw records do not contain every required interior truth")
    expected_arms = {A, B}
    differences = {}
    for truth in truths:
        for frac in fractions:
            reps = sorted({key[2] for key in indexed if key[0] == truth and key[1] == frac})
            for rep in reps:
                arm_rows = {arm: indexed.get((truth, frac, rep, arm)) for arm in expected_arms}
                if any(row is None for row in arm_rows.values()):
                    raise ValueError(
                        f"missing paired A/B record at H0={truth}, f={frac}, replicate={rep}"
                    )
                a, b = arm_rows[A], arm_rows[B]
                for column in ("host_redshift", "standard_normal", "true_distance_mpc",
                               "observed_distance_mpc", "generated_sigma_mpc"):
                    if float(a[column]) != float(b[column]):
                        raise ValueError(
                            f"A/B input pairing mismatch for {column} at {truth}, {frac}, {rep}"
                        )
                h = truth
                metric_values = {
                    "coverage68_diff": float(float(a["lower_68"]) <= h <= float(a["upper_68"])) -
                                        float(float(b["lower_68"]) <= h <= float(b["upper_68"])),
                    "absolute_median_error_diff": (
                        abs(float(a["median"]) - h) - abs(float(b["median"]) - h)
                    ),
                    "squared_median_error_diff": (
                        (float(a["median"]) - h) ** 2 - (float(b["median"]) - h) ** 2
                    ),
                    "width68_diff": (float(a["upper_68"]) - float(a["lower_68"])) -
                                    (float(b["upper_68"]) - float(b["lower_68"])),
                    "width95_diff": (float(a["upper_95"]) - float(a["lower_95"])) -
                                    (float(b["upper_95"]) - float(b["lower_95"])),
                }
                differences[(truth, frac, rep)] = metric_values

    results = []
    for frac in fractions:
        cell_rows = {
            truth: [
                differences[(truth, frac, rep)]
                for rep in sorted(rep for (t, f, rep) in differences if t == truth and f == frac)
            ]
            for truth in INTERIOR_TRUTHS
        }
        if any(len(cell_rows[t]) == 0 for t in INTERIOR_TRUTHS):
            raise ValueError(f"missing interior truth-cell records for fractional error {frac}")
        rng = np.random.default_rng(seed + int(round(frac * 1000)))
        for metric in METRICS:
            per_cell = [np.asarray([row[metric] for row in cell_rows[t]], dtype=float)
                        for t in INTERIOR_TRUTHS]
            point = float(np.mean([values.mean() for values in per_cell]))
            boot = np.empty(bootstrap_resamples, dtype=float)
            for index in range(bootstrap_resamples):
                boot[index] = np.mean([
                    values[rng.integers(0, len(values), size=len(values))].mean()
                    for values in per_cell
                ])
            results.append({
                "fractional_error": frac,
                "truth_cells": list(INTERIOR_TRUTHS),
                "replicates_per_cell": len(per_cell[0]),
                "estimand": (
                    "equal-weight mean of paired A-minus-B replicate differences over "
                    "the seven fixed interior H0 cells"
                ),
                "contrast": metric,
                "direction": "A_minus_B",
                "mean_difference": point,
                "stratified_bootstrap_95_low": float(np.quantile(boot, 0.025)),
                "stratified_bootstrap_95_high": float(np.quantile(boot, 0.975)),
                "bootstrap_resamples": bootstrap_resamples,
                "bootstrap_seed": seed + int(round(frac * 1000)),
            })
    return {
        "status": "PASS",
        "scope": "synthetic known-host low-redshift sensitivity diagnostic only",
        "raw_record_count": len(rows),
        "interior_truth_cells": list(INTERIOR_TRUTHS),
        "all_a_b_pairs_checked": len(differences),
        "interior_a_b_pairs_checked": sum(1 for key in differences if key[0] in INTERIOR_TRUTHS),
        "pairing_verified": True,
        "results": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("replicates_csv", type=Path)
    parser.add_argument("--seed", type=int, default=4182026)
    parser.add_argument("--bootstrap-resamples", type=int, default=10000)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = summarize(args.replicates_csv, args.seed, args.bootstrap_resamples)
    output = args.output or args.replicates_csv.with_name("paired_contrast_stratified_summary.json")
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "results"}, indent=2))
    for item in report["results"]:
        print(f"f={item['fractional_error']:.2f} {item['contrast']}: "
              f"{item['mean_difference']:.6g} "
              f"[{item['stratified_bootstrap_95_low']:.6g}, "
              f"{item['stratified_bootstrap_95_high']:.6g}]")


if __name__ == "__main__":
    main()
