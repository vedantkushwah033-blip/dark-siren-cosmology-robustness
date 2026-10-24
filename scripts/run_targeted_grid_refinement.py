"""Run high-resolution posterior summaries for cases flagged by grid accuracy audit.

Input is the evaluator's case_accuracy_assessment.csv. Each unique simulated
case flagged for any posterior summary is recomputed at 8,001--64,001 points,
holding observation, redshifts, noise model, prior, and weighting fixed.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary

Z1 = 0.02
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
GRID_SIZES = (8001, 16001, 32001, 64001)
IDENTITY = ("H0_true", "z2", "noise_fraction", "replicate", "inference_weighting")
SUMMARY_KEYS = {
    "posterior_median": "median",
    "lower_68": "lower_68",
    "upper_68": "upper_68",
    "lower_95": "lower_95",
    "upper_95": "upper_95",
}


def run(flagged_csv: Path, per_run_csv: Path):
    with flagged_csv.open(newline="", encoding="utf-8") as f:
        flagged = list(csv.DictReader(f))
    wanted = {tuple(row[k] for k in IDENTITY) for row in flagged if row["flag_for_targeted_refinement"].lower() == "true"}
    with per_run_csv.open(newline="", encoding="utf-8") as f:
        source = list(csv.DictReader(f))
    baseline = {}
    for row in source:
        key = tuple(row[k] for k in IDENTITY)
        if key in wanted and int(row["grid_size"]) == 4001:
            baseline[key] = row
    missing = wanted - set(baseline)
    if missing:
        raise ValueError(f"Flagged cases missing 4001-point baseline rows: {len(missing)}")
    output = []
    for key in sorted(wanted):
        row = baseline[key]
        h0, z2, noise = float(row["H0_true"]), float(row["z2"]), float(row["noise_fraction"])
        d_obs = float(row["observed_distance_mpc"])
        redshifts = np.array([Z1, z2], dtype=float)
        sigma = noise * C_KM_S * 0.03 / h0
        weights = np.array([0.7, 0.3]) if row["inference_weighting"] == "correct" else np.array([0.3, 0.7])
        for size in GRID_SIZES:
            grid = np.linspace(25.0, 125.0, size)
            density = posterior_h0(grid, redshifts, d_obs, sigma, weights=weights,
                                   prior_min=PRIOR_MIN, prior_max=PRIOR_MAX)
            summary = posterior_summary(grid, density)
            prior_support = (grid >= PRIOR_MIN) & (grid <= PRIOR_MAX)
            independently_integrated_area = float(trapezoid(density[prior_support], grid[prior_support]))
            record = dict(zip(IDENTITY, key))
            record.update({
                "observed_distance_mpc": d_obs,
                "grid_size": size,
                "posterior_area_on_prior_support": independently_integrated_area,
                "posterior_median": summary["median"],
                "lower_68": summary["lower_68"],
                "upper_68": summary["upper_68"],
                "lower_95": summary["lower_95"],
                "upper_95": summary["upper_95"],
            })
            output.append(record)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flagged-cases", required=True)
    parser.add_argument("--per-run", required=True)
    parser.add_argument("--output", default="results/grid-refinement-audit-full/targeted_refinement.csv")
    args = parser.parse_args()
    rows = run(Path(args.flagged_cases), Path(args.per_run))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else list(IDENTITY) + ["observed_distance_mpc", "grid_size", *SUMMARY_KEYS]
    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Targeted refinement: {len(rows)} rows for {len(rows)//len(GRID_SIZES) if rows else 0} flagged cases; output={output}")


if __name__ == "__main__":
    main()
