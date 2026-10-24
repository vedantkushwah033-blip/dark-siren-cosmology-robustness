"""Compare toy selection-normalized inference across lower detection thresholds.

This is a deterministic sensitivity audit, not a realistic GW selection model.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np
from dark_siren_robustness.hardened_inference import posterior_h0, posterior_summary
from dark_siren_robustness.selection import detection_efficiency_h0, posterior_h0_selected

HOST_REDSHIFTS = np.array([0.02, 0.04])
HOST_WEIGHTS = np.array([0.70, 0.30])
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
SIGMA_MPC = 20.0
DETECTION_MAX = 180.0
LOWER_THRESHOLDS = (0.0, 10.0, 20.0, 40.0)
OBSERVED_DISTANCES = (50.0, 90.0, 140.0)

def run_audit(grid_size=4001):
    if grid_size < 101:
        raise ValueError("grid_size must be >= 101.")
    grid = np.linspace(25.0, 125.0, grid_size)
    prior_support = (grid >= PRIOR_MIN) & (grid <= PRIOR_MAX)
    rows = []
    for lower in LOWER_THRESHOLDS:
        beta = detection_efficiency_h0(
            grid, HOST_REDSHIFTS, SIGMA_MPC, DETECTION_MAX,
            weights=HOST_WEIGHTS, detection_distance_min=lower,
        )
        for distance in OBSERVED_DISTANCES:
            uncorrected = posterior_h0(
                grid, HOST_REDSHIFTS, distance, SIGMA_MPC,
                weights=HOST_WEIGHTS, prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
            )
            corrected = posterior_h0_selected(
                grid, HOST_REDSHIFTS, distance, SIGMA_MPC, DETECTION_MAX,
                weights=HOST_WEIGHTS, prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
                detection_distance_min=lower,
            )
            for label, posterior in (("uncorrected", uncorrected), ("selection_normalized", corrected)):
                row = {
                    "observed_distance_mpc": distance, "method": label,
                    "grid_size": grid_size, "detection_min_mpc": lower,
                    "detection_max_mpc": DETECTION_MAX,
                    "beta_min_prior_support": float(beta[prior_support].min()),
                    "beta_max_prior_support": float(beta[prior_support].max()),
                }
                row.update(posterior_summary(grid, posterior))
                rows.append(row)
    return rows

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-size", type=int, default=4001)
    parser.add_argument("--output-dir", default="results/lower-threshold-sensitivity")
    args = parser.parse_args()
    rows = run_audit(args.grid_size)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    with (out / "lower_threshold_sensitivity.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"PASS: {len(rows)} posterior summaries across {len(LOWER_THRESHOLDS)} lower thresholds.")
    for row in rows:
        print(row)
if __name__ == "__main__":
    main()
