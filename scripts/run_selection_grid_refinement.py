"""Deterministic grid-refinement audit for the toy selection-conditioned posterior."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path
import numpy as np
from dark_siren_robustness.hardened_inference import posterior_summary
from dark_siren_robustness.selection import posterior_h0_selected

HOST_REDSHIFTS = np.array([0.02, 0.04])
HOST_WEIGHTS = np.array([0.70, 0.30])
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
SIGMA_MPC, DETECTION_MAX = 20.0, 180.0
OBSERVED_DISTANCES = (70.0, 100.0, 140.0)

def run_audit(grid_sizes=(1001, 4001, 8001)):
    if len(grid_sizes) < 2 or any(int(n) < 101 for n in grid_sizes):
        raise ValueError("Provide at least two grid sizes, each >= 101.")
    if sorted(set(grid_sizes)) != list(grid_sizes):
        raise ValueError("Grid sizes must be unique and increasing.")
    rows = []
    for distance in OBSERVED_DISTANCES:
        summaries = {}
        for n in grid_sizes:
            grid = np.linspace(25.0, 125.0, int(n))
            posterior = posterior_h0_selected(
                grid, HOST_REDSHIFTS, distance, SIGMA_MPC, DETECTION_MAX,
                weights=HOST_WEIGHTS, prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
            )
            summaries[n] = posterior_summary(grid, posterior)
        reference = summaries[grid_sizes[-1]]
        for n in grid_sizes:
            row = {"observed_distance_mpc": distance, "grid_size": n}
            row.update(summaries[n])
            row["abs_median_diff_vs_finest"] = abs(row["median"] - reference["median"])
            row["max_abs_quantile_diff_vs_finest"] = max(
                abs(row[k] - reference[k]) for k in ("lower_68", "upper_68", "lower_95", "upper_95")
            )
            rows.append(row)
    return rows

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="results/selection-grid-refinement")
    args = parser.parse_args()
    rows = run_audit()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    with (out / "grid_refinement.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"PASS: {len(rows)} posterior summaries across {len(OBSERVED_DISTANCES)} fixed observations.")
    for row in rows:
        print(row)
if __name__ == "__main__":
    main()
