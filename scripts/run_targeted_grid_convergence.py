"""Reproduce the targeted upper-95% posterior-grid convergence case.

This is a deterministic numerical diagnostic for one synthetic observation. It is
not an observational H0 measurement and does not replace the full-design audit.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import (
    C_KM_S,
    posterior_h0,
    posterior_summary,
)

GRID_SIZES = (1001, 2001, 4001, 8001, 16001, 32001, 64001)
H0_TRUE = 52.0
Z1, Z2 = 0.02, 0.04
NOISE_FRACTION = 0.03
OBSERVED_DISTANCE_MPC = 119.47415214556634
WEIGHTS = np.array([0.3, 0.7])
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0


def calculate_convergence():
    """Return posterior summaries and successive-grid changes for the fixed case."""
    sigma = NOISE_FRACTION * C_KM_S * 0.03 / H0_TRUE
    rows = []
    previous = None
    for size in GRID_SIZES:
        grid = np.linspace(25.0, 125.0, size)
        density = posterior_h0(
            grid,
            np.array([Z1, Z2]),
            OBSERVED_DISTANCE_MPC,
            sigma,
            weights=WEIGHTS,
            prior_min=PRIOR_MIN,
            prior_max=PRIOR_MAX,
        )
        summary = posterior_summary(grid, density)
        prior_support = (grid >= PRIOR_MIN) & (grid <= PRIOR_MAX)
        independently_integrated_area = float(trapezoid(density[prior_support], grid[prior_support]))
        row = {
            "grid_size": size,
            "H0_true": H0_TRUE,
            "z1": Z1,
            "z2": Z2,
            "noise_fraction": NOISE_FRACTION,
            "observed_distance_mpc": OBSERVED_DISTANCE_MPC,
            "inference_weighting": "reversed_0.3_0.7",
            "posterior_area_on_prior_support": independently_integrated_area,
            "posterior_median": summary["median"],
            "lower_68": summary["lower_68"],
            "upper_68": summary["upper_68"],
            "lower_95": summary["lower_95"],
            "upper_95": summary["upper_95"],
        }
        if previous is None:
            row["abs_upper95_change_from_previous_grid"] = ""
        else:
            row["abs_upper95_change_from_previous_grid"] = abs(
                row["upper_95"] - previous["upper_95"]
            )
        rows.append(row)
        previous = row
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="results/targeted-grid-convergence.csv",
        help="Output CSV path (default: results/targeted-grid-convergence.csv)",
    )
    args = parser.parse_args()
    rows = calculate_convergence()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print("Targeted synthetic-case grid convergence:")
    for row in rows:
        print(
            f"grid={row['grid_size']:>5} "
            f"median={row['posterior_median']:.8f} "
            f"upper95={row['upper_95']:.8f} "
            f"successive_upper95_change="
            f"{row['abs_upper95_change_from_previous_grid']}"
        )
    print(f"Wrote {len(rows)} rows to {output}")


if __name__ == "__main__":
    main()
