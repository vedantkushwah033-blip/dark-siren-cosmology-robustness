"""Deterministic posterior-grid stability diagnostic for the hardened synthetic model.

This is a numerical diagnostic, not a scientific result from the factorial experiment.
It separates grid discretization/support checks from changing the prior bounds.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

from dark_siren_robustness.hardened_inference import posterior_h0, posterior_summary


def evaluate_case(name, grid, prior_min, prior_max):
    posterior = posterior_h0(
        grid,
        redshifts=[0.012, 0.018, 0.026, 0.037, 0.051],
        gw_distance=165.0,
        gw_sigma=22.0,
        weights=[1.0, 1.5, 0.7, 2.0, 0.8],
        missing_weight=0.12,
        missing_redshifts=[0.009, 0.021, 0.044],
        missing_weights=[0.3, 0.5, 0.2],
        prior_min=prior_min,
        prior_max=prior_max,
    )
    summary = posterior_summary(grid, posterior)
    area = float(trapezoid(posterior, grid))
    # Measure posterior mass near the prior boundaries, not the numerical grid
    # boundaries. For a wider grid with an unchanged prior, grid-edge mass is
    # mechanically zero and would conceal possible prior-boundary concentration.
    band_width = 0.02 * (prior_max - prior_min)
    lower_mask = (grid >= prior_min) & (grid <= prior_min + band_width)
    upper_mask = (grid <= prior_max) & (grid >= prior_max - band_width)
    return {
        "case": name,
        "grid_min": float(grid[0]),
        "grid_max": float(grid[-1]),
        "grid_points": int(grid.size),
        "grid_spacing": float(grid[1] - grid[0]),
        "prior_min": float(prior_min),
        "prior_max": float(prior_max),
        "posterior_area": area,
        "lower_prior_boundary_band_mass": float(trapezoid(posterior[lower_mask], grid[lower_mask])),
        "upper_prior_boundary_band_mass": float(trapezoid(posterior[upper_mask], grid[upper_mask])),
        **summary,
    }


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="artifacts/posterior_grid_stability.csv")
    args = parser.parse_args()

    cases = [
        evaluate_case("baseline", np.linspace(50.0, 90.0, 1601), 50.0, 90.0),
        evaluate_case("baseline_finer", np.linspace(50.0, 90.0, 3201), 50.0, 90.0),
        evaluate_case("wider_grid_same_prior", np.linspace(30.0, 110.0, 6401), 50.0, 90.0),
        evaluate_case("wider_grid_wider_prior", np.linspace(30.0, 110.0, 6401), 30.0, 110.0),
    ]
    reference = next(row for row in cases if row["case"] == "baseline_finer")
    for row in cases:
        row["median_abs_diff_vs_finer_baseline"] = abs(row["median"] - reference["median"])
        row["lower68_abs_diff_vs_finer_baseline"] = abs(row["lower_68"] - reference["lower_68"])
        row["upper68_abs_diff_vs_finer_baseline"] = abs(row["upper_68"] - reference["upper_68"])
        row["lower95_abs_diff_vs_finer_baseline"] = abs(row["lower_95"] - reference["lower_95"])
        row["upper95_abs_diff_vs_finer_baseline"] = abs(row["upper_95"] - reference["upper_95"])

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)
    print(json.dumps({"output": str(output), "cases": cases}, indent=2))


if __name__ == "__main__":
    run()
