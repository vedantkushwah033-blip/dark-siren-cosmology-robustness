#!/usr/bin/env python3
"""Independently recompute raw distance-uncertainty posteriors on a finer H0 grid."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

C_KM_S = 299_792.458
SUMMARY_COLUMNS = ("median", "lower_68", "upper_68", "lower_95", "upper_95")
QUANTILES = (0.50, 0.16, 0.84, 0.025, 0.975)


def _posterior_quantiles(row: dict[str, str], grid: np.ndarray) -> dict[str, float]:
    z = float(row["host_redshift"])
    observed = float(row["observed_distance_mpc"])
    fractional_error = float(row["fractional_error"])
    sigma_recorded = float(row["generated_sigma_mpc"])
    mean = C_KM_S * z / grid
    if row["arm"] == "B_fractional_sigma_H0":
        sigma = fractional_error * mean
    elif row["arm"] in ("A_fixed_realized_sigma", "C_fixed_reference_sigma"):
        sigma = np.full_like(grid, sigma_recorded)
    else:
        raise ValueError(f"unknown likelihood arm: {row['arm']}")
    if np.any(~np.isfinite(sigma)) or np.any(sigma <= 0):
        raise ValueError("non-finite or non-positive likelihood uncertainty")
    log_likelihood = (
        -0.5 * ((observed - mean) / sigma) ** 2
        - np.log(sigma)
        - 0.5 * math.log(2.0 * math.pi)
    )
    density = np.exp(log_likelihood - np.max(log_likelihood))
    area = float(np.trapezoid(density, grid))
    if not math.isfinite(area) or area <= 0:
        raise ValueError("invalid fine-grid posterior normalization")
    density /= area
    cdf = np.r_[0.0, np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(grid))]
    cdf /= cdf[-1]
    values = [float(np.interp(q, cdf, grid)) for q in QUANTILES]
    # Flag probability near either edge using one spacing of the reference 1,601-point grid.
    # This is a diagnostic only: it does not fail the audit by itself.
    coarse_spacing = 40.0 / 1600.0
    edge_mask = (grid <= grid[0] + coarse_spacing) | (grid >= grid[-1] - coarse_spacing)
    boundary_mass = float(np.trapezoid(density[edge_mask], grid[edge_mask]))
    result = dict(zip(SUMMARY_COLUMNS, values))
    result["boundary_mass_within_one_coarse_spacing"] = boundary_mass
    return result


def check_grid_stability(raw_csv: Path, fine_grid_points: int = 3201,
                         max_allowed_change: float = 1e-3) -> dict:
    if fine_grid_points < 1601 or fine_grid_points % 2 == 0:
        raise ValueError("fine_grid_points must be an odd integer >= 1601")
    if not math.isfinite(max_allowed_change) or max_allowed_change <= 0:
        raise ValueError("max_allowed_change must be finite and positive")
    grid = np.linspace(50.0, 90.0, fine_grid_points)
    changes = {column: [] for column in SUMMARY_COLUMNS}
    boundary_masses = []
    boundary_mass_threshold = 0.05
    coverage = {}
    seen_keys = set()
    boundary_flags = []
    arms_by_replicate = {}
    replicates_by_cell = {}
    expected_arms = {
        "A_fixed_realized_sigma",
        "B_fractional_sigma_H0",
        "C_fixed_reference_sigma",
    }
    row_count = 0
    with raw_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"true_h0", "arm", "replicate", "lower_68", "upper_68", *SUMMARY_COLUMNS,
                    "host_redshift", "observed_distance_mpc", "fractional_error",
                    "generated_sigma_mpc"}
        if not required.issubset(reader.fieldnames or []):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"raw CSV missing columns: {missing}")
        for row in reader:
            row_count += 1
            truth = float(row["true_h0"])
            fractional_error = float(row["fractional_error"])
            replicate = int(row["replicate"])
            arm = row["arm"]
            key = (truth, fractional_error, replicate, arm)
            if key in seen_keys:
                raise ValueError(f"duplicate raw record key: {key}")
            seen_keys.add(key)
            replicate_key = (truth, fractional_error, replicate)
            arms_by_replicate.setdefault(replicate_key, set()).add(arm)
            cell_key = (truth, fractional_error)
            replicates_by_cell.setdefault(cell_key, set()).add(replicate)
            if arm not in coverage:
                coverage[arm] = {"coarse_successes": 0, "fine_successes": 0, "n": 0}
            fine = _posterior_quantiles(row, grid)
            boundary_mass = fine["boundary_mass_within_one_coarse_spacing"]
            boundary_masses.append(boundary_mass)
            if boundary_mass > boundary_mass_threshold:
                boundary_flags.append({
                    "true_h0": truth,
                    "fractional_error": fractional_error,
                    "replicate": replicate,
                    "arm": arm,
                    "boundary_mass": boundary_mass,
                })
            for column in SUMMARY_COLUMNS:
                coarse = float(row[column])
                delta = abs(coarse - fine[column])
                changes[column].append(delta)
                if delta > max_allowed_change:
                    raise ValueError(
                        f"grid change {delta} exceeds tolerance {max_allowed_change} "
                        f"for row {row_count}, column {column}"
                    )
            coarse_covered = float(row["lower_68"]) <= truth <= float(row["upper_68"])
            fine_covered = fine["lower_68"] <= truth <= fine["upper_68"]
            coverage[arm]["coarse_successes"] += int(coarse_covered)
            coverage[arm]["fine_successes"] += int(fine_covered)
            coverage[arm]["n"] += 1
    if row_count == 0:
        raise ValueError("raw CSV is empty")
    for key, observed_arms in arms_by_replicate.items():
        if observed_arms != expected_arms:
            missing = sorted(expected_arms - observed_arms)
            unexpected = sorted(observed_arms - expected_arms)
            raise ValueError(
                f"incomplete or unexpected likelihood arms for replicate {key}: "
                f"missing={missing}, unexpected={unexpected}"
            )
    replicate_sets = list(replicates_by_cell.values())
    if replicate_sets:
        reference_replicates = replicate_sets[0]
        if reference_replicates != set(range(max(reference_replicates) + 1)):
            raise ValueError("replicate IDs must be contiguous from zero within each truth/error cell")
        for cell, replicates in replicates_by_cell.items():
            if replicates != reference_replicates:
                raise ValueError(
                    f"inconsistent replicate IDs across truth/error cells at {cell}: "
                    f"expected={sorted(reference_replicates)}, observed={sorted(replicates)}"
                )
    for arm, item in coverage.items():
        item["coarse_coverage_68"] = item["coarse_successes"] / item["n"]
        item["fine_coverage_68"] = item["fine_successes"] / item["n"]
        item["coverage_unchanged"] = item["coarse_successes"] == item["fine_successes"]
        if not item["coverage_unchanged"]:
            raise ValueError(f"68% coverage changed under grid refinement for arm {arm}")
    quantile_changes = {}
    for column, values in changes.items():
        arr = np.asarray(values, dtype=float)
        quantile_changes[column] = {
            "median_absolute_change": float(np.median(arr)),
            "p95_absolute_change": float(np.quantile(arr, 0.95)),
            "maximum_absolute_change": float(np.max(arr)),
        }
    return {
        "status": "PASS",
        "scope": "synthetic known-host low-redshift sensitivity diagnostic only",
        "raw_record_count": row_count,
        "unique_record_key_count": len(seen_keys),
        "fine_grid_points": fine_grid_points,
        "fine_grid_support": [50.0, 90.0],
        "max_allowed_quantile_change": max_allowed_change,
        "quantile_changes": quantile_changes,
        "boundary_mass_diagnostic": {
            "definition": "posterior probability within one 1,601-point-grid spacing (0.025 H0 units) of either endpoint",
            "threshold": boundary_mass_threshold,
            "flagged_record_count": len(boundary_flags),
            "maximum_mass": float(max(boundary_masses)),
            "median_mass": float(np.median(boundary_masses)),
            "flags_are_diagnostic_not_automatic_failures": True,
            "flagged_records": boundary_flags,
        },
        "coverage_68_by_arm": coverage,
        "method": (
            "independent Gaussian-likelihood and trapezoidal-grid reconstruction "
            "from raw CSV; does not import experiment code"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("replicates_csv", type=Path)
    parser.add_argument("--fine-grid-points", type=int, default=3201)
    parser.add_argument("--max-allowed-change", type=float, default=1e-3)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = check_grid_stability(
        args.replicates_csv, args.fine_grid_points, args.max_allowed_change
    )
    output = args.output or args.replicates_csv.with_name("grid_stability_report.json")
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
