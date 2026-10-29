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
    left_mask = grid <= grid[0] + coarse_spacing
    right_mask = grid >= grid[-1] - coarse_spacing
    boundary_mass = float(
        np.trapezoid(density[left_mask], grid[left_mask])
        + np.trapezoid(density[right_mask], grid[right_mask])
    )
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
    grid_violations = []
    coverage_changes = []
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
                    grid_violations.append({
                        "row_number": row_count,
                        "record_key": {
                            "true_h0": truth,
                            "fractional_error": fractional_error,
                            "replicate": replicate,
                            "arm": arm,
                        },
                        "summary": column,
                        "absolute_change": delta,
                        "tolerance": max_allowed_change,
                    })
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
            coverage_changes.append({
                "arm": arm,
                "coarse_successes": item["coarse_successes"],
                "fine_successes": item["fine_successes"],
                "n": item["n"],
            })
    quantile_changes = {}
    for column, values in changes.items():
        arr = np.asarray(values, dtype=float)
        quantile_changes[column] = {
            "median_absolute_change": float(np.median(arr)),
            "p95_absolute_change": float(np.quantile(arr, 0.95)),
            "maximum_absolute_change": float(np.max(arr)),
        }
    return {
        "status": "FAIL" if grid_violations or coverage_changes else "PASS",
        "scope": "synthetic known-host low-redshift sensitivity diagnostic only",
        "raw_record_count": row_count,
        "unique_record_key_count": len(seen_keys),
        "fine_grid_points": fine_grid_points,
        "fine_grid_support": [50.0, 90.0],
        "max_allowed_quantile_change": max_allowed_change,
        "quantile_changes": quantile_changes,
        "grid_violation_count": len(grid_violations),
        "grid_violations": grid_violations,
        "coverage_changed_arms": coverage_changes,
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



def compare_grid_resolutions(
    raw_csv: Path,
    grid_points: tuple[int, ...] = (801, 1601, 3201, 6401),
    median_abs_tolerance: float = 0.02,
    median_relative_tolerance: float = 5e-4,
    interval_relative_tolerance: float = 1e-3,
) -> dict:
    """Compare independently reconstructed posterior summaries across grid resolutions.

    Tolerances are reporting thresholds, not universal scientific guarantees. The
    report preserves all adjacent-resolution comparisons so users can assess
    whether apparent convergence is monotone and whether interval membership flips.
    """
    points = tuple(int(n) for n in grid_points)
    if len(points) < 2 or any(n < 3 or n % 2 == 0 for n in points):
        raise ValueError("grid_points must contain at least two odd integers >= 3")
    if tuple(sorted(set(points))) != points:
        raise ValueError("grid_points must be strictly increasing and unique")
    tolerances = (median_abs_tolerance, median_relative_tolerance, interval_relative_tolerance)
    if any(not math.isfinite(x) or x <= 0 for x in tolerances):
        raise ValueError("all convergence tolerances must be finite and positive")

    grids = {n: np.linspace(50.0, 90.0, n) for n in points}
    summaries = {n: [] for n in points}
    row_metadata = []
    required = {
        "true_h0", "arm", "replicate", "fractional_error",
        "host_redshift", "observed_distance_mpc", "generated_sigma_mpc",
        "median", "lower_68", "upper_68", "lower_95", "upper_95",
    }
    with raw_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"raw CSV missing columns: {missing}")
        for row in reader:
            truth = float(row["true_h0"])
            if not math.isfinite(truth):
                raise ValueError("true_h0 contains a non-finite value")
            row_metadata.append({
                "truth": truth,
                "arm": row["arm"],
                "key": (truth, float(row["fractional_error"]), int(row["replicate"]), row["arm"]),
            })
            for n in points:
                summaries[n].append(_posterior_quantiles(row, grids[n]))
    if not row_metadata:
        raise ValueError("raw CSV is empty")

    comparisons = []
    summary_columns = ("median", "lower_68", "upper_68", "lower_95", "upper_95")
    for coarse_n, fine_n in zip(points[:-1], points[1:]):
        by_summary = {}
        for column in summary_columns:
            absolute = np.asarray([
                abs(a[column] - b[column])
                for a, b in zip(summaries[coarse_n], summaries[fine_n])
            ], dtype=float)
            relative = absolute / np.maximum(
                np.asarray([abs(b[column]) for b in summaries[fine_n]], dtype=float),
                np.finfo(float).eps,
            )
            by_summary[column] = {
                "median_absolute_change": float(np.median(absolute)),
                "p95_absolute_change": float(np.quantile(absolute, 0.95)),
                "maximum_absolute_change": float(np.max(absolute)),
                "median_relative_change": float(np.median(relative)),
                "p95_relative_change": float(np.quantile(relative, 0.95)),
                "maximum_relative_change": float(np.max(relative)),
            }
        flips = {}
        for arm in sorted({m["arm"] for m in row_metadata}):
            indices = [i for i, m in enumerate(row_metadata) if m["arm"] == arm]
            changed = sum(
                (summaries[coarse_n][i]["lower_68"] <= row_metadata[i]["truth"] <= summaries[coarse_n][i]["upper_68"])
                != (summaries[fine_n][i]["lower_68"] <= row_metadata[i]["truth"] <= summaries[fine_n][i]["upper_68"])
                for i in indices
            )
            flips[arm] = {"changed_68_coverage_classifications": int(changed), "n": len(indices)}
        median_abs = by_summary["median"]["maximum_absolute_change"]
        median_rel = by_summary["median"]["maximum_relative_change"]
        interval_rel = max(
            by_summary["lower_68"]["maximum_relative_change"],
            by_summary["upper_68"]["maximum_relative_change"],
            by_summary["lower_95"]["maximum_relative_change"],
            by_summary["upper_95"]["maximum_relative_change"],
        )
        comparisons.append({
            "coarse_grid_points": coarse_n,
            "fine_grid_points": fine_n,
            "summary_changes": by_summary,
            "coverage_classification_changes_by_arm": flips,
            "meets_configured_thresholds": bool(
                median_abs <= median_abs_tolerance
                and median_rel <= median_relative_tolerance
                and interval_rel <= interval_relative_tolerance
            ),
            "thresholds": {
                "maximum_absolute_median_change": median_abs_tolerance,
                "maximum_relative_median_change": median_relative_tolerance,
                "maximum_relative_interval_endpoint_change": interval_relative_tolerance,
            },
        })
    return {
        "status": "PASS" if comparisons[-1]["meets_configured_thresholds"] else "REVIEW",
        "interpretation": "The final adjacent-grid comparison meets the configured numerical thresholds; this is not evidence of model adequacy or universal convergence.",
        "scope": "synthetic known-host low-redshift sensitivity diagnostic only",
        "raw_record_count": len(row_metadata),
        "grid_points": list(points),
        "support": [50.0, 90.0],
        "comparisons": comparisons,
        "method": "independent Gaussian-likelihood and trapezoidal-grid reconstruction from raw CSV",
    }

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("replicates_csv", type=Path)
    parser.add_argument("--fine-grid-points", type=int, default=3201)
    parser.add_argument("--max-allowed-change", type=float, default=1e-3)
    parser.add_argument("--resolution-sweep", action="store_true", help="also compare 801/1601/3201/6401-point grids")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = check_grid_stability(
        args.replicates_csv, args.fine_grid_points, args.max_allowed_change
    )
    if args.resolution_sweep:
        report["multi_resolution_convergence"] = compare_grid_resolutions(args.replicates_csv)
    output = args.output or args.replicates_csv.with_name("grid_stability_report.json")
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
