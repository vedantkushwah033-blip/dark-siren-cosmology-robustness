"""Independent audit of paired host-density H0-grid sensitivity artifacts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from audit_host_density_benchmark import audit_artifacts


METRICS = ("posterior_mean", "posterior_median", "ci68_low", "ci68_high", "ci68_width")


def audit_grid_sensitivity(folder: Path, tolerance: float = 0.01) -> dict:
    manifest_path = folder / "manifest.json"
    if not manifest_path.is_file():
        raise ValueError("manifest.json is missing")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("artifact_schema_version") != 1 or manifest.get("status") != (
        "synthetic_exploratory_host_density_h0_grid_sensitivity_not_real_H0_measurement"
    ):
        raise ValueError("unsupported schema or scientific-scope status")
    points = manifest.get("grid_points")
    if not isinstance(points, list) or len(points) < 2 or points != sorted(set(points)):
        raise ValueError("invalid grid_points manifest")
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive")
    expected_artifacts = {item["grid_points"]: item for item in manifest.get("artifacts", [])}
    if set(expected_artifacts) != set(points):
        raise ValueError("manifest artifact list does not match grid_points")

    by_grid = {}
    storage = {}
    for grid_points in points:
        item = expected_artifacts[grid_points]
        raw_path = folder / item["raw_csv"]
        summary_path = folder / item["summary_json"]
        raw_bytes = raw_path.read_bytes()
        if hashlib.sha256(raw_bytes).hexdigest() != item["raw_csv_sha256"]:
            raise ValueError(f"raw CSV checksum mismatch for grid {grid_points}")
        if len(raw_bytes) != item["raw_csv_bytes"]:
            raise ValueError(f"raw CSV byte count mismatch for grid {grid_points}")
        report = audit_artifacts(raw_path.parent)
        if report.get("status") != "PASS":
            raise ValueError(f"posterior artifact audit failed for grid {grid_points}")
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if summary.get("h0_grid", [None, None, None])[2] != grid_points:
            raise ValueError(f"summary H0 grid mismatch for grid {grid_points}")
        with raw_path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != item["record_count"]:
            raise ValueError(f"record count mismatch for grid {grid_points}")
        records = {}
        for row in rows:
            key = (int(row["replicate"]), row["arm"])
            if key in records:
                raise ValueError(f"duplicate replicate/arm row at grid {grid_points}: {key}")
            records[key] = row
        by_grid[grid_points] = records
        storage[grid_points] = len(raw_bytes)

    finest = points[-1]
    reference_keys = set(by_grid[finest])
    if len(reference_keys) != manifest["n_replicates"] * 3:
        raise ValueError("finest-grid record count does not match expected paired arms")
    per_grid = []
    for grid_points in points:
        records = by_grid[grid_points]
        if set(records) != reference_keys:
            raise ValueError(f"replicate/arm keys differ at grid {grid_points}")
        changes = {metric: [] for metric in METRICS}
        for key, row in records.items():
            ref = by_grid[finest][key]
            if int(row["seed"]) != int(ref["seed"]):
                raise ValueError(f"paired seed mismatch for {key}")
            for metric in METRICS:
                delta = abs(float(row[metric]) - float(ref[metric]))
                if not math.isfinite(delta):
                    raise ValueError(f"non-finite grid difference for {key}, {metric}")
                changes[metric].append(delta)
        maxima = {metric: max(values) for metric, values in changes.items()}
        per_grid.append({
            "grid_points": grid_points,
            "raw_csv_bytes": storage[grid_points],
            "max_abs_difference_vs_finest": maxima,
            "within_tolerance": all(value <= tolerance for value in maxima.values()),
        })

    # A convergence claim requires the two finest distinct resolutions to agree,
    # not merely a coarse grid to agree with itself.
    finest_pair = per_grid[-2:]
    converged = all(item["within_tolerance"] for item in finest_pair)
    coarsest_acceptable = next(
        (item["grid_points"] for item in per_grid if item["within_tolerance"]), None
    )
    comparison_path = folder / manifest["comparison_csv"]
    comparison_bytes = comparison_path.read_bytes()
    if hashlib.sha256(comparison_bytes).hexdigest() != manifest["comparison_csv_sha256"]:
        raise ValueError("comparison CSV checksum mismatch")
    with comparison_path.open(newline="", encoding="utf-8") as handle:
        comparison_rows = list(csv.DictReader(handle))
    expected_comparison_count = len(points) * len(reference_keys)
    if len(comparison_rows) != expected_comparison_count:
        raise ValueError("comparison CSV row count mismatch")
    comparison_seen = set()
    for row in comparison_rows:
        grid = int(row["grid_points"])
        key = (int(row["replicate"]), row["arm"])
        combined_key = (grid, *key)
        if combined_key in comparison_seen or grid not in by_grid or key not in reference_keys:
            raise ValueError(f"duplicate or unexpected comparison key: {combined_key}")
        comparison_seen.add(combined_key)
        raw = by_grid[grid][key]
        ref = by_grid[finest][key]
        if int(row["seed"]) != int(raw["seed"]):
            raise ValueError(f"comparison seed mismatch: {combined_key}")
        for metric in METRICS:
            if not math.isclose(float(row[metric]), float(raw[metric]), rel_tol=1e-10, abs_tol=1e-10):
                raise ValueError(f"comparison metric mismatch for {combined_key}, {metric}")
            expected_delta = abs(float(raw[metric]) - float(ref[metric]))
            recorded_delta = float(row[f"abs_diff_vs_{finest}_{metric}"])
            if not math.isclose(recorded_delta, expected_delta, rel_tol=1e-8, abs_tol=1e-10):
                raise ValueError(f"comparison difference mismatch for {combined_key}, {metric}")
    return {
        "status": "PASS" if converged else "UNDERRESOLVED",
        "scope": "synthetic host-density H0-grid numerical sensitivity only",
        "tolerance_h0": tolerance,
        "finest_grid_points": finest,
        "finest_pair_converged": converged,
        "coarsest_grid_with_all_metrics_within_tolerance": coarsest_acceptable,
        "replicate_arm_pairs": len(reference_keys),
        "per_grid": per_grid,
        "comparison_rows": len(comparison_rows),
        "comparison_csv_bytes": len(comparison_bytes),
        "scientific_limit": "Numerical grid stability does not establish posterior calibration, model correctness, or real-data validity.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--tolerance", type=float, default=0.01)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = audit_grid_sensitivity(args.artifact_dir, args.tolerance)
    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    if report["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
