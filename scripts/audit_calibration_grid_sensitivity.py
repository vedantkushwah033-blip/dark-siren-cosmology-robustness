"""Independently audit paired PIT grid-sensitivity comparison tables."""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np

from scripts.audit_selected_likelihood_sbc import audit


RAW_FIELDS = (
    "truth_h0", "host_index", "host_redshift",
    "observed_distance_mpc", "rejection_attempts",
)
FLOAT_FIELDS = (
    "pit_mean", "pit_mean_monte_carlo_se", "pit_variance",
    "ks_statistic", "ks_p_value_descriptive_only",
    "max_abs_decile_frequency_deviation",
    "mean_pit_delta_vs_finest", "mean_pit_delta_monte_carlo_se",
    "mean_pit_delta_normal_95_low", "mean_pit_delta_normal_95_high",
    "rms_pit_delta_vs_finest", "max_abs_pit_delta_vs_finest",
)


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _integer(value: str, field: str) -> int:
    """Parse integer-valued CSV metadata without silently truncating decimals."""
    if not isinstance(value, str) or re.fullmatch(r"-?[0-9]+", value.strip()) is None:
        raise ValueError(f"comparison field {field} must be an integer")
    return int(value.strip())


def _close(observed: str, expected: float, field: str, grid: int) -> None:
    value = float(observed)
    if not math.isfinite(value) or not math.isclose(value, expected, rel_tol=2e-9, abs_tol=2e-11):
        raise ValueError(f"grid {grid}: comparison field {field} does not match independent reconstruction")


def audit_grid_sensitivity(root: str | Path) -> dict:
    root = Path(root)
    comparison_path = root / "grid_sensitivity.csv"
    rows = _read_csv(comparison_path)
    if len(rows) < 2:
        raise ValueError("comparison must contain at least two grid resolutions")
    grids = [_integer(row.get("grid_points"), "grid_points") for row in rows]
    if len(set(grids)) != len(grids) or any(g < 2 for g in grids):
        raise ValueError("grid resolutions must be unique integers >= 2")
    reference_grid = max(grids)
    try:
        replicate_counts = [_integer(row.get("n_replicates"), "n_replicates") for row in rows]
        seeds = [_integer(row.get("seed"), "seed") for row in rows]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("comparison rows must contain integer replicate counts and seeds") from exc
    if len(set(replicate_counts)) != 1:
        raise ValueError("all grid resolutions must use the same replicate count")
    if len(set(seeds)) != 1:
        raise ValueError("all grid resolutions must use the same random seed for paired comparison")
    # Enforce the same predeclared minimum as the experiment runner. An audit
    # must not certify a hand-edited or externally produced tiny sample merely
    # because its arithmetic is internally consistent.
    if replicate_counts[0] < 100 or seeds[0] < 0:
        raise ValueError("replicate count must be >= 100 and seed must be non-negative")
    records: dict[int, list[dict]] = {}
    summaries: dict[int, dict] = {}

    for row, grid in zip(rows, grids):
        if row.get("independent_audit") != "PASS":
            raise ValueError(f"grid {grid}: comparison table does not record a passing independent audit")
        case = root / f"grid_{grid}"
        raw_path, summary_path = case / "pit_replicates.csv", case / "pit_summary.json"
        report = audit(raw_path, summary_path)
        if report.get("status") != "PASS":
            raise ValueError(f"grid {grid}: raw-record audit did not pass")
        raw = _read_csv(raw_path)
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if len(raw) != _integer(row.get("n_replicates"), "n_replicates") or int(summary["n_replicates"]) != len(raw):
            raise ValueError(f"grid {grid}: replicate-count mismatch")
        if int(summary["seed"]) != _integer(row.get("seed"), "seed"):
            raise ValueError(f"grid {grid}: seed mismatch")
        if int(summary["posterior_grid"]["points"]) != grid:
            raise ValueError(f"grid {grid}: grid metadata mismatch")
        records[grid], summaries[grid] = raw, summary
        for field in ("pit_mean", "pit_mean_monte_carlo_se", "pit_variance", "ks_statistic",
                      "ks_p_value_descriptive_only", "max_abs_decile_frequency_deviation"):
            _close(row[field], float(summary[field]), field, grid)

    base = records[grids[0]]
    for grid in grids[1:]:
        if len(records[grid]) != len(base):
            raise ValueError("paired grids have different replicate counts")
        for field in RAW_FIELDS:
            if [r[field] for r in records[grid]] != [r[field] for r in base]:
                raise ValueError(f"paired simulation records differ across grids in {field}")

    ref = np.asarray([float(r["pit"]) for r in records[reference_grid]], dtype=float)
    if not np.all(np.isfinite(ref)) or np.any((ref < 0) | (ref > 1)):
        raise ValueError("reference PIT values must be finite and in [0, 1]")
    for row, grid in zip(rows, grids):
        pit = np.asarray([float(r["pit"]) for r in records[grid]], dtype=float)
        if not np.all(np.isfinite(pit)) or np.any((pit < 0) | (pit > 1)):
            raise ValueError(f"grid {grid}: PIT values must be finite and in [0, 1]")
        delta = pit - ref
        mean = float(np.mean(delta))
        se = float(np.std(delta, ddof=1) / np.sqrt(len(delta))) if len(delta) > 1 else 0.0
        expected = {
            "reference_grid_points": reference_grid,
            "mean_pit_delta_vs_finest": mean,
            "mean_pit_delta_monte_carlo_se": se,
            "mean_pit_delta_normal_95_low": mean - 1.96 * se,
            "mean_pit_delta_normal_95_high": mean + 1.96 * se,
            "rms_pit_delta_vs_finest": float(np.sqrt(np.mean(delta**2))),
            "max_abs_pit_delta_vs_finest": float(np.max(np.abs(delta))),
        }
        for field, value in expected.items():
            if field == "reference_grid_points":
                if _integer(row.get(field), field) != value:
                    raise ValueError(f"grid {grid}: incorrect reference grid")
            else:
                _close(row[field], value, field, grid)

    return {
        "status": "PASS",
        "n_grid_resolutions": len(grids),
        "grid_points": grids,
        "reference_grid_points": reference_grid,
        "n_replicates_per_grid": len(base),
        "paired_raw_records_identical": True,
        "comparison_statistics_reconstructed": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", default="artifacts/calibration_grid_sensitivity")
    args = parser.parse_args()
    result = audit_grid_sensitivity(args.directory)
    out = Path(args.directory) / "grid_sensitivity_audit.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
