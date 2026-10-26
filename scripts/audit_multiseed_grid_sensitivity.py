"""Independently audit across-seed grid-sensitivity aggregation."""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

# Support both direct execution (python scripts/audit_multiseed_grid_sensitivity.py)
# and package import from tests or other repository scripts.
_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(_REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPOSITORY_ROOT))

from scripts.audit_calibration_grid_sensitivity import audit_grid_sensitivity


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _integer(value: str, field: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must be an integer") from exc
    if isinstance(value, str) and value.strip() != str(parsed):
        raise ValueError(f"{field} must be an integer")
    return parsed


def _sample_sd(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean = sum(values) / len(values)
    return math.sqrt(sum((value - mean) ** 2 for value in values) / (len(values) - 1))


def _close(observed: str, expected: float, field: str, grid: int) -> None:
    try:
        value = float(observed)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"grid {grid}: aggregate field {field} must be numeric") from exc
    if not math.isfinite(value) or not math.isclose(value, expected, rel_tol=2e-9, abs_tol=2e-11):
        raise ValueError(f"grid {grid}: aggregate field {field} does not match independent reconstruction")


def audit_multiseed_sensitivity(root: str | Path) -> dict:
    root = Path(root)
    per_seed_rows = _read_csv(root / "per_seed_paired_contrasts.csv")
    aggregate_rows = _read_csv(root / "multiseed_grid_sensitivity.csv")
    if not per_seed_rows or len(aggregate_rows) < 2:
        raise ValueError("multi-seed audit requires per-seed rows and at least two grid resolutions")

    seed_values = [_integer(row.get("seed"), "seed") for row in per_seed_rows]
    seeds = list(dict.fromkeys(seed_values))
    grids = [_integer(row.get("grid_points"), "grid_points") for row in aggregate_rows]
    if len(seeds) < 2 or any(seed < 0 for seed in seeds):
        raise ValueError("multi-seed audit requires at least two distinct non-negative seeds")
    if len(set(grids)) != len(grids) or any(grid < 2 for grid in grids):
        raise ValueError("aggregate grid resolutions must be unique integers >= 2")
    expected_pairs = {(seed, grid) for seed in seeds for grid in grids}
    observed_pairs = {
        (_integer(row.get("seed"), "seed"), _integer(row.get("grid_points"), "grid_points"))
        for row in per_seed_rows
    }
    if len(per_seed_rows) != len(expected_pairs) or observed_pairs != expected_pairs:
        raise ValueError("per-seed table must contain exactly one row per seed and grid")

    replicate_counts = {_integer(row.get("n_replicates"), "n_replicates") for row in per_seed_rows}
    if len(replicate_counts) != 1 or min(replicate_counts) < 100:
        raise ValueError("per-seed rows must have a common replicate count >= 100")
    n_replicates = next(iter(replicate_counts))

    # When final run metadata exists, verify that it describes these exact
    # artifacts. During the runner's pre-metadata audit this file is absent;
    # the dedicated post-run CI audit checks it after the PASS marker is written.
    metadata_path = root / "multiseed_metadata.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("status") != "PASS":
            raise ValueError("multi-seed metadata status must be PASS")
        if metadata.get("grid_points") != grids:
            raise ValueError("metadata grid_points do not match aggregate rows")
        if metadata.get("seeds") != seeds:
            raise ValueError("metadata seeds do not match per-seed rows")
        if metadata.get("replicates_per_seed_per_grid") != n_replicates:
            raise ValueError("metadata replicate count does not match per-seed rows")
        if metadata.get("per_seed_raw_and_comparison_audits") != "PASS":
            raise ValueError("metadata does not confirm passing per-seed audits")
        if metadata.get("independent_aggregate_audit") != "PASS":
            raise ValueError("metadata does not confirm passing aggregate audit")

    # Re-run each seed's paired-comparison audit, and verify that the
    # aggregation input rows are exactly those independently audited rows.
    for seed in seeds:
        seed_dir = root / f"seed_{seed}"
        report = audit_grid_sensitivity(seed_dir)
        if report.get("status") != "PASS":
            raise ValueError(f"seed {seed}: paired-comparison audit failed")
        saved_report = json.loads((seed_dir / "grid_sensitivity_audit.json").read_text(encoding="utf-8"))
        if saved_report != report:
            raise ValueError(f"seed {seed}: saved paired-comparison audit does not match reconstruction")
        comparison_rows = _read_csv(seed_dir / "grid_sensitivity.csv")
        aggregate_inputs = [
            {key: value for key, value in row.items() if key != "seed"}
            for row in per_seed_rows
            if _integer(row.get("seed"), "seed") == seed
        ]
        if comparison_rows != aggregate_inputs:
            raise ValueError(f"seed {seed}: aggregate input rows differ from audited comparison table")

    finest = max(grids)
    for row, grid in zip(aggregate_rows, grids):
        if _integer(row.get("reference_grid_points"), "reference_grid_points") != finest:
            raise ValueError(f"grid {grid}: incorrect aggregate reference grid")
        if _integer(row.get("n_seeds"), "n_seeds") != len(seeds):
            raise ValueError(f"grid {grid}: incorrect aggregate seed count")
        if _integer(row.get("replicates_per_seed"), "replicates_per_seed") != n_replicates:
            raise ValueError(f"grid {grid}: incorrect aggregate replicate count")
        if _integer(row.get("total_replicates_across_seeds"), "total_replicates_across_seeds") != len(seeds) * n_replicates:
            raise ValueError(f"grid {grid}: incorrect total replicate count")
        if row.get("seed_list") != ";".join(str(seed) for seed in seeds):
            raise ValueError(f"grid {grid}: seed list does not match per-seed records")
        if row.get("all_seed_grid_audits_passed") != "PASS":
            raise ValueError(f"grid {grid}: aggregate does not record passing seed audits")
        selected = [
            item for item in per_seed_rows
            if _integer(item.get("grid_points"), "grid_points") == grid
        ]
        pit_means = [float(item["pit_mean"]) for item in selected]
        deltas = [float(item["mean_pit_delta_vs_finest"]) for item in selected]
        if not all(math.isfinite(value) for value in pit_means + deltas):
            raise ValueError(f"grid {grid}: per-seed statistics must be finite")
        pit_sd = _sample_sd(pit_means)
        delta_sd = _sample_sd(deltas)
        expected = {
            "mean_pit_mean_across_seeds": sum(pit_means) / len(seeds),
            "sd_pit_mean_across_seeds": pit_sd,
            "se_pit_mean_across_seeds": pit_sd / math.sqrt(len(seeds)),
            "mean_paired_pit_delta_vs_finest_across_seeds": sum(deltas) / len(seeds),
            "sd_paired_pit_delta_across_seeds": delta_sd,
            "se_paired_pit_delta_across_seeds": delta_sd / math.sqrt(len(seeds)),
        }
        for field, value in expected.items():
            _close(row.get(field), value, field, grid)

    return {
        "status": "PASS",
        "n_seeds": len(seeds),
        "seeds": seeds,
        "n_grid_resolutions": len(grids),
        "grid_points": grids,
        "replicates_per_seed_per_grid": n_replicates,
        "per_seed_pairs_complete": True,
        "seed_level_audits_reconstructed": True,
        "aggregate_statistics_reconstructed": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", default="artifacts/multiseed_grid_sensitivity")
    args = parser.parse_args()
    result = audit_multiseed_sensitivity(args.directory)
    out = Path(args.directory) / "multiseed_audit_report.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
