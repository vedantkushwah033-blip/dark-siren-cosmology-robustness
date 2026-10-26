"""Run independently audited paired grid-sensitivity experiments across seeds.

This extends the single-seed screening experiment to expose seed-to-seed
variation. It remains a synthetic numerical diagnostic, not astrophysical
validation. Each seed's raw records and grid comparisons are audited before
its aggregate statistics are included.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(_REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPOSITORY_ROOT))

from scripts.audit_calibration_grid_sensitivity import audit_grid_sensitivity
from scripts.audit_multiseed_grid_sensitivity import audit_multiseed_sensitivity
from scripts.run_calibration_grid_sensitivity import run_sensitivity
from scripts.run_selected_likelihood_sbc import N_REPLICATES, SEED


def run_multiseed_sensitivity(
    grid_points: list[int],
    seeds: list[int],
    n_replicates: int = 300,
    output_dir: str | Path = "artifacts/multiseed_grid_sensitivity",
) -> list[dict]:
    if isinstance(n_replicates, bool) or not isinstance(n_replicates, int) or n_replicates < 100:
        raise ValueError("n_replicates must be an integer of at least 100.")
    if len(seeds) < 2:
        raise ValueError("At least two distinct seeds are required for a multi-seed study.")
    if any(isinstance(seed, bool) or not isinstance(seed, int) or seed < 0 for seed in seeds):
        raise ValueError("Seeds must be non-negative integers.")
    if len(set(seeds)) != len(seeds):
        raise ValueError("Seeds must be unique.")
    if len(grid_points) < 2:
        raise ValueError("Provide at least two grid resolutions for comparison.")
    if any(isinstance(points, bool) or not isinstance(points, int) or points < 2 for points in grid_points):
        raise ValueError("Grid resolutions must be integers of at least two points.")
    if len(set(grid_points)) != len(grid_points):
        raise ValueError("Grid resolutions must be unique.")

    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    aggregate_path = root / "multiseed_grid_sensitivity.csv"
    aggregate_path.unlink(missing_ok=True)
    per_seed_path = root / "per_seed_paired_contrasts.csv"
    per_seed_path.unlink(missing_ok=True)
    # Do not let a failed rerun leave an old metadata PASS marker behind.
    (root / "multiseed_metadata.json").unlink(missing_ok=True)
    (root / "multiseed_audit_report.json").unlink(missing_ok=True)

    per_seed_rows = []
    for seed in seeds:
        seed_dir = root / f"seed_{seed}"
        seed_results = run_sensitivity(
            grid_points=grid_points,
            n_replicates=n_replicates,
            seed=seed,
            output_dir=seed_dir,
        )
        # Auditing each grid's raw records is necessary but not sufficient:
        # independently reconstruct the paired comparison table before using
        # any of its rows in the across-seed aggregate.
        comparison_audit = audit_grid_sensitivity(seed_dir)
        (seed_dir / "grid_sensitivity_audit.json").write_text(
            json.dumps(comparison_audit, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        if comparison_audit.get("status") != "PASS":
            raise RuntimeError(f"Seed {seed}: paired comparison audit did not pass.")
        for row in seed_results:
            if row.get("independent_audit") != "PASS":
                raise RuntimeError(f"Seed {seed}, grid {row.get('grid_points')}: independent audit did not pass.")
            per_seed_rows.append({"seed": seed, **row})

    # Treat seeds, not individual replicates, as the units for uncertainty
    # about seed-to-seed variation. Keep the within-seed paired contrast SE
    # distinct from the across-seed standard error reported below.
    aggregate_rows = []
    for points in grid_points:
        selected = [row for row in per_seed_rows if row["grid_points"] == points]
        if len(selected) != len(seeds):
            raise RuntimeError(f"Missing seed-level results for grid {points}.")
        pit_means = [float(row["pit_mean"]) for row in selected]
        paired_deltas = [float(row["mean_pit_delta_vs_finest"]) for row in selected]
        if not all(math.isfinite(value) for value in pit_means + paired_deltas):
            raise ValueError(f"Non-finite seed-level statistic for grid {points}.")
        pit_sd = _sample_sd(pit_means)
        delta_sd = _sample_sd(paired_deltas)
        aggregate_rows.append({
            "grid_points": points,
            "reference_grid_points": max(grid_points),
            "n_seeds": len(seeds),
            "replicates_per_seed": n_replicates,
            "total_replicates_across_seeds": len(seeds) * n_replicates,
            "seed_list": ";".join(str(seed) for seed in seeds),
            "mean_pit_mean_across_seeds": sum(pit_means) / len(pit_means),
            "sd_pit_mean_across_seeds": pit_sd,
            "se_pit_mean_across_seeds": pit_sd / math.sqrt(len(seeds)),
            "mean_paired_pit_delta_vs_finest_across_seeds": sum(paired_deltas) / len(paired_deltas),
            "sd_paired_pit_delta_across_seeds": delta_sd,
            "se_paired_pit_delta_across_seeds": delta_sd / math.sqrt(len(seeds)),
            "all_seed_grid_audits_passed": "PASS",
        })

    with per_seed_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(per_seed_rows[0]))
        writer.writeheader()
        writer.writerows(per_seed_rows)
    with aggregate_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(aggregate_rows[0]))
        writer.writeheader()
        writer.writerows(aggregate_rows)

    # Reconstruct the across-seed aggregation from saved per-seed records.
    # Only publish the success metadata after this independent audit passes.
    aggregate_audit = audit_multiseed_sensitivity(root)
    (root / "multiseed_audit_report.json").write_text(
        json.dumps(aggregate_audit, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if aggregate_audit.get("status") != "PASS":
        raise RuntimeError("Independent multi-seed aggregate audit did not pass.")
    metadata = {
        "status": "PASS",
        "interpretation": "Synthetic numerical sensitivity only; not astrophysical validation.",
        "grid_points": grid_points,
        "seeds": seeds,
        "replicates_per_seed_per_grid": n_replicates,
        "seed_is_uncertainty_unit_for_aggregate_se": True,
        "per_seed_raw_and_comparison_audits": "PASS",
        "independent_aggregate_audit": "PASS",
        "aggregate_csv": aggregate_path.name,
        "per_seed_csv": per_seed_path.name,
    }
    (root / "multiseed_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return aggregate_rows


def _sample_sd(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean = sum(values) / len(values)
    return math.sqrt(sum((value - mean) ** 2 for value in values) / (len(values) - 1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=300)
    parser.add_argument("--seeds", default=f"{SEED},20261026,20261027")
    parser.add_argument("--grid-points", default="1201,2401,4801")
    parser.add_argument("--output-dir", default="artifacts/multiseed_grid_sensitivity")
    args = parser.parse_args()
    seeds = [int(item.strip()) for item in args.seeds.split(",") if item.strip()]
    points = [int(item.strip()) for item in args.grid_points.split(",") if item.strip()]
    results = run_multiseed_sensitivity(points, seeds, args.replicates, args.output_dir)
    print(json.dumps(results, indent=2))
    print(f"wrote {Path(args.output_dir) / 'multiseed_grid_sensitivity.csv'}")


if __name__ == "__main__":
    main()
