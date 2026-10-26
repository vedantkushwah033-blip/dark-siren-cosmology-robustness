"""Compare selected-likelihood PIT summaries across posterior-grid resolutions.

This is a synthetic numerical-sensitivity diagnostic. Every grid's raw records
and summary are independently audited before the comparison CSV is written.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

# Direct execution sets sys.path to scripts/, whereas pytest/import execution
# normally starts at the repository root. Normalize both entry points.
_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(_REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPOSITORY_ROOT))

from scripts.audit_selected_likelihood_sbc import audit
from scripts.run_selected_likelihood_sbc import N_REPLICATES, SEED, run_calibration


def run_sensitivity(
    grid_points: list[int],
    n_replicates: int = N_REPLICATES,
    seed: int = SEED,
    output_dir: str | Path = "artifacts/calibration_grid_sensitivity",
) -> list[dict]:
    if len(grid_points) < 2:
        raise ValueError("Provide at least two grid resolutions for comparison.")
    if len(set(grid_points)) != len(grid_points):
        raise ValueError("Grid resolutions must be unique.")
    if any(points < 2 for points in grid_points):
        raise ValueError("Each grid resolution must contain at least two points.")

    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    results = []
    for points in grid_points:
        case_dir = root / f"grid_{points}"
        case_dir.mkdir(parents=True, exist_ok=True)
        rows, summary = run_calibration(
            n_replicates=n_replicates,
            seed=seed,
            grid_points=points,
        )
        raw_path = case_dir / "pit_replicates.csv"
        with raw_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        summary_path = case_dir / "pit_summary.json"
        summary_path.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        audit_report = audit(raw_path, summary_path)
        (case_dir / "audit_report.json").write_text(
            json.dumps(audit_report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        results.append({
            "grid_points": points,
            "n_replicates": n_replicates,
            "seed": seed,
            "pit_mean": summary["pit_mean"],
            "pit_mean_monte_carlo_se": summary["pit_mean_monte_carlo_se"],
            "pit_variance": summary["pit_variance"],
            "ks_statistic": summary["ks_statistic"],
            "ks_p_value_descriptive_only": summary["ks_p_value_descriptive_only"],
            "max_abs_decile_frequency_deviation": summary["max_abs_decile_frequency_deviation"],
            "independent_audit": audit_report["status"],
        })

    comparison_path = root / "grid_sensitivity.csv"
    with comparison_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=N_REPLICATES)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--grid-points", default="1201,2401,4801")
    parser.add_argument("--output-dir", default="artifacts/calibration_grid_sensitivity")
    args = parser.parse_args()
    points = [int(item.strip()) for item in args.grid_points.split(",") if item.strip()]
    results = run_sensitivity(points, args.replicates, args.seed, args.output_dir)
    print(json.dumps(results, indent=2))
    print(f"wrote {Path(args.output_dir) / 'grid_sensitivity.csv'}")


if __name__ == "__main__":
    main()
