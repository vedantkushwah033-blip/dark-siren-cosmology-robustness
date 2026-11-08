"""Paired H0-grid convergence experiment for the synthetic host-density benchmark."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from dark_siren_robustness.host_density_benchmark import run_benchmark


METRICS = ("posterior_mean", "posterior_median", "ci68_low", "ci68_high", "ci68_width")


def run_grid_sensitivity(*, output_dir: Path, grid_points: list[int], n_replicates: int = 10,
                         n_parent: int = 4000, n_detected: int = 30,
                         seed: int = 20261108, z_points: int = 401) -> dict:
    if len(grid_points) < 2 or grid_points != sorted(set(grid_points)):
        raise ValueError("grid_points must contain at least two unique increasing values")
    if any(isinstance(p, bool) or p < 31 for p in grid_points):
        raise ValueError("every H0 grid must contain at least 31 points")
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    artifacts = []
    for points in grid_points:
        folder = output_dir / f"grid_{points}"
        folder.mkdir(parents=True, exist_ok=True)
        rows, summary = run_benchmark(
            n_replicates=n_replicates, n_parent=n_parent, n_detected=n_detected,
            seed=seed, h0_points=points, z_points=z_points,
        )
        raw_path = folder / "replicate_metrics.csv"
        summary_path = folder / "summary.json"
        rows.to_csv(raw_path, index=False)
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        raw_bytes = raw_path.read_bytes()
        artifacts.append({
            "grid_points": points,
            "raw_csv": str(raw_path.relative_to(output_dir)),
            "summary_json": str(summary_path.relative_to(output_dir)),
            "raw_csv_bytes": len(raw_bytes),
            "raw_csv_sha256": hashlib.sha256(raw_bytes).hexdigest(),
            "record_count": len(rows),
        })
        for row in rows.to_dict(orient="records"):
            records.append({key: row[key] for key in ("replicate", "seed", "arm", *METRICS)})
            records[-1]["grid_points"] = points

    long = pd.DataFrame(records)
    finest = max(grid_points)
    reference = long[long.grid_points == finest].set_index(["replicate", "arm"])
    comparison = []
    for row in long.to_dict(orient="records"):
        ref = reference.loc[(row["replicate"], row["arm"])]
        item = {key: row[key] for key in ("grid_points", "replicate", "seed", "arm", *METRICS)}
        for metric in METRICS:
            item[f"abs_diff_vs_{finest}_{metric}"] = abs(float(row[metric]) - float(ref[metric]))
        comparison.append(item)
    comparison_path = output_dir / "grid_comparison.csv"
    pd.DataFrame(comparison).to_csv(comparison_path, index=False)
    manifest = {
        "artifact_schema_version": 1,
        "status": "synthetic_exploratory_host_density_h0_grid_sensitivity_not_real_H0_measurement",
        "method": "paired simulation seeds and identical parent/event draws across H0 grid resolutions",
        "grid_points": grid_points,
        "finest_reference_grid_points": finest,
        "n_replicates": n_replicates,
        "n_parent": n_parent,
        "n_detected": n_detected,
        "seed_start": seed,
        "z_points": z_points,
        "comparison_csv": "grid_comparison.csv",
        "comparison_csv_sha256": hashlib.sha256(comparison_path.read_bytes()).hexdigest(),
        "artifacts": artifacts,
        "limitations": [
            "Synthetic low-redshift toy benchmark only; no LVK data or measured H0.",
            "Grid convergence is numerical stability, not posterior calibration or physical-model validation.",
            "The oracle arm is a diagnostic reference and not an operational estimator.",
            "This experiment varies H0 grid resolution while holding the redshift grid and simulation seeds fixed.",
        ],
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-points", default="161,321,641,1281,3201")
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--n-parent", type=int, default=4000)
    parser.add_argument("--n-detected", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20261108)
    parser.add_argument("--z-points", type=int, default=401)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/host_density_grid_sensitivity"))
    args = parser.parse_args()
    points = [int(value) for value in args.grid_points.split(",")]
    print(json.dumps(run_grid_sensitivity(
        output_dir=args.output_dir, grid_points=points, n_replicates=args.replicates,
        n_parent=args.n_parent, n_detected=args.n_detected, seed=args.seed,
        z_points=args.z_points,
    ), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
