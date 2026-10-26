"""Small, paired synthetic sensitivity runner for the three distance-noise models."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dark_siren_robustness.uncertainty_sensitivity import run_experiment


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261026)
    parser.add_argument("--grid-points", type=int, default=1601)
    parser.add_argument("--output-dir", default="artifacts/distance_uncertainty_sensitivity")
    args = parser.parse_args()
    raw, summaries = run_experiment(args.replicates, args.seed, args.grid_points)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    raw_path = output / "replicates.csv"
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(raw[0]))
        writer.writeheader()
        writer.writerows(raw)
    summary_path = output / "summary.json"
    metadata = {
        "status": "synthetic sensitivity diagnostic; not an observational GW analysis",
        "seed": args.seed,
        "replicates_per_cell": args.replicates,
        "grid_points": args.grid_points,
        "prior": [50.0, 90.0],
        "truths": [50, 55, 60, 65, 70, 75, 80, 85, 90],
        "fractional_errors": [0.05, 0.10, 0.15, 0.20],
        "arms": ["A_fixed_realized_sigma", "B_fractional_sigma_H0", "C_fixed_reference_sigma"],
        "raw_record_count": len(raw),
        "summary_cell_count": len(summaries),
        "boundary_truths_must_be_interpreted_separately": True,
        "summary": summaries,
    }
    summary_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {len(raw)} raw records to {raw_path}")
    print(f"Wrote {len(summaries)} cell summaries to {summary_path}")


if __name__ == "__main__":
    main()
