"""Small, paired synthetic sensitivity runner for the three distance-noise models."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import platform
import sys
from importlib.metadata import version

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from dark_siren_robustness.uncertainty_sensitivity import paired_ab_contrasts, run_experiment


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261026)
    parser.add_argument("--grid-points", type=int, default=1601)
    parser.add_argument("--bootstrap-resamples", type=int, default=2000)
    parser.add_argument("--output-dir", default="artifacts/distance_uncertainty_sensitivity")
    args = parser.parse_args()
    raw, summaries = run_experiment(args.replicates, args.seed, args.grid_points)
    paired_contrasts = paired_ab_contrasts(raw, seed=args.seed, bootstrap_resamples=args.bootstrap_resamples)
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
        "code_revision": os.environ.get("GITHUB_SHA", "not_recorded_outside_ci"),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": version("scipy"),
        "failed_run_count": 0,
        "likelihood_equations": {
            "A_fixed_realized_sigma": "logL=-0.5*((d_obs-c*z/H0)/sigma_realized)^2-log(sigma_realized)-0.5*log(2*pi)",
            "B_fractional_sigma_H0": "sigma(H0)=f*c*z/H0; logL=-0.5*((d_obs-c*z/H0)/sigma(H0))^2-log(sigma(H0))-0.5*log(2*pi)",
            "C_fixed_reference_sigma": "sigma_ref=f*c*z/70; logL=-0.5*((d_obs-c*z/H0)/sigma_ref)^2-log(sigma_ref)-0.5*log(2*pi)",
        },
        "seed": args.seed,
        "replicates_per_cell": args.replicates,
        "grid_points": args.grid_points,
        "bootstrap_resamples": args.bootstrap_resamples,
        "prior": [50.0, 90.0],
        "truths": [50, 55, 60, 65, 70, 75, 80, 85, 90],
        "fractional_errors": [0.05, 0.10, 0.15, 0.20],
        "arms": ["A_fixed_realized_sigma", "B_fractional_sigma_H0", "C_fixed_reference_sigma"],
        "raw_record_count": len(raw),
        "summary_cell_count": len(summaries),
        "boundary_truths_must_be_interpreted_separately": True,
        "summary": summaries,
        "paired_ab_contrasts": paired_contrasts,
        "paired_ab_contrast_count": len(paired_contrasts),
    }
    summary_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {len(raw)} raw records to {raw_path}")
    print(f"Wrote {len(summaries)} cell summaries to {summary_path}")


if __name__ == "__main__":
    main()
