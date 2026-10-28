"""Run a reproducible, controlled audit of redshift/distance-dependent catalogue thinning.

This is a software/statistical calibration experiment, not a real-survey selection model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dark_siren_robustness.completeness import (
    distance_dependent_completeness,
    thin_catalogue,
)

DISTANCE_MAX_MPC = 1200.0
SCALE_MPC = 450.0
FLOOR = 0.05
BIN_EDGES_MPC = [0.0, 300.0, 600.0, 900.0, 1200.0000001]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(n_parent: int, seed: int, output_dir: Path) -> None:
    if n_parent < 1000:
        raise ValueError("n_parent must be at least 1000 for a useful calibration audit")
    output_dir.mkdir(parents=True, exist_ok=True)

    # A fixed, deliberately simple parent catalogue isolates thinning from
    # population-sampling variability. This is conditional on these distances.
    distance = np.linspace(0.0, DISTANCE_MAX_MPC, n_parent, endpoint=False)
    probability = distance_dependent_completeness(distance, SCALE_MPC, FLOOR)
    retained = thin_catalogue(probability, np.random.default_rng(seed))
    records = pd.DataFrame(
        {
            "object_id": np.arange(n_parent, dtype=np.int64),
            "distance_mpc": distance,
            "inclusion_probability": probability,
            "retained": retained.astype(np.int8),
        }
    )
    records_path = output_dir / "catalogue_thinning_records.csv"
    records.to_csv(records_path, index=False, float_format="%.12g")

    records["distance_bin"] = pd.cut(
        records["distance_mpc"], BIN_EDGES_MPC, right=False, include_lowest=True
    )
    rows = []
    for label, group in records.groupby("distance_bin", observed=True, sort=True):
        p = group["inclusion_probability"].to_numpy()
        observed = int(group["retained"].sum())
        expected = float(p.sum())
        variance = float(np.sum(p * (1.0 - p)))
        rows.append(
            {
                "distance_bin": str(label),
                "n_parent": int(len(group)),
                "expected_retained": expected,
                "observed_retained": observed,
                "binomial_sd": float(np.sqrt(variance)),
                "standardized_residual": float((observed - expected) / np.sqrt(variance)),
                "mean_inclusion_probability": float(p.mean()),
                "observed_retention_fraction": float(observed / len(group)),
            }
        )

    summary = {
        "experiment": "controlled_distance_dependent_catalogue_thinning",
        "status": "synthetic_calibration_only_not_a_real_survey_model",
        "seed": int(seed),
        "n_parent": int(n_parent),
        "distance_domain_mpc": [0.0, DISTANCE_MAX_MPC],
        "completeness_model": {
            "formula": "p(d) = floor + (1 - floor) * exp(-d / scale_mpc)",
            "scale_mpc": SCALE_MPC,
            "floor": FLOOR,
        },
        "n_retained": int(retained.sum()),
        "overall_expected_retained": float(probability.sum()),
        "overall_observed_retention_fraction": float(retained.mean()),
        "strata": rows,
        "records_sha256": sha256_file(records_path),
        "interpretation": (
            "Conditional Bernoulli-thinning calibration on a fixed parent distance grid. "
            "It does not validate a real survey selection function, a galaxy redshift "
            "intensity model, or any H0 inference."
        ),
    }
    summary_path = output_dir / "catalogue_thinning_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"records": str(records_path), "summary": str(summary_path), "n_retained": int(retained.sum())}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-parent", type=int, default=40000)
    parser.add_argument("--seed", type=int, default=20261028)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/catalogue_thinning"))
    args = parser.parse_args()
    run(args.n_parent, args.seed, args.output_dir)


if __name__ == "__main__":
    main()
