"""Independently reconstruct catalogue-thinning calibration summaries from raw records."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

BIN_EDGES_MPC = [0.0, 300.0, 600.0, 900.0, 1200.0000001]
EXPECTED_DOMAIN = [0.0, 1200.0]
EXPECTED_SCALE_MPC = 450.0
EXPECTED_FLOOR = 0.05


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit(output_dir: Path) -> dict:
    records_path = output_dir / "catalogue_thinning_records.csv"
    summary_path = output_dir / "catalogue_thinning_summary.json"
    data = pd.read_csv(records_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))

    if summary.get("experiment") != "controlled_distance_dependent_catalogue_thinning":
        raise AssertionError("Unexpected experiment identifier")
    if summary.get("status") != "synthetic_calibration_only_not_a_real_survey_model":
        raise AssertionError("Summary must retain the synthetic-only status")
    if summary.get("distance_domain_mpc") != EXPECTED_DOMAIN:
        raise AssertionError("Declared distance domain differs from the audited design")
    model = summary.get("completeness_model", {})
    if model.get("formula") != "p(d) = floor + (1 - floor) * exp(-d / scale_mpc)":
        raise AssertionError("Unexpected completeness formula declaration")
    if model.get("scale_mpc") != EXPECTED_SCALE_MPC or model.get("floor") != EXPECTED_FLOOR:
        raise AssertionError("Declared completeness parameters differ from audited constants")
    if not isinstance(summary.get("seed"), int) or summary["seed"] < 0:
        raise AssertionError("Seed must be a non-negative integer")

    required = {"object_id", "distance_mpc", "inclusion_probability", "retained"}
    if set(data.columns) != required:
        raise AssertionError(f"Unexpected raw-record columns: {list(data.columns)}")
    if len(data) != summary["n_parent"] or len(data) < 1000:
        raise AssertionError("Raw-record count does not match summary or is too small")
    if data["object_id"].duplicated().any() or not np.array_equal(
        data["object_id"].to_numpy(), np.arange(len(data))
    ):
        raise AssertionError("object_id must be unique and consecutive from zero")
    distance = data["distance_mpc"].to_numpy(dtype=float)
    p = data["inclusion_probability"].to_numpy(dtype=float)
    retained = data["retained"].to_numpy(dtype=int)
    if not np.all(np.isfinite(distance)) or not np.all(np.isfinite(p)):
        raise AssertionError("Non-finite distance or probability in raw records")
    if np.any(distance < EXPECTED_DOMAIN[0]) or np.any(distance >= EXPECTED_DOMAIN[1]):
        raise AssertionError("Distances fall outside the declared half-open domain")
    if np.any(np.diff(distance) <= 0):
        raise AssertionError("Parent-catalogue distances must be strictly increasing")
    if np.any((p < 0) | (p > 1)) or not np.all(np.isin(retained, [0, 1])):
        raise AssertionError("Invalid inclusion probabilities or retained flags")
    # Recompute the declared model without importing the production function.
    expected_p = EXPECTED_FLOOR + (1.0 - EXPECTED_FLOOR) * np.exp(
        -distance / EXPECTED_SCALE_MPC
    )
    if not np.allclose(p, expected_p, rtol=1e-10, atol=1e-12):
        raise AssertionError("Raw inclusion probabilities do not match the declared formula")
    if sha256_file(records_path) != summary["records_sha256"]:
        raise AssertionError("Raw-record SHA-256 does not match the summary")

    if int(retained.sum()) != summary["n_retained"]:
        raise AssertionError("Total retained count does not match summary")
    if not np.isclose(p.sum(), summary["overall_expected_retained"], rtol=1e-10):
        raise AssertionError("Overall expected retained count does not match raw records")
    if not np.isclose(
        retained.mean(), summary["overall_observed_retention_fraction"], rtol=1e-12
    ):
        raise AssertionError("Overall observed retention fraction does not match raw records")

    data = data.assign(
        distance_bin=pd.cut(
            distance, BIN_EDGES_MPC, right=False, include_lowest=True
        )
    )
    rebuilt = []
    for label, group in data.groupby("distance_bin", observed=True, sort=True):
        gp = group["inclusion_probability"].to_numpy(dtype=float)
        gr = group["retained"].to_numpy(dtype=int)
        expected = float(gp.sum())
        variance = float(np.sum(gp * (1.0 - gp)))
        observed = int(gr.sum())
        rebuilt.append(
            {
                "distance_bin": str(label),
                "n_parent": int(len(group)),
                "expected_retained": expected,
                "observed_retained": observed,
                "standardized_residual": float((observed - expected) / np.sqrt(variance)),
                "mean_inclusion_probability": float(gp.mean()),
                "observed_retention_fraction": float(observed / len(group)),
            }
        )
        # A six-sigma guard is a gross-error detector, not a hypothesis test
        # or proof that the selection model is astrophysically correct.
        if abs(rebuilt[-1]["standardized_residual"]) > 6.0:
            raise AssertionError(f"Stratum {label} deviates by more than six binomial SDs")
    stored = summary["strata"]
    if len(rebuilt) != len(stored):
        raise AssertionError("Stratum count differs between raw records and summary")
    for actual, reported in zip(rebuilt, stored, strict=True):
        for key, value in actual.items():
            if isinstance(value, float):
                if not np.isclose(value, reported[key], rtol=1e-9, atol=1e-9):
                    raise AssertionError(f"Summary mismatch for {actual['distance_bin']}:{key}")
            elif value != reported[key]:
                raise AssertionError(f"Summary mismatch for {actual['distance_bin']}:{key}")
    report = {
        "status": "PASS",
        "audit_scope": "raw-record integrity and independent arithmetic reconstruction only",
        "n_parent": int(len(data)),
        "n_retained": int(retained.sum()),
        "records_sha256": sha256_file(records_path),
        "strata_reconstructed": len(rebuilt),
        "scientific_limit": (
            "Passing this audit does not validate the toy selection law as a real survey model "
            "and does not validate any cosmological inference."
        ),
    }
    (output_dir / "catalogue_thinning_audit.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path, nargs="?", default=Path("artifacts/catalogue_thinning"))
    args = parser.parse_args()
    audit(args.output_dir)


if __name__ == "__main__":
    main()
