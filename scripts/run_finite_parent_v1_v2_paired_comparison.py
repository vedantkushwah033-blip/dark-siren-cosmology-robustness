"""Run a paired V1/V2 finite-parent benchmark on identical simulated inputs.

The two arms are generated from the same seed and must have identical parent,
catalogue, GW-distance, and validation-sample provenance. This is a solver
comparison only; it does not compare model families or establish superiority.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


IDENTITY_COLUMNS = (
    "scenario", "draw", "validation_size", "validation_seed",
    "simulation_seed", "thinning_seed", "true_intercept", "true_slope",
    "true_quadratic", "h0_true", "sigma_mpc", "detection_threshold_mpc",
    "validation_indices", "catalogue_validation_overlap_indices",
    "catalogue_validation_overlap_count", "n_validation_included",
    "n_catalogue", "n_detected", "parent_size", "parent_frame_sha256", "catalogue_sha256",
    "gw_distance_sha256",
)


def _run_arm(args: argparse.Namespace, version: str, output: Path) -> None:
    command = [
        sys.executable, "scripts/run_finite_parent_validation_benchmark.py",
        "--scenarios", str(args.scenarios),
        "--draws", str(args.draws),
        "--n-parent", str(args.n_parent),
        "--n-detected", str(args.n_detected),
        "--validation-sizes", args.validation_sizes,
        "--seed", str(args.seed),
        "--solver-version", version,
        "--output-dir", str(output),
    ]
    subprocess.run(command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Paired V1/V2 finite-parent solver comparison on identical samples."
    )
    parser.add_argument("--scenarios", type=int, default=2)
    parser.add_argument("--draws", type=int, default=3)
    parser.add_argument("--n-parent", type=int, default=1200)
    parser.add_argument("--n-detected", type=int, default=8)
    parser.add_argument("--validation-sizes", default="100,300")
    parser.add_argument("--seed", type=int, default=20261107)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    out = args.output_dir
    v1_dir, v2_dir = out / "v1", out / "v2"
    v1_dir.mkdir(parents=True, exist_ok=True)
    v2_dir.mkdir(parents=True, exist_ok=True)
    _run_arm(args, "v1", v1_dir)
    _run_arm(args, "v2", v2_dir)

    v1 = pd.read_csv(v1_dir / "records.csv")
    v2 = pd.read_csv(v2_dir / "records.csv")
    manifests = {
        "v1": json.loads((v1_dir / "summary.json").read_text(encoding="utf-8")),
        "v2": json.loads((v2_dir / "summary.json").read_text(encoding="utf-8")),
    }
    common_manifest_fields = (
        "schema", "status", "n_attempted", "validation_sizes", "scenarios",
        "draws_per_scenario_and_size", "seed", "h0_points", "z_points",
        "h0_grid_bounds", "redshift_grid_bounds", "model_family",
    )
    for field in common_manifest_fields:
        if manifests["v1"].get(field) != manifests["v2"].get(field):
            raise AssertionError(f"paired arm manifests disagree on {field!r}")
    if manifests["v1"].get("schema") != "finite_parent_srswor_smoke_v3":
        raise AssertionError("paired arm manifest has an unsupported schema")
    if manifests["v1"].get("status") != "SMOKE_ONLY_NOT_CONFIRMATORY":
        raise AssertionError("paired arms must remain explicitly non-confirmatory")
    requested_sizes = [int(item) for item in args.validation_sizes.split(",") if item.strip()]
    if manifests["v1"].get("seed") != args.seed:
        raise AssertionError("paired arm manifest seed does not match requested seed")
    if manifests["v1"].get("scenarios") != args.scenarios:
        raise AssertionError("paired arm manifest scenario count does not match request")
    if manifests["v1"].get("draws_per_scenario_and_size") != args.draws:
        raise AssertionError("paired arm manifest draw count does not match request")
    if manifests["v1"].get("validation_sizes") != requested_sizes:
        raise AssertionError("paired arm manifest validation sizes do not match request")
    expected_attempts = args.scenarios * args.draws * len(requested_sizes)
    for label, frame, expected in (("V1", v1, "v1"), ("V2", v2, "v2")):
        manifest = manifests[expected]
        if manifest.get("solver_version") != expected:
            raise AssertionError(f"{label} manifest has an unexpected solver version")
        if manifest.get("solver_revision") != (
            "2.1-armijo-noise-slack" if expected == "v2" else "v1-lbfgsb-kkt"
        ):
            raise AssertionError(f"{label} manifest has an unexpected solver revision")
        if len(frame) != expected_attempts or len(frame) != manifest.get("n_attempted"):
            raise AssertionError(f"{label} artifact attempt count disagrees with requested design")
        n_success = int(frame["status"].eq("success").sum())
        n_failed = int(frame["status"].eq("failed").sum())
        if n_success != manifest.get("n_successful") or n_failed != manifest.get("n_failed"):
            raise AssertionError(f"{label} manifest success/failure counts disagree with records")
        if n_success + n_failed != expected_attempts:
            raise AssertionError(f"{label} attempt accounting does not reconcile")
        if not frame["status"].isin(["success", "failed"]).all():
            raise AssertionError(f"{label} artifact contains an unknown attempt status")
        if not frame["solver_version"].eq(expected).all():
            raise AssertionError(f"{label} artifact contains an unexpected solver version")
    if len(v1) != len(v2):
        raise AssertionError("paired V1/V2 arms have different attempt counts")
    for column in IDENTITY_COLUMNS:
        left = v1[column].fillna("<NA>").astype(str).to_numpy()
        right = v2[column].fillna("<NA>").astype(str).to_numpy()
        if not np.array_equal(left, right):
            mismatch = int(np.flatnonzero(left != right)[0])
            raise AssertionError(
                f"paired arms differ in fixed input/provenance column {column!r} "
                f"at record {mismatch}"
            )

    keys = ["scenario", "draw", "validation_size"]
    if v1.duplicated(keys).any() or v2.duplicated(keys).any():
        raise AssertionError("paired arm contains duplicate scenario/draw/size keys")
    pair = v1[keys + ["status", "posterior_mean", "density_ise", "covered_68"]].merge(
        v2[keys + ["status", "posterior_mean", "density_ise", "covered_68"]],
        on=keys, how="outer", suffixes=("_v1", "_v2"), validate="one_to_one",
        indicator=True,
    )
    if not pair["_merge"].eq("both").all():
        raise AssertionError("V1/V2 arms do not contain identical attempt keys")

    both = pair["status_v1"].eq("success") & pair["status_v2"].eq("success")
    v1_only = pair["status_v1"].eq("success") & pair["status_v2"].eq("failed")
    v2_only = pair["status_v1"].eq("failed") & pair["status_v2"].eq("success")
    both_failed = pair["status_v1"].eq("failed") & pair["status_v2"].eq("failed")
    summary = {
        "schema": "finite_parent_v1_v2_paired_comparison_v1",
        "status": "SMOKE_ONLY_NOT_CONFIRMATORY",
        "interpretation": (
            "Paired numerical solver comparison on identical synthetic frames and "
            "SRSWOR validation samples; no estimator-superiority claim."
        ),
        "seed": int(args.seed),
        "n_attempted": int(len(pair)),
        "n_both_success": int(both.sum()),
        "n_v1_only_success": int(v1_only.sum()),
        "n_v2_only_success": int(v2_only.sum()),
        "n_both_failed": int(both_failed.sum()),
        "posterior_mean_delta_v2_minus_v1_mean": (
            float((pair.loc[both, "posterior_mean_v2"] -
                   pair.loc[both, "posterior_mean_v1"]).mean()) if both.any() else None
        ),
        "posterior_mean_delta_v2_minus_v1_max_abs": (
            float((pair.loc[both, "posterior_mean_v2"] -
                   pair.loc[both, "posterior_mean_v1"]).abs().max()) if both.any() else None
        ),
        "density_ise_delta_v2_minus_v1_mean": (
            float((pair.loc[both, "density_ise_v2"] -
                   pair.loc[both, "density_ise_v1"]).mean()) if both.any() else None
        ),
        "coverage_disagreements_among_both_success": (
            int((pair.loc[both, "covered_68_v1"] !=
                 pair.loc[both, "covered_68_v2"]).sum()) if both.any() else 0
        ),
        "identity_columns_checked": list(IDENTITY_COLUMNS),
        "v1_artifact": "v1/",
        "v2_artifact": "v2/",
        "limitations": [
            "Both arms use the same linear logistic inclusion model; only the optimizer differs.",
            "The smoke sample is not a confirmatory calibration or coverage study.",
            "Paired successes do not establish that either solver is statistically superior.",
            "Each arm must also pass its own independent artifact audit.",
        ],
    }
    pair.drop(columns=["_merge"]).to_csv(out / "paired_records.csv", index=False)
    (out / "paired_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
