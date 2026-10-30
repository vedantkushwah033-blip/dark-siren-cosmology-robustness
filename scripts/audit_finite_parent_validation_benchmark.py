"""Independent structural audit for finite-parent SRSWOR benchmark outputs.

This audit verifies design bookkeeping and replay of validation samples. It does
not independently reconstruct the statistical estimator; that is a separate
validation requirement before confirmatory claims.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit

from dark_siren_robustness.generative_selection import generate_synthetic_population


def audit(artifact_dir: Path) -> dict:
    records = pd.read_csv(artifact_dir / "records.csv")
    summary = json.loads((artifact_dir / "summary.json").read_text(encoding="utf-8"))
    required = {
        "scenario", "draw", "validation_size", "status", "validation_seed",
        "simulation_seed", "thinning_seed", "true_intercept", "true_slope",
        "true_quadratic", "sigma_mpc", "detection_threshold_mpc",
        "validation_indices", "parent_size", "parent_frame_sha256",
        "catalogue_sha256", "gw_distance_sha256",
    }
    missing = required.difference(records.columns)
    if missing:
        raise AssertionError(f"records missing columns: {sorted(missing)}")
    if summary.get("schema") != "finite_parent_srswor_smoke_v1":
        raise AssertionError("unsupported finite-parent artifact schema")
    if summary.get("status") != "SMOKE_ONLY_NOT_CONFIRMATORY":
        raise AssertionError("artifact must retain non-confirmatory status")
    attempted = len(records)
    successful = int(records["status"].eq("success").sum())
    failed = int(records["status"].eq("failed").sum())
    if attempted != int(summary["n_attempted"]) or successful != int(summary["n_successful"]):
        raise AssertionError("summary attempt/success counts disagree with records")
    if failed != int(summary["n_failed"]) or attempted != successful + failed:
        raise AssertionError("failed-attempt accounting does not reconcile")
    if not records["status"].isin(["success", "failed"]).all():
        raise AssertionError("unknown attempt status")

    sizes = [int(v) for v in summary["validation_sizes"]]
    n_scenarios = int(summary["scenarios"])
    draws = int(summary["draws_per_scenario_and_size"])
    expected = n_scenarios * draws * len(sizes)
    if attempted != expected:
        raise AssertionError(f"expected {expected} attempts, found {attempted}")

    hash_cols = ["parent_frame_sha256", "catalogue_sha256", "gw_distance_sha256"]
    for (scenario, size), group in records.groupby(["scenario", "validation_size"]):
        if len(group) != draws:
            raise AssertionError("scenario/size group has incorrect draw count")
        if group["draw"].astype(int).tolist() != list(range(draws)):
            raise AssertionError("draw identifiers are missing or out of order")
    for scenario, group in records.groupby("scenario"):
        for col in hash_cols:
            values = group[col].dropna().astype(str).unique()
            if len(values) != 1 or len(values[0]) != 64:
                raise AssertionError(f"scenario {scenario} has inconsistent {col}")
            if any(c not in "0123456789abcdef" for c in values[0]):
                raise AssertionError(f"scenario {scenario} has malformed {col}")
        if group["validation_seed"].dropna().nunique() != 1:
            raise AssertionError("validation seed must be fixed within scenario")

        seed_values = group["validation_seed"].dropna().unique()
        if len(seed_values) != 1:
            raise AssertionError("missing validation seed")
        first = group.iloc[0]
        simulation = generate_synthetic_population(
            seed=int(first["simulation_seed"]),
            n_parent=int(first["parent_size"]),
            n_detected=int(first["n_detected"]),
            h0_true=float(first["h0_true"]),
            sigma_mpc=float(first["sigma_mpc"]),
            detection_threshold_mpc=float(first["detection_threshold_mpc"]),
            z_min=0.01, z_max=0.40,
        )
        parent_z = simulation["parent_truth"]["redshift"].to_numpy(dtype=float)
        distances = simulation["gw_events_observed"]["distance_mpc"].to_numpy(dtype=float)
        midpoint, half_range = 0.205, 0.195
        x_parent = (parent_z - midpoint) / half_range
        p_inclusion = expit(
            float(first["true_intercept"])
            + float(first["true_slope"]) * x_parent
            + float(first["true_quadratic"]) * (x_parent*x_parent - 1.0/3.0)
        )
        retained = np.random.default_rng(int(first["thinning_seed"])).random(len(parent_z)) < p_inclusion
        observed_z = parent_z[retained]
        hashes = {
            "parent_frame_sha256": hashlib.sha256(np.ascontiguousarray(parent_z).view(np.uint8)).hexdigest(),
            "catalogue_sha256": hashlib.sha256(np.ascontiguousarray(observed_z).view(np.uint8)).hexdigest(),
            "gw_distance_sha256": hashlib.sha256(np.ascontiguousarray(distances).view(np.uint8)).hexdigest(),
        }
        for col, digest in hashes.items():
            if str(first[col]) != digest:
                raise AssertionError(f"independent regeneration does not match {col}")
        rng = np.random.default_rng(int(seed_values[0]))
        ordered = group.sort_values(["validation_size", "draw"])
        for row in ordered.itertuples(index=False):
            if row.validation_indices is None or pd.isna(row.validation_indices) or not str(row.validation_indices):
                if row.status == "success":
                    raise AssertionError("successful attempt is missing sample indices")
                continue
            indices = np.fromstring(str(row.validation_indices), sep=",", dtype=int)
            size = int(row.validation_size)
            if indices.size != size or np.unique(indices).size != size:
                raise AssertionError("stored validation sample is not exact SRSWOR")
            if np.any((indices < 0) | (indices >= int(row.parent_size))):
                raise AssertionError("validation index outside parent frame")
            expected_indices = rng.choice(int(row.parent_size), size=size, replace=False)
            if not np.array_equal(indices, expected_indices):
                raise AssertionError("validation sample does not replay from declared seed")
            if int(row.n_validation_included) != int(np.sum(retained[indices])):
                raise AssertionError("validation inclusion count does not match regenerated parent frame")
    return {
        "status": "STRUCTURAL_AUDIT_PASS",
        "n_attempted": attempted,
        "n_successful": successful,
        "n_failed": failed,
        "n_scenarios": n_scenarios,
        "validation_sizes": sizes,
        "scope": "independent synthetic-input regeneration, SRSWOR bookkeeping, provenance, and seed replay; not estimator reconstruction",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    args = parser.parse_args()
    result = audit(args.artifact_dir)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
