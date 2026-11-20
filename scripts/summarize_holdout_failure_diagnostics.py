#!/usr/bin/env python3
"""Summarize frozen holdout failures and inspect quadrature-ESS audit discrepancies."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _cell_name(replicate: int, gamma: float) -> str:
    return f"rep{replicate:03d}_g{gamma:+g}".replace("+", "p").replace("-", "m")


def inspect_quadrature_ess(root: Path, failed: pd.DataFrame) -> dict:
    """Compare the saved E-arm quadrature ESS with its archived normalized weights."""
    children = failed[["replicate", "gamma"]].drop_duplicates().sort_values(["gamma", "replicate"])
    records = []
    missing = []
    for row in children.itertuples(index=False):
        cell = _cell_name(int(row.replicate), float(row.gamma))
        child_dir = root / "replicates" / cell
        diag_path = child_dir / "results" / "weight_diagnostics.csv"
        archive_path = child_dir / "results" / "log_likelihoods.npz"
        if not diag_path.exists() or not archive_path.exists():
            missing.append(cell)
            continue
        diagnostics = pd.read_csv(diag_path)
        with np.load(archive_path, allow_pickle=False) as archive:
            for n in sorted(int(v) for v in diagnostics["validation_size"].unique()):
                rows = diagnostics.loc[
                    (diagnostics["validation_size"] == n)
                    & (diagnostics["arm"] == "E_uncertainty_integrated_ipw")
                ]
                key = f"posterior_weights_n{n}"
                if len(rows) != 1 or key not in archive:
                    missing.append(f"{cell}:n{n}:missing-E-row-or-{key}")
                    continue
                weights = np.asarray(archive[key], dtype=float)
                expected = float(1.0 / np.sum(weights * weights))
                observed = float(rows.iloc[0]["quadrature_ess"])
                records.append({
                    "cell": cell,
                    "replicate": int(row.replicate),
                    "gamma": float(row.gamma),
                    "validation_size": n,
                    "stored_quadrature_ess": observed,
                    "archive_weight_reconstruction": expected,
                    "absolute_difference": abs(observed - expected),
                    "weight_count": int(weights.size),
                    "weight_sum": float(np.sum(weights)),
                })
    frame = pd.DataFrame(records)
    if len(frame):
        summary = {
            "comparisons": int(len(frame)),
            "maximum_absolute_difference": float(frame["absolute_difference"].max()),
            "median_absolute_difference": float(frame["absolute_difference"].median()),
            "comparisons_over_audit_tolerance": int(
                (frame["absolute_difference"] > (2e-12 + 2e-10 * np.abs(frame["archive_weight_reconstruction"]))).sum()
            ),
            "example_largest_differences": frame.sort_values(
                "absolute_difference", ascending=False
            ).head(10).to_dict(orient="records"),
        }
    else:
        summary = {"comparisons": 0}
    summary["missing_child_files"] = missing
    return summary


def summarize_ledger(ledger: pd.DataFrame, config: dict, root: Path) -> tuple[dict, pd.DataFrame]:
    required = {"cell_id", "replicate", "gamma", "validation_size", "arm", "status", "error"}
    missing = required - set(ledger.columns)
    if missing:
        raise ValueError(f"replicate ledger is missing columns: {sorted(missing)}")
    if ledger["cell_id"].duplicated().any():
        raise ValueError("replicate ledger contains duplicate cell_id values")

    failed = ledger.loc[ledger["status"] != "SUCCESS"].copy()
    failed["error"] = failed["error"].fillna("").astype(str)
    failed["error_signature"] = failed["error"].map(
        lambda value: value.splitlines()[0][:300] if value.strip() else "UNKNOWN_ERROR"
    )
    child = (
        ledger.groupby(["replicate", "gamma"], sort=True)
        .agg(
            planned_cells=("cell_id", "size"),
            failed_cells=("status", lambda values: int((values != "SUCCESS").sum())),
            error_signatures=("error", lambda values: int(values.fillna("").astype(str).nunique())),
        )
        .reset_index()
    )
    affected = child.loc[child["failed_cells"] > 0].copy()
    signatures = (
        failed.groupby(["gamma", "error_signature"], dropna=False)
        .size().rename("failed_cells").reset_index()
        .sort_values(["failed_cells", "gamma"], ascending=[False, True])
    )
    examples = []
    for _, row in signatures.iterrows():
        match = failed.loc[
            (failed["gamma"] == row["gamma"])
            & (failed["error_signature"] == row["error_signature"])
        ].iloc[0]
        examples.append({
            "gamma": float(row["gamma"]),
            "error_signature": str(row["error_signature"]),
            "failed_cells": int(row["failed_cells"]),
            "example_replicate": int(match["replicate"]),
            "example_validation_size": int(match["validation_size"]),
            "example_arm": str(match["arm"]),
            "error_excerpt": str(match["error"])[:2000],
        })

    expected_cells_per_child = len(config["validation_sizes"]) * 5
    bad_child_denominators = child.loc[child["planned_cells"] != expected_cells_per_child]
    quadrature_ess = inspect_quadrature_ess(root, failed)
    report = {
        "status": "PASS_DIAGNOSTIC" if len(bad_child_denominators) == 0 else "FAIL_LEDGER_STRUCTURE",
        "master_seed": int(config["master_seed"]),
        "truth_mode": str(config.get("truth_mode", "fixed_truth")),
        "planned_cells": int(len(ledger)),
        "successful_cells": int((ledger["status"] == "SUCCESS").sum()),
        "failed_cells": int(len(failed)),
        "affected_replicate_gamma_runs": int(len(affected)),
        "planned_replicate_gamma_runs": int(len(child)),
        "failed_cells_by_gamma": {str(float(k)): int(v) for k, v in failed.groupby("gamma").size().to_dict().items()},
        "affected_child_runs_by_gamma": {str(float(k)): int(v) for k, v in affected.groupby("gamma").size().to_dict().items()},
        "error_signature_count": int(len(signatures)),
        "error_signatures": examples,
        "child_denominator_mismatches": int(len(bad_child_denominators)),
        "quadrature_ess_reconstruction": quadrature_ess,
        "interpretation": (
            "Frozen-artifact diagnosis only. This does not convert an incomplete holdout into a valid "
            "predictive test and does not authorize replacing the predeclared seed."
        ),
    }
    return report, child


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/holdout_failure_diagnostics"))
    args = parser.parse_args()
    config = json.loads((args.artifact_dir / "config.json").read_text(encoding="utf-8"))
    ledger = pd.read_csv(args.artifact_dir / "replicate_ledger.csv")
    report, child = summarize_ledger(ledger, config, args.artifact_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "failure_diagnostics.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    child.to_csv(args.output_dir / "replicate_gamma_failure_counts.csv", index=False)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
