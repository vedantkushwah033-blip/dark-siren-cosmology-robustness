"""Independently audit multi-replicate five-arm lifecycle and aggregate outputs.

This auditor intentionally does not import the multi-replicate runner. It
reconstructs the seed schedule, expected replicate/arm cells, paired contrasts,
and failure-inclusive metric bounds from saved raw posterior summaries.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ARMS = [
    "A_catalogue_only",
    "B_parent_oracle",
    "C_known_probability_ipw",
    "D_map_plugin_ipw",
    "E_uncertainty_integrated_ipw",
]


def _seed_schedule(master_seed: int, replicates: int) -> list[tuple[int, int]]:
    children = np.random.SeedSequence(master_seed).spawn(2 * replicates)
    values = [int(child.generate_state(1, dtype=np.uint32)[0]) for child in children]
    pairs = list(zip(values[0::2], values[1::2]))
    if len({value for pair in pairs for value in pair}) != 2 * replicates:
        raise AssertionError("seed schedule contains a collision")
    return pairs


def _expected_contrasts(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (replicate, n_validation), group in summary.groupby(
        ["replicate", "validation_size"], sort=True
    ):
        by_arm = {str(row.arm): row for row in group.itertuples(index=False)}
        if not all(arm in by_arm and by_arm[arm].status == "SUCCESS" for arm in ARMS):
            continue
        for arm in ARMS[1:]:
            for reference in ("A_catalogue_only", "C_known_probability_ipw"):
                left, right = by_arm[arm], by_arm[reference]
                rows.append({
                    "replicate": int(replicate),
                    "validation_size": int(n_validation),
                    "arm": arm,
                    "reference_arm": reference,
                    "median_error_difference": float(left.median_error - right.median_error),
                    "absolute_error_difference": float(abs(left.median_error) - abs(right.median_error)),
                    "squared_error_difference": float(left.median_error**2 - right.median_error**2),
                    "coverage_68_difference": int(left.q16 <= left.h0_true <= left.q84)
                    - int(right.q16 <= right.h0_true <= right.q84),
                })
    return pd.DataFrame(rows, columns=[
        "replicate", "validation_size", "arm", "reference_arm",
        "median_error_difference", "absolute_error_difference",
        "squared_error_difference", "coverage_68_difference",
    ])


def _expected_aggregate(summary: pd.DataFrame, ledger: pd.DataFrame, replicates: int) -> pd.DataFrame:
    rows = []
    for n_validation in sorted(ledger.validation_size.unique()):
        for arm in ARMS:
            observed = summary.loc[
                (summary.validation_size == n_validation)
                & (summary.arm == arm)
                & (summary.status == "SUCCESS")
            ]
            errors = observed.median_error.to_numpy(float)
            widths = (observed.q84 - observed.q16).to_numpy(float)
            covered = (
                (observed.q16.to_numpy(float) <= observed.h0_true.to_numpy(float))
                & (observed.h0_true.to_numpy(float) <= observed.q84.to_numpy(float))
            )
            failures = replicates - len(observed)
            squared_sum = float(np.sum(errors**2))
            absolute_sum = float(np.sum(np.abs(errors)))
            rows.append({
                "validation_size": int(n_validation),
                "arm": arm,
                "planned_replicates": replicates,
                "successful_replicates": len(observed),
                "failed_or_non_evaluable": failures,
                "success_only_bias": float(np.mean(errors)) if len(errors) else np.nan,
                "success_only_rmse": float(np.sqrt(np.mean(errors**2))) if len(errors) else np.nan,
                "success_only_mean_abs_error": float(np.mean(np.abs(errors))) if len(errors) else np.nan,
                "success_only_68pct_coverage": float(np.mean(covered)) if len(covered) else np.nan,
                "success_only_mean_68pct_width": float(np.mean(widths)) if len(widths) else np.nan,
                "failure_inclusive_mean_68pct_width_lower": float(np.sum(widths) / replicates),
                "failure_inclusive_mean_68pct_width_upper": float((np.sum(widths) + failures * 40.0) / replicates),
                "failure_inclusive_coverage_lower": float(np.sum(covered) / replicates),
                "failure_inclusive_coverage_upper": float((np.sum(covered) + failures) / replicates),
                "failure_inclusive_rmse_lower": float(np.sqrt(squared_sum / replicates)),
                "failure_inclusive_rmse_upper": float(np.sqrt((squared_sum + failures * 400.0) / replicates)),
                "failure_inclusive_abs_error_lower": float(absolute_sum / replicates),
                "failure_inclusive_abs_error_upper": float((absolute_sum + failures * 20.0) / replicates),
                "failure_inclusive_bias_lower": float((np.sum(errors) - failures * 20.0) / replicates),
                "failure_inclusive_bias_upper": float((np.sum(errors) + failures * 20.0) / replicates),
            })
    return pd.DataFrame(rows)


def _assert_frame_equal(actual: pd.DataFrame, expected: pd.DataFrame, label: str) -> None:
    if list(actual.columns) != list(expected.columns):
        raise AssertionError(f"{label}: column schema mismatch")
    if len(actual) != len(expected):
        raise AssertionError(f"{label}: row count {len(actual)} != expected {len(expected)}")
    for column in expected.columns:
        left, right = actual[column], expected[column]
        if pd.api.types.is_numeric_dtype(right):
            if not np.allclose(
                left.to_numpy(float), right.to_numpy(float),
                rtol=2e-10, atol=2e-10, equal_nan=True,
            ):
                raise AssertionError(f"{label}: numeric mismatch in {column}")
        elif not left.astype(str).equals(right.astype(str)):
            raise AssertionError(f"{label}: mismatch in {column}")


def audit(root: Path) -> dict[str, object]:
    manifest_path = root / "manifest.json"
    if not manifest_path.exists():
        raise AssertionError("top-level SHA-256 manifest is missing")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS":
        raise AssertionError("top-level manifest status is not PASS")
    listed = manifest.get("files", {})
    actual = {str(path.relative_to(root)) for path in root.rglob("*")
              if path.is_file() and path.name not in {"manifest.json", "independent_audit.json"}}
    if set(listed) != actual or int(manifest.get("file_count", -1)) != len(actual):
        raise AssertionError("top-level manifest file set/count mismatch")
    for relative, metadata in listed.items():
        path = root / relative
        if path.stat().st_size != int(metadata["bytes"]):
            raise AssertionError(f"top-level manifest size mismatch: {relative}")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != metadata["sha256"]:
            raise AssertionError(f"top-level manifest SHA-256 mismatch: {relative}")

    config = json.loads((root / "config.json").read_text(encoding="utf-8"))
    if config.get("status") != "engineering_pilot_not_confirmatory_not_observational":
        raise AssertionError("unexpected experiment status")
    replicates = int(config["replicates"])
    sizes = [int(value) for value in config["validation_sizes"]]
    expected_seeds = _seed_schedule(int(config["master_seed"]), replicates)
    declared = config["seed_schedule"]
    if len(declared) != replicates:
        raise AssertionError("seed schedule length mismatch")
    for index, (expected, declared_item) in enumerate(zip(expected_seeds, declared), start=1):
        if (
            int(declared_item["replicate"]) != index
            or int(declared_item["analysis_seed"]) != expected[0]
            or int(declared_item["validation_seed"]) != expected[1]
        ):
            raise AssertionError(f"seed schedule mismatch at replicate {index}")

    ledger = pd.read_csv(root / "replicate_ledger.csv")
    summary = pd.read_csv(root / "replicate_posterior_summaries.csv")
    audits = pd.read_csv(root / "replicate_audits.csv")
    expected_cells = replicates * len(sizes) * len(ARMS)
    if len(ledger) != expected_cells:
        raise AssertionError(f"ledger has {len(ledger)} rows; expected {expected_cells}")
    if ledger.duplicated(["replicate", "validation_size", "arm"]).any():
        raise AssertionError("duplicate lifecycle ledger cells")
    if set(ledger.status) - {"SUCCESS", "FAILED", "NON_EVALUABLE"}:
        raise AssertionError("unknown lifecycle status")
    if set(ledger.audit_status) - {"PASS", "NOT_RUN", "NOT_PASS"}:
        raise AssertionError("unknown independent-audit status")

    successful_replicates = 0
    for index, (analysis_seed, validation_seed) in enumerate(expected_seeds, start=1):
        rows = ledger.loc[ledger.replicate == index]
        if len(rows) != len(sizes) * len(ARMS):
            raise AssertionError(f"replicate {index} ledger is incomplete")
        if set(rows.analysis_seed) != {analysis_seed} or set(rows.validation_seed) != {validation_seed}:
            raise AssertionError(f"replicate {index} seed mismatch")
        if set(rows.status) == {"SUCCESS"} and set(rows.audit_status) == {"PASS"}:
            successful_replicates += 1
            rep_summary = summary.loc[summary.replicate == index]
            if len(rep_summary) != len(sizes) * len(ARMS):
                raise AssertionError(f"successful replicate {index} summary rows incomplete")
            if set(rep_summary.validation_size.astype(int)) != set(sizes):
                raise AssertionError(f"successful replicate {index} validation-size set mismatch")
            for size in sizes:
                if len(rep_summary.loc[rep_summary.validation_size == size]) != len(ARMS):
                    raise AssertionError(f"successful replicate {index}, n={size}: arm rows incomplete")
            if set(rep_summary.audit_status) != {"PASS"}:
                raise AssertionError(f"successful replicate {index} missing audit PASS")
            report_path = root / "replicates" / f"replicate_{index:04d}" / "audit.json"
            if not report_path.exists():
                raise AssertionError(f"replicate {index} raw-artifact audit missing")
            report = json.loads(report_path.read_text())
            audit_row = audits.loc[audits.replicate == index]
            if report.get("status") != "PASS" or len(audit_row) != 1:
                raise AssertionError(f"replicate {index} raw-artifact audit missing or failed")
            checks = {
                "manifest_files_verified": int(report["manifest_files_verified"]),
                "posterior_summary_rows_reconstructed": int(report["posterior_summary_rows_reconstructed"]),
            }
            for field, expected_value in checks.items():
                if int(audit_row.iloc[0][field]) != expected_value:
                    raise AssertionError(f"replicate {index} audit ledger mismatch in {field}")
            if not np.isclose(
                float(audit_row.iloc[0]["maximum_log_likelihood_absolute_difference"]),
                float(report["maximum_log_likelihood_absolute_difference"]),
                rtol=1e-12, atol=1e-15,
            ):
                raise AssertionError(f"replicate {index} audit likelihood discrepancy mismatch")
        else:
            if set(rows.status) not in ({"FAILED"}, {"NON_EVALUABLE"}):
                raise AssertionError(f"replicate {index} has inconsistent partial lifecycle statuses")
            if len(summary.loc[summary.replicate == index]):
                raise AssertionError(f"failed replicate {index} leaked partial posterior summaries")

    if len(audits) != replicates or audits.replicate.tolist() != list(range(1, replicates + 1)):
        raise AssertionError("replicate audit ledger is incomplete or unordered")
    expected_audit_status = [
        "PASS" if set(ledger.loc[ledger.replicate == i, "status"]) == {"SUCCESS"} else "NOT_PASS"
        for i in range(1, replicates + 1)
    ]
    if audits.status.tolist() != expected_audit_status:
        raise AssertionError("replicate audit statuses disagree with lifecycle ledger")

    expected_arms = set(ARMS)
    if len(summary):
        if summary.duplicated(["replicate", "validation_size", "arm"]).any():
            raise AssertionError("duplicate posterior summary cells")
        if set(summary.arm) != expected_arms:
            raise AssertionError("posterior summary arm set mismatch")
        if not set(summary.status).issubset({"SUCCESS"}):
            raise AssertionError("non-success posterior rows must not be included")
        if not np.all(np.isfinite(summary[["median_error", "q16", "q84"]].to_numpy(float))):
            raise AssertionError("non-finite posterior summary value")

    expected_aggregate = _expected_aggregate(summary, ledger, replicates)
    actual_aggregate = pd.read_csv(root / "aggregate_metrics.csv")
    _assert_frame_equal(actual_aggregate, expected_aggregate, "aggregate metrics")
    expected_contrasts = _expected_contrasts(summary)
    actual_contrasts = pd.read_csv(root / "paired_contrasts.csv")
    _assert_frame_equal(actual_contrasts, expected_contrasts, "paired contrasts")
    return {
        "status": "PASS",
        "planned_replicates": replicates,
        "successful_replicates": successful_replicates,
        "failed_replicates": replicates - successful_replicates,
        "planned_lifecycle_cells": expected_cells,
        "posterior_summary_rows": len(summary),
        "aggregate_metric_rows_reconstructed": len(expected_aggregate),
        "paired_contrast_rows_reconstructed": len(expected_contrasts),
        "claim_limit": "Audits pilot orchestration and arithmetic; does not establish scientific calibration or estimator superiority.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = audit(args.artifact_dir)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
