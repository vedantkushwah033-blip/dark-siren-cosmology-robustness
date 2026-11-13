"""Run a paired multi-replicate five-arm engineering pilot.

Every replicate has one shared parent/catalogue/event realization and one
independent nested validation frame. All attempts, including failures, remain
in the lifecycle ledger. This is not a confirmatory Monte Carlo study.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from run_five_arm_smoke import run as run_single_replicate


ARMS = [
    "A_catalogue_only",
    "B_parent_oracle",
    "C_known_probability_ipw",
    "D_map_plugin_ipw",
    "E_uncertainty_integrated_ipw",
]


def _parse_sizes(raw: str) -> list[int]:
    try:
        values = [int(part.strip()) for part in raw.split(",") if part.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("validation sizes must be comma-separated integers") from exc
    if not values or sorted(set(values)) != values or values[0] < 5:
        raise argparse.ArgumentTypeError("validation sizes must be unique, increasing integers >= 5")
    return values


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _seed_schedule(master_seed: int, replicates: int) -> list[tuple[int, int]]:
    if master_seed < 0 or replicates < 1:
        raise ValueError("master_seed must be nonnegative and replicates positive")
    children = np.random.SeedSequence(master_seed).spawn(2 * replicates)
    values = [int(child.generate_state(1, dtype=np.uint32)[0]) for child in children]
    pairs = list(zip(values[0::2], values[1::2]))
    flattened = [value for pair in pairs for value in pair]
    if len(set(flattened)) != len(flattened):
        raise ArithmeticError("seed derivation produced a collision")
    return pairs


def _paired_contrasts(summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
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
                    "coverage_68_difference": int(left.q16 <= left.h0_true <= left.q84) - int(right.q16 <= right.h0_true <= right.q84),
                })
    return pd.DataFrame(rows, columns=[
        "replicate", "validation_size", "arm", "reference_arm",
        "median_error_difference", "absolute_error_difference",
        "squared_error_difference", "coverage_68_difference",
    ])


def _aggregate(summary: pd.DataFrame, ledger: pd.DataFrame, replicates: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    planned_sizes = sorted(int(value) for value in ledger["validation_size"].unique())
    for n_validation in planned_sizes:
        for arm in ARMS:
            observed = summary.loc[
                (summary.validation_size == n_validation)
                & (summary.arm == arm)
                & (summary.status == "SUCCESS")
            ]
            errors = observed["median_error"].to_numpy(float)
            widths = (observed["q84"] - observed["q16"]).to_numpy(float)
            covered = (
                (observed["q16"].to_numpy(float) <= observed["h0_true"].to_numpy(float))
                & (observed["h0_true"].to_numpy(float) <= observed["q84"].to_numpy(float))
            )
            failures = replicates - len(observed)
            squared_sum = float(np.sum(errors**2))
            absolute_sum = float(np.sum(np.abs(errors)))
            rows.append({
                "validation_size": n_validation,
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
    aggregate = pd.DataFrame(rows)
    contrasts = _paired_contrasts(summary)
    return aggregate, contrasts


def run(args: argparse.Namespace) -> dict[str, object]:
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError("output directory must be empty; choose a fresh path to preserve prior artifacts")
    sizes = args.validation_sizes
    if sizes[-1] > args.validation_frame:
        raise ValueError("largest validation size cannot exceed validation frame")
    if not 50.0 <= float(args.h0_true) <= 90.0:
        raise ValueError("failure-inclusive bounds require h0_true within the declared [50, 90] H0 prior support")
    seed_pairs = _seed_schedule(args.master_seed, args.replicates)
    ledger_rows: list[dict[str, object]] = []
    summary_parts: list[pd.DataFrame] = []
    audit_rows: list[dict[str, object]] = []
    for index, (analysis_seed, validation_seed) in enumerate(seed_pairs, start=1):
        replicate_dir = output / "replicates" / f"replicate_{index:04d}"
        replicate_dir.mkdir(parents=True, exist_ok=True)
        status = "SUCCESS"
        error_message = ""
        audit_status = "NOT_RUN"
        try:
            single_args = SimpleNamespace(
                seed=analysis_seed,
                validation_seed=validation_seed,
                n_parent=args.n_parent,
                n_detected=args.n_detected,
                validation_frame=args.validation_frame,
                validation_sizes=sizes,
                catalogue_floor=args.catalogue_floor,
                h0_true=args.h0_true,
                h0_points=args.h0_points,
                scott_multiplier=args.scott_multiplier,
                quadrature_order=args.quadrature_order,
                sigma_mpc=args.sigma_mpc,
                detection_threshold_mpc=args.detection_threshold_mpc,
                output_dir=replicate_dir,
            )
            run_single_replicate(single_args)
            audit_path = replicate_dir / "audit.json"
            command = [
                sys.executable,
                str(Path(__file__).with_name("audit_five_arm_smoke.py")),
                str(replicate_dir),
                "--output",
                str(audit_path),
            ]
            completed = subprocess.run(command, check=False, capture_output=True, text=True)
            if completed.returncode != 0:
                audit_status = "NOT_PASS"
                raise RuntimeError(
                    "independent artifact audit failed: "
                    + (completed.stderr[-3000:] or completed.stdout[-3000:])
                )
            audit_report = json.loads(audit_path.read_text(encoding="utf-8"))
            if audit_report.get("status") != "PASS":
                audit_status = "NOT_PASS"
                raise RuntimeError("independent artifact auditor did not return PASS")
            audit_status = "PASS"
            part = pd.read_csv(replicate_dir / "results" / "posterior_summaries.csv")
            part.insert(0, "replicate", index)
            part.insert(1, "analysis_seed", analysis_seed)
            part.insert(2, "validation_seed", validation_seed)
            part["audit_status"] = audit_status
            summary_parts.append(part)
            audit_rows.append({
                "replicate": index,
                "status": "PASS",
                "manifest_files_verified": audit_report["manifest_files_verified"],
                "posterior_summary_rows_reconstructed": audit_report["posterior_summary_rows_reconstructed"],
                "maximum_log_likelihood_absolute_difference": audit_report["maximum_log_likelihood_absolute_difference"],
            })
        except Exception as exc:  # Preserve the attempt; never silently drop it.
            status = "FAILED"
            error_message = f"{type(exc).__name__}: {exc}"[:4000]
            audit_rows.append({
                "replicate": index,
                "status": "NOT_PASS",
                "manifest_files_verified": 0,
                "posterior_summary_rows_reconstructed": 0,
                "maximum_log_likelihood_absolute_difference": np.nan,
            })
        for n_validation in sizes:
            for arm in ARMS:
                ledger_rows.append({
                    "replicate": index,
                    "analysis_seed": analysis_seed,
                    "validation_seed": validation_seed,
                    "validation_size": n_validation,
                    "arm": arm,
                    "status": status,
                    "audit_status": audit_status,
                    "error": error_message,
                })
    ledger = pd.DataFrame(ledger_rows)
    summaries = pd.concat(summary_parts, ignore_index=True) if summary_parts else pd.DataFrame(
        columns=[
            "replicate", "analysis_seed", "validation_seed", "validation_size", "arm",
            "status", "audit_status", "h0_true", "median_error", "posterior_median",
            "q025", "q05", "q16", "q84", "q95", "q975",
        ]
    )
    aggregate, contrasts = _aggregate(summaries, ledger, args.replicates)
    ledger.to_csv(output / "replicate_ledger.csv", index=False)
    summaries.to_csv(output / "replicate_posterior_summaries.csv", index=False)
    aggregate.to_csv(output / "aggregate_metrics.csv", index=False)
    contrasts.to_csv(output / "paired_contrasts.csv", index=False)
    pd.DataFrame(audit_rows).to_csv(output / "replicate_audits.csv", index=False)
    config = {
        "experiment": "matched_five_arm_multi_replicate_engineering_pilot",
        "status": "engineering_pilot_not_confirmatory_not_observational",
        "master_seed": args.master_seed,
        "seed_schedule": [
            {"replicate": i, "analysis_seed": pair[0], "validation_seed": pair[1]}
            for i, pair in enumerate(seed_pairs, start=1)
        ],
        "replicates": args.replicates,
        "n_parent": args.n_parent,
        "n_detected": args.n_detected,
        "validation_frame": args.validation_frame,
        "validation_sizes": sizes,
        "catalogue_floor": args.catalogue_floor,
        "h0_true": args.h0_true,
        "h0_points": args.h0_points,
        "quadrature_order": args.quadrature_order,
        "sigma_mpc": args.sigma_mpc,
        "detection_threshold_mpc": args.detection_threshold_mpc,
        "scott_multiplier": args.scott_multiplier,
        "bandwidth_selection_status": "engineering_fixture_not_training_holdout_selected",
        "failure_bound_assumptions": {
            "h0_prior_support": [50.0, 90.0],
            "absolute_median_error_max": 20.0,
            "squared_median_error_max": 400.0,
            "coverage_failure_indicator_bounds": [0, 1],
            "interval_width_bounds": [0.0, 40.0],
        },
        "claim_limit": "Engineering pilot only; no estimator ranking, calibration, or observational inference.",
    }
    (output / "config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    files = sorted(
        path for path in output.rglob("*")
        if path.is_file() and path.name != "manifest.json"
    )
    manifest = {
        "status": "PASS",
        "file_count": len(files),
        "files": {
            str(path.relative_to(output)): {
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for path in files
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return {
        "status": "PASS" if (ledger.status == "SUCCESS").all() and (ledger.audit_status == "PASS").all() else "INCOMPLETE",
        "replicates": args.replicates,
        "successful_replicates": int((pd.DataFrame(audit_rows).status == "PASS").sum()),
        "failed_replicates": int((ledger.groupby("replicate").status.first() != "SUCCESS").sum()),
        "posterior_summary_rows": len(summaries),
        "output_dir": str(output),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master-seed", type=int, default=20261115)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--n-parent", type=int, default=1200)
    parser.add_argument("--n-detected", type=int, default=4)
    parser.add_argument("--validation-frame", type=int, default=1500)
    parser.add_argument("--validation-sizes", type=_parse_sizes, default=_parse_sizes("100,300,1000"))
    parser.add_argument("--catalogue-floor", type=float, default=0.10)
    parser.add_argument("--h0-true", type=float, default=70.0)
    parser.add_argument("--h0-points", type=int, default=41)
    parser.add_argument("--scott-multiplier", type=float, default=0.7)
    parser.add_argument("--quadrature-order", type=int, default=5)
    parser.add_argument("--sigma-mpc", type=float, default=140.0)
    parser.add_argument("--detection-threshold-mpc", type=float, default=250.0)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.replicates < 1 or args.master_seed < 0:
        parser.error("replicates must be positive and master seed nonnegative")
    print(json.dumps(run(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
