"""Paired multi-replicate degree-one versus degree-two inclusion-model pilot.

Both degrees receive identical parent/event/validation seeds. Each degree's
aggregate is independently audited by the existing replicate auditor, then
cross-degree raw inputs and posterior summaries are paired. Engineering pilot
only: this is not a calibration or observational analysis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

INPUTS = [
    "observable/catalogue_observed.csv",
    "observable/gw_events_observed.csv",
    "truth/parent_truth.csv",
    "truth/gw_event_truth.csv",
    "truth/validation_frame_truth.csv",
]


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _paired_child_inputs(degree1: Path, degree2: Path, sizes: list[int]) -> int:
    checked = 0
    paths = INPUTS + [f"observable/validation_n{n}.csv" for n in sizes]
    for rel in paths:
        left, right = pd.read_csv(degree1 / rel), pd.read_csv(degree2 / rel)
        # The catalogue and validation outcomes must be identical across fitted
        # inclusion-model degrees; only the fitted model changes.
        pd.testing.assert_frame_equal(left, right, check_dtype=False, check_exact=True)
        checked += 1
    return checked


def _degree_contrasts(one: pd.DataFrame, two: pd.DataFrame) -> pd.DataFrame:
    keys = ["replicate", "gamma", "validation_size", "arm"]
    fields = ["status", "h0_true", "median_error", "posterior_median",
              "q025", "q16", "q84", "q975"]
    a = one[keys + fields].rename(columns={c: f"{c}_d1" for c in fields})
    b = two[keys + fields].rename(columns={c: f"{c}_d2" for c in fields})
    merged = a.merge(b, on=keys, how="outer", validate="one_to_one", indicator=True)
    merged = merged.loc[merged["_merge"] == "both"].copy()
    merged = merged.loc[
        (merged["status_d1"] == "SUCCESS") & (merged["status_d2"] == "SUCCESS")
    ].copy()
    if merged.empty:
        return merged.drop(columns="_merge")
    if not np.allclose(merged["h0_true_d1"], merged["h0_true_d2"], rtol=0, atol=0):
        raise AssertionError("H0 truth differs across model degrees")
    merged["degree2_minus_degree1_median"] = (
        merged["posterior_median_d2"] - merged["posterior_median_d1"]
    )
    merged["degree2_minus_degree1_median_error"] = (
        merged["median_error_d2"] - merged["median_error_d1"]
    )
    merged["absolute_error_change_d2_minus_d1"] = (
        np.abs(merged["median_error_d2"]) - np.abs(merged["median_error_d1"])
    )
    merged["squared_error_change_d2_minus_d1"] = (
        merged["median_error_d2"] ** 2 - merged["median_error_d1"] ** 2
    )
    merged["coverage68_d1"] = (
        (merged["q16_d1"] <= merged["h0_true_d1"]) &
        (merged["h0_true_d1"] <= merged["q84_d1"])
    )
    merged["coverage68_d2"] = (
        (merged["q16_d2"] <= merged["h0_true_d2"]) &
        (merged["h0_true_d2"] <= merged["q84_d2"])
    )
    merged["coverage68_change_d2_minus_d1"] = (
        merged["coverage68_d2"].astype(int) - merged["coverage68_d1"].astype(int)
    )
    return merged.drop(columns="_merge")


def run(args: argparse.Namespace) -> dict:
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("output directory must be empty")
    if args.replicates < 2 or args.n_detected < 1 or args.n_parent < 100:
        raise ValueError("invalid replicate or simulation dimensions")
    if args.quadrature_order < 3 or args.quadrature_order > 9:
        raise ValueError("quadrature order must be between 3 and 9")
    if not args.gammas or sorted(set(args.gammas)) != args.gammas:
        raise ValueError("gammas must be unique and increasing")
    if not args.validation_sizes or sorted(set(args.validation_sizes)) != args.validation_sizes:
        raise ValueError("validation sizes must be unique and increasing")
    if args.validation_sizes[0] < 5 or args.validation_sizes[-1] > args.validation_frame:
        raise ValueError("validation sizes must fit inside the validation frame")
    if not args.h0_prior_min <= args.h0_true <= args.h0_prior_max:
        raise ValueError("H0 truth must lie inside prior support")

    runner = Path(__file__).with_name("run_quadratic_five_arm_replicates.py")
    auditor = Path(__file__).with_name("audit_quadratic_five_arm_replicates.py")
    degree_roots = {}
    for degree in (1, 2):
        root = out / f"degree{degree}"
        cmd = [
            sys.executable, str(runner),
            "--master-seed", str(args.master_seed),
            "--replicates", str(args.replicates),
            "--gammas", ",".join(map(str, args.gammas)),
            "--n-parent", str(args.n_parent),
            "--n-detected", str(args.n_detected),
            "--validation-frame", str(args.validation_frame),
            "--validation-sizes", ",".join(map(str, args.validation_sizes)),
            "--curvature-beta", str(args.curvature_beta),
            "--inclusion-degree", str(degree),
            "--h0-true", str(args.h0_true),
            "--h0-prior-min", str(args.h0_prior_min),
            "--h0-prior-max", str(args.h0_prior_max),
            "--h0-points", str(args.h0_points),
            "--scott-multiplier", str(args.scott_multiplier),
            "--quadrature-order", str(args.quadrature_order),
            "--sigma-mpc", str(args.sigma_mpc),
            "--detection-threshold-mpc", str(args.detection_threshold_mpc),
            "--output-dir", str(root),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        (out / f"degree{degree}.stdout.txt").write_text(proc.stdout)
        (out / f"degree{degree}.stderr.txt").write_text(proc.stderr)
        if proc.returncode:
            raise RuntimeError(f"degree-{degree} aggregate runner failed: {proc.stderr[-4000:]}")
        audit_path = root / "independent_aggregate_audit.json"
        audited = subprocess.run(
            [sys.executable, str(auditor), str(root), "--output", str(audit_path)],
            capture_output=True, text=True,
        )
        (out / f"degree{degree}.audit.stdout.txt").write_text(audited.stdout)
        (out / f"degree{degree}.audit.stderr.txt").write_text(audited.stderr)
        if audited.returncode:
            raise RuntimeError(f"degree-{degree} independent aggregate audit failed: {audited.stderr[-4000:]}")
        report = json.loads(audit_path.read_text())
        if report.get("status") != "PASS":
            raise AssertionError(f"degree-{degree} aggregate auditor did not PASS")
        degree_roots[degree] = root

    summaries = {
        d: pd.read_csv(degree_roots[d] / "replicate_posterior_summaries.csv")
        for d in (1, 2)
    }
    diagnostics = {
        d: pd.read_csv(degree_roots[d] / "replicate_weight_diagnostics.csv")
        for d in (1, 2)
    }
    ledgers = {
        d: pd.read_csv(degree_roots[d] / "replicate_ledger.csv")
        for d in (1, 2)
    }
    pairing_rows = []
    for rep in range(args.replicates):
        for gamma in args.gammas:
            cell = f"rep{rep:03d}_g{gamma:+g}".replace("+", "p").replace("-", "m")
            left = degree_roots[1] / "replicates" / cell
            right = degree_roots[2] / "replicates" / cell
            if not left.is_dir() or not right.is_dir():
                pairing_rows.append({
                    "replicate": rep, "gamma": gamma,
                    "paired_input_files_checked": 0, "status": "NOT_EVALUABLE_CHILD_FAILURE",
                })
                continue
            checked = _paired_child_inputs(left, right, args.validation_sizes)
            pairing_rows.append({
                "replicate": rep, "gamma": gamma,
                "paired_input_files_checked": checked, "status": "PASS",
            })

    contrasts = _degree_contrasts(summaries[1], summaries[2])
    contrasts.to_csv(out / "paired_degree_posterior_contrasts.csv", index=False)
    pd.DataFrame(pairing_rows).to_csv(out / "cross_degree_input_pairing.csv", index=False)
    pd.concat([
        summaries[1].assign(inclusion_degree=1),
        summaries[2].assign(inclusion_degree=2),
    ], ignore_index=True).to_csv(out / "all_degree_posterior_summaries.csv", index=False)
    pd.concat([
        diagnostics[1].assign(inclusion_degree=1),
        diagnostics[2].assign(inclusion_degree=2),
    ], ignore_index=True).to_csv(out / "all_degree_weight_diagnostics.csv", index=False)
    pd.concat([
        ledgers[1].assign(inclusion_degree=1),
        ledgers[2].assign(inclusion_degree=2),
    ], ignore_index=True).to_csv(out / "combined_failure_ledger.csv", index=False)

    aggregate_rows = []
    for degree in (1, 2):
        frame = pd.read_csv(degree_roots[degree] / "aggregate_metrics.csv")
        frame.insert(0, "inclusion_degree", degree)
        aggregate_rows.append(frame)
    pd.concat(aggregate_rows, ignore_index=True).to_csv(out / "aggregate_metrics_by_degree.csv", index=False)

    config = {
        "experiment": "paired_degree1_degree2_multireplicate_quadratic_logit_inclusion_pilot",
        "status": "engineering_pilot_not_confirmatory_not_observational",
        "master_seed": args.master_seed, "replicates": args.replicates,
        "gammas": args.gammas, "n_parent": args.n_parent,
        "n_detected": args.n_detected, "validation_frame": args.validation_frame,
        "validation_sizes": args.validation_sizes, "quadrature_order": args.quadrature_order,
        "h0_true": args.h0_true, "h0_prior_support": [args.h0_prior_min, args.h0_prior_max],
        "degree1_aggregate_audit": "PASS", "degree2_aggregate_audit": "PASS",
        "cross_degree_pairing_blocks": len(pairing_rows),
        "paired_input_files_checked": int(sum(x["paired_input_files_checked"] for x in pairing_rows)),
        "paired_posterior_cells": len(contrasts),
        "claim_limit": "Paired synthetic sensitivity pilot only; not calibration proof, estimator superiority, or observational inference.",
    }
    (out / "config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "manifest.json")
    manifest = {"status": "PASS", "file_count": len(files), "files": {
        str(p.relative_to(out)): {"bytes": p.stat().st_size, "sha256": _sha(p)}
        for p in files
    }}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    statuses = pd.concat([ledgers[1], ledgers[2]]).status.value_counts().to_dict()
    return {
        "status": "PASS" if all(v == "SUCCESS" for v in statuses) else "INCOMPLETE",
        "degree1_degree2_child_status_counts": statuses,
        "cross_degree_pairing_blocks": len(pairing_rows),
        "paired_input_files_checked": config["paired_input_files_checked"],
        "paired_posterior_cells": len(contrasts),
        "degree1_aggregate_audit": "PASS", "degree2_aggregate_audit": "PASS",
        "output_dir": str(out),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--master-seed", type=int, default=202611201)
    p.add_argument("--replicates", type=int, default=3)
    p.add_argument("--gammas", type=lambda s: [float(x) for x in s.split(",")], default=[-1.0, 0.0, 1.0])
    p.add_argument("--n-parent", type=int, default=1200)
    p.add_argument("--n-detected", type=int, default=30)
    p.add_argument("--validation-frame", type=int, default=500)
    p.add_argument("--validation-sizes", type=lambda s: [int(x) for x in s.split(",")], default=[100, 300, 500])
    p.add_argument("--curvature-beta", type=float, default=-3.9)
    p.add_argument("--h0-true", type=float, default=70.0)
    p.add_argument("--h0-prior-min", type=float, default=50.0)
    p.add_argument("--h0-prior-max", type=float, default=90.0)
    p.add_argument("--h0-points", type=int, default=81)
    p.add_argument("--scott-multiplier", type=float, default=0.7)
    p.add_argument("--quadrature-order", type=int, default=5)
    p.add_argument("--sigma-mpc", type=float, default=140.0)
    p.add_argument("--detection-threshold-mpc", type=float, default=250.0)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    print(json.dumps(run(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
