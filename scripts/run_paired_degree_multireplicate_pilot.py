"""Paired multi-replicate degree-one versus degree-two inclusion-model pilot.

Runs the existing independently audited five-arm replicate workflow twice with
identical master seeds. It verifies identical synthetic inputs across model
degrees and exact invariance of A/B/C likelihood curves, then computes paired
degree contrasts. This is a finite synthetic engineering pilot, not
confirmatory calibration or an observational analysis.
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
from scipy.stats import t

ARMS = [
    "A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
    "D_map_plugin_ipw", "E_uncertainty_integrated_ipw",
]
COMMON_ARMS = ARMS[:3]
INPUT_FILES = [
    "observable/catalogue_observed.csv",
    "observable/gw_events_observed.csv",
    "truth/parent_truth.csv",
    "truth/gw_event_truth.csv",
    "truth/validation_frame_truth.csv",
]


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _write_manifest(out: Path, status: str) -> dict:
    """Write an integrity manifest carrying the actual experiment lifecycle status."""
    files = sorted(
        p for p in out.rglob("*")
        if p.is_file() and p.name != "manifest.json"
    )
    manifest = {
        "status": status,
        "file_count": len(files),
        "files": {
            str(p.relative_to(out)): {
                "bytes": p.stat().st_size,
                "sha256": _sha(p),
            }
            for p in files
        },
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return manifest


def _paired_summary(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    """Summarize degree-2 minus degree-1 replicate-paired posterior errors."""
    keys = ["replicate", "gamma", "validation_size", "arm"]
    a = left[keys + ["median_error", "q16", "q84", "h0_true"]].copy()
    b = right[keys + ["median_error", "q16", "q84", "h0_true"]].copy()
    paired = a.merge(b, on=keys, suffixes=("_d1", "_d2"), validate="one_to_one")
    if len(paired) != len(a) or len(paired) != len(b):
        raise AssertionError("degree-one and degree-two summary cells are not paired")
    if not np.array_equal(paired.h0_true_d1.to_numpy(), paired.h0_true_d2.to_numpy()):
        raise AssertionError("H0 truths differ across model degrees")
    paired["median_error_difference_d2_minus_d1"] = (
        paired.median_error_d2 - paired.median_error_d1
    )
    paired["absolute_error_difference_d2_minus_d1"] = (
        paired.median_error_d2.abs() - paired.median_error_d1.abs()
    )
    paired["squared_error_difference_d2_minus_d1"] = (
        paired.median_error_d2.pow(2) - paired.median_error_d1.pow(2)
    )
    paired["coverage68_d1"] = (
        (paired.q16_d1 <= paired.h0_true_d1) & (paired.h0_true_d1 <= paired.q84_d1)
    )
    paired["coverage68_d2"] = (
        (paired.q16_d2 <= paired.h0_true_d2) & (paired.h0_true_d2 <= paired.q84_d2)
    )
    return paired


def _aggregate_paired(paired: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = [
        "median_error_difference_d2_minus_d1",
        "absolute_error_difference_d2_minus_d1",
        "squared_error_difference_d2_minus_d1",
    ]
    for (gamma, nval, arm), group in paired.groupby(
        ["gamma", "validation_size", "arm"], sort=True
    ):
        row = {
            "gamma": float(gamma), "validation_size": int(nval), "arm": str(arm),
            "paired_replicates": int(len(group)),
            "degree1_bias": float(group.median_error_d1.mean()),
            "degree2_bias": float(group.median_error_d2.mean()),
            "degree1_rmse": float(np.sqrt(np.mean(group.median_error_d1.to_numpy() ** 2))),
            "degree2_rmse": float(np.sqrt(np.mean(group.median_error_d2.to_numpy() ** 2))),
            "degree1_68pct_coverage": float(group.coverage68_d1.mean()),
            "degree2_68pct_coverage": float(group.coverage68_d2.mean()),
        }
        for metric in metrics:
            values = group[metric].to_numpy(float)
            mean = float(values.mean())
            se = float(values.std(ddof=1) / np.sqrt(len(values))) if len(values) > 1 else np.nan
            critical = float(t.ppf(0.975, len(values) - 1)) if len(values) > 1 else np.nan
            row[metric + "_mean"] = mean
            row[metric + "_mcse"] = se
            row[metric + "_ci95_low"] = mean - critical * se if len(values) > 1 else np.nan
            row[metric + "_ci95_high"] = mean + critical * se if len(values) > 1 else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _verify_child_pair(d1: Path, d2: Path) -> dict:
    checked = 0
    for rel in INPUT_FILES:
        a, b = pd.read_csv(d1 / rel), pd.read_csv(d2 / rel)
        pd.testing.assert_frame_equal(a, b, check_dtype=False, check_exact=True)
        checked += 1
    for f in sorted((d1 / "observable").glob("validation_n*.csv")):
        rel = f.relative_to(d1)
        pd.testing.assert_frame_equal(
            pd.read_csv(f), pd.read_csv(d2 / rel), check_dtype=False, check_exact=True
        )
        checked += 1
    with np.load(d1 / "results" / "log_likelihoods.npz") as a, np.load(
        d2 / "results" / "log_likelihoods.npz"
    ) as b:
        if not np.array_equal(a["H0_grid"], b["H0_grid"]):
            raise AssertionError("H0 grids differ across inclusion-model degrees")
        maximum_common_arm_delta = 0.0
        for nval in sorted(int(x.split("_n")[-1]) for x in a.files
                           if x.startswith("A_catalogue_only_n")):
            for arm in COMMON_ARMS:
                key = f"{arm}_n{nval}"
                delta = float(np.max(np.abs(a[key] - b[key])))
                maximum_common_arm_delta = max(maximum_common_arm_delta, delta)
                if not np.allclose(a[key], b[key], rtol=1e-12, atol=1e-12):
                    raise AssertionError(
                        f"common arm {key} changed beyond floating-point tolerance; "
                        f"max abs delta={delta}"
                    )
    return {
        "input_files_compared": checked,
        "common_arm_likelihoods_within_tolerance": True,
        "maximum_common_arm_log_likelihood_difference": maximum_common_arm_delta,
    }


def run(args: argparse.Namespace) -> dict:
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("output directory must be empty")
    runner = Path(__file__).with_name("run_quadratic_five_arm_replicates.py")
    auditor = Path(__file__).with_name("audit_quadratic_five_arm_replicates.py")
    child_roots = {}
    runs = []
    for degree in (1, 2):
        child = out / f"degree{degree}"
        child.mkdir(parents=True, exist_ok=True)
        command = [
            sys.executable, str(runner),
            "--master-seed", str(args.master_seed),
            "--replicates", str(args.replicates),
            "--gammas=" + ",".join(map(str, args.gammas)),
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
            "--output-dir", str(child),
        ]
        proc = subprocess.run(command, capture_output=True, text=True)
        (out / f"degree{degree}.stdout.txt").write_text(proc.stdout)
        (out / f"degree{degree}.stderr.txt").write_text(proc.stderr)
        if not child.exists():
            raise RuntimeError(f"degree-{degree} runner did not create its output directory")
        audit = subprocess.run(
            [sys.executable, str(auditor), str(child), "--output", str(child / "audit_report.json")],
            capture_output=True, text=True,
        )
        (out / f"degree{degree}.audit.stdout.txt").write_text(audit.stdout)
        (out / f"degree{degree}.audit.stderr.txt").write_text(audit.stderr)
        runner_result = None
        if proc.stdout.strip():
            try:
                runner_result = json.loads(proc.stdout)
            except json.JSONDecodeError:
                runner_result = None
        audit_report = None
        if (child / "audit_report.json").is_file():
            audit_report = json.loads((child / "audit_report.json").read_text())
        config_path = child / "config.json"
        child_config = json.loads(config_path.read_text()) if config_path.exists() else {}
        child_status = "PASS" if (
            proc.returncode == 0 and runner_result
            and runner_result.get("status") == "PASS"
            and audit.returncode == 0
            and audit_report and audit_report.get("status") == "PASS"
        ) else "INCOMPLETE"
        runs.append({
            "degree": degree, "runner_exit_code": proc.returncode,
            "auditor_exit_code": audit.returncode, "status": child_status,
            "audit_report": audit_report, "config": child_config,
        })
        child_roots[degree] = child

    paired_checks = []
    complete = all(x["status"] == "PASS" for x in runs)
    paired_frames = []
    aggregate_frames = []
    if complete:
        for key in ("master_seed", "replicates", "gammas", "n_parent", "n_detected",
                    "validation_frame", "validation_sizes", "curvature_beta",
                    "h0_true", "h0_prior_support", "h0_points", "scott_multiplier",
                    "quadrature_order", "sigma_mpc", "detection_threshold_mpc"):
            if runs[0]["config"].get(key) != runs[1]["config"].get(key):
                raise AssertionError(f"degree runs differ in frozen setting {key}")
        for rep in range(args.replicates):
            for gamma in args.gammas:
                label = f"rep{rep:03d}_g{float(gamma):+g}".replace("+", "p").replace("-", "m")
                left = child_roots[1] / "replicates" / label
                right = child_roots[2] / "replicates" / label
                paired_checks.append({
                    "replicate": rep, "gamma": gamma,
                    **_verify_child_pair(left, right),
                })
        left_summary = pd.read_csv(child_roots[1] / "replicate_posterior_summaries.csv")
        right_summary = pd.read_csv(child_roots[2] / "replicate_posterior_summaries.csv")
        paired = _paired_summary(left_summary, right_summary)
        paired.to_csv(out / "paired_degree_contrasts.csv", index=False)
        aggregate = _aggregate_paired(paired)
        aggregate.to_csv(out / "paired_degree_aggregate_metrics.csv", index=False)
        complete = len(paired) == args.replicates * len(args.gammas) * len(args.validation_sizes) * len(ARMS)
        if not complete:
            raise AssertionError("paired summary rows are incomplete")
    pd.DataFrame(paired_checks).to_csv(out / "paired_input_and_likelihood_checks.csv", index=False)

    config = {
        "experiment": "paired_degree1_degree2_quadratic_logit_multireplicate_pilot",
        "status": "PASS" if complete else "INCOMPLETE",
        "claim_limit": "Finite paired synthetic pilot only; no estimator superiority, calibration, or observational inference.",
        "degree_prior_parameterization": (
            "Degree one uses [1,x] with independent Gaussian priors on (alpha,beta); "
            "degree two uses the centered basis [1,x,x^2-c_g] with the same marginal "
            "prior on (alpha,beta) and an independent Gaussian prior on gamma."
        ),
        "master_seed": args.master_seed, "replicates": args.replicates,
        "gammas": args.gammas, "n_parent": args.n_parent, "n_detected": args.n_detected,
        "validation_frame": args.validation_frame, "validation_sizes": args.validation_sizes,
        "curvature_beta": args.curvature_beta, "h0_true": args.h0_true,
        "h0_prior_support": [args.h0_prior_min, args.h0_prior_max],
        "h0_points": args.h0_points, "scott_multiplier": args.scott_multiplier,
        "quadrature_order": args.quadrature_order, "sigma_mpc": args.sigma_mpc,
        "detection_threshold_mpc": args.detection_threshold_mpc,
        "degree_runs": runs, "paired_child_checks": len(paired_checks),
    }
    (out / "config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    _write_manifest(out, config["status"])
    return {
        "status": config["status"], "degree_runs": [
            {k: v for k, v in x.items() if k != "config"} for x in runs
        ], "paired_child_checks": len(paired_checks),
        "paired_rows": sum(len(pd.read_csv(p)) for p in [
            out / "paired_degree_contrasts.csv"
        ] if p.exists()),
        "output_dir": str(out),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--master-seed", type=int, default=20261120)
    p.add_argument("--replicates", type=int, default=4)
    p.add_argument("--gammas", type=lambda x: [float(v) for v in x.split(",")],
                   default=[-1.0, 0.0, 1.0])
    p.add_argument("--n-parent", type=int, default=1200)
    p.add_argument("--n-detected", type=int, default=30)
    p.add_argument("--validation-frame", type=int, default=500)
    p.add_argument("--validation-sizes", type=lambda x: [int(v) for v in x.split(",")],
                   default=[100, 300, 500])
    p.add_argument("--curvature-beta", type=float, default=-3.9)
    p.add_argument("--h0-true", type=float, default=70.0)
    p.add_argument("--h0-prior-min", type=float, default=30.0)
    p.add_argument("--h0-prior-max", type=float, default=110.0)
    p.add_argument("--h0-points", type=int, default=281)
    p.add_argument("--scott-multiplier", type=float, default=0.7)
    p.add_argument("--quadrature-order", type=int, default=5)
    p.add_argument("--sigma-mpc", type=float, default=140.0)
    p.add_argument("--detection-threshold-mpc", type=float, default=250.0)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()
    result = run(args)
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
