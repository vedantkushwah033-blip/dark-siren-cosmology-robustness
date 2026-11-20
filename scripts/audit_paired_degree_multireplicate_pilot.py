"""Independent audit for the paired degree-one/degree-two pilot.

Does not import the paired runner. Verifies the outer and child manifests,
child audit attestations, frozen configuration, paired raw inputs, common-arm
likelihood invariance, and independently reconstructs paired summary contrasts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


INPUT_FILES = [
    "observable/catalogue_observed.csv",
    "observable/gw_events_observed.csv",
    "truth/parent_truth.csv",
    "truth/gw_event_truth.csv",
    "truth/validation_frame_truth.csv",
]
COMMON_ARMS = [
    "A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
]


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _verify_manifest(root: Path, manifest_name: str = "manifest.json") -> int:
    manifest = json.loads((root / manifest_name).read_text())
    if manifest.get("status") != "PASS":
        raise AssertionError(f"{root}: manifest status is not PASS")
    count = 0
    for rel, meta in manifest["files"].items():
        path = (root / rel).resolve()
        if root.resolve() not in path.parents:
            raise AssertionError(f"manifest path escapes artifact root: {rel}")
        if not path.is_file():
            raise AssertionError(f"manifest file missing: {rel}")
        if path.stat().st_size != int(meta["bytes"]) or _sha(path) != meta["sha256"]:
            raise AssertionError(f"manifest hash/size mismatch: {rel}")
        count += 1
    if count != int(manifest["file_count"]):
        raise AssertionError("manifest file_count mismatch")
    return count


def _compare_frames(left: Path, right: Path, rel: str) -> int:
    pd.testing.assert_frame_equal(
        pd.read_csv(left / rel), pd.read_csv(right / rel),
        check_dtype=False, check_exact=True,
    )
    return 1


def _verify_pair(left: Path, right: Path) -> dict:
    checks = 0
    for rel in INPUT_FILES:
        checks += _compare_frames(left, right, rel)
    for f in sorted((left / "observable").glob("validation_n*.csv")):
        checks += _compare_frames(left, right, str(f.relative_to(left)))
    with np.load(left / "results" / "log_likelihoods.npz") as a, np.load(
        right / "results" / "log_likelihoods.npz"
    ) as b:
        if not np.array_equal(a["H0_grid"], b["H0_grid"]):
            raise AssertionError("H0 grid mismatch across degrees")
        keys = [k for k in a.files if any(k.startswith(f"{arm}_n") for arm in COMMON_ARMS)]
        if not keys:
            raise AssertionError("no common-arm likelihood arrays found")
        maximum_common_arm_delta = 0.0
        for key in keys:
            if key not in b.files:
                raise AssertionError(f"common-arm likelihood missing in degree two: {key}")
            delta = float(np.max(np.abs(a[key] - b[key])))
            maximum_common_arm_delta = max(maximum_common_arm_delta, delta)
            if not np.allclose(a[key], b[key], rtol=1e-12, atol=1e-12):
                raise AssertionError(
                    f"common-arm likelihood mismatch beyond tolerance: {key}; max abs delta={delta}"
                )
    return {
        # Match the orchestrator record schema; values are independently recomputed.
        "input_files_compared": checks,
        "common_arm_likelihoods_within_tolerance": True,
        "maximum_common_arm_log_likelihood_difference": maximum_common_arm_delta,
    }


def _reconstruct_paired(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
    keys = ["replicate", "gamma", "validation_size", "arm"]
    cols = keys + ["median_error", "q16", "q84", "h0_true"]
    joined = left[cols].merge(right[cols], on=keys, suffixes=("_d1", "_d2"),
                              validate="one_to_one")
    if len(joined) != len(left) or len(joined) != len(right):
        raise AssertionError("unpaired degree summary rows")
    if not np.array_equal(joined.h0_true_d1, joined.h0_true_d2):
        raise AssertionError("degree truths differ")
    joined["median_error_difference_d2_minus_d1"] = joined.median_error_d2 - joined.median_error_d1
    joined["absolute_error_difference_d2_minus_d1"] = (
        joined.median_error_d2.abs() - joined.median_error_d1.abs()
    )
    joined["squared_error_difference_d2_minus_d1"] = (
        joined.median_error_d2.pow(2) - joined.median_error_d1.pow(2)
    )
    joined["coverage68_d1"] = (
        (joined.q16_d1 <= joined.h0_true_d1) & (joined.h0_true_d1 <= joined.q84_d1)
    )
    joined["coverage68_d2"] = (
        (joined.q16_d2 <= joined.h0_true_d2) & (joined.h0_true_d2 <= joined.q84_d2)
    )
    return joined


def _reconstruct_aggregate(paired: pd.DataFrame) -> pd.DataFrame:
    """Independently rebuild the runner's grouped paired-effect summaries."""
    from scipy.stats import t

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


def _verify_aggregate_metrics(paired: pd.DataFrame, observed: pd.DataFrame) -> int:
    """Fail closed if the stored aggregate table differs from paired rows."""
    expected = _reconstruct_aggregate(paired)
    sort_keys = ["gamma", "validation_size", "arm"]
    if not set(sort_keys).issubset(observed.columns):
        raise AssertionError("paired aggregate table is missing grouping keys")
    pd.testing.assert_frame_equal(
        observed.sort_values(sort_keys).reset_index(drop=True),
        expected.sort_values(sort_keys).reset_index(drop=True),
        check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12,
    )
    return len(expected)


def audit(root: Path) -> dict:
    root = root.resolve()
    cfg = json.loads((root / "config.json").read_text())
    if cfg.get("status") != "PASS":
        raise AssertionError("paired experiment is incomplete; cannot audit as PASS")
    outer_count = _verify_manifest(root)
    runs = cfg["degree_runs"]
    if len(runs) != 2 or [int(x["degree"]) for x in runs] != [1, 2]:
        raise AssertionError("expected exactly degree-one and degree-two runs")
    for run in runs:
        if run.get("status") != "PASS" or not run.get("audit_report"):
            raise AssertionError("a degree run did not pass lifecycle and independent audit")
        if run["audit_report"].get("status") != "PASS":
            raise AssertionError("nested degree audit did not PASS")
    left_root, right_root = root / "degree1", root / "degree2"
    child_manifest_count = 0
    for degree_root in (left_root, right_root):
        child_manifest_count += _verify_manifest(degree_root)
    for key in ("master_seed", "replicates", "gammas", "n_parent", "n_detected",
                "validation_frame", "validation_sizes", "h0_true", "h0_prior_support",
                "h0_points", "quadrature_order"):
        if runs[0]["config"].get(key) != runs[1]["config"].get(key):
            raise AssertionError(f"degree config mismatch: {key}")

    paired_checks = []
    for rep in range(int(cfg["replicates"])):
        for gamma in cfg["gammas"]:
            label = f"rep{rep:03d}_g{float(gamma):+g}".replace("+", "p").replace("-", "m")
            d1 = left_root / "replicates" / label
            d2 = right_root / "replicates" / label
            for child in (d1, d2):
                report = json.loads((child / "audit_report.json").read_text())
                if report.get("status") != "PASS":
                    raise AssertionError(f"child audit failed: {child}")
                child_manifest_count += _verify_manifest(child)
            paired_checks.append({"replicate": rep, "gamma": gamma, **_verify_pair(d1, d2)})

    recorded_checks = pd.read_csv(root / "paired_input_and_likelihood_checks.csv")
    recomputed_checks = pd.DataFrame(paired_checks)
    pd.testing.assert_frame_equal(
        recorded_checks.sort_values(["replicate", "gamma"]).reset_index(drop=True),
        recomputed_checks.sort_values(["replicate", "gamma"]).reset_index(drop=True),
        check_dtype=False, check_exact=True,
    )
    left = pd.read_csv(left_root / "replicate_posterior_summaries.csv")
    right = pd.read_csv(right_root / "replicate_posterior_summaries.csv")
    expected = _reconstruct_paired(left, right)
    observed = pd.read_csv(root / "paired_degree_contrasts.csv")
    pd.testing.assert_frame_equal(
        observed.sort_values(["replicate", "gamma", "validation_size", "arm"]).reset_index(drop=True),
        expected.sort_values(["replicate", "gamma", "validation_size", "arm"]).reset_index(drop=True),
        check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12,
    )
    observed_aggregate = pd.read_csv(root / "paired_degree_aggregate_metrics.csv")
    aggregate_rows = _verify_aggregate_metrics(expected, observed_aggregate)
    return {
        "status": "PASS", "outer_manifest_entries_verified": outer_count,
        "nested_manifest_entries_verified": child_manifest_count,
        "paired_child_blocks_verified": len(paired_checks),
        "paired_summary_rows_reconstructed": len(expected),\n        "paired_aggregate_rows_reconstructed": aggregate_rows,
        "common_arm_likelihoods_invariant": True,
        "interpretation": "Artifact integrity and paired arithmetic only; not calibration evidence.",
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("artifact_dir", type=Path)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    result = audit(args.artifact_dir)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered)


if __name__ == "__main__":
    main()
