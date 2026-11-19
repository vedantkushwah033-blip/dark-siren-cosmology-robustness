"""Pool two independently seeded 600-cell H0 calibration replications.

Verifies both archived manifests, orchestration audits, planned-cell ledgers,
seed independence, and pairing before computing pooled and per-replication
fixed-truth summaries. Non-successes remain in all-attempt denominators.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t


KEYS = ["cell_id", "replicate", "gamma", "n_validation",
        "analysis_seed", "validation_seed", "h0_true"]
STATUSES = {"SUCCESS", "NON_EVALUABLE", "FAILED", "NOT_RUN", "RUNNING"}


def _verify_manifest(root: Path) -> int:
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise ValueError(f"missing aggregate manifest: {root}")
    manifest = json.loads(manifest_path.read_text())
    bad = []
    for name, expected in manifest.items():
        path = root / name
        actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        if actual != expected:
            bad.append(name)
    if bad:
        raise AssertionError(f"manifest mismatch in {root}: {bad[:10]}")
    return len(manifest)


def _load_replication(root: Path, label: str, expected_cells: int) -> tuple[pd.DataFrame, dict]:
    root = root.resolve()
    manifest_entries = _verify_manifest(root)
    audit_path = root / "independent_orchestration_audit.json"
    if not audit_path.is_file():
        raise ValueError(f"missing independent orchestration audit in {root}")
    audit = json.loads(audit_path.read_text())
    if audit.get("status") != "PASS" or int(audit.get("planned_cells", -1)) != expected_cells:
        raise AssertionError(f"{label} orchestration audit is not a complete PASS: {audit}")
    planned = pd.read_csv(root / "planned_cells.csv")
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    if len(planned) != expected_cells or len(ledger) != expected_cells:
        raise AssertionError(f"{label} planned/ledger row count differs from {expected_cells}")
    for name, frame in (("planned", planned), ("ledger", ledger)):
        missing = set(KEYS) - set(frame.columns)
        if missing:
            raise AssertionError(f"{label} {name} missing columns: {sorted(missing)}")
        if frame.cell_id.duplicated().any() or frame.duplicated(
            ["replicate", "gamma", "n_validation"]
        ).any():
            raise AssertionError(f"{label} {name} has duplicate cell keys")
    left = planned[KEYS].sort_values(KEYS[:4]).reset_index(drop=True)
    right = ledger[KEYS].sort_values(KEYS[:4]).reset_index(drop=True)
    pd.testing.assert_frame_equal(left, right, check_dtype=False, check_exact=True)
    if not set(ledger.status).issubset(STATUSES):
        raise AssertionError(f"{label} has unknown execution statuses")
    if ledger.status.isin(["RUNNING", "NOT_RUN"]).any():
        raise AssertionError(f"{label} has unfinished planned cells")
    observed_counts = {str(k): int(v) for k,v in ledger.status.value_counts().items()}
    audit_counts = {str(k): int(v) for k,v in audit.get("status_counts", {}).items()}
    if observed_counts != audit_counts:
        raise AssertionError(f"{label} ledger/audit status counts disagree")
    if not set(ledger.gamma.astype(float)) == {-1.0, 0.0, 1.0}:
        raise AssertionError(f"{label} gamma factorial is incomplete")
    if not set(ledger.n_validation.astype(int)) == {100, 500}:
        raise AssertionError(f"{label} validation-size factorial is incomplete")
    if set(ledger.replicate.astype(int)) != set(range(100)):
        raise AssertionError(f"{label} replicate IDs are not 0..99")
    if ledger.groupby("replicate").analysis_seed.nunique().ne(1).any():
        raise AssertionError(f"{label} analysis seed is not shared within replicate")
    if ledger.groupby("replicate").validation_seed.nunique().ne(1).any():
        raise AssertionError(f"{label} validation seed is not shared within replicate")
    if ledger.groupby(["gamma", "n_validation"]).size().ne(100).any():
        raise AssertionError(f"{label} has an incomplete gamma × validation-size cell")
    first = ledger[ledger.status == "SUCCESS"].sort_values(
        ["replicate", "gamma", "n_validation"]
    ).iloc[0]
    config_path = root / "replicates" / str(first.cell_id) / "config.json"
    if not config_path.is_file():
        raise AssertionError(f"{label} representative child config is missing")
    cfg = json.loads(config_path.read_text())
    analysis = cfg["analysis_generator_config"]
    validation = cfg["validation_generator_config"]
    inference = cfg["inference"]
    analysis_keys = (
        "experiment", "n_parent", "n_detected_events", "h0_true_km_s_mpc",
        "sigma_mpc", "detection_threshold_mpc", "redshift_domain",
        "catalogue_selection", "merger_weight", "distance_model",
    )
    validation_keys = (
        "independent_frame", "n_validation_frame", "n_validation",
        "sampling_design", "catalogue_selection",
    )
    signature = {
        "analysis": {key: analysis[key] for key in analysis_keys},
        "validation": {key: validation[key] for key in validation_keys},
        "inference": inference,
    }
    if not np.all(ledger.h0_true.to_numpy(float) == 70.0):
        raise AssertionError(f"{label} has an unexpected H0 truth")
    ledger = ledger.copy()
    ledger["source_run"] = label
    ledger["global_replicate"] = ledger.replicate.astype(int) + (0 if label == "replication_1" else 100)
    ledger["manifest_entries_verified"] = manifest_entries
    return ledger, audit, signature


def _wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n <= 0:
        return float("nan"), float("nan")
    p = k/n
    den = 1+z*z/n
    center = (p+z*z/(2*n))/den
    half = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return float(max(0, center-half)), float(min(1, center+half))


def _mean_ci(x: np.ndarray) -> tuple[float, float, float, float]:
    v = np.asarray(x, float)
    v = v[np.isfinite(v)]
    n = len(v)
    if n == 0:
        return float("nan"), float("nan"), float("nan"), float("nan")
    mean = float(v.mean())
    sd = float(v.std(ddof=1)) if n > 1 else float("nan")
    mcse = sd/np.sqrt(n) if n > 1 else float("nan")
    critical = float(t.ppf(.975, n-1)) if n > 1 else float("nan")
    return mean, mcse, mean-critical*mcse, mean+critical*mcse


def _paired_contrast(left: pd.DataFrame, right: pd.DataFrame,
                     left_name: str, right_name: str, label: str) -> list[dict]:
    merged = left.merge(right, on=["source_run", "replicate"], suffixes=(f"_{left_name}", f"_{right_name}"),
                        validate="one_to_one")
    rows = []
    if merged.empty:
        return rows
    for key in ("analysis_seed", "validation_seed", "h0_true"):
        if not np.array_equal(merged[f"{key}_{left_name}"], merged[f"{key}_{right_name}"]):
            raise AssertionError(f"{label} paired {key} mismatch")
    metrics = {
        "signed_median_bias": merged[f"median_bias_{left_name}"]-merged[f"median_bias_{right_name}"],
        "absolute_median_error": (
            np.abs(merged[f"median_bias_{left_name}"]) -
            np.abs(merged[f"median_bias_{right_name}"])
        ),
        "width68": (
            (merged[f"h0_q84_{left_name}"]-merged[f"h0_q16_{left_name}"]) -
            (merged[f"h0_q84_{right_name}"]-merged[f"h0_q16_{right_name}"])
        ),
    }
    for metric, values in metrics.items():
        x = np.asarray(values, float)
        x = x[np.isfinite(x)]
        mean, mcse, low, high = _mean_ci(x)
        rows.append({
            "contrast": label, "metric": metric, "paired_replicates": len(x),
            "mean_difference": mean, "mcse": mcse, "t95_low": low, "t95_high": high,
        })
    return rows


def analyze(replication_1: Path, replication_2: Path, output: Path,
            expected_cells: int = 600) -> dict:
    a, audit_a, signature_a = _load_replication(replication_1, "replication_1", expected_cells)
    b, audit_b, signature_b = _load_replication(replication_2, "replication_2", expected_cells)
    if signature_a != signature_b:
        raise AssertionError("replication scientific configurations differ; pooling is invalid")
    seeds_a = set(a.analysis_seed.astype(int)) | set(a.validation_seed.astype(int))
    seeds_b = set(b.analysis_seed.astype(int)) | set(b.validation_seed.astype(int))
    if seeds_a.intersection(seeds_b):
        raise AssertionError("replication seed streams overlap")
    ledger = pd.concat([a, b], ignore_index=True)
    if len(ledger) != 2*expected_cells:
        raise AssertionError("pooled ledger has wrong planned denominator")
    success = ledger[ledger.status == "SUCCESS"].copy()
    numeric = ["h0_true", "h0_median", "h0_q16", "h0_q84", "median_bias",
               "covered68", "covered90", "covered95",
               "min_catalogue_kish_ess", "max_normalized_weight"]
    if not np.isfinite(success[numeric].to_numpy(float)).all():
        raise AssertionError("SUCCESS rows contain non-finite inference metrics")

    summaries = []
    replication_summaries = []
    for group_keys, group in ledger.groupby(["gamma", "n_validation"], sort=True):
        gamma, nval = group_keys
        ok = group[group.status == "SUCCESS"]
        planned, n = len(group), len(ok)
        nonsuccess = planned-n
        bias = ok.median_bias.to_numpy(float)
        mean, mcse, low, high = _mean_ci(bias)
        row = {
            "gamma": float(gamma), "n_validation": int(nval),
            "planned": planned, "success": n, "non_success": nonsuccess,
            "non_evaluable": int((group.status == "NON_EVALUABLE").sum()),
            "failed": int((group.status == "FAILED").sum()),
            "mean_median_bias": mean, "bias_mcse": mcse,
            "bias_t95_low": low, "bias_t95_high": high,
            "rmse": float(np.sqrt(np.mean(bias*bias))) if n else float("nan"),
            "mean_width68": float((ok.h0_q84-ok.h0_q16).mean()) if n else float("nan"),
            "median_min_kish_ess": float(ok.min_catalogue_kish_ess.median()) if n else float("nan"),
            "median_max_normalized_weight": float(ok.max_normalized_weight.median()) if n else float("nan"),
        }
        for level in (68, 90, 95):
            k = int(ok[f"covered{level}"].sum()) if n else 0
            wilson_low, wilson_high = _wilson(k, n)
            row.update({
                f"coverage{level}": k/n if n else float("nan"),
                f"coverage{level}_wilson95_low": wilson_low,
                f"coverage{level}_wilson95_high": wilson_high,
                f"coverage{level}_all_attempt_lower": k/planned,
                f"coverage{level}_all_attempt_upper": (k+nonsuccess)/planned,
            })
        summaries.append(row)
        for source, sub in group.groupby("source_run"):
            good = sub[sub.status == "SUCCESS"]
            replication_summaries.append({
                "source_run": source, "gamma": float(gamma),
                "n_validation": int(nval), "planned": len(sub),
                "success": len(good), "non_success": len(sub)-len(good),
                "mean_median_bias": float(good.median_bias.mean()) if len(good) else float("nan"),
                "rmse": float(np.sqrt(np.mean(good.median_bias.to_numpy(float)**2))) if len(good) else float("nan"),
                "coverage68": float(good.covered68.mean()) if len(good) else float("nan"),
                "coverage90": float(good.covered90.mean()) if len(good) else float("nan"),
                "coverage95": float(good.covered95.mean()) if len(good) else float("nan"),
            })
    pooled = pd.DataFrame(summaries)
    per_rep = pd.DataFrame(replication_summaries)

    contrasts = []
    for nval in (100, 500):
        subset = success[success.n_validation == nval]
        left = subset[subset.gamma == 1.0]
        right = subset[subset.gamma == -1.0]
        contrasts += _paired_contrast(left, right, "plus", "minus",
                                      f"gamma +1 minus -1 at n={nval}")
    for gamma in (-1.0, 0.0, 1.0):
        subset = success[success.gamma == gamma]
        left = subset[subset.n_validation == 500]
        right = subset[subset.n_validation == 100]
        contrasts += _paired_contrast(left, right, "large", "small",
                                      f"validation 500 minus 100 at gamma={gamma:g}")
    contrast_df = pd.DataFrame(contrasts)

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    pooled.to_csv(output/"pooled_cell_summary.csv", index=False)
    per_rep.to_csv(output/"replication_cell_summary.csv", index=False)
    contrast_df.to_csv(output/"pooled_paired_contrasts.csv", index=False)
    report = {
        "status": "PASS",
        "planned_cells_per_replication": expected_cells,
        "pooled_planned_cells": len(ledger),
        "pooled_success": int((ledger.status == "SUCCESS").sum()),
        "pooled_non_success": int((ledger.status != "SUCCESS").sum()),
        "status_counts": {str(k): int(v) for k,v in ledger.status.value_counts().items()},
        "replication_1_manifest_entries_verified": int(a.manifest_entries_verified.iloc[0]),
        "replication_2_manifest_entries_verified": int(b.manifest_entries_verified.iloc[0]),
        "replication_1_audit": audit_a,
        "replication_2_audit": audit_b,
        "seed_streams_disjoint": True,
        "cell_summaries": len(pooled),
        "paired_contrast_rows": len(contrast_df),
        "interpretation": "Pooled fixed-truth repeated-sampling summaries across two independent seeds. Success-only metrics are conditional on successful posterior production; all-attempt bounds retain non-successes. Not joint SBC or observational H0 inference.",
    }
    (output/"pooling_audit.json").write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("replication_1",type=Path)
    p.add_argument("replication_2",type=Path)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--expected-cells",type=int,default=600)
    a=p.parse_args()
    print(json.dumps(analyze(a.replication_1,a.replication_2,a.output_dir,a.expected_cells),indent=2))


if __name__=="__main__":
    main()
