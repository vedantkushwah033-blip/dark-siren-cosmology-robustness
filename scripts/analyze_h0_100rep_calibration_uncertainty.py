"""Failure-inclusive post-hoc uncertainty analysis for the 100-replicate H0 pilot.

Validates the archived plan, full ledger, and aggregate SHA-256 manifest before
computing success-conditional calibration metrics and all-attempt coverage
bounds. This is descriptive fixed-truth repeated-sampling analysis, not SBC.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t


IDENTITY = ["replicate", "gamma", "n_validation", "analysis_seed",
            "validation_seed", "h0_true"]
ALLOWED = {"SUCCESS", "NON_EVALUABLE", "FAILED", "RUNNING", "NOT_RUN"}


def verify_manifest(root: Path) -> int:
    path = root / "manifest.json"
    if not path.is_file():
        raise ValueError("aggregate manifest.json is missing")
    manifest = json.loads(path.read_text())
    bad = []
    for name, digest in manifest.items():
        target = root / name
        if not target.is_file() or hashlib.sha256(target.read_bytes()).hexdigest() != digest:
            bad.append(name)
    if bad:
        raise AssertionError(f"aggregate manifest mismatch: {bad[:10]}")
    return len(manifest)


def verify_plan(root: Path, ledger: pd.DataFrame) -> int:
    path = root / "planned_cells.csv"
    if not path.is_file():
        raise ValueError("planned_cells.csv is missing")
    planned = pd.read_csv(path)
    keys = ["cell_id", *IDENTITY]
    for label, frame in (("planned", planned), ("ledger", ledger)):
        missing = set(keys) - set(frame.columns)
        if missing:
            raise ValueError(f"{label} ledger missing columns: {sorted(missing)}")
        if frame.cell_id.duplicated().any():
            raise ValueError(f"{label} ledger has duplicate cell_id")
        if frame.duplicated(["replicate", "gamma", "n_validation"]).any():
            raise ValueError(f"{label} ledger has duplicate replicate/gamma/size")
    left = planned[keys].sort_values(keys[:4]).reset_index(drop=True)
    right = ledger[keys].sort_values(keys[:4]).reset_index(drop=True)
    try:
        pd.testing.assert_frame_equal(left, right, check_dtype=False, check_exact=True)
    except AssertionError as exc:
        raise AssertionError("planned cell inventory, seeds, or H0 truths do not match result ledger") from exc
    return len(planned)


def mean_t_interval(values: np.ndarray) -> tuple[float, float, float]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if not len(x):
        return float("nan"), float("nan"), float("nan")
    mean = float(x.mean())
    if len(x) < 2:
        return mean, float("nan"), float("nan")
    half = float(t.ppf(.975, len(x)-1) * x.std(ddof=1) / np.sqrt(len(x)))
    return mean, mean-half, mean+half


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    p = k/n
    den = 1 + z*z/n
    center = (p + z*z/(2*n))/den
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return float(max(0, center-half)), float(min(1, center+half))


def bootstrap_mean_interval(values: np.ndarray, rng: np.random.Generator,
                            draws: int = 10000) -> tuple[float, float]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return float("nan"), float("nan")
    indices = rng.integers(0, len(x), size=(draws, len(x)))
    means = x[indices].mean(axis=1)
    return tuple(map(float, np.quantile(means, [.025, .975])))


def _assert_pairing_keys(paired: pd.DataFrame, left: str, right: str) -> None:
    for key in ("analysis_seed", "validation_seed", "h0_true"):
        a = paired[f"{key}_{left}"].to_numpy()
        b = paired[f"{key}_{right}"].to_numpy()
        if not np.array_equal(a, b):
            raise ValueError(f"paired replicate mismatch in {key}")


def _paired_rows(left: pd.DataFrame, right: pd.DataFrame,
                 left_suffix: str, right_suffix: str) -> pd.DataFrame:
    merged = left.merge(right, on="replicate", suffixes=(f"_{left_suffix}", f"_{right_suffix}"),
                        validate="one_to_one")
    if not merged.empty:
        _assert_pairing_keys(merged, left_suffix, right_suffix)
    return merged


def analyze(root: Path, output: Path, seed: int = 20261119) -> dict:
    root = root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest_entries = verify_manifest(root)
    ledger_path = root / "replicate_ledger.csv"
    if not ledger_path.is_file():
        raise ValueError("replicate_ledger.csv is missing")
    ledger = pd.read_csv(ledger_path)
    required = {"cell_id", *IDENTITY, "status", "h0_median", "h0_q16", "h0_q84",
                "covered68", "covered90", "covered95",
                "min_catalogue_kish_ess", "max_normalized_weight"}
    missing = required - set(ledger.columns)
    if missing:
        raise ValueError(f"ledger missing columns: {sorted(missing)}")
    planned_count = verify_plan(root, ledger)
    if not set(ledger.status).issubset(ALLOWED):
        raise ValueError("ledger contains an unknown execution status")
    if ledger.empty:
        raise ValueError("replicate ledger is empty")
    success = ledger[ledger.status == "SUCCESS"].copy()
    numeric = ["h0_true", "h0_median", "h0_q16", "h0_q84",
               "covered68", "covered90", "covered95",
               "min_catalogue_kish_ess", "max_normalized_weight"]
    if not np.isfinite(success[numeric].to_numpy(float)).all():
        raise ValueError("SUCCESS rows contain non-finite inference metrics")
    if (success.h0_q16 > success.h0_median).any() or (success.h0_median > success.h0_q84).any():
        raise ValueError("SUCCESS posterior quantiles are not ordered")
    if not success.empty and not success[["covered68", "covered90", "covered95"]].isin([0, 1]).all().all():
        raise ValueError("SUCCESS coverage indicators must be binary")

    rng = np.random.default_rng(seed)
    rows = []
    for (gamma, nval), group in ledger.groupby(["gamma", "n_validation"], sort=True):
        ok = group[group.status == "SUCCESS"]
        n_success = len(ok)
        planned = len(group)
        non_success = planned - n_success
        bias = (ok.h0_median - ok.h0_true).to_numpy(float)
        mean_bias, bias_low, bias_high = mean_t_interval(bias)
        rmse = float(np.sqrt(np.mean(bias*bias))) if n_success else float("nan")
        rmse_low, rmse_high = bootstrap_mean_interval(bias*bias, rng)
        if np.isfinite([rmse_low, rmse_high]).all():
            rmse_low, rmse_high = float(np.sqrt(max(0, rmse_low))), float(np.sqrt(max(0, rmse_high)))
        else:
            rmse_low = rmse_high = float("nan")
        row = {
            "gamma": float(gamma), "n_validation": int(nval),
            "planned": planned, "success": n_success, "non_success": non_success,
            "non_evaluable": int((group.status == "NON_EVALUABLE").sum()),
            "failed": int((group.status == "FAILED").sum()),
            "running": int((group.status == "RUNNING").sum()),
            "not_run": int((group.status == "NOT_RUN").sum()),
            "mean_median_bias": mean_bias, "bias_t95_low": bias_low,
            "bias_t95_high": bias_high,
            "median_bias": float(np.median(bias)) if n_success else float("nan"),
            "rmse_posterior_median": rmse,
            "rmse_bootstrap95_low": rmse_low, "rmse_bootstrap95_high": rmse_high,
            "mean_width68": float((ok.h0_q84-ok.h0_q16).mean()) if n_success else float("nan"),
            "median_min_kish_ess": float(ok.min_catalogue_kish_ess.median()) if n_success else float("nan"),
            "median_max_normalized_weight": float(ok.max_normalized_weight.median()) if n_success else float("nan"),
        }
        for level in (68, 90, 95):
            col = f"covered{level}"
            k = int(ok[col].sum()) if n_success else 0
            low, high = wilson(k, n_success)
            row.update({
                f"coverage{level}": k/n_success if n_success else float("nan"),
                f"coverage{level}_wilson95_low": low,
                f"coverage{level}_wilson95_high": high,
                f"coverage{level}_all_attempt_lower": k/planned,
                f"coverage{level}_all_attempt_upper": (k+non_success)/planned,
            })
        rows.append(row)
    cell_summary = pd.DataFrame(rows)
    cell_summary.to_csv(output / "cell_uncertainty_summary.csv", index=False)

    paired_gamma = []
    gammas = sorted(ledger.gamma.unique())
    if 0.0 in gammas:
        for gamma in gammas:
            if float(gamma) == 0.0:
                continue
            for nval in sorted(ledger.n_validation.unique()):
                left = success[(success.gamma == gamma) & (success.n_validation == nval)]
                right = success[(success.gamma == 0.0) & (success.n_validation == nval)]
                paired = _paired_rows(left, right, "gamma", "zero")
                if paired.empty:
                    continue
                metrics = {
                    "signed_median_bias": paired.h0_median_gamma - paired.h0_median_zero,
                    "absolute_median_error": abs(paired.h0_median_gamma-paired.h0_true_gamma)
                                             - abs(paired.h0_median_zero-paired.h0_true_zero),
                    "width68": (paired.h0_q84_gamma-paired.h0_q16_gamma)
                               - (paired.h0_q84_zero-paired.h0_q16_zero),
                }
                for metric, values in metrics.items():
                    x = np.asarray(values, float)
                    avg, low, high = mean_t_interval(x)
                    boot_low, boot_high = bootstrap_mean_interval(x, rng)
                    paired_gamma.append({
                        "contrast": f"gamma_{gamma:g}_minus_gamma_0",
                        "gamma": float(gamma), "n_validation": int(nval),
                        "metric": metric, "paired_replicates": len(x),
                        "mean_difference": avg, "paired_t95_low": low,
                        "paired_t95_high": high, "paired_bootstrap95_low": boot_low,
                        "paired_bootstrap95_high": boot_high,
                    })
    pd.DataFrame(paired_gamma).to_csv(output / "paired_gamma_contrasts.csv", index=False)

    paired_sizes = []
    for gamma, group in success.groupby("gamma"):
        sizes = sorted(group.n_validation.unique())
        if len(sizes) < 2:
            continue
        low = group[group.n_validation == sizes[0]]
        high = group[group.n_validation == sizes[-1]]
        paired = _paired_rows(high, low, "high", "low")
        if paired.empty:
            continue
        metrics = {
            "signed_median_bias": paired.h0_median_high - paired.h0_median_low,
            "absolute_median_error": abs(paired.h0_median_high-paired.h0_true_high)
                                     - abs(paired.h0_median_low-paired.h0_true_low),
            "width68": (paired.h0_q84_high-paired.h0_q16_high)
                       - (paired.h0_q84_low-paired.h0_q16_low),
        }
        for metric, values in metrics.items():
            x = np.asarray(values, float)
            avg, low_ci, high_ci = mean_t_interval(x)
            boot_low, boot_high = bootstrap_mean_interval(x, rng)
            paired_sizes.append({
                "gamma": float(gamma), "comparison": f"{sizes[-1]}_minus_{sizes[0]}",
                "metric": metric, "paired_replicates": len(x),
                "mean_difference": avg, "paired_t95_low": low_ci,
                "paired_t95_high": high_ci, "paired_bootstrap95_low": boot_low,
                "paired_bootstrap95_high": boot_high,
            })
    pd.DataFrame(paired_sizes).to_csv(output / "paired_validation_size_contrasts.csv", index=False)

    status_counts = {str(k): int(v) for k, v in ledger.status.value_counts().items()}
    complete = len(ledger) == planned_count and not (ledger.status.isin(["RUNNING", "NOT_RUN"])).any()
    report = {
        "status": "PASS" if complete else "INCOMPLETE",
        "planned_cells": planned_count, "ledger_rows": len(ledger),
        "status_counts": status_counts, "manifest_entries_verified": manifest_entries,
        "cell_summaries": len(cell_summary),
        "paired_gamma_contrasts": len(paired_gamma),
        "paired_validation_size_contrasts": len(paired_sizes),
        "interpretation": "Success-only calibration metrics are conditional on evaluable posterior production. All-attempt coverage bounds preserve every planned cell by treating non-successes as uncovered or covered. Fixed-truth repeated-sampling only; not joint SBC or observational H0 inference.",
    }
    (output / "posthoc_report.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("h0_pilot_posthoc"))
    parser.add_argument("--seed", type=int, default=20261119)
    args = parser.parse_args()
    print(json.dumps(analyze(args.artifact_dir, args.output_dir, args.seed), indent=2))


if __name__ == "__main__":
    main()
