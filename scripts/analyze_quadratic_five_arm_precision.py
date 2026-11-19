"""Precision readout for paired five-arm 30-event replicate artifacts.

Produces Wilson intervals for empirical coverage and replicate-paired t intervals
for the prespecified contrasts E-D, E-C, D-C, and B-C. Failed cells remain in
the planned denominator; success-only estimates and failure-inclusive bounds
are both retained. This is descriptive Monte Carlo analysis, not calibration
certification or observational inference.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, t

ARMS = ["A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
        "D_map_plugin_ipw", "E_uncertainty_integrated_ipw"]
CONTRASTS = [("E_uncertainty_integrated_ipw", "D_map_plugin_ipw", "E_minus_D"),
             ("E_uncertainty_integrated_ipw", "C_known_probability_ipw", "E_minus_C"),
             ("D_map_plugin_ipw", "C_known_probability_ipw", "D_minus_C"),
             ("B_parent_oracle", "C_known_probability_ipw", "B_minus_C")]


def wilson_interval(successes: int, n: int, confidence: float = 0.95) -> tuple[float, float]:
    if n <= 0:
        return float("nan"), float("nan")
    z = float(norm.ppf(0.5 + confidence / 2))
    p = successes / n
    den = 1 + z*z/n
    center = (p + z*z/(2*n)) / den
    half = z * np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return float(center-half), float(center+half)


def paired_t_interval(values: np.ndarray, confidence: float = 0.95) -> tuple[float, float, float, float]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    n = int(x.size)
    if n == 0:
        return float("nan"), float("nan"), float("nan"), float("nan")
    mean = float(x.mean())
    if n == 1:
        return mean, float("nan"), float("nan"), float("nan")
    se = float(x.std(ddof=1) / np.sqrt(n))
    q = float(t.ppf(0.5 + confidence/2, n-1))
    return mean, se, mean-q*se, mean+q*se


def analyze(root: Path, output_prefix: Path) -> dict:
    root = root.resolve()
    cfg = json.loads((root / "config.json").read_text())
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    summaries = pd.read_csv(root / "replicate_posterior_summaries.csv")
    reps = int(cfg["replicates"])
    truth = float(cfg["h0_true"])
    prior_lo, prior_hi = map(float, cfg["h0_prior_support"])
    max_error = max(abs(truth-prior_lo), abs(prior_hi-truth))
    planned = reps * len(cfg["gammas"]) * len(cfg["validation_sizes"]) * len(ARMS)
    if len(ledger) != planned:
        raise AssertionError(f"ledger has {len(ledger)} rows; expected {planned}")
    if ledger.cell_id.duplicated().any():
        raise AssertionError("duplicate lifecycle cell IDs")

    cell_rows = []
    for (gamma, nval, arm), planned_rows in ledger.groupby(
            ["gamma", "validation_size", "arm"], sort=True):
        g = summaries[(summaries.gamma == gamma) &
                      (summaries.validation_size == nval) &
                      (summaries.arm == arm) &
                      (summaries.status == "SUCCESS")]
        errors = g.median_error.to_numpy(float)
        widths = (g.q84-g.q16).to_numpy(float)
        covered = ((g.q16.to_numpy(float) <= truth) &
                   (truth <= g.q84.to_numpy(float)))
        n_ok = len(g)
        n_fail = reps - n_ok
        wilson_lo, wilson_hi = wilson_interval(int(covered.sum()), n_ok)
        cell_rows.append({
            "gamma": float(gamma), "validation_size": int(nval), "arm": arm,
            "planned_n": reps, "successful_n": n_ok, "failed_or_non_evaluable_n": n_fail,
            "success_only_bias": float(errors.mean()) if n_ok else np.nan,
            "success_only_rmse": float(np.sqrt(np.mean(errors**2))) if n_ok else np.nan,
            "success_only_68_coverage": float(covered.mean()) if n_ok else np.nan,
            "coverage_wilson95_low": wilson_lo, "coverage_wilson95_high": wilson_hi,
            "success_only_mean_68_width": float(widths.mean()) if n_ok else np.nan,
            "failure_inclusive_rmse_low": float(np.sqrt(np.sum(errors**2)/reps)),
            "failure_inclusive_rmse_high": float(np.sqrt((np.sum(errors**2)+n_fail*max_error**2)/reps),
            "failure_inclusive_coverage_low": float(covered.sum()/reps),
            "failure_inclusive_coverage_high": float((covered.sum()+n_fail)/reps),
        })
    cells = pd.DataFrame(cell_rows)
    cells.to_csv(output_prefix.with_name(output_prefix.name + "_cells.csv"), index=False)

    paired_rows = []
    for (gamma, nval), group in summaries.groupby(["gamma", "validation_size"], sort=True):
        by_rep = {arm: frame.set_index("replicate") for arm, frame in group.groupby("arm")}
        for left, right, label in CONTRASTS:
            if left not in by_rep or right not in by_rep:
                continue
            common = by_rep[left].join(by_rep[right], how="inner", lsuffix="_left",
                                       rsuffix="_right")
            common = common[(common.status_left == "SUCCESS") &
                            (common.status_right == "SUCCESS")]
            metrics = {
                "median_error": common.median_error_left-common.median_error_right,
                "absolute_error": common.median_error_left.abs()-common.median_error_right.abs(),
                "squared_error": common.median_error_left**2-common.median_error_right**2,
                "coverage68": ((common.q16_left <= truth) & (truth <= common.q84_left)).astype(float) -
                              ((common.q16_right <= truth) & (truth <= common.q84_right)).astype(float),
                "interval_width68": (common.q84_left-common.q16_left) -
                                    (common.q84_right-common.q16_right),
            }
            for metric, values in metrics.items():
                mean, se, lo, hi = paired_t_interval(values.to_numpy(float))
                paired_rows.append({
                    "gamma": float(gamma), "validation_size": int(nval),
                    "contrast": label, "left_arm": left, "right_arm": right,
                    "metric": metric, "paired_success_n": int(values.notna().sum()),
                    "mean_paired_difference": mean, "paired_se": se,
                    "t95_low": lo, "t95_high": hi,
                })
    paired = pd.DataFrame(paired_rows)
    paired.to_csv(output_prefix.with_name(output_prefix.name + "_paired_contrasts.csv"),
                  index=False)

    result = {
        "status": "PASS",
        "planned_cells": planned,
        "lifecycle_status_counts": ledger.status.value_counts().to_dict(),
        "cell_summary_rows": len(cells),
        "paired_contrast_summary_rows": len(paired),
        "contrasts": [label for _, _, label in CONTRASTS],
        "coverage_interval": "Wilson 95% interval among successful replicates",
        "paired_intervals": "two-sided Student-t 95% intervals on replicate-paired differences",
        "failure_handling": "success-only estimates shown separately; failure-inclusive RMSE and coverage bounds retain all planned replicates",
        "interpretation": "Descriptive synthetic engineering-pilot readout only; not calibration certification or observational inference.",
    }
    output_prefix.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact_dir", type=Path)
    p.add_argument("--output-prefix", type=Path, default=Path("five_arm_precision_readout"))
    args = p.parse_args()
    print(json.dumps(analyze(args.artifact_dir, args.output_prefix), indent=2))


if __name__ == "__main__":
    main()
