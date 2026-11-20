#!/usr/bin/env python3
"""Independent descriptive readout for prior-predictive five-arm SBC pilots.

Consumes the audited aggregate CSVs. It does not rerun or modify inference.
All arm comparisons are paired within replicate/gamma/validation-size.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t

ARMS = [
    "A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
    "D_map_plugin_ipw", "E_uncertainty_integrated_ipw",
]
COMPARISONS = [
    ("B_parent_oracle", "A_catalogue_only"),
    ("C_known_probability_ipw", "A_catalogue_only"),
    ("D_map_plugin_ipw", "C_known_probability_ipw"),
    ("E_uncertainty_integrated_ipw", "C_known_probability_ipw"),
    ("E_uncertainty_integrated_ipw", "D_map_plugin_ipw"),
]

def _mean_ci(values: np.ndarray) -> tuple[float, float, float, float]:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or x.size == 0 or not np.isfinite(x).all():
        raise AssertionError("paired Monte Carlo summaries require a non-empty finite vector")
    n = len(x)
    if not n:
        return (float("nan"),) * 4
    mean = float(x.mean())
    if n < 2:
        return mean, float("nan"), float("nan"), float("nan")
    se = float(x.std(ddof=1) / np.sqrt(n))
    critical = float(t.ppf(0.975, n - 1))
    return mean, se, mean - critical * se, mean + critical * se

def analyze(root: Path, output: Path) -> dict:
    root = root.resolve()
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    summary = pd.read_csv(root / "replicate_posterior_summaries.csv")
    diagnostics = pd.read_csv(root / "replicate_weight_diagnostics.csv")
    key = ["replicate", "gamma", "validation_size", "arm"]
    if ledger.duplicated(key).any() or summary.duplicated(key).any() or diagnostics.duplicated(key).any():
        raise AssertionError("duplicate replicate/gamma/size/arm key")
    if not (ledger.status == "SUCCESS").all() or not (ledger.audit_status == "PASS").all():
        raise AssertionError("readout requires every planned cell to be successful and audited")
    if set(map(tuple, ledger[key].itertuples(index=False, name=None))) != set(map(tuple, summary[key].itertuples(index=False, name=None))):
        raise AssertionError("summary rows do not exactly cover successful lifecycle cells")
    if set(map(tuple, summary[key].itertuples(index=False, name=None))) != set(map(tuple, diagnostics[key].itertuples(index=False, name=None))):
        raise AssertionError("weight diagnostics do not exactly cover summary cells")

    rows = []
    for (gamma, nval, arm), g in summary.groupby(["gamma", "validation_size", "arm"], sort=True):
        ranks = g.posterior_rank.to_numpy(float)
        errors = g.median_error.to_numpy(float)
        if not np.isfinite(ranks).all() or np.any((ranks < 0) | (ranks > 1)):
            raise AssertionError("invalid posterior ranks")
        required_numeric = ["h0_true", "median_error", "q16", "q84"]
        if not np.isfinite(g[required_numeric].to_numpy(dtype=float)).all():
            raise AssertionError("non-finite posterior summary field")
        if np.any(g.q16.to_numpy(float) > g.q84.to_numpy(float)):
            raise AssertionError("posterior interval endpoints are reversed")
        rank_mean, rank_se, rank_lo, rank_hi = _mean_ci(ranks)
        bias, bias_se, bias_lo, bias_hi = _mean_ci(errors)
        sq = errors**2
        rmse = float(np.sqrt(sq.mean()))
        coverage = float(((g.q16 <= g.h0_true) & (g.h0_true <= g.q84)).mean())
        wd = diagnostics[(diagnostics.gamma == gamma) & (diagnostics.validation_size == nval) & (diagnostics.arm == arm)]
        ess_col = "minimum_node_kish_ess" if arm == "E_uncertainty_integrated_ipw" else "kish_ess"
        ess = wd[ess_col].to_numpy(float)
        rows.append({
            "gamma": float(gamma), "validation_size": int(nval), "arm": arm,
            "replicates": len(g), "mean_rank": rank_mean, "rank_mcse": rank_se,
            "rank_mean_t95_low": rank_lo, "rank_mean_t95_high": rank_hi,
            "rank_mean_abs_deviation_from_half": float(np.mean(np.abs(ranks - 0.5))),
            "median_error_bias": bias, "bias_mcse": bias_se,
            "bias_t95_low": bias_lo, "bias_t95_high": bias_hi,
            "median_error_rmse": rmse, "median_error_mae": float(np.mean(np.abs(errors))),
            "coverage_68": coverage,
            "mean_68_interval_width": float(np.mean(g.q84.to_numpy(float)-g.q16.to_numpy(float))),
            "mean_minimum_kish_ess": float(np.mean(ess)),
            "mean_abs_rank_deviation_spearman_ess": (
                float(spearmanr(ess, np.abs(ranks-0.5)).statistic)
                if len(ess) >= 3 and np.std(ess) > 0 and np.std(np.abs(ranks-0.5)) > 0 else float("nan")
            ),
            "mean_abs_error_spearman_ess": (
                float(spearmanr(ess, np.abs(errors)).statistic)
                if len(ess) >= 3 and np.std(ess) > 0 and np.std(np.abs(errors)) > 0 else float("nan")
            ),
        })
    # Each replicate represents one independently drawn truth, shared across
    # all curvature settings, validation sizes, and inference arms.
    truth_counts = summary.groupby("replicate").h0_true.nunique()
    if not (truth_counts == 1).all():
        raise AssertionError("H0 truth must be identical across arms, curvatures, and validation sizes within replicate")

    cell_counts = ledger.groupby(["replicate", "gamma", "validation_size"]).arm.agg(list)
    expected_arms = set(ARMS)
    if not all(len(arms) == len(ARMS) and set(arms) == expected_arms for arms in cell_counts):
        raise AssertionError("each replicate/gamma/validation-size cell must contain exactly the five distinct arms")

    non_e_rows = diagnostics[diagnostics.arm != "E_uncertainty_integrated_ipw"]
    non_e_ess = non_e_rows["kish_ess"].to_numpy(float)
    if not np.isfinite(non_e_ess).all() or not (non_e_ess > 0).all():
        raise AssertionError("non-integrated-arm Kish ESS diagnostics must be finite and positive")
    e_rows = diagnostics[diagnostics.arm == "E_uncertainty_integrated_ipw"]
    if not np.isfinite(e_rows["minimum_node_kish_ess"].to_numpy(float)).all() or not (e_rows["minimum_node_kish_ess"].to_numpy(float) > 0).all():
        raise AssertionError("integrated-arm node ESS must be finite and positive")

    strata = pd.DataFrame(rows)

    paired_rows = []
    for (gamma, nval), group in summary.groupby(["gamma", "validation_size"], sort=True):
        wide = group.pivot(index="replicate", columns="arm", values=["median_error", "posterior_rank", "q16", "q84", "h0_true"])
        if wide.isna().any().any() or set(wide.columns.get_level_values(1)) != set(ARMS):
            raise AssertionError("incomplete paired arm block")
        for arm, reference in COMPARISONS:
            for metric in ("median_error", "posterior_rank"):
                delta = (wide[(metric, arm)] - wide[(metric, reference)]).to_numpy(float)
                mean, se, lo, hi = _mean_ci(delta)
                paired_rows.append({
                    "gamma": float(gamma), "validation_size": int(nval),
                    "arm": arm, "reference_arm": reference, "metric": metric,
                    "mean_paired_difference": mean, "paired_mcse": se,
                    "paired_t95_low": lo, "paired_t95_high": hi,
                    "fraction_negative_difference": float(np.mean(delta < 0)),
                    "replicates": len(delta),
                })
            # Compare squared median error; positive means arm is worse.
            delta_sq = (wide[("median_error", arm)]**2 - wide[("median_error", reference)]**2).to_numpy(float)
            mean, se, lo, hi = _mean_ci(delta_sq)
            paired_rows.append({
                "gamma": float(gamma), "validation_size": int(nval),
                "arm": arm, "reference_arm": reference, "metric": "squared_median_error",
                "mean_paired_difference": mean, "paired_mcse": se,
                "paired_t95_low": lo, "paired_t95_high": hi,
                "fraction_negative_difference": float(np.mean(delta_sq < 0)),
                "replicates": len(delta_sq),
            })
    paired = pd.DataFrame(paired_rows)
    output.mkdir(parents=True, exist_ok=True)
    strata.to_csv(output / "prior_predictive_detailed_strata.csv", index=False)
    paired.to_csv(output / "prior_predictive_paired_contrasts.csv", index=False)
    result = {
        "status": "PASS",
        "planned_cells": int(len(ledger)),
        "successful_cells": int((ledger.status == "SUCCESS").sum()),
        "strata_rows": int(len(strata)),
        "paired_contrast_rows": int(len(paired)),
        "interpretation": "Descriptive Monte Carlo summaries. Paired t intervals quantify finite-replicate uncertainty under the simulated design; small pilot strata and correlated repeated validation sizes limit confirmatory interpretation. ESS correlations are exploratory, not causal.",
    }
    (output / "detailed_readout_report.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact_dir", type=Path)
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args()
    result = analyze(a.artifact_dir, a.output_dir)
    print(json.dumps(result, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()
