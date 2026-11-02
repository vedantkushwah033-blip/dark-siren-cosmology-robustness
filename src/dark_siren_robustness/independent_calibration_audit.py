"""Auditor-owned reconstruction of nested-bootstrap calibration summaries.

This module deliberately does not import the production calibration summarizer.
It rebuilds scenario-level and scenario-balanced tables from raw calibration
records so the replay auditor can compare stored summaries against a separately
implemented calculation.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


_REQUIRED = {
    "scenario", "validation_size", "outer_draw", "method",
    "outer_status", "bootstrap_status", "outer_posterior_mean",
    "bootstrap_se", "interval_low", "interval_high", "finite_parent_target",
}


def _mean_t_interval(values: pd.Series, confidence: float = 0.95):
    values = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    values = values[np.isfinite(values)]
    count = int(values.size)
    if count == 0:
        return float("nan"), float("nan"), float("nan"), 0
    mean = float(np.mean(values))
    if count == 1:
        return mean, float("nan"), float("nan"), count
    standard_error = float(np.std(values, ddof=1) / np.sqrt(count))
    critical = float(student_t.ppf(0.5 + confidence / 2.0, count - 1))
    return (
        mean,
        mean - critical * standard_error,
        mean + critical * standard_error,
        count,
    )


def reconstruct_calibration_summaries(records: pd.DataFrame):
    """Independently rebuild both published summary tables from raw rows."""
    missing = sorted(_REQUIRED.difference(records.columns))
    if missing:
        raise ValueError(f"audit records missing columns: {missing}")
    if records.empty:
        raise ValueError("audit records are empty")

    frame = records.copy()
    keys = ["scenario", "validation_size", "outer_draw", "method"]
    if frame[keys].isna().any().any() or frame.duplicated(keys).any():
        raise ValueError("audit design keys must be nonmissing and unique")
    if not frame["outer_status"].isin(["success", "failed"]).all():
        raise ValueError("audit outer status is invalid")
    if not frame["bootstrap_status"].isin(["success", "failed", "not_run"]).all():
        raise ValueError("audit bootstrap status is invalid")

    for name in ("outer_posterior_mean", "bootstrap_se", "interval_low",
                 "interval_high", "finite_parent_target"):
        frame[name] = pd.to_numeric(frame[name], errors="coerce")

    if not np.isfinite(frame["finite_parent_target"].to_numpy(dtype=float)).all():
        raise ValueError("audit finite-parent target must be finite")
    outer_success = frame["outer_status"].eq("success")
    boot_success = frame["bootstrap_status"].eq("success")
    if frame.loc[outer_success, "outer_posterior_mean"].isna().any():
        raise ValueError("audit successful outer rows require posterior means")
    if frame.loc[~outer_success, "outer_posterior_mean"].notna().any():
        raise ValueError("audit failed outer rows cannot carry posterior means")
    if frame.loc[boot_success, ["bootstrap_se", "interval_low", "interval_high"]].isna().any().any():
        raise ValueError("audit successful bootstrap rows require interval metrics")
    if frame.loc[~boot_success, ["bootstrap_se", "interval_low", "interval_high"]].notna().any().any():
        raise ValueError("audit unsuccessful bootstrap rows cannot carry interval metrics")
    if not np.isfinite(frame.loc[outer_success, "outer_posterior_mean"]).all():
        raise ValueError("audit successful outer means must be finite")
    metrics = frame.loc[boot_success, ["bootstrap_se", "interval_low", "interval_high"]]
    if not np.isfinite(metrics.to_numpy(dtype=float)).all():
        raise ValueError("audit successful bootstrap metrics must be finite")
    if (metrics["bootstrap_se"] < 0).any():
        raise ValueError("audit bootstrap standard errors cannot be negative")
    if (metrics["interval_low"] > metrics["interval_high"]).any():
        raise ValueError("audit interval endpoints are reversed")
    if (frame.loc[~outer_success, "bootstrap_status"] != "not_run").any():
        raise ValueError("audit bootstrap must not run after outer failure")
    if (frame.loc[outer_success, "bootstrap_status"] == "not_run").any():
        raise ValueError("audit bootstrap status missing after successful outer fit")

    for _, cell in frame.groupby(["scenario", "validation_size"], sort=True):
        method_set = set(cell["method"].astype(str))
        per_draw = cell.groupby("outer_draw")["method"].apply(lambda x: set(x.astype(str)))
        if any(observed != method_set for observed in per_draw):
            raise ValueError("audit method set differs across outer draws")

    consistency = frame.groupby(
        ["scenario", "validation_size", "outer_draw"], dropna=False
    ).agg(
        status_count=("outer_status", "nunique"),
        mean_count=("outer_posterior_mean", "nunique"),
        target_count=("finite_parent_target", "nunique"),
    )
    if (consistency["status_count"] > 1).any() or (consistency["mean_count"] > 1).any():
        raise ValueError("audit outer result differs across methods")
    if (consistency["target_count"] > 1).any():
        raise ValueError("audit finite-parent target differs within an outer draw")
    if frame.groupby("scenario")["finite_parent_target"].nunique().gt(1).any():
        raise ValueError("audit finite-parent target changes within a scenario")

    # Outer variability is computed once per shared outer draw, not once per
    # method. Failed outer draws remain in the attempted denominator.
    outer = frame.drop_duplicates(["scenario", "validation_size", "outer_draw"])
    outer_rows = []
    for (scenario, size), group in outer.groupby(
        ["scenario", "validation_size"], sort=True
    ):
        estimates = group.loc[
            group["outer_status"].eq("success"), "outer_posterior_mean"
        ].to_numpy(dtype=float)
        outer_rows.append({
            "scenario": scenario,
            "validation_size": int(size),
            "n_outer_attempted": int(len(group)),
            "n_outer_successful": int(estimates.size),
            "outer_failure_fraction": float(1.0 - estimates.size / len(group)),
            "empirical_srswor_sd": (
                float(np.std(estimates, ddof=1)) if estimates.size >= 2 else float("nan")
            ),
        })
    outer_summary = pd.DataFrame(outer_rows)

    scenario_rows = []
    for (scenario, size, method), group in frame.groupby(
        ["scenario", "validation_size", "method"], sort=True, dropna=False
    ):
        successful = group.loc[group["bootstrap_status"].eq("success")]
        outer_info = outer_summary.loc[
            outer_summary["scenario"].eq(scenario)
            & outer_summary["validation_size"].eq(size)
        ].iloc[0]
        empirical_sd = float(outer_info["empirical_srswor_sd"])
        mean_se = (
            float(successful["bootstrap_se"].mean())
            if not successful.empty else float("nan")
        )
        ratio = (
            mean_se / empirical_sd
            if np.isfinite(empirical_sd) and empirical_sd > 0 else float("nan")
        )
        covered = (
            (successful["interval_low"] <= successful["finite_parent_target"])
            & (successful["finite_parent_target"] <= successful["interval_high"])
        )
        n_outer = int(outer_info["n_outer_successful"])
        scenario_rows.append({
            "scenario": scenario,
            "validation_size": int(size),
            "method": str(method),
            "n_outer_attempted": int(outer_info["n_outer_attempted"]),
            "n_outer_successful": n_outer,
            "outer_failure_fraction": float(outer_info["outer_failure_fraction"]),
            "n_bootstrap_attempted": n_outer,
            "n_bootstrap_successful": int(len(successful)),
            "bootstrap_failure_fraction": (
                float(1.0 - len(successful) / n_outer) if n_outer else float("nan")
            ),
            "empirical_srswor_sd": empirical_sd,
            "mean_inner_bootstrap_se": mean_se,
            "mean_se_to_empirical_sd_ratio": ratio,
            "finite_parent_target_interval_coverage": (
                float(covered.mean()) if len(successful) else float("nan")
            ),
            "mean_interval_width": (
                float((successful["interval_high"] - successful["interval_low"]).mean())
                if len(successful) else float("nan")
            ),
            "n_intervals_evaluated": int(len(successful)),
        })
    scenario_level = pd.DataFrame(scenario_rows)

    aggregate_rows = []
    for (size, method), group in scenario_level.groupby(
        ["validation_size", "method"], sort=True
    ):
        ratio = group["mean_se_to_empirical_sd_ratio"].dropna()
        coverage = group["finite_parent_target_interval_coverage"].dropna()
        widths = group["mean_interval_width"].dropna()
        ratio_mean, ratio_low, ratio_high, ratio_n = _mean_t_interval(ratio)
        coverage_mean, coverage_low, coverage_high, coverage_n = _mean_t_interval(coverage)
        width_mean, width_low, width_high, width_n = _mean_t_interval(widths)
        outer_mean, outer_low, outer_high, _ = _mean_t_interval(
            group["outer_failure_fraction"]
        )
        inner_mean, inner_low, inner_high, _ = _mean_t_interval(
            group["bootstrap_failure_fraction"]
        )
        aggregate_rows.append({
            "validation_size": int(size),
            "method": str(method),
            "n_scenarios": int(group["scenario"].nunique()),
            "n_scenarios_with_estimable_sd_ratio": int(ratio_n),
            "scenario_balanced_mean_se_ratio": ratio_mean,
            "scenario_balanced_mean_se_ratio_ci95_low": ratio_low,
            "scenario_balanced_mean_se_ratio_ci95_high": ratio_high,
            "scenario_balanced_sd_se_ratio": (
                float(ratio.std(ddof=1)) if ratio.size > 1 else float("nan")
            ),
            "scenario_balanced_mean_coverage": coverage_mean,
            "scenario_balanced_mean_coverage_ci95_low": coverage_low,
            "scenario_balanced_mean_coverage_ci95_high": coverage_high,
            "n_scenarios_with_coverage": int(coverage_n),
            "scenario_balanced_sd_coverage": (
                float(coverage.std(ddof=1)) if coverage.size > 1 else float("nan")
            ),
            "scenario_balanced_mean_interval_width": width_mean,
            "scenario_balanced_mean_interval_width_ci95_low": width_low,
            "scenario_balanced_mean_interval_width_ci95_high": width_high,
            "n_scenarios_with_interval_width": int(width_n),
            "scenario_balanced_outer_failure_fraction": outer_mean,
            "scenario_balanced_outer_failure_fraction_ci95_low": outer_low,
            "scenario_balanced_outer_failure_fraction_ci95_high": outer_high,
            "scenario_balanced_bootstrap_failure_fraction": inner_mean,
            "scenario_balanced_bootstrap_failure_fraction_ci95_low": inner_low,
            "scenario_balanced_bootstrap_failure_fraction_ci95_high": inner_high,
        })
    return scenario_level, pd.DataFrame(aggregate_rows)
