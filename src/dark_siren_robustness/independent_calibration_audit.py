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
    if frame[keys].isna().any().any():
        raise ValueError("audit design keys must be nonmissing")
    if frame.duplicated(keys).any():
        raise ValueError("audit design keys must be unique")
    size_values = pd.to_numeric(frame["validation_size"], errors="coerce")
    if (
        size_values.isna().any() or (size_values < 1).any()
        or not np.equal(size_values, np.floor(size_values)).all()
    ):
        raise ValueError("audit validation_size must contain positive integer values")
    frame["validation_size"] = size_values.astype(int)
    draw_values = pd.to_numeric(frame["outer_draw"], errors="coerce")
    if (
        draw_values.isna().any() or (draw_values < 0).any()
        or not np.equal(draw_values, np.floor(draw_values)).all()
    ):
        raise ValueError("audit outer_draw must contain nonnegative integer identifiers")
    if frame["method"].astype(str).str.strip().str.len().eq(0).any():
        raise ValueError("audit method labels must be nonempty or whitespace-only")
    if not frame["outer_status"].isin(["success", "failed"]).all():
        raise ValueError("audit outer status is invalid")
    if not frame["bootstrap_status"].isin(["success", "non_evaluable", "failed", "not_run"]).all():
        raise ValueError("audit bootstrap status is invalid")

    count_columns = (
        "n_bootstrap_attempted", "n_bootstrap_successful", "n_bootstrap_failed"
    )
    has_count_ledger = any(column in frame.columns for column in count_columns)
    if has_count_ledger and not all(column in frame.columns for column in count_columns):
        raise ValueError("audit bootstrap count ledger is incomplete")
    for column in count_columns:
        if column not in frame:
            frame[column] = np.nan
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    if has_count_ledger:
        for _, row in frame.iterrows():
            status = str(row["bootstrap_status"])
            values = [row[column] for column in count_columns]
            if status in ("success", "non_evaluable"):
                if any(pd.isna(value) or not np.isfinite(value) or value < 0 or value != int(value) for value in values):
                    raise ValueError("audit evaluated procedures require finite nonnegative integer counts")
                attempted, successful, failed = map(int, values)
                if successful + failed != attempted:
                    raise ValueError("audit inner bootstrap counts do not reconcile")
                if (status == "success" and successful < 2) or (status == "non_evaluable" and successful >= 2):
                    raise ValueError("audit bootstrap status disagrees with successful inner count")
            elif status == "failed":
                if any(pd.notna(value) for value in values):
                    raise ValueError("audit wrapper failures must have unknown inner counts")
            elif status == "not_run":
                if any(pd.isna(value) or int(value) != 0 for value in values):
                    raise ValueError("audit not-run procedures must have zero inner counts")

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
            "n_bootstrap_attempted": int(group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].ne("not_run").sum()) if "bootstrap_is_derived" in group else int(group["bootstrap_status"].ne("not_run").sum()),
            "n_bootstrap_successful": int(group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].eq("success").sum()) if "bootstrap_is_derived" in group else int(group["bootstrap_status"].eq("success").sum()),
            "n_bootstrap_procedures_attempted": int(group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].ne("not_run").sum()) if "bootstrap_is_derived" in group else int(group["bootstrap_status"].ne("not_run").sum()),
            "n_bootstrap_procedures_successful": int(group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].eq("success").sum()) if "bootstrap_is_derived" in group else int(group["bootstrap_status"].eq("success").sum()),
            "n_bootstrap_procedures_non_evaluable": int(group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].eq("non_evaluable").sum()) if "bootstrap_is_derived" in group else int(group["bootstrap_status"].eq("non_evaluable").sum()),
            "n_bootstrap_procedures_failed": int(group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].eq("failed").sum()) if "bootstrap_is_derived" in group else int(group["bootstrap_status"].eq("failed").sum()),
            "n_bootstrap_procedures_not_run": int(group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].eq("not_run").sum()) if "bootstrap_is_derived" in group else int(group["bootstrap_status"].eq("not_run").sum()),
            "n_bootstrap_procedures_with_unknown_inner_count": int(group.loc[~group["bootstrap_is_derived"].fillna(False), "n_bootstrap_attempted"].isna().where(group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].ne("not_run"), False).sum()) if "bootstrap_is_derived" in group else int(group["n_bootstrap_attempted"].isna().where(group["bootstrap_status"].ne("not_run"), False).sum()),
            "n_inner_resamples_attempted": int(group.loc[~group["bootstrap_is_derived"].fillna(False), "n_bootstrap_attempted"].dropna().sum()) if "bootstrap_is_derived" in group else int(group["n_bootstrap_attempted"].dropna().sum()),
            "n_inner_resamples_successful": int(group.loc[~group["bootstrap_is_derived"].fillna(False), "n_bootstrap_successful"].dropna().sum()) if "bootstrap_is_derived" in group else int(group["n_bootstrap_successful"].dropna().sum()),
            "n_inner_resamples_failed": int(group.loc[~group["bootstrap_is_derived"].fillna(False), "n_bootstrap_failed"].dropna().sum()) if "bootstrap_is_derived" in group else int(group["n_bootstrap_failed"].dropna().sum()),
            "inner_resample_failure_fraction": (
                float(group.loc[~group["bootstrap_is_derived"].fillna(False), "n_bootstrap_failed"].dropna().sum()
                      / group.loc[~group["bootstrap_is_derived"].fillna(False), "n_bootstrap_attempted"].dropna().sum())
                if "bootstrap_is_derived" in group and group.loc[~group["bootstrap_is_derived"].fillna(False), "n_bootstrap_attempted"].dropna().sum() > 0
                else (float(group["n_bootstrap_failed"].dropna().sum() / group["n_bootstrap_attempted"].dropna().sum())
                      if "bootstrap_is_derived" not in group and group["n_bootstrap_attempted"].dropna().sum() > 0 else float("nan"))
            ),
            "bootstrap_wrapper_failure_fraction": (
                float(group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].eq("failed").sum()
                      / group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].ne("not_run").sum())
                if "bootstrap_is_derived" in group and group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].ne("not_run").sum() > 0
                else (float(group["bootstrap_status"].eq("failed").sum() / group["bootstrap_status"].ne("not_run").sum())
                      if "bootstrap_is_derived" not in group and group["bootstrap_status"].ne("not_run").sum() > 0 else float("nan"))
            ),
            "bootstrap_failure_fraction": (
                float(group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].isin(["non_evaluable", "failed"]).sum()
                      / group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].ne("not_run").sum())
                if "bootstrap_is_derived" in group and group.loc[~group["bootstrap_is_derived"].fillna(False), "bootstrap_status"].ne("not_run").sum() > 0
                else (float(group["bootstrap_status"].isin(["non_evaluable", "failed"]).sum() / group["bootstrap_status"].ne("not_run").sum())
                      if "bootstrap_is_derived" not in group and group["bootstrap_status"].ne("not_run").sum() > 0 else float("nan"))
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
            "scenario_balanced_inner_resample_failure_fraction": (
                float(group["inner_resample_failure_fraction"].dropna().mean())
                if group["inner_resample_failure_fraction"].notna().any() else float("nan")
            ),
            "scenario_balanced_wrapper_failure_fraction": (
                float(group["bootstrap_wrapper_failure_fraction"].dropna().mean())
                if group["bootstrap_wrapper_failure_fraction"].notna().any() else float("nan")
            ),
        })
    return scenario_level, pd.DataFrame(aggregate_rows)
