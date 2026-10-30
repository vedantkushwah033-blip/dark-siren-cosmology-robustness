"""Scenario-balanced calibration summaries for nested SRSWOR/bootstrap experiments.

This module summarizes already generated outer repeated-SRSWOR and inner
bootstrap results. It does not generate data or establish validity by itself.
The empirical outer SRSWOR standard deviation is the reference for uncertainty
calibration, while interval coverage is evaluated against the fixed
scenario-specific finite-parent posterior-mean target.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


_REQUIRED = {
    "scenario", "validation_size", "outer_draw", "method", "status",
    "posterior_mean", "bootstrap_se", "interval_low", "interval_high",
    "finite_parent_target",
}


def summarize_bootstrap_calibration(records: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return scenario-level and scenario-balanced method/size summaries.

    Required input grain is one row per (scenario, validation_size, outer_draw,
    method). Failed rows remain in the input and denominator. The empirical
    repeated-SRSWOR SD is computed separately for each scenario and validation
    size from successful outer estimates, then compared with the mean reported
    inner-bootstrap SE. No outer draws are treated as independent scenarios.
    """
    missing = sorted(_REQUIRED.difference(records.columns))
    if missing:
        raise ValueError(f"records missing required columns: {missing}")
    frame = records.copy()
    if frame.empty:
        raise ValueError("records must not be empty")
    key = ["scenario", "validation_size", "outer_draw", "method"]
    if frame[key].isna().any().any():
        raise ValueError("design identifiers and method must be nonmissing")
    if frame.duplicated(key).any():
        raise ValueError("records must be unique at scenario/size/draw/method grain")
    if not frame["status"].isin(["success", "failed"]).all():
        raise ValueError("status must be success or failed")
    numeric = ["posterior_mean", "bootstrap_se", "interval_low", "interval_high", "finite_parent_target"]
    for column in numeric:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    successful = frame["status"].eq("success")
    if frame.loc[successful, numeric].isna().any().any():
        raise ValueError("successful records require all finite calibration metrics")
    if not np.isfinite(frame.loc[successful, numeric].to_numpy(dtype=float)).all():
        raise ValueError("successful calibration metrics must be finite")
    if (frame.loc[successful, "bootstrap_se"] < 0).any():
        raise ValueError("bootstrap_se must be nonnegative")
    if (frame.loc[successful, "interval_low"] > frame.loc[successful, "interval_high"]).any():
        raise ValueError("interval_low cannot exceed interval_high")
    if frame.groupby(["scenario", "validation_size"])["finite_parent_target"].nunique(dropna=True).gt(1).any():
        raise ValueError("finite_parent_target must be fixed within each scenario and validation size")

    rows: list[dict] = []
    for (scenario, size, method), group in frame.groupby(
        ["scenario", "validation_size", "method"], sort=True, dropna=False
    ):
        ok = group.loc[group["status"].eq("success")]
        estimates = ok["posterior_mean"].to_numpy(dtype=float)
        empirical_sd = float(np.std(estimates, ddof=1)) if estimates.size >= 2 else float("nan")
        mean_se = float(ok["bootstrap_se"].mean()) if not ok.empty else float("nan")
        ratios = mean_se / empirical_sd if np.isfinite(empirical_sd) and empirical_sd > 0 else float("nan")
        covered = (
            (ok["interval_low"] <= ok["finite_parent_target"])
            & (ok["finite_parent_target"] <= ok["interval_high"])
        )
        rows.append({
            "scenario": scenario,
            "validation_size": int(size),
            "method": str(method),
            "n_outer_attempted": int(len(group)),
            "n_outer_successful": int(len(ok)),
            "outer_failure_fraction": float(1.0 - len(ok) / len(group)),
            "empirical_srswor_sd": empirical_sd,
            "mean_inner_bootstrap_se": mean_se,
            "mean_se_to_empirical_sd_ratio": ratios,
            "finite_parent_target_interval_coverage": float(covered.mean()) if len(ok) else float("nan"),
            "mean_interval_width": float((ok["interval_high"] - ok["interval_low"]).mean()) if len(ok) else float("nan"),
            "n_intervals_evaluated": int(len(ok)),
        })
    scenario_level = pd.DataFrame(rows)

    aggregate_rows: list[dict] = []
    for (size, method), group in scenario_level.groupby(["validation_size", "method"], sort=True):
        ratio = group["mean_se_to_empirical_sd_ratio"].dropna()
        coverage = group["finite_parent_target_interval_coverage"].dropna()
        widths = group["mean_interval_width"].dropna()
        aggregate_rows.append({
            "validation_size": int(size),
            "method": str(method),
            "n_scenarios": int(group["scenario"].nunique()),
            "n_scenarios_with_estimable_sd_ratio": int(ratio.size),
            "scenario_balanced_mean_se_ratio": float(ratio.mean()) if ratio.size else float("nan"),
            "scenario_balanced_sd_se_ratio": float(ratio.std(ddof=1)) if ratio.size > 1 else float("nan"),
            "scenario_balanced_mean_coverage": float(coverage.mean()) if coverage.size else float("nan"),
            "scenario_balanced_sd_coverage": float(coverage.std(ddof=1)) if coverage.size > 1 else float("nan"),
            "scenario_balanced_mean_interval_width": float(widths.mean()) if widths.size else float("nan"),
            "scenario_balanced_outer_failure_fraction": float(group["outer_failure_fraction"].mean()),
        })
    return scenario_level, pd.DataFrame(aggregate_rows)
