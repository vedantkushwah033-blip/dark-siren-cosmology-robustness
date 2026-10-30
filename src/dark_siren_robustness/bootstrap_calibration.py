"""Scenario-balanced calibration summaries for nested SRSWOR/bootstrap experiments.

The outer SRSWOR estimator and inner bootstrap are separate layers. Empirical
outer-sampling SD is calculated once from outer estimates, independently of
whether a particular inner bootstrap method succeeds on each outer draw.
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


def _scenario_mean_t_interval(values: pd.Series, confidence: float = 0.95) -> tuple[float, float, float, int]:
    """Student-t interval across independent parent scenarios, never pooled draws."""
    x = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    x = x[np.isfinite(x)]
    n = int(x.size)
    if n == 0:
        return float("nan"), float("nan"), float("nan"), 0
    mean = float(np.mean(x))
    if n < 2:
        return mean, float("nan"), float("nan"), n
    se = float(np.std(x, ddof=1) / np.sqrt(n))
    critical = float(student_t.ppf(0.5 + confidence / 2.0, df=n - 1))
    return mean, mean - critical * se, mean + critical * se, n


def summarize_bootstrap_calibration(records: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return scenario-level and scenario-balanced method/size summaries.

    Input grain is one row per (scenario, validation_size, outer_draw, method).
    Outer estimates are shared across methods; inner bootstrap success/failure
    is method-specific. Failed outer draws remain in the outer denominator,
    while inner failures remain in the bootstrap denominator.
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
    size_values = pd.to_numeric(frame["validation_size"], errors="coerce")
    if (
        size_values.isna().any() or (size_values < 1).any()
        or not np.equal(size_values, np.floor(size_values)).all()
    ):
        raise ValueError("validation_size must contain positive integer values")
    frame["validation_size"] = size_values.astype(int)
    draw_values = pd.to_numeric(frame["outer_draw"], errors="coerce")
    if draw_values.isna().any() or (draw_values < 0).any() or not np.equal(draw_values, np.floor(draw_values)).all():
        raise ValueError("outer_draw must contain nonnegative integer identifiers")
    if frame["method"].astype(str).str.strip().str.len().eq(0).any():
        raise ValueError("method labels must be nonempty or whitespace-only")
    if frame.duplicated(key).any():
        raise ValueError("records must be unique at scenario/size/draw/method grain")
    if not frame["outer_status"].isin(["success", "failed"]).all():
        raise ValueError("outer_status must be success or failed")
    if not frame["bootstrap_status"].isin(["success", "failed", "not_run"]).all():
        raise ValueError("bootstrap_status must be success, failed, or not_run")

    numeric = [
        "outer_posterior_mean", "bootstrap_se", "interval_low",
        "interval_high", "finite_parent_target",
    ]
    for column in numeric:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    if frame["finite_parent_target"].isna().any() or not np.isfinite(
        frame["finite_parent_target"].to_numpy(dtype=float)
    ).all():
        raise ValueError("finite_parent_target must be finite for every attempted outer draw")
    outer_ok = frame["outer_status"].eq("success")
    boot_ok = frame["bootstrap_status"].eq("success")
    if frame.loc[outer_ok, ["outer_posterior_mean", "finite_parent_target"]].isna().any().any():
        raise ValueError("successful outer draws require a posterior mean and finite-parent target")
    if frame.loc[~outer_ok, "outer_posterior_mean"].notna().any():
        raise ValueError("failed outer draws must not carry an outer posterior mean")
    if frame.loc[boot_ok, ["bootstrap_se", "interval_low", "interval_high"]].isna().any().any():
        raise ValueError("successful bootstrap records require SE and interval endpoints")
    if frame.loc[~boot_ok, ["bootstrap_se", "interval_low", "interval_high"]].notna().any().any():
        raise ValueError("failed or not-run bootstrap records must not carry successful interval metrics")
    finite_columns = ["outer_posterior_mean", "bootstrap_se", "interval_low", "interval_high", "finite_parent_target"]
    if not np.isfinite(frame.loc[outer_ok, ["outer_posterior_mean", "finite_parent_target"]].to_numpy(dtype=float)).all():
        raise ValueError("successful outer metrics must be finite")
    if not np.isfinite(frame.loc[boot_ok, ["bootstrap_se", "interval_low", "interval_high"]].to_numpy(dtype=float)).all():
        raise ValueError("successful bootstrap metrics must be finite")
    if (frame.loc[boot_ok, "bootstrap_se"] < 0).any():
        raise ValueError("bootstrap_se must be nonnegative")
    if (frame.loc[boot_ok, "interval_low"] > frame.loc[boot_ok, "interval_high"]).any():
        raise ValueError("interval_low cannot exceed interval_high")
    if (frame.loc[~outer_ok, "bootstrap_status"] != "not_run").any():
        raise ValueError("bootstrap must be not_run when the outer estimator failed")
    if (frame.loc[outer_ok, "bootstrap_status"] == "not_run").any():
        raise ValueError("bootstrap status must be success or failed after a successful outer fit")

    # Every scenario/size cell must contain the same declared method set
    # on every outer draw. Otherwise method-specific missing rows could silently
    # change the empirical comparison population.
    for (scenario, size), cell in frame.groupby(["scenario", "validation_size"], sort=True):
        expected_methods = set(cell["method"].astype(str))
        method_sets = cell.groupby("outer_draw")["method"].apply(lambda values: set(values.astype(str)))
        if any(method_set != expected_methods for method_set in method_sets):
            raise ValueError(
                f"every outer draw must contain the complete method set for scenario={scenario}, size={size}"
            )

    # Every method must report the same outer result/status for a given draw.
    outer_consistency = frame.groupby(
        ["scenario", "validation_size", "outer_draw"], dropna=False
    ).agg(
        outer_status_count=("outer_status", "nunique"),
        outer_mean_count=("outer_posterior_mean", "nunique"),
        target_count=("finite_parent_target", "nunique"),
    )
    if (outer_consistency["outer_status_count"] > 1).any() or (outer_consistency["outer_mean_count"] > 1).any():
        raise ValueError("outer status and posterior mean must be identical across bootstrap methods")
    if (outer_consistency["target_count"] > 1).any():
        raise ValueError("finite-parent target must be fixed within each outer draw")

    # The finite-frame target is computed from the full parent frame and
    # therefore cannot depend on validation sample size or outer draw.
    if frame.groupby("scenario")["finite_parent_target"].nunique(dropna=True).gt(1).any():
        raise ValueError("finite_parent_target must be fixed within each scenario across validation sizes")

    # Deduplicate outer draws before computing empirical SRSWOR variability.
    outer = frame.drop_duplicates(["scenario", "validation_size", "outer_draw"])
    outer_rows = []
    for (scenario, size), group in outer.groupby(["scenario", "validation_size"], sort=True):
        estimates = group.loc[group["outer_status"].eq("success"), "outer_posterior_mean"].to_numpy(dtype=float)
        outer_rows.append({
            "scenario": scenario,
            "validation_size": int(size),
            "n_outer_attempted": int(len(group)),
            "n_outer_successful": int(estimates.size),
            "outer_failure_fraction": float(1.0 - estimates.size / len(group)),
            "empirical_srswor_sd": float(np.std(estimates, ddof=1)) if estimates.size >= 2 else float("nan"),
        })
    outer_summary = pd.DataFrame(outer_rows)

    rows: list[dict] = []
    for (scenario, size, method), group in frame.groupby(
        ["scenario", "validation_size", "method"], sort=True, dropna=False
    ):
        boot = group.loc[group["bootstrap_status"].eq("success")]
        outer_info = outer_summary.loc[
            outer_summary["scenario"].eq(scenario)
            & outer_summary["validation_size"].eq(size)
        ].iloc[0]
        empirical_sd = float(outer_info["empirical_srswor_sd"])
        mean_se = float(boot["bootstrap_se"].mean()) if not boot.empty else float("nan")
        ratio = mean_se / empirical_sd if np.isfinite(empirical_sd) and empirical_sd > 0 else float("nan")
        covered = (
            (boot["interval_low"] <= boot["finite_parent_target"])
            & (boot["finite_parent_target"] <= boot["interval_high"])
        )
        n_outer_success = int(outer_info["n_outer_successful"])
        n_boot_attempted = n_outer_success
        rows.append({
            "scenario": scenario,
            "validation_size": int(size),
            "method": str(method),
            "n_outer_attempted": int(outer_info["n_outer_attempted"]),
            "n_outer_successful": n_outer_success,
            "outer_failure_fraction": float(outer_info["outer_failure_fraction"]),
            "n_bootstrap_attempted": n_boot_attempted,
            "n_bootstrap_successful": int(len(boot)),
            "bootstrap_failure_fraction": float(1.0 - len(boot) / n_boot_attempted) if n_boot_attempted else float("nan"),
            "empirical_srswor_sd": empirical_sd,
            "mean_inner_bootstrap_se": mean_se,
            "mean_se_to_empirical_sd_ratio": ratio,
            "finite_parent_target_interval_coverage": float(covered.mean()) if len(boot) else float("nan"),
            "mean_interval_width": float((boot["interval_high"] - boot["interval_low"]).mean()) if len(boot) else float("nan"),
            "n_intervals_evaluated": int(len(boot)),
        })
    scenario_level = pd.DataFrame(rows)

    aggregate_rows: list[dict] = []
    for (size, method), group in scenario_level.groupby(["validation_size", "method"], sort=True):
        ratio = group["mean_se_to_empirical_sd_ratio"].dropna()
        coverage = group["finite_parent_target_interval_coverage"].dropna()
        widths = group["mean_interval_width"].dropna()
        ratio_mean, ratio_low, ratio_high, ratio_n = _scenario_mean_t_interval(ratio)
        coverage_mean, coverage_low, coverage_high, coverage_n = _scenario_mean_t_interval(coverage)
        width_mean, width_low, width_high, width_n = _scenario_mean_t_interval(widths)
        outer_fail_mean, outer_fail_low, outer_fail_high, outer_fail_n = _scenario_mean_t_interval(
            group["outer_failure_fraction"]
        )
        inner_fail_mean, inner_fail_low, inner_fail_high, inner_fail_n = _scenario_mean_t_interval(
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
            "scenario_balanced_sd_se_ratio": float(ratio.std(ddof=1)) if ratio.size > 1 else float("nan"),
            "scenario_balanced_mean_coverage": coverage_mean,
            "scenario_balanced_mean_coverage_ci95_low": coverage_low,
            "scenario_balanced_mean_coverage_ci95_high": coverage_high,
            "n_scenarios_with_coverage": int(coverage_n),
            "scenario_balanced_sd_coverage": float(coverage.std(ddof=1)) if coverage.size > 1 else float("nan"),
            "scenario_balanced_mean_interval_width": width_mean,
            "scenario_balanced_mean_interval_width_ci95_low": width_low,
            "scenario_balanced_mean_interval_width_ci95_high": width_high,
            "n_scenarios_with_interval_width": int(width_n),
            "scenario_balanced_outer_failure_fraction": outer_fail_mean,
            "scenario_balanced_outer_failure_fraction_ci95_low": outer_fail_low,
            "scenario_balanced_outer_failure_fraction_ci95_high": outer_fail_high,
            "scenario_balanced_bootstrap_failure_fraction": inner_fail_mean,
            "scenario_balanced_bootstrap_failure_fraction_ci95_low": inner_fail_low,
            "scenario_balanced_bootstrap_failure_fraction_ci95_high": inner_fail_high,
        })
    return scenario_level, pd.DataFrame(aggregate_rows)
