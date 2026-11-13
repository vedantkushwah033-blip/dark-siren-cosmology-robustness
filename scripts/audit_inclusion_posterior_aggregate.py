"""Independently reconstruct aggregate calibration summaries from the replicate ledger.

This audit deliberately does not import the calibration runner. It recomputes
failure-inclusive bounds and paired validation-size contrasts from the ledger,
then compares them with the stored summary. It complements (but does not replace)
the per-replicate raw-artifact auditor.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _close(actual: Any, expected: Any, *, atol: float = 1e-12) -> bool:
    try:
        return bool(np.isclose(float(actual), float(expected), rtol=1e-10, atol=atol, equal_nan=True))
    except (TypeError, ValueError):
        return False


def _assert_close(actual: Any, expected: Any, label: str) -> None:
    if not _close(actual, expected):
        raise AssertionError(f"{label}: stored={actual!r}, reconstructed={expected!r}")


def _attempt_bounds(frame: pd.DataFrame, indicator: str) -> list[float]:
    attempted = frame.loc[frame["status"] != "NOT_RUN"]
    if attempted.empty:
        return [float("nan"), float("nan")]
    evaluable = attempted.loc[attempted["status"] == "SUCCESS"]
    observed = int(evaluable[indicator].astype(int).sum()) if not evaluable.empty else 0
    unresolved = len(attempted) - len(evaluable)
    return [observed / len(attempted), (observed + unresolved) / len(attempted)]


def _validate_all_attempt_bounds(ledger: pd.DataFrame, summary: dict) -> int:
    checked = 0
    for key, stored in summary["summary_by_cell"].items():
        subset = ledger.loc[
            (ledger["experiment"] == stored["experiment"])
            & np.isclose(ledger["catalogue_floor"].astype(float), float(stored["catalogue_floor"]))
            & (ledger["n_validation"].astype(int) == int(stored["n_validation"]))
        ]
        attempted = subset.loc[subset["status"] != "NOT_RUN"]
        successful = subset.loc[subset["status"] == "SUCCESS"]
        if len(attempted) != int(stored["attempted"]):
            raise AssertionError(f"{key}: attempted denominator disagrees with ledger")
        if len(successful) != int(stored["success"]):
            raise AssertionError(f"{key}: evaluable denominator disagrees with ledger")

        experiment = str(stored["experiment"])
        if experiment == "H0_PIT_FIXED_INCLUSION_TRUTH":
            columns = {
                level: f"h0_pit_interval_included_{level}" for level in (68, 90, 95)
            }
            prefix = "prior_predictive_interval_inclusion"
        elif experiment == "FIXED_TRUTH":
            columns = {level: f"covered_{level}" for level in (68, 90, 95)}
            prefix = "coverage"
        else:
            continue

        for level, column in columns.items():
            field = f"{prefix}_{level}_all_attempt_bounds"
            if field not in stored:
                raise AssertionError(f"{key}: missing stored failure-inclusive field {field}")
            reconstructed = _attempt_bounds(subset, column)
            if len(stored[field]) != 2:
                raise AssertionError(f"{key}: {field} must contain two endpoints")
            _assert_close(stored[field][0], reconstructed[0], f"{key} {field} lower")
            _assert_close(stored[field][1], reconstructed[1], f"{key} {field} upper")
            checked += 1
    return checked


def _paired_metrics(experiment: str) -> dict[str, tuple[str, str | None]]:
    metrics: dict[str, tuple[str, str | None]] = {
        "posterior_median": ("h0_median", None),
        "signed_median_error": ("h0_median", "signed_error"),
        "absolute_median_error": ("h0_median", "absolute_error"),
        "width_68": ("h0_q84", "h0_q16"),
        "width_90": ("h0_q95", "h0_q05"),
        "width_95": ("h0_q975", "h0_q025"),
    }
    if experiment == "H0_PIT_FIXED_INCLUSION_TRUTH":
        metrics["pit"] = ("pit", None)
        metrics.update({
            f"prior_predictive_interval_inclusion_{level}":
                (f"h0_pit_interval_included_{level}", None)
            for level in (68, 90, 95)
        })
    elif experiment == "FIXED_TRUTH":
        metrics.update({
            f"coverage_{level}": (f"covered_{level}", None)
            for level in (68, 90, 95)
        })
    return metrics


def _metric_values(row: pd.Series, column: str, kind: str | None) -> float:
    value = float(row[column])
    if kind == "signed_error":
        value -= float(row["h0_true"])
    elif kind == "absolute_error":
        value = abs(value - float(row["h0_true"]))
    elif kind is not None:
        value -= float(row[kind])
    return value


def _reconstruct_paired_contrasts(ledger: pd.DataFrame) -> dict[str, dict]:
    output: dict[str, dict] = {}
    grouping = ledger.groupby(["experiment", "catalogue_floor"], dropna=False)
    for (experiment_raw, floor_raw), group in grouping:
        experiment = str(experiment_raw)
        floor = float(floor_raw)
        sizes = sorted(int(value) for value in group["n_validation"].dropna().unique())
        metrics = _paired_metrics(experiment)
        for i, low_size in enumerate(sizes):
            for high_size in sizes[i + 1:]:
                low = group.loc[group["n_validation"].astype(int) == low_size]
                high = group.loc[group["n_validation"].astype(int) == high_size]
                if low["replicate"].duplicated().any() or high["replicate"].duplicated().any():
                    raise AssertionError(
                        f"duplicate replicate IDs within {experiment}, floor={floor}, sizes={low_size}/{high_size}"
                    )
                low_by_id = low.set_index("replicate", drop=False)
                high_by_id = high.set_index("replicate", drop=False)
                common = sorted(set(low_by_id.index).intersection(high_by_id.index))
                valid_ids = [
                    rid for rid in common
                    if low_by_id.loc[rid, "status"] == "SUCCESS"
                    and high_by_id.loc[rid, "status"] == "SUCCESS"
                ]
                key = f"{experiment}|floor={floor:g}|{low_size}->{high_size}"
                item = {
                    "experiment": experiment,
                    "catalogue_floor": floor,
                    "smaller_validation_size": low_size,
                    "larger_validation_size": high_size,
                    "attempted_common_replicates": len(common),
                    "evaluable_pairs": len(valid_ids),
                    "excluded_pairs": len(common) - len(valid_ids),
                    "difference_convention": "larger_validation_size_minus_smaller_validation_size",
                    "metrics": {},
                }
                for metric, (column, kind) in metrics.items():
                    deltas = []
                    for rid in valid_ids:
                        lo = low_by_id.loc[rid]
                        hi = high_by_id.loc[rid]
                        delta = _metric_values(hi, column, kind) - _metric_values(lo, column, kind)
                        if np.isfinite(delta):
                            deltas.append(float(delta))
                    n = len(deltas)
                    mean = float(np.mean(deltas)) if n else None
                    se = float(np.std(deltas, ddof=1) / math.sqrt(n)) if n > 1 else None
                    interval = (
                        [mean - 1.959963984540054 * se, mean + 1.959963984540054 * se]
                        if mean is not None and se is not None else None
                    )
                    item["metrics"][metric] = {
                        "n_pairs": n,
                        "mean_delta": mean,
                        "paired_mc_se": se,
                        "normal_approx_95_interval": interval,
                    }
                output[key] = item
    return output


def _validate_paired_contrasts(ledger: pd.DataFrame, summary: dict) -> int:
    stored = summary.get("paired_validation_size_contrasts")
    if not isinstance(stored, dict):
        raise AssertionError("summary lacks paired_validation_size_contrasts")
    reconstructed = _reconstruct_paired_contrasts(ledger)
    if set(stored) != set(reconstructed):
        missing = sorted(set(reconstructed) - set(stored))
        extra = sorted(set(stored) - set(reconstructed))
        raise AssertionError(f"paired contrast keys differ; missing={missing}, extra={extra}")

    checked = 0
    for key, expected in reconstructed.items():
        observed = stored[key]
        for field in (
            "experiment", "smaller_validation_size", "larger_validation_size",
            "attempted_common_replicates", "evaluable_pairs", "excluded_pairs",
            "difference_convention",
        ):
            if observed.get(field) != expected[field]:
                raise AssertionError(
                    f"{key}: paired contrast {field} differs: "
                    f"stored={observed.get(field)!r}, reconstructed={expected[field]!r}"
                )
        _assert_close(observed.get("catalogue_floor"), expected["catalogue_floor"], f"{key} floor")
        if set(observed.get("metrics", {})) != set(expected["metrics"]):
            raise AssertionError(f"{key}: paired metric names differ")
        for metric, values in expected["metrics"].items():
            actual = observed["metrics"][metric]
            for field in ("n_pairs",):
                if actual.get(field) != values[field]:
                    raise AssertionError(f"{key} {metric}: {field} differs")
            for field in ("mean_delta", "paired_mc_se"):
                left, right = actual.get(field), values[field]
                if left is None or right is None:
                    if left is not None or right is not None:
                        raise AssertionError(f"{key} {metric}: {field} nullability differs")
                else:
                    _assert_close(left, right, f"{key} {metric} {field}")
            left_interval = actual.get("normal_approx_95_interval")
            right_interval = values["normal_approx_95_interval"]
            if left_interval is None or right_interval is None:
                if left_interval is not None or right_interval is not None:
                    raise AssertionError(f"{key} {metric}: interval nullability differs")
            else:
                for endpoint, (left, right) in enumerate(zip(left_interval, right_interval)):
                    _assert_close(left, right, f"{key} {metric} interval endpoint {endpoint}")
            checked += 1
    return checked


def audit(root: Path) -> dict:
    root = root.resolve()
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    summary = json.loads((root / "calibration_summary.json").read_text(encoding="utf-8"))
    required = {"cell_id", "replicate", "experiment", "catalogue_floor", "n_validation", "status"}
    missing = required - set(ledger.columns)
    if missing:
        raise AssertionError(f"replicate ledger missing required columns: {sorted(missing)}")
    if ledger["cell_id"].duplicated().any():
        raise AssertionError("replicate ledger contains duplicate cell IDs")
    if ledger["replicate"].isna().any():
        raise AssertionError("replicate IDs must be present")
    if not set(ledger["status"].dropna()).issubset(
        {"SUCCESS", "NON_EVALUABLE", "FAILED", "NOT_RUN"}
    ):
        raise AssertionError("ledger contains an unknown lifecycle status")
    bounds = _validate_all_attempt_bounds(ledger, summary)
    contrasts = _validate_paired_contrasts(ledger, summary)
    return {
        "status": "PASS",
        "ledger_rows": int(len(ledger)),
        "failure_inclusive_bounds_reconstructed": bounds,
        "paired_metric_summaries_reconstructed": contrasts,
        "scope": (
            "Aggregate arithmetic and pairing only. This does not independently "
            "reconstruct replicate-level scientific artifacts or establish calibration."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = audit(args.artifact_dir)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
