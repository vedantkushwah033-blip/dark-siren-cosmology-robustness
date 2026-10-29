"""Independent consistency audit for estimated-inclusion pilot output files.

This auditor reads only persisted CSV/JSON artifacts and independently
reconstructs arm summaries and paired absolute-error contrasts. It does not
import the simulation or aggregation implementation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _wilson(successes: int, trials: int) -> tuple[float, float]:
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("invalid binomial counts")
    z = 1.959963984540054
    p = successes / trials
    denominator = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / denominator
    half = z * np.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / denominator
    return max(0.0, float(center - half)), min(1.0, float(center + half))


def _close(actual: object, expected: object, label: str, *, atol: float = 1e-9) -> None:
    if expected is None:
        if actual is not None:
            raise AssertionError(f"{label}: expected null, got {actual!r}")
        return
    if actual is None or not np.isfinite(float(actual)) or not np.isclose(
        float(actual), float(expected), rtol=1e-9, atol=atol
    ):
        raise AssertionError(f"{label}: stored={actual!r}, reconstructed={expected!r}")



def _strict_bool(value: object, label: str) -> bool:
    """Parse persisted Boolean values without Python's truthiness trap."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1"}:
            return True
        if normalized in {"false", "0"}:
            return False
    raise AssertionError(f"{label}: expected an explicit Boolean/0/1 value, got {value!r}")


def _required_finite_number(value: object, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AssertionError(f"{label}: expected a finite number") from exc
    if not np.isfinite(number):
        raise AssertionError(f"{label}: expected a finite number")
    return number

def _integer_values(values: pd.Series, label: str) -> np.ndarray:
    """Parse integer-valued CSV columns without silently truncating decimals."""
    numeric = pd.to_numeric(values, errors="raise").to_numpy(dtype=float)
    if not np.all(np.isfinite(numeric)) or not np.all(numeric == np.floor(numeric)):
        raise AssertionError(f"{label}: expected finite integer-valued entries")
    return numeric.astype(np.int64)

def audit_pilot_artifacts(output_dir: str | Path) -> dict:
    root = Path(output_dir)
    csv_path, json_path = root / "replicates.csv", root / "summary.json"
    if not csv_path.is_file() or not json_path.is_file():
        raise FileNotFoundError("pilot output must contain replicates.csv and summary.json")
    frame = pd.read_csv(csv_path)
    summary = json.loads(json_path.read_text(encoding="utf-8"))
    required = {
        "replicate", "replicate_seed", "catalogue_seed", "audit_seed", "arm",
        "validation_size", "status", "error", "fit_status", "fit_boundary_limited",
        "n_catalogue", "n_detected",
    }
    missing = required - set(frame.columns)
    if missing:
        raise AssertionError(f"replicate table missing columns: {sorted(missing)}")
    if frame.empty or frame.duplicated(["replicate", "arm"]).any():
        raise AssertionError("replicate table is empty or has duplicate replicate/arm rows")
    if not set(frame["status"]).issubset({"success", "failed"}):
        raise AssertionError("unknown fit status in replicate table")
    success = frame["status"].eq("success")
    failed_errors = frame.loc[~success, "error"].fillna("").astype(str)
    if not failed_errors.str.len().gt(0).all():
        raise AssertionError("failed rows must preserve a nonempty error reason")
    numeric = [
        "posterior_mean", "posterior_median", "ci68_low", "ci68_high", "ci68_width",
        "absolute_error", "density_ise", "weight_ess", "weight_ess_fraction",
        "max_inverse_weight",
    ]
    if success.any():
        for column in [*numeric, "covered_68", "inclusion_probability_mae"]:
            if column not in frame.columns:
                raise AssertionError(f"replicate table missing column {column}")
    for column in numeric:
        if column in frame.columns:
            values = frame.loc[success, column].to_numpy(dtype=float)
            if values.size and not np.all(np.isfinite(values)):
                raise AssertionError(f"successful rows contain non-finite values in {column}")
    if success.any():
        if (frame.loc[success, "ci68_low"] > frame.loc[success, "ci68_high"]).any():
            raise AssertionError("credible interval endpoints are reversed")
        if (frame.loc[success, "absolute_error"] < 0).any():
            raise AssertionError("absolute posterior-mean errors cannot be negative")
    boundary_flags = np.array([_strict_bool(v, "fit_boundary_limited") for v in frame.loc[success, "fit_boundary_limited"]], dtype=bool)
    boundary_status = frame.loc[success, "fit_status"].eq("boundary_limited").to_numpy(dtype=bool)
    if not np.array_equal(boundary_flags, boundary_status):
        raise AssertionError("fit boundary flags disagree with fit_status")

    h0_true = _required_finite_number(summary.get("h0_true"), "summary.h0_true")
    requested_replicates = summary.get("n_replicates_requested")
    if isinstance(requested_replicates, bool) or not isinstance(requested_replicates, int) or requested_replicates < 1:
        raise AssertionError("summary.n_replicates_requested must be a positive integer")
    validation_sizes = summary.get("validation_sizes")
    if not isinstance(validation_sizes, list) or not validation_sizes:
        raise AssertionError("summary.validation_sizes must be a nonempty list")
    if any(isinstance(v, bool) or not isinstance(v, int) or v < 1 for v in validation_sizes):
        raise AssertionError("summary.validation_sizes must contain positive integers")
    expected_arms = {"catalogue_only", "oracle"} | {f"estimated_{v}" for v in validation_sizes}
    arms_in_rows = set(frame["arm"])
    if arms_in_rows != expected_arms:
        raise AssertionError(f"replicate arm set differs from declared design: expected={sorted(expected_arms)}, observed={sorted(arms_in_rows)}")
    expected_replicates = set(range(requested_replicates))
    replicate_values = _integer_values(frame["replicate"], "replicate")
    observed_replicates = set(np.unique(replicate_values).tolist())
    if observed_replicates != expected_replicates:
        raise AssertionError(f"replicate IDs differ from declared design: expected={sorted(expected_replicates)}, observed={sorted(observed_replicates)}")
    expected_pairs = {(replicate, arm) for replicate in expected_replicates for arm in expected_arms}
    observed_pairs = set(zip(replicate_values.tolist(), frame["arm"]))
    if observed_pairs != expected_pairs:
        missing_pairs = sorted(expected_pairs - observed_pairs)[:5]
        extra_pairs = sorted(observed_pairs - expected_pairs)[:5]
        raise AssertionError(f"incomplete replicate-arm matrix; missing={missing_pairs}, extra={extra_pairs}")
    expected_validation_size = {"catalogue_only": 0, "oracle": 0, **{f"estimated_{v}": v for v in validation_sizes}}
    for arm, expected_size in expected_validation_size.items():
        sizes = set(_integer_values(frame.loc[frame["arm"] == arm, "validation_size"], f"{arm}.validation_size").tolist())
        if sizes != {expected_size}:
            raise AssertionError(f"{arm}: validation_size does not match declared design")
    if not isinstance(summary.get("arm_summaries"), list):
        raise AssertionError("summary.arm_summaries must be a list")
    if any(not isinstance(item, dict) or "arm" not in item for item in summary["arm_summaries"]):
        raise AssertionError("each summary arm entry must be an object with an arm name")
    summary_arm_names = [item["arm"] for item in summary["arm_summaries"]]
    if len(summary_arm_names) != len(set(summary_arm_names)):
        raise AssertionError("summary contains duplicate arm summaries")
    arms_in_summary = set(summary_arm_names)
    if arms_in_rows != arms_in_summary:
        raise AssertionError("summary arm set does not match replicate table")
    for replicate, group in frame.groupby("replicate"):
        for seed_column in ("replicate_seed", "catalogue_seed", "audit_seed"):
            if group[seed_column].nunique(dropna=False) != 1:
                raise AssertionError(f"{seed_column} varies across arms in replicate {replicate}")
        if group["n_catalogue"].nunique(dropna=False) != 1 or group["n_detected"].nunique(dropna=False) != 1:
            raise AssertionError(f"catalogue/event counts differ across arms in replicate {replicate}")

    summary_by_arm = {item["arm"]: item for item in summary["arm_summaries"]}
    for arm in arms_in_rows:
        all_rows = frame.loc[frame["arm"] == arm]
        group = all_rows.loc[all_rows["status"] == "success"]
        stored = summary_by_arm[arm]
        _close(stored["n_attempted"], len(all_rows), f"{arm}.n_attempted")
        _close(stored["n_success"], len(group), f"{arm}.n_success")
        _close(stored["n_failed"], len(all_rows) - len(group), f"{arm}.n_failed")
        if group.empty:
            if stored.get("status") != "no_successful_fits":
                raise AssertionError(f"{arm}: zero-success arm must be marked no_successful_fits")
            continue
        if stored.get("status") == "no_successful_fits":
            raise AssertionError(f"{arm}: summary marks an arm non-evaluable despite successful rows")
        posterior_mean = group["posterior_mean"].to_numpy(dtype=float)
        ci_low = group["ci68_low"].to_numpy(dtype=float)
        ci_high = group["ci68_high"].to_numpy(dtype=float)
        ci_width = group["ci68_width"].to_numpy(dtype=float)
        row_abs_error = group["absolute_error"].to_numpy(dtype=float)
        error = posterior_mean - h0_true
        if not np.allclose(ci_width, ci_high - ci_low, rtol=1e-9, atol=1e-9):
            raise AssertionError(f"{arm}: ci68_width disagrees with interval endpoints")
        if not np.allclose(row_abs_error, np.abs(error), rtol=1e-9, atol=1e-9):
            raise AssertionError(f"{arm}: absolute_error disagrees with posterior_mean and h0_true")
        if (group["density_ise"].to_numpy(dtype=float) < 0).any():
            raise AssertionError(f"{arm}: density_ise cannot be negative")
        ess = group["weight_ess"].to_numpy(dtype=float)
        ess_fraction = group["weight_ess_fraction"].to_numpy(dtype=float)
        n_catalogue = group["n_catalogue"].to_numpy(dtype=float)
        if (ess <= 0).any() or (ess > n_catalogue + 1e-9).any():
            raise AssertionError(f"{arm}: weight_ess must be positive and no greater than n_catalogue")
        if not np.allclose(ess_fraction, ess / n_catalogue, rtol=1e-9, atol=1e-9):
            raise AssertionError(f"{arm}: weight_ess_fraction disagrees with weight_ess / n_catalogue")
        if (group["max_inverse_weight"].to_numpy(dtype=float) < 1.0 - 1e-9).any():
            raise AssertionError(f"{arm}: max_inverse_weight cannot be below one")
        mae = group["inclusion_probability_mae"].to_numpy(dtype=float)
        if arm != "catalogue_only" and not np.all(np.isfinite(mae)):
            raise AssertionError(f"{arm}: inclusion_probability_mae must be finite for successful rows")
        if np.any(np.isfinite(mae) & (mae < 0)):
            raise AssertionError(f"{arm}: inclusion_probability_mae cannot be negative")
        covered = np.array([_strict_bool(v, f"{arm}.covered_68") for v in group["covered_68"]], dtype=bool)
        reconstructed_coverage = ((ci_low <= h0_true) & (h0_true <= ci_high))
        if not np.array_equal(covered, reconstructed_coverage):
            raise AssertionError(f"{arm}: covered_68 flags disagree with interval endpoints and h0_true")
        successes, trials = int(np.sum(covered)), int(len(covered))
        low, high = _wilson(successes, trials)
        reconstructed = {
            "mean_bias": float(np.mean(error)),
            "rmse": float(np.sqrt(np.mean(error**2))),
            "coverage_68": float(np.mean(covered)),
            "coverage_numerator": successes,
            "coverage_denominator": trials,
            "coverage_wilson_95_low": low,
            "coverage_wilson_95_high": high,
            "mean_interval_width_68": float(group["ci68_width"].mean()),
            "mean_density_ise": float(group["density_ise"].mean()),
            "mean_weight_ess_fraction": float(group["weight_ess_fraction"].mean()),
            "mean_max_inverse_weight": float(group["max_inverse_weight"].mean()),
        }
        for key, value in reconstructed.items():
            _close(stored.get(key), value, f"{arm}.{key}")
        mae_values = group["inclusion_probability_mae"].dropna().to_numpy(dtype=float)
        reconstructed_mae = None if mae_values.size == 0 else float(np.mean(mae_values))
        _close(stored.get("mean_inclusion_probability_mae"), reconstructed_mae,
               f"{arm}.mean_inclusion_probability_mae")

    successful = frame.loc[success]
    raw_contrasts = summary.get("contrasts", [])
    if not isinstance(raw_contrasts, list):
        raise AssertionError("summary.contrasts must be a list")
    if any(not isinstance(item, dict) or "comparison_arm" not in item for item in raw_contrasts):
        raise AssertionError("each contrast entry must be an object with comparison_arm")
    contrast_names = [item["comparison_arm"] for item in raw_contrasts]
    if len(contrast_names) != len(set(contrast_names)):
        raise AssertionError("summary contains duplicate comparison contrasts")
    stored_contrasts = {item["comparison_arm"]: item for item in raw_contrasts}
    expected_comparisons = [arm for arm in arms_in_rows if arm != "catalogue_only"]
    unexpected_contrasts = set(stored_contrasts) - set(expected_comparisons)
    if unexpected_contrasts:
        raise AssertionError(f"unexpected stored contrasts: {sorted(unexpected_contrasts)}")
    if successful.empty and stored_contrasts:
        raise AssertionError("all-failed artifact must not contain paired contrasts")
    baseline = (
        successful.loc[successful["arm"] == "catalogue_only"].set_index("replicate")
        if not successful.empty else None
    )
    for arm in expected_comparisons:
        if successful.empty:
            continue
        comparison = successful.loc[successful["arm"] == arm].set_index("replicate")
        paired = baseline[["absolute_error"]].join(
            comparison[["absolute_error"]], how="inner", lsuffix="_baseline", rsuffix="_comparison"
        )
        if paired.empty:
            if arm in stored_contrasts:
                raise AssertionError(f"{arm}: contrast exists but has no paired successful replicates")
            continue
        values = (
            paired["absolute_error_baseline"] - paired["absolute_error_comparison"]
        ).to_numpy(dtype=float)
        stored = stored_contrasts.get(arm)
        if stored is None:
            raise AssertionError(f"{arm}: missing paired contrast")
        if stored.get("reference_arm") != "catalogue_only":
            raise AssertionError(f"{arm}: contrast reference_arm must be catalogue_only")
        _close(stored.get("n_paired_replicates"), values.size, f"{arm}.contrast.n")
        mean = float(np.mean(values))
        _close(stored.get("mean_absolute_error_reduction"), mean, f"{arm}.contrast.mean")
        if values.size < 2:
            se = low = high = None
        else:
            se = float(np.std(values, ddof=1) / np.sqrt(values.size))
            low, high = mean - 1.96 * se, mean + 1.96 * se
        _close(stored.get("standard_error"), se, f"{arm}.contrast.se")
        _close(stored.get("normal_95_low"), low, f"{arm}.contrast.low")
        _close(stored.get("normal_95_high"), high, f"{arm}.contrast.high")

    boundary_success = np.array([
        _strict_bool(v, "fit_boundary_limited")
        for v in frame.loc[success, "fit_boundary_limited"]
    ], dtype=bool)
    per_arm_counts = {}
    for arm in sorted(arms_in_rows):
        arm_rows = frame.loc[frame["arm"] == arm]
        arm_success = arm_rows["status"].eq("success")
        per_arm_counts[arm] = {
            "n_attempted": int(len(arm_rows)),
            "n_success": int(arm_success.sum()),
            "n_failed": int((~arm_success).sum()),
            "failure_fraction": float((~arm_success).mean()),
            "status": "evaluable" if arm_success.any() else "non_evaluable_no_successful_fits",
        }
    n_boundary = int(boundary_success.sum())
    comparison_diagnostics = {}
    for arm in sorted(arms_in_rows - {"catalogue_only"}):
        baseline_all = frame.loc[frame["arm"] == "catalogue_only"].set_index("replicate")
        comparison_all = frame.loc[frame["arm"] == arm].set_index("replicate")
        paired_status = baseline_all[["status"]].join(
            comparison_all[["status"]], how="inner", lsuffix="_baseline", rsuffix="_comparison"
        )
        baseline_ok = paired_status["status_baseline"].eq("success")
        comparison_ok = paired_status["status_comparison"].eq("success")
        comparison_diagnostics[arm] = {
            "n_attempted_common_replicates": int(len(paired_status)),
            "n_paired_successful_replicates": int((baseline_ok & comparison_ok).sum()),
            "n_baseline_success_comparison_failure": int((baseline_ok & ~comparison_ok).sum()),
            "n_baseline_failure_comparison_success": int((~baseline_ok & comparison_ok).sum()),
            "n_both_failed": int((~baseline_ok & ~comparison_ok).sum()),
            "successful_pair_contrast_conditions_on_fit_success": True,
        }
    return {
        "status": "PASS",
        "n_rows": int(len(frame)),
        "n_success": int(success.sum()),
        "n_failed": int((~success).sum()),
        "n_arms": int(len(arms_in_rows)),
        "n_replicates_expected": int(requested_replicates),
        "n_replicate_arm_pairs_expected": int(len(expected_pairs)),
        "n_replicate_arm_pairs_observed": int(len(observed_pairs)),
        "failure_fraction": float((~success).mean()),
        "n_boundary_limited_successful_fits": n_boundary,
        "boundary_limited_fraction_of_successful_fits": (
            None if not success.any() else float(n_boundary / int(success.sum()))
        ),
        "per_arm_counts": per_arm_counts,
        "comparison_diagnostics": comparison_diagnostics,
        "audit_scope": "persisted-table integrity and independent summary reconstruction only",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Independently audit estimated-inclusion pilot artifacts.")
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    report = audit_pilot_artifacts(args.output_dir)
    (args.output_dir / "independent_audit.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
