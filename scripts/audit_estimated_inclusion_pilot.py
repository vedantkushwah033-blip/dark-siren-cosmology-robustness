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


def audit_pilot_artifacts(output_dir: str | Path) -> dict:
    root = Path(output_dir)
    csv_path, json_path = root / "replicates.csv", root / "summary.json"
    if not csv_path.is_file() or not json_path.is_file():
        raise FileNotFoundError("pilot output must contain replicates.csv and summary.json")
    frame = pd.read_csv(csv_path)
    summary = json.loads(json_path.read_text(encoding="utf-8"))
    required = {
        "replicate", "replicate_seed", "catalogue_seed", "audit_seed", "arm",
        "validation_size", "status", "error", "n_catalogue", "n_detected",
    }
    missing = required - set(frame.columns)
    if missing:
        raise AssertionError(f"replicate table missing columns: {sorted(missing)}")
    if frame.empty or frame.duplicated(["replicate", "arm"]).any():
        raise AssertionError("replicate table is empty or has duplicate replicate/arm rows")
    if not set(frame["status"]).issubset({"success", "failed"}):
        raise AssertionError("unknown fit status in replicate table")
    success = frame["status"].eq("success")
    if not frame.loc[~success, "error"].fillna("").str.len().gt(0).all():
        raise AssertionError("failed rows must preserve a nonempty error reason")
    numeric = [
        "posterior_mean", "ci68_low", "ci68_high", "ci68_width",
        "absolute_error", "density_ise", "weight_ess", "weight_ess_fraction",
        "max_inverse_weight",
    ]
    for column in numeric:
        if column not in frame.columns:
            raise AssertionError(f"replicate table missing numerical column {column}")
        values = frame.loc[success, column].to_numpy(dtype=float)
        if values.size and not np.all(np.isfinite(values)):
            raise AssertionError(f"successful rows contain non-finite values in {column}")
    if (frame.loc[success, "ci68_low"] > frame.loc[success, "ci68_high"]).any():
        raise AssertionError("credible interval endpoints are reversed")
    if (frame.loc[success, "absolute_error"] < 0).any():
        raise AssertionError("absolute posterior-mean errors cannot be negative")

    arms_in_rows = set(frame["arm"])
    arms_in_summary = {item["arm"] for item in summary.get("arm_summaries", [])}
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
            continue
        h0_true = float(summary.get("h0_true", 70.0))
        error = group["posterior_mean"].to_numpy(dtype=float) - h0_true
        covered = group["covered_68"].astype(bool).to_numpy()
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
    baseline = successful.loc[successful["arm"] == "catalogue_only"].set_index("replicate")
    stored_contrasts = {
        item["comparison_arm"]: item for item in summary.get("contrasts", [])
    }
    expected_comparisons = [arm for arm in arms_in_rows if arm != "catalogue_only"]
    for arm in expected_comparisons:
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

    return {
        "status": "PASS",
        "n_rows": int(len(frame)),
        "n_success": int(success.sum()),
        "n_failed": int((~success).sum()),
        "n_arms": int(len(arms_in_rows)),
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
