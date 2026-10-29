"""Paired grid and KDE-bandwidth sensitivity summary for the selection pilot."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from audit_selection_misspecification import audit_artifacts

ARMS = (
    "analytic_truth", "oracle", "catalogue_only", "correct_weighting",
    "shifted_transition", "shallow_slope",
)
KEYS = ["replicate", "arm"]


def _load_run(folder: Path) -> pd.DataFrame:
    report = audit_artifacts(folder)
    if report.get("status") != "PASS":
        raise ValueError(f"Artifact audit did not pass for {folder}: {report}")
    frame = pd.read_csv(folder / "replicate_metrics.csv")
    if frame.duplicated(KEYS).any():
        raise ValueError(f"Duplicate replicate/arm rows in {folder}")
    counts = frame.groupby("replicate")["arm"].agg(list)
    if not all(set(arms) == set(ARMS) and len(arms) == len(ARMS) for arms in counts):
        raise ValueError(f"Incomplete six-arm pairing in {folder}")
    return frame


def _bootstrap_mean_interval(values: np.ndarray, *, seed: int, resamples: int) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size < 2 or not np.all(np.isfinite(values)):
        raise ValueError("Bootstrap input must be a finite one-dimensional paired vector")
    rng = np.random.default_rng(seed)
    draws = np.empty(resamples, dtype=float)
    chunk_size = 1000
    for start in range(0, resamples, chunk_size):
        stop = min(start + chunk_size, resamples)
        indices = rng.integers(0, values.size, size=(stop - start, values.size))
        draws[start:stop] = values[indices].mean(axis=1)
    low, high = np.quantile(draws, [0.025, 0.975])
    return float(low), float(high)


def summarize_paired_frames(
    baseline: pd.DataFrame,
    comparison: pd.DataFrame,
    *,
    name: str,
    seed: int = 20261031,
    resamples: int = 10000,
) -> tuple[list[dict], dict]:
    """Summarize paired changes without treating events or grid points as replicates."""
    for label, frame in (("baseline", baseline), ("comparison", comparison)):
        if frame.duplicated(KEYS).any():
            raise ValueError(f"{label} contains duplicate replicate/arm rows")
        if set(frame["arm"]) != set(ARMS):
            raise ValueError(f"{label} does not contain all six arms")
    merged = baseline.merge(
        comparison, on=KEYS, how="outer", suffixes=("_baseline", "_comparison"),
        indicator=True, validate="one_to_one",
    )
    if not merged["_merge"].eq("both").all():
        raise ValueError(f"Replicate/arm pairing differs for sensitivity configuration {name}")
    for column in ("seed", "h0_true", "n_parent", "n_catalogue_observed", "n_detected"):
        left, right = f"{column}_baseline", f"{column}_comparison"
        if not np.array_equal(merged[left].to_numpy(), merged[right].to_numpy()):
            raise ValueError(f"Paired design mismatch in {column} for {name}")

    arm_rows = []
    for arm in ARMS:
        group = merged.loc[merged["arm"] == arm].sort_values("replicate")
        mean_delta = group["posterior_mean_comparison"] - group["posterior_mean_baseline"]
        low_delta = group["ci68_low_comparison"] - group["ci68_low_baseline"]
        high_delta = group["ci68_high_comparison"] - group["ci68_high_baseline"]
        width_delta = group["ci68_width_comparison"] - group["ci68_width_baseline"]
        ise_delta = group["density_ise_comparison"] - group["density_ise_baseline"]
        truth = group["h0_true_baseline"].to_numpy(dtype=float)
        base_error = group["posterior_mean_baseline"].to_numpy(dtype=float) - truth
        comp_error = group["posterior_mean_comparison"].to_numpy(dtype=float) - truth
        base_covered = group["covered_68_baseline"].astype(bool).to_numpy()
        comp_covered = group["covered_68_comparison"].astype(bool).to_numpy()
        arm_rows.append({
            "configuration": name,
            "arm": arm,
            "n_paired_replicates": int(len(group)),
            "baseline_rmse": float(np.sqrt(np.mean(base_error ** 2))),
            "comparison_rmse": float(np.sqrt(np.mean(comp_error ** 2))),
            "rmse_change": float(np.sqrt(np.mean(comp_error ** 2)) - np.sqrt(np.mean(base_error ** 2))),
            "baseline_coverage": float(np.mean(base_covered)),
            "comparison_coverage": float(np.mean(comp_covered)),
            "coverage_change_percentage_points": float(100 * (np.mean(comp_covered) - np.mean(base_covered))),
            "mean_posterior_mean_change": float(mean_delta.mean()),
            "p95_abs_posterior_mean_change": float(np.quantile(np.abs(mean_delta), 0.95)),
            "max_abs_posterior_mean_change": float(np.max(np.abs(mean_delta))),
            "p95_abs_ci_endpoint_change": float(np.quantile(np.abs(np.concatenate([low_delta, high_delta])), 0.95)),
            "max_abs_ci_endpoint_change": float(np.max(np.abs(np.concatenate([low_delta, high_delta])))),
            "mean_ci_width_change": float(width_delta.mean()),
            "mean_density_ise_change": float(ise_delta.mean()),
        })

    def primary(frame: pd.DataFrame) -> pd.Series:
        means = frame.pivot(index="replicate", columns="arm", values="posterior_mean").sort_index()
        truths = frame.groupby("replicate")["h0_true"].first().sort_index()
        return (means["catalogue_only"] - truths).abs() - (means["correct_weighting"] - truths).abs()

    base_d = primary(baseline)
    comp_d = primary(comparison)
    if not base_d.index.equals(comp_d.index):
        raise ValueError(f"Primary paired contrast has mismatched replicates for {name}")
    delta_d = comp_d.to_numpy() - base_d.to_numpy()
    base_ci = _bootstrap_mean_interval(base_d.to_numpy(), seed=seed, resamples=resamples)
    comp_ci = _bootstrap_mean_interval(comp_d.to_numpy(), seed=seed + 1, resamples=resamples)
    delta_ci = _bootstrap_mean_interval(delta_d, seed=seed + 2, resamples=resamples)
    primary_row = {
        "configuration": name,
        "estimand": "abs_error_catalogue_only_minus_abs_error_correct_weighting",
        "n_paired_replicates": int(len(base_d)),
        "baseline_mean_contrast": float(base_d.mean()),
        "baseline_bootstrap_95_low": base_ci[0],
        "baseline_bootstrap_95_high": base_ci[1],
        "comparison_mean_contrast": float(comp_d.mean()),
        "comparison_bootstrap_95_low": comp_ci[0],
        "comparison_bootstrap_95_high": comp_ci[1],
        "contrast_change_vs_baseline": float(delta_d.mean()),
        "contrast_change_bootstrap_95_low": delta_ci[0],
        "contrast_change_bootstrap_95_high": delta_ci[1],
        "interpretation": "positive primary contrast favors correct weighting over catalogue-only on absolute H0 error in this synthetic design",
    }
    return arm_rows, primary_row


def compare_runs(
    baseline_dir: Path,
    comparisons: dict[str, Path],
    output_dir: Path,
    *,
    seed: int = 20261031,
    resamples: int = 10000,
) -> dict:
    baseline = _load_run(baseline_dir)
    arm_rows: list[dict] = []
    primary_rows: list[dict] = []
    for index, (name, folder) in enumerate(comparisons.items()):
        comparison = _load_run(folder)
        rows, primary = summarize_paired_frames(
            baseline, comparison, name=name, seed=seed + index * 10, resamples=resamples
        )
        for row in rows:
            if name.startswith("h0_") or name == "z_1601":
                row["numerical_grid_gate"] = (
                    "PASS"
                    if row["max_abs_posterior_mean_change"] <= 0.1
                    and row["max_abs_ci_endpoint_change"] <= 0.2
                    else "NO_GO_NUMERICAL_SENSITIVITY_EXCEEDS_PRESET_TOLERANCE"
                )
            else:
                row["numerical_grid_gate"] = "NOT_APPLICABLE_BANDWIDTH_SENSITIVITY"
        arm_rows.extend(rows)
        primary_rows.append(primary)

    output_dir.mkdir(parents=True, exist_ok=True)
    arm_path = output_dir / "configuration_arm_comparisons.csv"
    primary_path = output_dir / "primary_contrast_comparisons.csv"
    pd.DataFrame(arm_rows).to_csv(arm_path, index=False)
    pd.DataFrame(primary_rows).to_csv(primary_path, index=False)
    result = {
        "status": "SYNTHETIC_SENSITIVITY_AUDIT_PASS",
        "baseline": str(baseline_dir),
        "comparisons": {name: str(path) for name, path in comparisons.items()},
        "replicates_per_configuration": int(baseline["replicate"].nunique()),
        "arms": list(ARMS),
        "bootstrap_resamples": int(resamples),
        "bootstrap_seed_start": int(seed),
        "grid_acceptance_thresholds": {
            "max_abs_posterior_mean_change_km_s_mpc": 0.1,
            "max_abs_68_percent_interval_endpoint_change_km_s_mpc": 0.2,
            "applies_to": ["h0_321", "h0_641", "z_1601"],
        },
        "files": {
            arm_path.name: {"sha256": hashlib.sha256(arm_path.read_bytes()).hexdigest()},
            primary_path.name: {"sha256": hashlib.sha256(primary_path.read_bytes()).hexdigest()},
        },
        "scientific_scope": "controlled synthetic toy only; artifact consistency and paired sensitivity, not real LVK validation",
    }
    summary_path = output_dir / "sensitivity_summary.json"
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def _parse_comparison(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("comparison must be NAME=PATH")
    name, path = value.split("=", 1)
    if not name or not path:
        raise argparse.ArgumentTypeError("comparison must be NAME=PATH")
    return name, Path(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--comparison", action="append", type=_parse_comparison, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--bootstrap-resamples", type=int, default=10000)
    args = parser.parse_args()
    if args.bootstrap_resamples < 100:
        parser.error("--bootstrap-resamples must be at least 100")
    comparisons = dict(args.comparison)
    if len(comparisons) != len(args.comparison):
        parser.error("comparison names must be unique")
    result = compare_runs(
        args.baseline, comparisons, args.output_dir,
        seed=args.seed, resamples=args.bootstrap_resamples,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
