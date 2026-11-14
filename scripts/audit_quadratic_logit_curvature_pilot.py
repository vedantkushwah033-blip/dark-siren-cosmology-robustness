"""Cross-gamma auditor and paired summary for quadratic-logit density pilots."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

from audit_absolute_bandwidth_density_replicates import audit as audit_replicates

DEFAULT_GAMMAS = (-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0)
COMMON_CONFIG_KEYS = (
    "master_seed", "seed_stride", "replicates", "n_parent",
    "validation_frame", "validation_sizes", "beta", "quadrature_order",
    "grid_points", "absolute_widths", "scott_multiplier",
)


def _close(a, b, label: str, atol: float = 2e-10) -> None:
    if not np.allclose(
        np.asarray(a, dtype=float), np.asarray(b, dtype=float),
        rtol=2e-10, atol=atol, equal_nan=False,
    ):
        diff = np.max(np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)))
        raise AssertionError(f"{label} mismatch; max abs difference={diff}")


def _same_frame(reference: pd.DataFrame, candidate: pd.DataFrame,
                exact_columns: tuple[str, ...], numeric_columns: tuple[str, ...],
                label: str) -> None:
    if len(reference) != len(candidate):
        raise AssertionError(f"{label} row count differs across gamma")
    for column in exact_columns:
        if not reference[column].equals(candidate[column]):
            raise AssertionError(f"{label} {column} differs across gamma")
    for column in numeric_columns:
        _close(reference[column].to_numpy(float), candidate[column].to_numpy(float),
               f"{label} {column}")


def _paired_contrasts(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    estimator_arms = (
        "A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
        "D_map_plugin_ipw", "E_posterior_mean_density",
        "E_posterior_weighted_node_metric",
    )
    references = (
        "A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
    )
    primary = metrics[metrics.node == -1]
    for (gamma, n, width), group in primary.groupby(
        ["gamma", "validation_size", "absolute_bandwidth_sd"], sort=True
    ):
        wide = group.pivot(index="replicate", columns="arm", values="iae")
        for arm in estimator_arms:
            if arm not in wide:
                continue
            for reference in references:
                if reference == arm or reference not in wide:
                    continue
                difference = (wide[arm] - wide[reference]).dropna().to_numpy(float)
                count = int(difference.size)
                if count < 2:
                    continue
                mean = float(np.mean(difference))
                sd = float(np.std(difference, ddof=1))
                se = sd / np.sqrt(count)
                critical = float(t.ppf(0.975, count - 1))
                rows.append({
                    "gamma": float(gamma),
                    "validation_size": int(n),
                    "absolute_bandwidth_sd": float(width),
                    "arm": arm,
                    "reference_arm": reference,
                    "contrast": "arm_minus_reference_iae",
                    "replicate_count": count,
                    "mean_paired_difference": mean,
                    "sd_paired_difference": sd,
                    "standard_error": se,
                    "ci95_lower": mean - critical * se,
                    "ci95_upper": mean + critical * se,
                    "interval_method": "paired t interval; exploratory pilot only",
                })
    return pd.DataFrame(rows)


def audit(root: Path, expected_gammas=DEFAULT_GAMMAS,
          output_dir: Path | None = None) -> dict:
    root = root.resolve()
    expected_gammas = tuple(float(x) for x in expected_gammas)
    if len(expected_gammas) < 2 or len(set(expected_gammas)) != len(expected_gammas):
        raise ValueError("expected_gammas must contain at least two unique values")

    artifact_dirs = [p for p in root.iterdir() if p.is_dir() and (p / "config.json").is_file()]
    by_gamma = {}
    configs = {}
    reports = {}
    for path in artifact_dirs:
        cfg = json.loads((path / "config.json").read_text())
        if cfg.get("catalogue_selection_model") != "quadratic_logit":
            continue
        gamma = float(cfg["gamma"])
        if gamma in by_gamma:
            raise AssertionError(f"duplicate artifact for gamma={gamma}")
        by_gamma[gamma] = path
        configs[gamma] = cfg
        reports[gamma] = audit_replicates(path)
        if reports[gamma]["status"] != "PASS":
            raise AssertionError(f"replicate audit failed for gamma={gamma}")

    if set(by_gamma) != set(expected_gammas):
        raise AssertionError(
            f"gamma artifact set mismatch: expected {sorted(expected_gammas)}, "
            f"found {sorted(by_gamma)}"
        )

    reference_gamma = expected_gammas[0]
    reference_cfg = configs[reference_gamma]
    for gamma in expected_gammas:
        cfg = configs[gamma]
        for key in COMMON_CONFIG_KEYS:
            if cfg.get(key) != reference_cfg.get(key):
                raise AssertionError(f"common design field {key} differs at gamma={gamma}")
        if not np.isclose(float(cfg["gamma"]), gamma):
            raise AssertionError(f"gamma config mismatch for artifact {by_gamma[gamma]}")
        if cfg.get("catalogue_floor") is not None:
            raise AssertionError("quadratic-logit pilot must not carry a floor-logistic parameter")
        if not np.isclose(float(cfg["beta"]), float(reference_cfg["beta"])):
            raise AssertionError("linear logit coefficient differs across gamma")
        # Each child auditor independently verifies the parent centering and
        # parent-average calibration from the analytic parent measure.
        child_calibrations = []
        for rid in range(int(cfg["replicates"])):
            child = by_gamma[gamma] / "replicates" / f"replicate_{rid:04d}"
            child_cfg = json.loads((child / "config.json").read_text())
            child_calibrations.append(child_cfg["quadratic_logit_calibration"])
            if not np.isclose(float(child_cfg["quadratic_logit_calibration"]["gamma"]), gamma):
                raise AssertionError(f"child gamma mismatch at gamma={gamma}, replicate={rid}")
        first_cal = child_calibrations[0]
        for cal in child_calibrations[1:]:
            for key in ("alpha", "beta", "gamma", "target_rate", "achieved_rate",
                        "x2_center", "coordinate_center", "coordinate_scale"):
                _close([float(cal[key])], [float(first_cal[key])],
                       f"within-gamma calibration {key} gamma={gamma}")
        if not np.isclose(float(first_cal["achieved_rate"]), float(first_cal["target_rate"]), atol=2e-10):
            raise AssertionError(f"parent-average inclusion rate not preserved at gamma={gamma}")

    # Verify that the same random-number streams induce paired parent and
    # validation-frame redshifts/IDs across curvature values. Retention labels
    # are expected to differ because the inclusion probability changes.
    reference_path = by_gamma[reference_gamma]
    n_replicates = int(reference_cfg["replicates"])
    for gamma in expected_gammas[1:]:
        candidate_path = by_gamma[gamma]
        for rid in range(n_replicates):
            ref_child = reference_path / "replicates" / f"replicate_{rid:04d}"
            cand_child = candidate_path / "replicates" / f"replicate_{rid:04d}"
            ref_parent = pd.read_csv(ref_child / "truth" / "parent_truth.csv")
            cand_parent = pd.read_csv(cand_child / "truth" / "parent_truth.csv")
            _same_frame(ref_parent, cand_parent, ("object_id",), ("redshift",),
                        f"parent paired sample replicate={rid} gamma={gamma}")
            ref_vf = pd.read_csv(ref_child / "truth" / "validation_frame_truth.csv")
            cand_vf = pd.read_csv(cand_child / "truth" / "validation_frame_truth.csv")
            _same_frame(ref_vf, cand_vf,
                        ("frame_id", "selected_for_validation"), ("redshift",),
                        f"validation frame paired sample replicate={rid} gamma={gamma}")
            ref_v = pd.read_csv(ref_child / "observable" / "validation_max.csv")
            cand_v = pd.read_csv(cand_child / "observable" / "validation_max.csv")
            _same_frame(ref_v, cand_v, ("validation_id",), ("redshift",),
                        f"validation observed paired sample replicate={rid} gamma={gamma}")

    metric_frames = []
    for gamma in expected_gammas:
        frame = pd.read_csv(by_gamma[gamma] / "density_metrics_all_replicates.csv")
        frame["gamma"] = gamma
        metric_frames.append(frame)
    metrics = pd.concat(metric_frames, ignore_index=True)
    summary = (
        metrics.groupby(
            ["gamma", "validation_size", "absolute_bandwidth_sd", "arm", "node"],
            as_index=False,
        )
        .agg(
            count=("iae", "count"),
            mean_iae=("iae", "mean"),
            sd_iae=("iae", "std"),
            median_iae=("iae", "median"),
            mean_ise=("ise", "mean"),
            sd_ise=("ise", "std"),
            mean_boundary_iae=("boundary_iae", "mean"),
        )
    )
    contrasts = _paired_contrasts(metrics)

    report = {
        "status": "PASS",
        "gammas": list(expected_gammas),
        "replicates_per_gamma": int(reference_cfg["replicates"]),
        "paired_seed_schedule_verified": True,
        "paired_parent_and_validation_redshifts_verified": True,
        "gamma_artifact_audits": {
            str(gamma): reports[gamma] for gamma in expected_gammas
        },
        "metric_rows": int(len(metrics)),
        "summary_rows": int(len(summary)),
        "paired_contrast_rows": int(len(contrasts)),
        "scope": (
            "Cross-gamma config, calibration, paired-seed, parent/validation "
            "sample, and per-gamma raw-input audits. Density-only synthetic "
            "stress pilot; no bandwidth selection or H0 inference."
        ),
    }
    if output_dir is not None:
        output_dir = output_dir.resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        summary.to_csv(output_dir / "curvature_density_iae_summary.csv", index=False)
        contrasts.to_csv(output_dir / "curvature_paired_iae_contrasts.csv", index=False)
        (output_dir / "audit.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_root", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--expected-gammas", type=lambda x: tuple(float(v) for v in x.split(",")),
                        default=DEFAULT_GAMMAS)
    args = parser.parse_args()
    print(json.dumps(audit(args.artifact_root, args.expected_gammas, args.output_dir),
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
