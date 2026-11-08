"""Independent reconstruction of finite-validation probability/IPW artifacts.

This auditor deliberately does not import the benchmark or estimator. It
verifies the complete parent-frame partition and recomputes all reported scores.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit

KEYS = ["scenario_seed", "replicate", "gamma", "n_validation", "model_degree"]
ROLES = {"validation", "analysis_catalogue", "nonretained_holdout"}
METRICS = [
    "heldout_brier_score", "heldout_log_loss", "heldout_probability_rmse",
    "kish_ess", "entropy_ess", "max_normalized_weight",
    "top_1pct_weight_mass", "top_5pct_weight_mass",
    "oracle_kish_ess", "oracle_entropy_ess", "oracle_max_normalized_weight",
]


def _weight_stats(probability: np.ndarray) -> dict[str, float]:
    p = np.asarray(probability, dtype=float)
    names = ("kish_ess", "entropy_ess", "max_normalized_weight",
             "top_1pct_weight_mass", "top_5pct_weight_mass")
    if p.size == 0:
        return {name: float("nan") for name in names}
    if not np.all(np.isfinite(p)) or np.any((p <= 0) | (p > 1)):
        raise ValueError("invalid inclusion probabilities in raw records")
    scaled_w = float(np.min(p)) / p
    q = scaled_w / np.sum(scaled_w)
    ordered = np.sort(q)[::-1]
    return {
        "kish_ess": float(np.sum(scaled_w) ** 2 / np.sum(scaled_w ** 2)),
        "entropy_ess": float(np.exp(-np.sum(q * np.log(q)))),
        "max_normalized_weight": float(np.max(q)),
        "top_1pct_weight_mass": float(np.sum(ordered[:max(1, int(np.ceil(.01 * len(q))))])),
        "top_5pct_weight_mass": float(np.sum(ordered[:max(1, int(np.ceil(.05 * len(q))))])),
    }


def audit(summary_path: Path, raw_path: Path, tolerance: float = 1e-10) -> dict:
    summary = pd.read_csv(summary_path)
    raw = pd.read_csv(raw_path)
    required_summary = set(KEYS + [
        "fit_status", "n_parent", "n_validation", "n_analysis_catalogue",
        "coefficients_json", "covariance_json", "fit_center", "fit_scale",
    ] + METRICS)
    required_raw = set(KEYS + [
        "split_role", "object_id", "redshift", "retained", "p_true", "p_hat"
    ])
    if not required_summary.issubset(summary.columns):
        raise ValueError(f"summary missing columns: {sorted(required_summary - set(summary.columns))}")
    if not required_raw.issubset(raw.columns):
        raise ValueError(f"raw file missing columns: {sorted(required_raw - set(raw.columns))}")
    if summary.duplicated(KEYS).any():
        raise ValueError("duplicate summary experiment keys")
    if raw.duplicated(KEYS + ["object_id"]).any():
        raise ValueError("duplicate object within an experiment")
    if not set(raw["split_role"].dropna().unique()).issubset(ROLES):
        raise ValueError("unexpected split_role in raw records")
    if not raw["retained"].isin([0, 1]).all():
        raise ValueError("raw retained labels must be binary")
    if not np.isfinite(raw[["redshift", "p_true"]].to_numpy()).all():
        raise ValueError("raw redshift/p_true contains non-finite values")
    if ((raw["p_true"] <= 0) | (raw["p_true"] >= 1)).any():
        raise ValueError("raw true inclusion probabilities must be in (0,1)")
    if (raw["redshift"] <= 0).any():
        raise ValueError("raw redshifts must be positive")

    reconstructed = 0
    for row in summary.itertuples(index=False):
        key = dict(zip(KEYS, [getattr(row, k) for k in KEYS]))
        mask = np.ones(len(raw), dtype=bool)
        for k, value in key.items():
            mask &= raw[k].to_numpy() == value
        units = raw.loc[mask].sort_values("object_id")
        n_parent = int(row.n_parent)
        if len(units) != n_parent:
            raise ValueError(f"{key}: expected {n_parent} parent records, got {len(units)}")
        ids = units["object_id"].to_numpy(dtype=int)
        if not np.array_equal(np.sort(ids), np.arange(n_parent)):
            raise ValueError(f"{key}: object IDs do not form the complete parent frame")
        if units["object_id"].nunique() != n_parent:
            raise ValueError(f"{key}: object IDs are not unique")

        validation = units.loc[units["split_role"] == "validation"]
        analysis = units.loc[units["split_role"] == "analysis_catalogue"]
        nonretained = units.loc[units["split_role"] == "nonretained_holdout"]
        if len(validation) != int(row.n_validation):
            raise ValueError(f"{key}: validation count does not match design")
        if len(analysis) != int(row.n_analysis_catalogue):
            raise ValueError(f"{key}: analysis-catalogue count does not match summary")
        if not (validation["retained"].isin([0, 1]).all()
                and (analysis["retained"] == 1).all()
                and (nonretained["retained"] == 0).all()):
            raise ValueError(f"{key}: split roles disagree with retention labels")
        if len(validation) + len(analysis) + len(nonretained) != n_parent:
            raise ValueError(f"{key}: split roles do not partition parent frame")
        if len(analysis) != int(units.loc[units["split_role"] != "validation", "retained"].sum()):
            raise ValueError(f"{key}: retained non-validation units disagree with analysis set")

        oracle_stats = _weight_stats(analysis["p_true"].to_numpy(dtype=float))
        oracle_estimates = {
            "oracle_" + k: oracle_stats[k]
            for k in ("kish_ess", "entropy_ess", "max_normalized_weight")
        }
        if row.fit_status == "SUCCESS":
            p_all = units["p_hat"].to_numpy(dtype=float)
            if not np.all(np.isfinite(p_all)) or np.any((p_all <= 0) | (p_all >= 1)):
                raise ValueError(f"{key}: successful fit has invalid predictions")
            holdout = units.loc[units["split_role"] != "validation"]
            p = holdout["p_hat"].to_numpy(dtype=float)
            y = holdout["retained"].to_numpy(dtype=float)
            estimates = {
                "heldout_brier_score": float(np.mean((p - y) ** 2)),
                "heldout_log_loss": float(-np.mean(y * np.log(p) + (1-y) * np.log1p(-p))),
                "heldout_probability_rmse": float(np.sqrt(np.mean(
                    (p - holdout["p_true"].to_numpy(dtype=float)) ** 2
                ))),
            }
            estimates.update(_weight_stats(analysis["p_hat"].to_numpy(dtype=float)))
            if pd.isna(row.coefficients_json) or pd.isna(row.covariance_json):
                raise ValueError(f"{key}: successful fit missing coefficients/covariance")
            beta = np.asarray(json.loads(row.coefficients_json), dtype=float)
            covariance = np.asarray(json.loads(row.covariance_json), dtype=float)
            expected_dimension = int(row.model_degree) + 1
            if beta.shape != (expected_dimension,) or covariance.shape != (
                expected_dimension, expected_dimension
            ):
                raise ValueError(f"{key}: coefficient/covariance dimensions invalid")
            if not np.all(np.isfinite(beta)) or not np.all(np.isfinite(covariance)):
                raise ValueError(f"{key}: non-finite coefficient/covariance")
            if not np.allclose(covariance, covariance.T, rtol=0, atol=1e-10):
                raise ValueError(f"{key}: covariance is not symmetric")
            if np.min(np.linalg.eigvalsh(covariance)) <= 0:
                raise ValueError(f"{key}: covariance is not positive definite")
            if not np.isfinite(row.fit_center) or not np.isfinite(row.fit_scale) or row.fit_scale <= 0:
                raise ValueError(f"{key}: invalid standardization metadata")
            x = (validation["redshift"].to_numpy(dtype=float) - float(row.fit_center)) / float(row.fit_scale)
            design = np.column_stack([x ** power for power in range(expected_dimension)])
            y_validation = validation["retained"].to_numpy(dtype=float)
            p_validation = expit(design @ beta)
            score = design.T @ (p_validation - y_validation)
            if np.linalg.norm(score, ord=np.inf) > 1e-6:
                raise ValueError(f"{key}: stored coefficients fail independent score check")
            information = design.T @ ((p_validation * (1.0 - p_validation))[:, None] * design)
            if np.min(np.linalg.eigvalsh(information)) <= 0:
                raise ValueError(f"{key}: reconstructed observed information is not positive definite")
            covariance_rebuilt = np.linalg.inv(information)
            if not np.allclose(covariance, covariance_rebuilt, rtol=1e-7, atol=1e-9):
                raise ValueError(f"{key}: covariance does not match reconstructed information")
            prediction_rebuilt = expit(
                ((units["redshift"].to_numpy(dtype=float) - float(row.fit_center))
                 / float(row.fit_scale))[:, None]
                ** np.arange(expected_dimension) @ beta
            )
            if not np.allclose(units["p_hat"].to_numpy(dtype=float), prediction_rebuilt,
                               rtol=1e-10, atol=1e-12):
                raise ValueError(f"{key}: raw predicted probabilities do not match fitted coefficients")
        else:
            if units["p_hat"].notna().any():
                raise ValueError(f"{key}: unsuccessful fit has predicted probabilities")
            if not pd.isna(row.coefficients_json) or not pd.isna(row.covariance_json):
                raise ValueError(f"{key}: unsuccessful fit must not report coefficients")
            estimates = {name: float("nan") for name in METRICS if not name.startswith("oracle_")}
        estimates.update(oracle_estimates)
        for name, estimate in estimates.items():
            reported = getattr(row, name)
            if np.isnan(estimate):
                if not pd.isna(reported):
                    raise ValueError(f"{key}: {name} should be NaN, got {reported}")
            elif not np.isfinite(reported) or abs(float(reported) - estimate) > tolerance:
                raise ValueError(
                    f"{key}: {name} mismatch; summary={reported}, raw_recomputed={estimate}"
                )
        reconstructed += 1

    raw_groups = len(raw.groupby(KEYS, dropna=False))
    if len(summary) != raw_groups:
        raise ValueError("summary/raw experiment-group counts differ")
    return {
        "status": "PASS",
        "experiments_reconstructed": reconstructed,
        "raw_unit_records": int(len(raw)),
        "summary_rows": int(len(summary)),
        "complete_parent_partitions": reconstructed,
        "tolerance": tolerance,
        "claim_limit": "Numerical reconstruction only; no estimator superiority or H0 calibration claim.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = audit(
        args.directory / "replicates.csv",
        args.directory / "raw_unit_records.csv",
        tolerance=args.tolerance,
    )
    rendered = json.dumps(report, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\\n", encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
