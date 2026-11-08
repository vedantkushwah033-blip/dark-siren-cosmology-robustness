"""Independently reconstruct finite-validation probability and IPW summaries.

The auditor consumes the unit-level non-validation parent records and the
summary CSV. It does not import the benchmark or fitter implementation.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

KEYS = ["scenario_seed", "replicate", "gamma", "n_validation", "model_degree"]
METRICS = [
    "heldout_brier_score", "heldout_log_loss", "heldout_probability_rmse",
    "kish_ess", "entropy_ess", "max_normalized_weight",
    "top_1pct_weight_mass", "top_5pct_weight_mass",
    "oracle_kish_ess", "oracle_entropy_ess", "oracle_max_normalized_weight",
]


def _weight_stats(probability: np.ndarray) -> dict[str, float]:
    p = np.asarray(probability, dtype=float)
    if p.size == 0:
        return {name: float("nan") for name in (
            "kish_ess", "entropy_ess", "max_normalized_weight",
            "top_1pct_weight_mass", "top_5pct_weight_mass",
        )}
    if not np.all(np.isfinite(p)) or np.any((p <= 0) | (p > 1)):
        raise ValueError("invalid inclusion probabilities in raw records")
    w = 1.0 / p
    q = w / np.sum(w)
    ordered = np.sort(q)[::-1]
    return {
        "kish_ess": float(np.sum(w) ** 2 / np.sum(w ** 2)),
        "entropy_ess": float(np.exp(-np.sum(q * np.log(q)))),
        "max_normalized_weight": float(np.max(q)),
        "top_1pct_weight_mass": float(np.sum(ordered[:max(1, int(np.ceil(.01 * len(q))))])),
        "top_5pct_weight_mass": float(np.sum(ordered[:max(1, int(np.ceil(.05 * len(q))))])),
    }


def audit(summary_path: Path, raw_path: Path, tolerance: float = 1e-10) -> dict:
    summary = pd.read_csv(summary_path)
    raw = pd.read_csv(raw_path)
    required_summary = set(KEYS + ["fit_status", "n_parent"] + METRICS)
    required_raw = set(KEYS + [
        "record_type", "object_id", "redshift", "retained", "p_true", "p_hat"
    ])
    if not required_summary.issubset(summary.columns):
        raise ValueError(f"summary missing columns: {sorted(required_summary - set(summary.columns))}")
    if not required_raw.issubset(raw.columns):
        raise ValueError(f"raw file missing columns: {sorted(required_raw - set(raw.columns))}")
    if summary.duplicated(KEYS).any():
        raise ValueError("duplicate summary experiment keys")
    if raw["record_type"].ne("heldout_parent").any():
        raise ValueError("unexpected raw record_type")
    if raw.duplicated(KEYS + ["object_id"]).any():
        raise ValueError("duplicate object within an experiment")
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
        expected_n = int(row.n_parent) - int(row.n_validation)
        if len(units) != expected_n:
            raise ValueError(f"{key}: expected {expected_n} non-validation units, got {len(units)}")
        if units["object_id"].nunique() != expected_n:
            raise ValueError(f"{key}: object IDs are not unique")
        if not np.all(np.isfinite(units["p_true"])):
            raise ValueError(f"{key}: missing simulation truth probabilities")
        retained = units["retained"].to_numpy(dtype=int) == 1
        oracle_stats = _weight_stats(
            units.loc[retained, "p_true"].to_numpy(dtype=float)
        )
        oracle_estimates = {
            "oracle_" + k: oracle_stats[k]
            for k in ("kish_ess", "entropy_ess", "max_normalized_weight")
        }
        if row.fit_status == "SUCCESS":
            p = units["p_hat"].to_numpy(dtype=float)
            if not np.all(np.isfinite(p)) or np.any((p <= 0) | (p >= 1)):
                raise ValueError(f"{key}: successful fit has invalid predicted probabilities")
            y = units["retained"].to_numpy(dtype=float)
            estimates = {
                "heldout_brier_score": float(np.mean((p - y) ** 2)),
                "heldout_log_loss": float(-np.mean(y * np.log(p) + (1-y) * np.log1p(-p))),
                "heldout_probability_rmse": float(np.sqrt(np.mean(
                    (p - units["p_true"].to_numpy(dtype=float)) ** 2
                ))),
            }
            estimates.update(_weight_stats(p[retained]))
        else:
            if units["p_hat"].notna().any():
                raise ValueError(f"{key}: failed/non-evaluable fit has predicted probabilities")
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

    if len(summary) != len(raw.groupby(KEYS, dropna=False)):
        raise ValueError("summary/raw experiment-group counts differ")
    return {
        "status": "PASS",
        "experiments_reconstructed": reconstructed,
        "raw_unit_records": int(len(raw)),
        "summary_rows": int(len(summary)),
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
    import json
    rendered = json.dumps(report, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
