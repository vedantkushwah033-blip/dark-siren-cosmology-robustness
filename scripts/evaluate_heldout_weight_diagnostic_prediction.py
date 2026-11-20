"""Prospective held-out test of weight-tail diagnostics as H0-error predictors.

Training streams are fixed at master seeds 2026111929 and 2026111940.
The untouched holdout stream is fixed at master seed 2026111950.
Within each gamma × validation-size × estimated-IPW arm, fit a one-predictor
OLS model of absolute posterior-median error on the pre-posterior top-5% weight
mass. Compare against a training-mean-only baseline. The bootstrap resamples
both training and holdout replicate units and refits the model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

TRAIN_SEEDS = (2026111929, 2026111940)
HOLDOUT_SEED = 2026111950
GAMMAS = (-1.0, 0.0, 1.0)
VALIDATION_SIZES = (100, 300, 500)
ARMS = {
    "D_map_plugin_ipw": "top_5_percent_mass",
    "E_uncertainty_integrated_ipw": "maximum_node_top_5_percent_mass",
}
EXPECTED_CONFIG = {
    "replicates": 100,
    "gammas": [-1.0, 0.0, 1.0],
    "validation_sizes": [100, 300, 500],
    "n_parent": 1200,
    "n_detected": 30,
    "validation_frame": 500,
    "h0_true": 70.0,
    "h0_prior_support": [30.0, 110.0],
}


def _verify_manifest(root: Path) -> int:
    """Verify every file and byte count declared by the immutable source manifest."""
    root = root.resolve()
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise AssertionError(f"{root}: missing aggregate manifest")
    manifest = json.loads(manifest_path.read_text())
    files = manifest.get("files")
    if not isinstance(files, dict) or not files:
        raise AssertionError(f"{root}: aggregate manifest has no file map")
    checked = 0
    for relative, metadata in files.items():
        rel = Path(relative)
        if rel.is_absolute() or ".." in rel.parts:
            raise AssertionError(f"{root}: unsafe manifest path: {relative}")
        path = (root / rel).resolve()
        if root not in path.parents or not path.is_file():
            raise AssertionError(f"{root}: manifest file is missing or escapes artifact root: {relative}")
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        if path.stat().st_size != int(metadata["bytes"]) or digest.hexdigest() != metadata["sha256"]:
            raise AssertionError(f"{root}: manifest checksum/size mismatch: {relative}")
        checked += 1
    if checked != int(manifest.get("file_count", -1)):
        raise AssertionError(f"{root}: aggregate manifest file count mismatch")
    return checked


def _assert_exact_lifecycle_keys(
    ledger: pd.DataFrame,
    summaries: pd.DataFrame,
    diagnostics: pd.DataFrame,
    *,
    root: Path,
) -> None:
    """Fail closed unless all three artifact tables describe identical cells."""
    keys = ["replicate", "gamma", "validation_size", "arm"]
    for label, frame in [
        ("ledger", ledger), ("summaries", summaries), ("diagnostics", diagnostics)
    ]:
        missing = sorted(set(keys) - set(frame.columns))
        if missing:
            raise AssertionError(f"{root}: {label} missing lifecycle key columns: {missing}")
        if frame.duplicated(keys).any():
            raise AssertionError(f"{root}: duplicate lifecycle keys in {label}")

    def key_set(frame: pd.DataFrame) -> set[tuple]:
        return set(frame[keys].itertuples(index=False, name=None))

    ledger_keys = key_set(ledger)
    summary_keys = key_set(summaries)
    diagnostic_keys = key_set(diagnostics)
    if ledger_keys != summary_keys:
        raise AssertionError(
            f"{root}: ledger/summary lifecycle key mismatch "
            f"(ledger-only={len(ledger_keys - summary_keys)}, "
            f"summary-only={len(summary_keys - ledger_keys)})"
        )
    if ledger_keys != diagnostic_keys:
        raise AssertionError(
            f"{root}: ledger/diagnostic lifecycle key mismatch "
            f"(ledger-only={len(ledger_keys - diagnostic_keys)}, "
            f"diagnostic-only={len(diagnostic_keys - ledger_keys)})"
        )

    expected_per_cell = 100
    counts = ledger.groupby(["gamma", "validation_size", "arm"], dropna=False).replicate.nunique()
    if len(counts) != 45 or not (counts == expected_per_cell).all():
        raise AssertionError(
            f"{root}: expected exactly 100 unique replicate IDs in each of 45 gamma/size/arm cells"
        )


def _validate_artifact(root: Path, expected_seed: int) -> pd.DataFrame:
    root = root.resolve()
    _verify_manifest(root)
    cfg = json.loads((root / "config.json").read_text())
    if int(cfg.get("master_seed", -1)) != expected_seed:
        raise AssertionError(f"{root}: wrong master seed")
    for key, expected in EXPECTED_CONFIG.items():
        actual = cfg.get(key)
        if actual != expected:
            raise AssertionError(f"{root}: config mismatch for {key}: {actual!r} != {expected!r}")

    audit_path = root / "independent_orchestration_audit.json"
    if not audit_path.exists():
        raise AssertionError(f"{root}: missing independent orchestration audit")
    audit = json.loads(audit_path.read_text())
    if audit.get("status") != "PASS" or audit.get("planned_cells") != 4500:
        raise AssertionError(f"{root}: independent orchestration audit did not PASS")
    if audit.get("successful_child_runs") != 300:
        raise AssertionError(f"{root}: expected 300 independently audited child runs")

    ledger = pd.read_csv(root / "replicate_ledger.csv")
    summaries = pd.read_csv(root / "replicate_posterior_summaries.csv")
    diagnostics = pd.read_csv(root / "replicate_weight_diagnostics.csv")
    expected_rows = 100 * 3 * 3 * 5
    if len(ledger) != expected_rows or ledger.cell_id.duplicated().any():
        raise AssertionError(f"{root}: invalid planned lifecycle denominator")
    if set(ledger.status) != {"SUCCESS"} or set(ledger.audit_status) != {"PASS"}:
        raise AssertionError(f"{root}: failed or unaudited lifecycle cell")
    for label, frame in [("summaries", summaries), ("diagnostics", diagnostics)]:
        if len(frame) != expected_rows:
            raise AssertionError(f"{root}: invalid {label} denominator")
    if set(summaries.status) != {"SUCCESS"} or set(summaries.audit_status) != {"PASS"}:
        raise AssertionError(f"{root}: posterior summaries are not all audited")

    _assert_exact_lifecycle_keys(ledger, summaries, diagnostics, root=root)
    keys = ["replicate", "gamma", "validation_size", "arm"]
    joined = summaries.merge(diagnostics, on=keys, validate="one_to_one", suffixes=("", "_diag"))
    joined["master_seed"] = expected_seed
    joined["absolute_median_error"] = joined["median_error"].abs()
    if not np.isfinite(joined["absolute_median_error"]).all():
        raise AssertionError(f"{root}: non-finite posterior-median errors")
    return joined


def _fit_predict(x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Fit OLS with an intercept after training-only standardization."""
    x_train = np.asarray(x_train, dtype=float)
    y_train = np.asarray(y_train, dtype=float)
    x_test = np.asarray(x_test, dtype=float)
    center = float(np.mean(x_train))
    scale = float(np.std(x_train, ddof=0))
    baseline = float(np.mean(y_train))
    if not np.isfinite(scale) or scale <= 1e-12:
        return np.full_like(x_test, baseline), baseline, 0.0
    xs = (x_train - center) / scale
    xt = (x_test - center) / scale
    design = np.column_stack([np.ones(len(xs)), xs])
    coef, *_ = np.linalg.lstsq(design, y_train, rcond=None)
    pred = coef[0] + coef[1] * xt
    return pred, baseline, float(coef[1] / scale)


def _metrics(y: np.ndarray, pred: np.ndarray) -> tuple[float, float]:
    err = np.asarray(pred, float) - np.asarray(y, float)
    return float(np.mean(np.abs(err))), float(np.sqrt(np.mean(err**2)))


def _holm_adjust(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 0.0
    m = len(p)
    for rank, idx in enumerate(order):
        value = min(1.0, (m - rank) * p[idx])
        running = max(running, value)
        adjusted[idx] = running
    return adjusted.tolist()


def _one_cell(train: pd.DataFrame, holdout: pd.DataFrame, *, x_column: str,
              resamples: int, seed: int) -> dict:
    xtr = train[x_column].to_numpy(float)
    ytr = train["absolute_median_error"].to_numpy(float)
    xho = holdout[x_column].to_numpy(float)
    yho = holdout["absolute_median_error"].to_numpy(float)
    if len(train) != 200 or len(holdout) != 100:
        raise AssertionError(f"expected 200 training and 100 holdout replicates; got {len(train)}, {len(holdout)}")
    if not (np.isfinite(xtr).all() and np.isfinite(xho).all()):
        raise AssertionError(f"non-finite diagnostic values in {x_column}")

    pred, baseline, slope = _fit_predict(xtr, ytr, xho)
    model_mae, model_rmse = _metrics(yho, pred)
    baseline_mae, baseline_rmse = _metrics(yho, np.full(len(yho), baseline))
    delta_mae = model_mae - baseline_mae
    delta_rmse = model_rmse - baseline_rmse
    rho = float(spearmanr(xho, yho).statistic) if np.unique(xho).size > 1 and np.unique(yho).size > 1 else float("nan")

    rng = np.random.default_rng(seed)
    boot_delta = []
    ntr, nho = len(train), len(holdout)
    for _ in range(resamples):
        it = rng.integers(0, ntr, size=ntr)
        ih = rng.integers(0, nho, size=nho)
        bp, bb, _ = _fit_predict(xtr[it], ytr[it], xho[ih])
        bmae, _ = _metrics(yho[ih], bp)
        base_mae, _ = _metrics(yho[ih], np.full(nho, bb))
        if np.isfinite(bmae - base_mae):
            boot_delta.append(bmae - base_mae)
    if len(boot_delta) < int(0.9 * resamples):
        raise AssertionError("too few successful nested bootstrap resamples")
    lo, hi = np.quantile(boot_delta, [0.025, 0.975])
    # Center the bootstrap distribution at the null value before computing
    # the one-sided H1: delta_MAE < 0 tail probability. This is the
    # bootstrap-pivot test, not the uncentered fraction of bootstrap estimates
    # that happen to be positive.
    boot_null = np.asarray(boot_delta) - delta_mae
    p_improve = (1 + int(np.sum(boot_null <= delta_mae))) / (len(boot_null) + 1)
    return {
        "training_n": len(train), "holdout_n": len(holdout),
        "diagnostic": x_column, "training_slope_per_unit": slope,
        "holdout_spearman_rho": rho,
        "holdout_model_mae": model_mae, "holdout_baseline_mae": baseline_mae,
        "holdout_delta_mae_model_minus_baseline": delta_mae,
        "delta_mae_bootstrap_low": float(lo), "delta_mae_bootstrap_high": float(hi),
        "delta_mae_one_sided_bootstrap_p": float(p_improve),
        "holdout_model_rmse": model_rmse, "holdout_baseline_rmse": baseline_rmse,
        "holdout_delta_rmse_model_minus_baseline": delta_rmse,
        "bootstrap_resamples": int(resamples),
    }


def evaluate(training_dirs: list[Path], holdout_dir: Path, output_prefix: Path,
             *, resamples: int = 3000) -> dict:
    if resamples < 1000:
        raise ValueError("use at least 1,000 bootstrap resamples")
    training = [_validate_artifact(p, seed) for p, seed in zip(training_dirs, TRAIN_SEEDS, strict=True)]
    holdout = _validate_artifact(holdout_dir, HOLDOUT_SEED)
    train = pd.concat(training, ignore_index=True)
    if train["master_seed"].nunique() != 2 or set(train.master_seed.unique()) != set(TRAIN_SEEDS):
        raise AssertionError("training seed streams are not the frozen pair")
    if set(holdout.master_seed.unique()) != {HOLDOUT_SEED}:
        raise AssertionError("holdout seed does not match the frozen independent stream")
    if set(train.master_seed.unique()) & set(holdout.master_seed.unique()):
        raise AssertionError("training/holdout seed leakage")

    rows = []
    for gamma in GAMMAS:
        for nval in VALIDATION_SIZES:
            for arm, x_column in ARMS.items():
                tr = train.loc[
                    (train.gamma == gamma) & (train.validation_size == nval) & (train.arm == arm)
                ].copy()
                ho = holdout.loc[
                    (holdout.gamma == gamma) & (holdout.validation_size == nval) & (holdout.arm == arm)
                ].copy()
                # One row per replicate and arm-cell in each independent stream.
                if tr.duplicated(["master_seed", "replicate"]).any() or ho.duplicated(["master_seed", "replicate"]).any():
                    raise AssertionError("duplicate replicate within a predictive cell")
                result = _one_cell(
                    tr, ho, x_column=x_column, resamples=resamples,
                    seed=20261120 + int((gamma + 1) * 1000) + nval + (0 if arm.startswith("D_") else 100000),
                )
                rows.append({"gamma": gamma, "validation_size": nval, "arm": arm, **result})

    result = pd.DataFrame(rows)
    result["holm_adjusted_p"] = _holm_adjust(result.delta_mae_one_sided_bootstrap_p.tolist())
    result["predeclared_predictive_value"] = (
        (result.holdout_delta_mae_model_minus_baseline < 0)
        & (result.delta_mae_bootstrap_high < 0)
        & (result.holm_adjusted_p < 0.05)
    )
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_prefix.with_suffix(".csv"), index=False)
    report = {
        "status": "PASS" if len(result) == 18 else "FAIL",
        "training_master_seeds": list(TRAIN_SEEDS),
        "holdout_master_seed": HOLDOUT_SEED,
        "training_replicates_per_cell": 200,
        "holdout_replicates_per_cell": 100,
        "predictive_cells": len(result),
        "bootstrap_resamples_per_cell": int(resamples),
        "primary_metric": "holdout MAE difference: diagnostic-only OLS minus training-mean baseline; negative favors diagnostic",
        "familywise_adjustment": "Holm over all 18 prespecified gamma × validation-size × arm tests",
        "cells_meeting_predeclared_predictive_value": int(result.predeclared_predictive_value.sum()),
        "limitations": [
            "One fixed H0 truth and one synthetic DGP family only.",
            "The holdout assesses prediction of absolute posterior-median error, not coverage calibration or causal effects.",
            "Bootstrap intervals resample training and holdout replicates but do not validate transfer to other DGPs.",
            "No universal ESS threshold or general-purpose diagnostic rule is claimed.",
        ],
    }
    output_prefix.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(result.to_string(index=False, float_format=lambda x: f"{x:.5f}"))
    return report


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--train-seed-1929", type=Path, required=True)
    p.add_argument("--train-seed-1940", type=Path, required=True)
    p.add_argument("--holdout-seed-1950", type=Path, required=True)
    p.add_argument("--output-prefix", type=Path, default=Path("heldout_weight_diagnostic_prediction"))
    p.add_argument("--bootstrap-resamples", type=int, default=3000)
    args = p.parse_args()
    evaluate(
        [args.train_seed_1929, args.train_seed_1940], args.holdout_seed_1950,
        args.output_prefix, resamples=args.bootstrap_resamples,
    )


if __name__ == "__main__":
    main()
