"""Train/holdout validation of pre-posterior IPW weight diagnostics.

The training seed stream fixes (i) the direction of the diagnostic-risk
relationship, (ii) the 75th-percentile posterior-error threshold, and
(iii) the corresponding diagnostic alert threshold. All three are then
evaluated unchanged on a disjoint seed stream. This is a single-split
synthetic validation, not a calibrated universal failure rule.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ARMS = {
    "D_map_plugin_ipw": "kish_ess_fraction",
    "E_uncertainty_integrated_ipw": "minimum_node_kish_ess",
}
KEYS = ["replicate", "gamma", "validation_size", "arm"]
SUMMARY_COLUMNS = {
    "h0_true", "posterior_median", "median_error", "q025", "q05",
    "q16", "q84", "q95", "q975",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_audited_stream(root: Path):
    """Verify source manifest and complete audited summary/diagnostic tables."""
    root = root.resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    files = manifest.get("files", {})
    if len(files) != int(manifest.get("file_count", -1)):
        raise AssertionError(f"manifest file count mismatch: {root}")
    for relative, metadata in files.items():
        path = root / relative
        if (not path.is_file()
                or path.stat().st_size != int(metadata["bytes"])
                or _sha256(path) != metadata["sha256"]):
            raise AssertionError(f"manifest mismatch: {relative}")

    config = json.loads((root / "config.json").read_text())
    summaries = pd.read_csv(root / "replicate_posterior_summaries.csv")
    diagnostics = pd.read_csv(root / "replicate_weight_diagnostics.csv")
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    planned = pd.read_csv(root / "planned_cells.csv")
    expected = (
        int(config["replicates"]) * len(config["gammas"])
        * len(config["validation_sizes"]) * 5
    )
    for label, frame in [
        ("ledger", ledger), ("planned cells", planned),
        ("posterior summaries", summaries), ("weight diagnostics", diagnostics),
    ]:
        if len(frame) != expected:
            raise AssertionError(f"{label} denominator mismatch: {len(frame)} != {expected}")
    if not (ledger["status"] == "SUCCESS").all():
        raise AssertionError("source ledger contains non-success outcomes")
    if not (ledger["audit_status"] == "PASS").all():
        raise AssertionError("source ledger contains non-PASS audit outcomes")
    if not (summaries["status"] == "SUCCESS").all() or not (
        summaries["audit_status"] == "PASS"
    ).all():
        raise AssertionError("posterior summaries are not all child-audited SUCCESS")
    # The aggregate diagnostic table contains measured metrics and pairing keys,
    # not lifecycle status columns. Its provenance is covered by the audited
    # replicate ledger and each source artifact's manifest.
    for label, frame in [("ledger", ledger), ("planned", planned)]:
        if "cell_id" not in frame.columns or frame["cell_id"].duplicated().any():
            raise AssertionError(f"missing or duplicate cell_id in {label}")
    if set(ledger["cell_id"]) != set(planned["cell_id"]):
        raise AssertionError("planned-cell and terminal-ledger IDs differ")
    for frame, label in [(summaries, "summaries"), (diagnostics, "diagnostics")]:
        if frame.duplicated(KEYS).any():
            raise AssertionError(f"duplicate pairing keys in {label}")
    if not SUMMARY_COLUMNS.issubset(summaries.columns):
        raise AssertionError(
            f"missing posterior columns: {sorted(SUMMARY_COLUMNS-set(summaries.columns))}"
        )
    quantiles = summaries[
        ["q025", "q05", "q16", "posterior_median", "q84", "q95", "q975"]
    ].to_numpy(float)
    if not np.isfinite(quantiles).all() or np.any(np.diff(quantiles, axis=1) < -1e-10):
        raise AssertionError("non-finite or unordered posterior quantiles")
    expected_error = summaries["posterior_median"].to_numpy(float) - summaries[
        "h0_true"
    ].to_numpy(float)
    if not np.allclose(
        summaries["median_error"].to_numpy(float), expected_error, rtol=0, atol=1e-8
    ):
        raise AssertionError("median_error does not equal posterior_median - h0_true")
    if diagnostics[["replicate", "gamma", "validation_size", "arm"]].isna().any().any():
        raise AssertionError("missing diagnostic pairing keys")
    joined = summaries.merge(diagnostics, on=KEYS, suffixes=("", "_diagnostic"),
                             validate="one_to_one")
    joined["absolute_median_error"] = joined["median_error"].abs()
    return config, joined, len(files)


def _auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Mann-Whitney AUROC with average ranks for ties."""
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, int)
    positive = labels == 1
    n_pos = int(positive.sum())
    n_neg = int((~positive).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = rankdata(scores, method="average")
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def _stratified_auc_bootstrap_ci(
    scores: np.ndarray, labels: np.ndarray, *, seed: int, replicates: int = 2000
):
    """Percentile CI for AUROC, resampling positive and negative classes separately."""
    scores = np.asarray(scores, float)
    labels = np.asarray(labels, int)
    positive = np.flatnonzero(labels == 1)
    negative = np.flatnonzero(labels == 0)
    if len(positive) == 0 or len(negative) == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    estimates = np.empty(replicates, dtype=float)
    for i in range(replicates):
        sample = np.concatenate([
            rng.choice(positive, size=len(positive), replace=True),
            rng.choice(negative, size=len(negative), replace=True),
        ])
        estimates[i] = _auc(scores[sample], labels[sample])
    low, high = np.quantile(estimates, [0.025, 0.975])
    return float(low), float(high)


def _wilson(successes: int, total: int, z: float = 1.959963984540054):
    if total == 0:
        return float("nan"), float("nan")
    p = successes / total
    den = 1 + z * z / total
    center = (p + z * z / (2 * total)) / den
    half = z * np.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / den
    return float(center - half), float(center + half)


def evaluate_frozen_rule(train_x, train_error, test_x, test_error, *, min_abs_rho=0.10):
    """Freeze direction and thresholds on train; evaluate once on untouched holdout."""
    train_x = np.asarray(train_x, float)
    train_error = np.asarray(train_error, float)
    test_x = np.asarray(test_x, float)
    test_error = np.asarray(test_error, float)
    if (train_x.ndim != 1 or train_error.ndim != 1 or test_x.ndim != 1
            or test_error.ndim != 1 or len(train_x) != len(train_error)
            or len(test_x) != len(test_error)):
        raise ValueError("diagnostics and outcomes must be matching one-dimensional arrays")
    if not all(np.isfinite(v).all() for v in (train_x, train_error, test_x, test_error)):
        raise ValueError("non-finite values in train/holdout data")
    if min(len(train_x), len(test_x)) < 20:
        return {"status": "NON_EVALUABLE", "reason": "fewer than 20 rows"}
    if np.unique(train_x).size < 2 or np.unique(train_error).size < 2:
        return {"status": "NON_EVALUABLE", "reason": "constant training variable"}
    rho = float(spearmanr(train_x, train_error).statistic)
    if not np.isfinite(rho) or abs(rho) < min_abs_rho:
        return {"status": "NON_EVALUABLE", "reason": "training direction too weak",
                "training_rho": rho}
    direction = 1 if rho > 0 else -1
    oriented_train = direction * train_x
    oriented_test = direction * test_x
    diagnostic_cut = float(np.quantile(oriented_train, 0.75))
    error_cut = float(np.quantile(train_error, 0.75))
    train_alert = oriented_train >= diagnostic_cut
    test_alert = oriented_test >= diagnostic_cut
    test_high_error = test_error >= error_cut
    if np.unique(test_high_error).size < 2 or np.unique(oriented_test).size < 2:
        return {"status": "NON_EVALUABLE", "reason": "holdout class or diagnostic is constant",
                "training_rho": rho}
    tp = int(np.sum(test_alert & test_high_error))
    fp = int(np.sum(test_alert & ~test_high_error))
    fn = int(np.sum(~test_alert & test_high_error))
    tn = int(np.sum(~test_alert & ~test_high_error))
    positives = tp + fn
    negatives = tn + fp
    sensitivity = tp / positives
    specificity = tn / negatives
    sens_ci = _wilson(tp, positives)
    spec_ci = _wilson(tn, negatives)
    return {
        "status": "PASS",
        "training_rho": rho,
        "risk_direction": "higher diagnostic" if direction > 0 else "lower diagnostic",
        "diagnostic_threshold_train_q75_oriented": diagnostic_cut,
        "absolute_error_threshold_train_q75": error_cut,
        "train_alert_fraction": float(train_alert.mean()),
        "holdout_n": len(test_x),
        "holdout_high_error_prevalence": float(test_high_error.mean()),
        "holdout_alert_fraction": float(test_alert.mean()),
        "holdout_auroc": _auc(oriented_test, test_high_error),
        "holdout_auroc_bootstrap95_low": _stratified_auc_bootstrap_ci(
            oriented_test, test_high_error, seed=20261120 + int(round(abs(rho) * 1e6))
        )[0],
        "holdout_auroc_bootstrap95_high": _stratified_auc_bootstrap_ci(
            oriented_test, test_high_error, seed=20261120 + int(round(abs(rho) * 1e6))
        )[1],
        "holdout_sensitivity": sensitivity,
        "holdout_sensitivity_wilson95_low": sens_ci[0],
        "holdout_sensitivity_wilson95_high": sens_ci[1],
        "holdout_specificity": specificity,
        "holdout_specificity_wilson95_low": spec_ci[0],
        "holdout_specificity_wilson95_high": spec_ci[1],
        "holdout_precision": tp / (tp + fp) if tp + fp else float("nan"),
        "holdout_balanced_accuracy": 0.5 * (sensitivity + specificity),
        "holdout_confusion_tp_fp_fn_tn": [tp, fp, fn, tn],
    }


def analyze(train_root: Path, holdout_root: Path, output_prefix: Path):
    train_cfg, train, train_manifest_count = load_audited_stream(train_root)
    hold_cfg, holdout, hold_manifest_count = load_audited_stream(holdout_root)
    if int(train_cfg["master_seed"]) == int(hold_cfg["master_seed"]):
        raise AssertionError("training and holdout master seeds must be disjoint")
    for key in ["replicates", "gammas", "validation_sizes", "n_parent", "n_detected",
                "validation_frame", "h0_true", "h0_prior_support"]:
        if train_cfg[key] != hold_cfg[key]:
            raise AssertionError(f"training/holdout configuration mismatch: {key}")
    if list(map(float, train_cfg["h0_prior_support"])) != [30.0, 110.0]:
        raise AssertionError("expected the expanded [30,110] H0 prior in both streams")
    if set(train["arm"]) != set(holdout["arm"]):
        raise AssertionError("training and holdout arm sets differ")
    rows = []
    for (gamma, nval, arm), tr in train.groupby(["gamma", "validation_size", "arm"]):
        if arm not in ARMS:
            continue
        te = holdout[
            (holdout.gamma == gamma) & (holdout.validation_size == nval)
            & (holdout.arm == arm)
        ]
        metric = ARMS[arm]
        if metric not in tr or metric not in te:
            raise AssertionError(f"missing primary diagnostic {metric}")
        tr = tr.sort_values("replicate")
        te = te.sort_values("replicate")
        if tr["replicate"].duplicated().any() or te["replicate"].duplicated().any():
            raise AssertionError("duplicate replicate index within a stratum")
        result = evaluate_frozen_rule(
            tr[metric].to_numpy(float), tr["absolute_median_error"].to_numpy(float),
            te[metric].to_numpy(float), te["absolute_median_error"].to_numpy(float),
        )
        rows.append({
            "gamma": float(gamma), "validation_size": int(nval), "arm": arm,
            "primary_diagnostic": metric,
            "training_master_seed": int(train_cfg["master_seed"]),
            "holdout_master_seed": int(hold_cfg["master_seed"]),
            "training_n": len(tr), **result,
        })
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(output_prefix.with_suffix(".csv"), index=False)
    report = {
        "status": "PASS" if len(frame) == 18 else "INCOMPLETE",
        "training_master_seed": int(train_cfg["master_seed"]),
        "holdout_master_seed": int(hold_cfg["master_seed"]),
        "training_manifest_entries_verified": train_manifest_count,
        "holdout_manifest_entries_verified": hold_manifest_count,
        "training_rows": len(train),
        "holdout_rows": len(holdout),
        "strata_planned": 18,
        "strata_reported": len(frame),
        "strata_evaluable": int((frame.status == "PASS").sum()),
        "strata_non_evaluable": int((frame.status != "PASS").sum()),
        "rule": (
            "Within each gamma × validation-size × arm stratum, learn Spearman "
            "risk direction on training stream; freeze the 75th-percentile "
            "absolute-error cutoff and oriented diagnostic alert cutoff; evaluate "
            "unchanged on the disjoint holdout stream."
        ),
        "limitations": [
            "One fixed train/holdout split and one synthetic H0 truth.",
            "Only two 100-replicate seed streams; no repeated cross-validation.",
            "The 75th-percentile threshold is an operational exploratory label, not a scientific failure definition.",
            "Results are stratum-specific and not a universal diagnostic threshold.",
            "No observational H0 or real-catalogue inference.",
        ],
    }
    output_prefix.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(frame.to_string(index=False, float_format=lambda value: f"{value:.4f}"))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("training_artifact", type=Path)
    parser.add_argument("holdout_artifact", type=Path)
    parser.add_argument("--output-prefix", type=Path, default=Path("holdout_diagnostic_validation"))
    args = parser.parse_args()
    analyze(args.training_artifact, args.holdout_artifact, args.output_prefix)


if __name__ == "__main__":
    main()
