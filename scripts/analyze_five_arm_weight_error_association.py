"""Exploratory association of pre-posterior IPW diagnostics with H0 errors.

Uses audited synthetic five-arm replicate artifacts. Reports within-cell
Spearman associations and deterministic replicate bootstrap intervals. This is
an association screen, not a predictive model, causal analysis, or calibration
test. Arms A/B/C are evaluated at n_validation=500 only because their diagnostic
and posterior outputs are invariant across validation-size labels.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ARMS = [
    "A_catalogue_only",
    "B_parent_oracle",
    "C_known_probability_ipw",
    "D_map_plugin_ipw",
    "E_uncertainty_integrated_ipw",
]
BASE_DIAGNOSTICS = [
    "kish_ess_fraction",
    "entropy_ess",
    "maximum_normalized_weight",
    "top_1_percent_mass",
    "top_5_percent_mass",
]
INTEGRATED_DIAGNOSTICS = [
    "quadrature_ess",
    "minimum_node_kish_ess",
    "maximum_node_top_5_percent_mass",
]


def spearman_bootstrap(x, y, *, resamples=3000, seed=20261120):
    """Spearman rho and percentile bootstrap interval over paired replicates."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape:
        raise ValueError("x and y must be matching one-dimensional arrays")
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    n = len(x)
    if n < 10 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return {"n": n, "rho": np.nan, "bootstrap_low": np.nan,
                "bootstrap_high": np.nan, "status": "NON_EVALUABLE"}
    rho = float(spearmanr(x, y).statistic)
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(resamples):
        idx = rng.integers(0, n, size=n)
        if np.unique(x[idx]).size < 2 or np.unique(y[idx]).size < 2:
            continue
        value = float(spearmanr(x[idx], y[idx]).statistic)
        if np.isfinite(value):
            estimates.append(value)
    if len(estimates) < int(0.8 * resamples):
        return {"n": n, "rho": rho, "bootstrap_low": np.nan,
                "bootstrap_high": np.nan, "status": "BOOTSTRAP_UNSTABLE"}
    low, high = np.quantile(estimates, [0.025, 0.975])
    return {"n": n, "rho": rho, "bootstrap_low": float(low),
            "bootstrap_high": float(high), "bootstrap_successes": len(estimates),
            "status": "PASS"}


def load_artifact(root: Path):
    root = root.resolve()
    cfg = json.loads((root / "config.json").read_text())
    summaries = pd.read_csv(root / "replicate_posterior_summaries.csv")
    diagnostics = pd.read_csv(root / "replicate_weight_diagnostics.csv")
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    reps = int(cfg["replicates"])
    expected = reps * len(cfg["gammas"]) * len(cfg["validation_sizes"]) * len(ARMS)
    if len(ledger) != expected or ledger.cell_id.duplicated().any():
        raise AssertionError(f"lifecycle denominator or cell IDs invalid: {root}")
    if not (ledger.status == "SUCCESS").all() or not (ledger.audit_status == "PASS").all():
        raise AssertionError(f"artifact has failed/non-audited cells: {root}")
    if len(summaries) != expected or summaries.duplicated(
        ["replicate", "gamma", "validation_size", "arm"]
    ).any():
        raise AssertionError(f"posterior summary denominator/key invalid: {root}")
    if not (summaries.status == "SUCCESS").all() or not (summaries.audit_status == "PASS").all():
        raise AssertionError(f"posterior summaries are not fully audited: {root}")
    if diagnostics.duplicated(["replicate", "gamma", "validation_size", "arm"]).any():
        raise AssertionError(f"duplicate diagnostic rows: {root}")
    if len(diagnostics) != expected:
        raise AssertionError(f"diagnostic denominator mismatch: {root}")
    summaries = summaries.copy()
    diagnostics = diagnostics.copy()
    summaries["master_seed"] = int(cfg["master_seed"])
    diagnostics["master_seed"] = int(cfg["master_seed"])
    return cfg, summaries, diagnostics


def analyze(artifact_dirs, output_prefix: Path, *, resamples=3000):
    loaded = [load_artifact(Path(p)) for p in artifact_dirs]
    configs = [item[0] for item in loaded]
    seeds = [int(c["master_seed"]) for c in configs]
    if len(set(seeds)) != len(seeds):
        raise AssertionError("duplicate master seed streams")
    for cfg in configs[1:]:
        for key in ["replicates", "gammas", "validation_sizes", "n_parent",
                    "n_detected", "validation_frame", "h0_true", "h0_prior_support"]:
            if cfg[key] != configs[0][key]:
                raise AssertionError(f"artifact configuration mismatch: {key}")

    summaries = pd.concat([item[1] for item in loaded], ignore_index=True)
    diagnostics = pd.concat([item[2] for item in loaded], ignore_index=True)
    keys = ["master_seed", "replicate", "gamma", "validation_size", "arm"]
    joined = summaries.merge(diagnostics, on=keys, validate="one_to_one")
    joined["absolute_median_error"] = joined["median_error"].abs()
    joined["covered68"] = (
        (joined["q16"] <= joined["h0_true"]) &
        (joined["h0_true"] <= joined["q84"])
    ).astype(int)

    rows = []
    for (gamma, nval, arm), group in joined.groupby(
        ["gamma", "validation_size", "arm"], sort=True
    ):
        # A/B/C do not depend on validation n. Avoid counting their identical
        # repeated rows as separate samples by retaining only the largest n.
        if arm in ARMS[:3] and int(nval) != max(map(int, configs[0]["validation_sizes"])):
            continue
        candidates = BASE_DIAGNOSTICS.copy()
        if arm == "E_uncertainty_integrated_ipw":
            candidates = INTEGRATED_DIAGNOSTICS.copy()
        for metric_index, metric in enumerate(candidates):
            if metric not in group.columns:
                continue
            result = spearman_bootstrap(
                group[metric].to_numpy(float),
                group["absolute_median_error"].to_numpy(float),
                resamples=resamples,
                seed=20261120 + int(round(float(gamma) * 100)) + int(nval) + metric_index,
            )
            rows.append({
                "master_seeds": ";".join(map(str, seeds)),
                "gamma": float(gamma),
                "validation_size": int(nval),
                "arm": str(arm),
                "diagnostic": metric,
                "outcome": "absolute_posterior_median_error",
                **result,
            })
    result = pd.DataFrame(rows)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_prefix.with_suffix(".csv"), index=False)
    report = {
        "status": "PASS" if len(result) else "EMPTY",
        "master_seeds": seeds,
        "replicates_per_stream": [int(c["replicates"]) for c in configs],
        "pooled_replicates": int(len(summaries) // (
            len(configs[0]["gammas"]) * len(configs[0]["validation_sizes"]) * len(ARMS)
        )),
        "audited_posterior_rows": len(summaries),
        "audited_weight_diagnostic_rows": len(diagnostics),
        "association_rows": len(result),
        "bootstrap_resamples": int(resamples),
        "scope": "Exploratory within-cell association only; not out-of-sample prediction, causality, or calibration certification.",
        "handling": "A/B/C use only the largest validation-size label to avoid repeated identical rows; D/E retain each validation size.",
    }
    output_prefix.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(result.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact_dirs", nargs="+", type=Path)
    p.add_argument("--output-prefix", type=Path, default=Path("weight_error_association"))
    p.add_argument("--bootstrap-resamples", type=int, default=3000)
    args = p.parse_args()
    if args.bootstrap_resamples < 100:
        raise ValueError("use at least 100 bootstrap resamples")
    analyze(args.artifact_dirs, args.output_prefix, resamples=args.bootstrap_resamples)


if __name__ == "__main__":
    main()
