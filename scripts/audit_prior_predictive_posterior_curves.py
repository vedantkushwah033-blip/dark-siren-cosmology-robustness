#!/usr/bin/env python3
"""Independently reconstruct prior-predictive posterior ranks and quantiles."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

ARMS = ("A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
        "D_map_plugin_ipw", "E_uncertainty_integrated_ipw")
QUANTILES = {"q025": .025, "q05": .05, "q10": .10, "q16": .16,
             "q25": .25, "q75": .75, "q84": .84, "q90": .90,
             "q95": .95, "q975": .975}

def _gamma_token(gamma):
    return f"{float(gamma):+g}".replace("+", "p").replace("-", "m")

def _posterior_cdf_summary(grid, log_likelihood, truth):
    """Normalize exp(log_likelihood) by trapezoidal integration over dH0."""
    x, ll = np.asarray(grid, float), np.asarray(log_likelihood, float)
    if x.ndim != 1 or ll.shape != x.shape or x.size < 2:
        raise AssertionError("H0 grid and log-likelihood must be matching 1D arrays")
    if not (np.isfinite(x).all() and np.isfinite(ll).all()):
        raise AssertionError("non-finite H0 grid or log-likelihood")
    dx = np.diff(x)
    if np.any(dx <= 0):
        raise AssertionError("H0 grid must be strictly increasing")
    if not x[0] <= truth <= x[-1]:
        raise AssertionError("generating H0 truth lies outside posterior grid")
    density = np.exp(ll - np.max(ll))
    increments = .5 * (density[:-1] + density[1:]) * dx
    cumulative = np.concatenate(([0.0], np.cumsum(increments)))
    total = float(cumulative[-1])
    if not np.isfinite(total) or total <= 0:
        raise AssertionError("posterior normalization integral is not positive")
    cdf = cumulative / total
    if np.any(np.diff(cdf) < 0):
        raise AssertionError("posterior CDF is not monotone")
    out = {"posterior_rank": float(np.interp(truth, x, cdf))}
    out.update({k: float(np.interp(q, cdf, x)) for k, q in QUANTILES.items()})
    return out

def _verify_manifest(root):
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    files = manifest.get("files", {})
    if manifest.get("status") != "PASS" or manifest.get("file_count") != len(files):
        raise AssertionError("outer manifest status/count is invalid")
    checked = 0
    for relative, expected in files.items():
        path = root / relative
        if not path.is_file():
            raise AssertionError(f"manifest file is missing: {relative}")
        payload = path.read_bytes()
        if len(payload) != int(expected["bytes"]):
            raise AssertionError(f"manifest byte count mismatch: {relative}")
        if hashlib.sha256(payload).hexdigest() != expected["sha256"]:
            raise AssertionError(f"manifest SHA-256 mismatch: {relative}")
        checked += 1
    return checked

def audit(root, output_prefix=None, *, verify_manifest=True,
          rank_tolerance=5e-10, quantile_tolerance=5e-8):
    root = Path(root).resolve()
    cfg = json.loads((root / "config.json").read_text(encoding="utf-8"))
    if cfg.get("truth_mode") != "prior_predictive":
        raise ValueError("posterior-rank audit requires truth_mode=prior_predictive")
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    summaries = pd.read_csv(root / "replicate_posterior_summaries.csv")
    key = ["replicate", "gamma", "validation_size", "arm"]
    for name, frame in (("ledger", ledger), ("summaries", summaries)):
        if frame.duplicated(key).any():
            raise AssertionError(f"duplicate {name} keys")
    if not (ledger.status == "SUCCESS").all() or not (ledger.audit_status == "PASS").all():
        raise AssertionError("all planned lifecycle rows must be SUCCESS/PASS")
    if "status" in summaries and not (summaries.status == "SUCCESS").all():
        raise AssertionError("all posterior summary rows must be SUCCESS")
    if "audit_status" in summaries and not (summaries.audit_status == "PASS").all():
        raise AssertionError("all posterior summary rows must carry audit PASS")
    lk = set(map(tuple, ledger[key].itertuples(index=False, name=None)))
    sk = set(map(tuple, summaries[key].itertuples(index=False, name=None)))
    if lk != sk:
        raise AssertionError("posterior summaries do not exactly cover lifecycle rows")
    reps, gammas = int(cfg["replicates"]), list(map(float, cfg["gammas"]))
    sizes = list(map(int, cfg["validation_sizes"]))
    expected_n = reps * len(gammas) * len(sizes) * len(ARMS)
    if len(ledger) != expected_n:
        raise AssertionError(f"planned denominator mismatch: {len(ledger)} != {expected_n}")
    if set(summaries.arm) != set(ARMS):
        raise AssertionError("posterior summaries do not contain exactly the five arms")
    manifest_files = _verify_manifest(root) if verify_manifest else 0
    rows, child_audits = [], 0
    max_rank_error = max_quantile_error = 0.0
    worst_rank = worst_quantile = None
    for rep in range(reps):
        for gamma in gammas:
            name = f"rep{rep:03d}_g{_gamma_token(gamma)}"
            child = root / "replicates" / name
            child_audit = json.loads((child / "audit_report.json").read_text(encoding="utf-8"))
            if child_audit.get("status") != "PASS":
                raise AssertionError(f"child audit not PASS: {name}")
            child_cfg = json.loads((child / "config.json").read_text(encoding="utf-8"))
            if float(child_cfg["curvature_gamma"]) != gamma:
                raise AssertionError(f"child curvature mismatch: {name}")
            truth = float(child_cfg["h0_true"])
            child_audits += 1
            with np.load(child / "results" / "log_likelihoods.npz", allow_pickle=False) as ar:
                grid = np.asarray(ar["H0_grid"], float)
                if len(grid) != int(cfg["h0_points"]):
                    raise AssertionError(f"H0 grid length mismatch: {name}")
                support = list(map(float, cfg["h0_prior_support"]))
                if not np.isclose(grid[0], support[0]) or not np.isclose(grid[-1], support[1]):
                    raise AssertionError(f"H0 grid support mismatch: {name}")
                for nval in sizes:
                    block = summaries.loc[(summaries.replicate == rep)
                        & (summaries.gamma == gamma) & (summaries.validation_size == nval)]
                    if len(block) != len(ARMS):
                        raise AssertionError(f"wrong arm count: {name}, n={nval}")
                    for arm in ARMS:
                        saved = block.loc[block.arm == arm]
                        if len(saved) != 1:
                            raise AssertionError(f"missing/duplicate summary: {name}, {arm}, n={nval}")
                        row = saved.iloc[0]
                        curve = f"{arm}_n{nval}"
                        if curve not in ar:
                            raise AssertionError(f"likelihood curve missing: {name}/{curve}")
                        if not np.isclose(float(row.h0_true), truth, rtol=0, atol=1e-12):
                            raise AssertionError(f"summary truth differs from child config: {name}/{curve}")
                        rebuilt = _posterior_cdf_summary(grid, np.asarray(ar[curve], float), truth)
                        rank_err = abs(rebuilt["posterior_rank"] - float(row.posterior_rank))
                        if rank_err > rank_tolerance:
                            raise AssertionError(f"posterior rank mismatch {name}/{curve}: {rank_err}")
                        max_rank_error = max(max_rank_error, rank_err)
                        if rank_err == max_rank_error:
                            worst_rank = f"{name}/{curve}"
                        for field in QUANTILES:
                            err = abs(rebuilt[field] - float(row[field]))
                            if err > quantile_tolerance:
                                raise AssertionError(f"{field} mismatch {name}/{curve}: {err}")
                            max_quantile_error = max(max_quantile_error, err)
                            if err == max_quantile_error:
                                worst_quantile = f"{name}/{curve}:{field}"
                        rows.append({"replicate": rep, "gamma": gamma,
                            "validation_size": nval, "arm": arm,
                            "rebuilt_posterior_rank": rebuilt["posterior_rank"],
                            "saved_posterior_rank": float(row.posterior_rank),
                            "absolute_rank_difference": rank_err,
                            **{f"rebuilt_{k}": v for k, v in rebuilt.items()
                               if k != "posterior_rank"}})
    report = {
        "status": "PASS", "truth_mode": "prior_predictive",
        "planned_likelihood_curves": expected_n,
        "independently_reconstructed_posterior_ranks": len(rows),
        "successful_child_audits_rechecked": child_audits,
        "outer_manifest_files_verified": manifest_files,
        "maximum_absolute_rank_difference": max_rank_error,
        "worst_rank_curve": worst_rank,
        "maximum_absolute_quantile_difference": max_quantile_error,
        "worst_quantile": worst_quantile,
        "rank_tolerance": rank_tolerance, "quantile_tolerance": quantile_tolerance,
        "method": "Independent trapezoidal integration of exp(saved log-likelihood) over dH0, followed by CDF interpolation at the generating truth and inverse-CDF quantiles.",
        "limitations": [
            "Reconstructs summaries from saved H0-grid likelihood curves; does not recompute the likelihood on a finer H0 grid.",
            "Does not establish grid convergence or transfer to another synthetic DGP.",
            "Prior-predictive rank tests remain finite-sample descriptive at 50 replicates per curvature/arm stratum.",
        ],
    }
    if output_prefix is not None:
        output_prefix = Path(output_prefix)
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(output_prefix.with_suffix(".csv"), index=False)
        output_prefix.with_suffix(".json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact_dir", type=Path)
    p.add_argument("--output-prefix", type=Path, required=True)
    p.add_argument("--skip-manifest", action="store_true")
    a = p.parse_args()
    result = audit(a.artifact_dir, a.output_prefix, verify_manifest=not a.skip_manifest)
    print(json.dumps(result, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()
