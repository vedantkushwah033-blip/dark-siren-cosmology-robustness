"""Independent audit and pooled readout for two 100-replicate expanded-prior runs.

The two streams are independent (master seeds 2026111940 and 2026111929).
This reconstructs each orchestration audit, verifies both aggregate manifests,
recomputes every saved posterior median from its child likelihood curve, and
renormalizes the same likelihoods on [50,90] as a prior-support sensitivity.
It does not certify observational validity or estimator superiority.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.stats import t

ARMS = [
    "A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
    "D_map_plugin_ipw", "E_uncertainty_integrated_ipw",
]
SEEDS = (2026111940, 2026111929)


def posterior(grid: np.ndarray, logl: np.ndarray, lo: float, hi: float) -> dict:
    mask = (grid >= lo - 1e-12) & (grid <= hi + 1e-12)
    h = grid[mask]
    y = logl[mask]
    if len(h) < 3 or not np.isclose(h[0], lo) or not np.isclose(h[-1], hi):
        raise AssertionError(f"grid does not contain support endpoints {lo},{hi}")
    d = np.exp(y - np.max(y))
    area = float(trapezoid(d, h))
    if not np.isfinite(area) or area <= 0:
        raise AssertionError("invalid posterior normalization")
    d /= area
    cdf = np.r_[0.0, cumulative_trapezoid(d, h)]
    cdf /= cdf[-1]
    return {
        "median": float(np.interp(.5, cdf, h)),
        "q025": float(np.interp(.025, cdf, h)),
        "q16": float(np.interp(.16, cdf, h)),
        "q84": float(np.interp(.84, cdf, h)),
        "q975": float(np.interp(.975, cdf, h)),
        "density": d,
        "grid": h,
    }


def t_interval(x: np.ndarray) -> tuple[float, float, float, float]:
    x = np.asarray(x, dtype=float)
    n = len(x)
    mean = float(np.mean(x))
    se = float(np.std(x, ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
    crit = float(t.ppf(.975, n - 1)) if n > 1 else float("nan")
    return mean, se, mean - crit * se, mean + crit * se


def audit_one(root: Path, expected_seed: int, auditor: Path, out: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    root = root.resolve()
    # Always rerun the auditor; never trust a report embedded in the source ZIP.
    audit_path = out / f"independent_orchestration_audit_seed_{expected_seed}.json"
    subprocess.run([sys.executable, str(auditor), str(root), "--output", str(audit_path)], check=True)
    audit = json.loads(audit_path.read_text())
    cfg = json.loads((root / "config.json").read_text())
    if cfg.get("master_seed") != expected_seed or cfg.get("replicates") != 100:
        raise AssertionError(f"unexpected stream identity/config at {root}")
    if cfg.get("h0_prior_support") != [30.0, 110.0]:
        raise AssertionError("expanded support is not [30,110]")
    if audit.get("status") != "PASS" or audit.get("planned_cells") != 4500:
        raise AssertionError(f"orchestration audit failed: {audit}")
    if audit.get("successful_child_runs") != 300 or audit.get("aggregate_manifest_entries_verified", 0) < 5000:
        raise AssertionError(f"incomplete child or manifest audit: {audit}")
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    summaries = pd.read_csv(root / "replicate_posterior_summaries.csv")
    if len(ledger) != 4500 or len(summaries) != 4500:
        raise AssertionError("expected 4,500 lifecycle and posterior rows")
    if set(ledger.status) != {"SUCCESS"} or set(ledger.audit_status) != {"PASS"}:
        raise AssertionError("planned denominator contains a failed or unaudited cell")
    if set(summaries.status) != {"SUCCESS"}:
        raise AssertionError("posterior summaries contain non-success rows")

    rows = []
    for (rep, gamma), group in ledger.groupby(["replicate", "gamma"], sort=True):
        child = root / "replicates" / (
            f"rep{int(rep):03d}_g{float(gamma):+g}".replace("+", "p").replace("-", "m")
        )
        child_cfg = json.loads((child / "config.json").read_text())
        if child_cfg.get("h0_prior_support") != [30.0, 110.0]:
            raise AssertionError("child prior support mismatch")
        with np.load(child / "results" / "log_likelihoods.npz") as z:
            grid = np.asarray(z["H0_grid"], dtype=float)
            if len(grid) != int(cfg["h0_points"]):
                raise AssertionError("unexpected H0 grid size")
            for nval in cfg["validation_sizes"]:
                for arm in ARMS:
                    key = f"{arm}_n{int(nval)}"
                    if key not in z:
                        raise AssertionError(f"missing likelihood {key}")
                    logl = np.asarray(z[key], dtype=float)
                    full = posterior(grid, logl, 30.0, 110.0)
                    restricted = posterior(grid, logl, 50.0, 90.0)
                    saved = summaries.loc[
                        (summaries.replicate == rep) & (summaries.gamma == gamma)
                        & (summaries.validation_size == int(nval)) & (summaries.arm == arm)
                    ]
                    if len(saved) != 1:
                        raise AssertionError("posterior summary key is not unique")
                    row = saved.iloc[0]
                    if not np.isclose(full["median"], float(row.posterior_median), rtol=0, atol=2e-7):
                        raise AssertionError("saved posterior median failed likelihood reconstruction")
                    # Integrate only [30,50] and divide by full-support area.
                    low = posterior(grid, logl, 30.0, 110.0)
                    mask = (grid >= 30.0 - 1e-12) & (grid <= 50.0 + 1e-12)
                    low_area = float(trapezoid(np.exp(logl[mask] - np.max(logl)), grid[mask]))
                    full_area = float(trapezoid(np.exp(logl - np.max(logl)), grid))
                    rows.append({
                        "master_seed": expected_seed, "replicate": int(rep), "gamma": float(gamma),
                        "validation_size": int(nval), "arm": arm,
                        "full_median": full["median"], "restricted_median": restricted["median"],
                        "median_shift_30_110_minus_50_90": full["median"] - restricted["median"],
                        "full_q025": full["q025"], "restricted_q025": restricted["q025"],
                        "q025_shift_30_110_minus_50_90": full["q025"] - restricted["q025"],
                        "mass_below_50": low_area / full_area,
                    })
    support = pd.DataFrame(rows)
    if len(support) != 4500:
        raise AssertionError(f"expected 4,500 support reconstructions, got {len(support)}")
    support.to_csv(out / f"support_cells_seed_{expected_seed}.csv", index=False)
    return summaries.assign(master_seed=expected_seed), support


def summarize(summaries: pd.DataFrame, support: pd.DataFrame, out: Path) -> dict:
    summaries = summaries.copy()
    summaries["covered68"] = (summaries.q16 <= summaries.h0_true) & (summaries.h0_true <= summaries.q84)
    summaries["covered90"] = (summaries.q05 <= summaries.h0_true) & (summaries.h0_true <= summaries.q95)
    summaries["covered95"] = (summaries.q025 <= summaries.h0_true) & (summaries.h0_true <= summaries.q975)
    summaries["squared_error"] = summaries.median_error ** 2
    summaries["absolute_error"] = summaries.median_error.abs()
    pooled = summaries.assign(width68=summaries.q84-summaries.q16).groupby(
        ["gamma", "validation_size", "arm"], sort=True
    ).agg(
        n=("median_error", "size"), bias=("median_error", "mean"),
        rmse=("squared_error", lambda x: float(np.sqrt(np.mean(x)))),
        coverage68=("covered68", "mean"), coverage90=("covered90", "mean"),
        coverage95=("covered95", "mean"), mean_width68=("width68", "mean"),
    ).reset_index()
    pooled.to_csv(out / "pooled_metrics_200rep.csv", index=False)

    # Paired arm contrasts within each parent/event replicate.
    wide = summaries.pivot_table(
        index=["master_seed","replicate","gamma","validation_size"],
        columns="arm", values="median_error", aggfunc="first"
    ).reset_index()
    contrast_rows = []
    for gamma in sorted(wide.gamma.unique()):
        for nval in sorted(wide.validation_size.unique()):
            g = wide[(wide.gamma == gamma) & (wide.validation_size == nval)]
            for arm, ref in [("E_uncertainty_integrated_ipw","C_known_probability_ipw"),
                             ("D_map_plugin_ipw","C_known_probability_ipw"),
                             ("E_uncertainty_integrated_ipw","D_map_plugin_ipw")]:
                d = g[arm].to_numpy(float) ** 2 - g[ref].to_numpy(float) ** 2
                m,se,lo,hi = t_interval(d)
                contrast_rows.append({"gamma":gamma,"validation_size":nval,"contrast":f"{arm}_minus_{ref}",
                    "n_pairs":len(d),"mean_squared_error_difference":m,"mcse":se,
                    "ci95_low":lo,"ci95_high":hi})
    contrasts = pd.DataFrame(contrast_rows)
    contrasts.to_csv(out / "paired_squared_error_contrasts.csv", index=False)

    # Paired change in squared error from validation n=100 to n=500 for D/E.
    size_rows = []
    for gamma in sorted(wide.gamma.unique()):
        for arm in ["D_map_plugin_ipw","E_uncertainty_integrated_ipw"]:
            sub = wide[wide.gamma == gamma].pivot(
                index=["master_seed","replicate"], columns="validation_size", values=arm
            )
            d = sub[500].to_numpy(float)**2 - sub[100].to_numpy(float)**2
            m,se,lo,hi = t_interval(d)
            size_rows.append({"gamma":gamma,"arm":arm,"n_pairs":len(d),
                "mean_squared_error_change_n500_minus_n100":m,"mcse":se,"ci95_low":lo,"ci95_high":hi})
    pd.DataFrame(size_rows).to_csv(out / "validation_size_paired_changes.csv", index=False)

    supp = support.groupby(["gamma","validation_size","arm"], sort=True).agg(
        n=("replicate","size"), mean_mass_below_50=("mass_below_50","mean"),
        fraction_mass_below_50_gt_10pct=("mass_below_50",lambda x: float((x>.10).mean())),
        mean_median_shift=("median_shift_30_110_minus_50_90","mean"),
        mean_q025_shift=("q025_shift_30_110_minus_50_90","mean"),
    ).reset_index()
    supp.to_csv(out / "pooled_prior_support_sensitivity_200rep.csv", index=False)

    report = {
        "status":"PASS", "independent_streams":2, "replicates_per_stream":100,
        "pooled_replicates_per_gamma_validation_arm":200,
        "lifecycle_cells_per_stream":4500, "posterior_reconstructions":len(summaries),
        "same_likelihood_prior_support_reconstructions":len(support),
        "interpretation":"Synthetic fixed-truth sensitivity only; not observational H0 inference or universal estimator ranking."
    }
    (out/"audit_summary.json").write_text(json.dumps(report,indent=2)+"\\n")
    # Compact numerical tables are printed to the workflow log as well as saved to CSV.
    print("\\n=== POOLED 200-REPLICATE METRICS (n_validation=500) ===")
    print(pooled[pooled.validation_size == 500].to_string(index=False))
    print("\\n=== PAIRED SQUARED-ERROR CONTRASTS (n_validation=500) ===")
    print(contrasts[contrasts.validation_size == 500].to_string(index=False))
    print("\\n=== VALIDATION-SIZE CHANGE (n=500 minus n=100; paired squared-error) ===")
    print(pd.DataFrame(size_rows).to_string(index=False))
    print("\\n=== PRIOR-SUPPORT SENSITIVITY (n_validation=500) ===")
    print(supp[supp.validation_size == 500].to_string(index=False))
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--seed1940",type=Path,required=True)
    p.add_argument("--seed1929",type=Path,required=True)
    p.add_argument("--auditor",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    a.output_dir.mkdir(parents=True,exist_ok=True)
    s1,p1=audit_one(a.seed1940,SEEDS[0],a.auditor,a.output_dir)
    s2,p2=audit_one(a.seed1929,SEEDS[1],a.auditor,a.output_dir)
    summaries=pd.concat([s1,s2],ignore_index=True)
    support=pd.concat([p1,p2],ignore_index=True)
    summaries.to_csv(a.output_dir/"all_posterior_summaries_200rep.csv",index=False)
    support.to_csv(a.output_dir/"all_support_cells_200rep.csv",index=False)
    report=summarize(summaries,support,a.output_dir)
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
