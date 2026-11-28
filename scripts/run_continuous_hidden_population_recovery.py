"""Paired parameter-recovery pilot for the continuous hidden-host model.

This is a synthetic log-distance analogue, separate from the existing
positive-truncated Gaussian-distance experiments.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.stats import truncnorm

from dark_siren_robustness.hidden_population import (
    C_KM_S, posterior_grid, posterior_marginals,
)

TRUTH = {"H0": 70.0, "fraction": 0.30, "mu_h": 0.075, "sigma_h": 0.012}
LISTED_Z = (0.02, 0.04, 0.06)
LISTED_W = (1/3, 1/3, 1/3)
Z_BOUNDS = (0.005, 0.10)
SIGMA_U = 0.10


def generate_dataset(n_events, rng):
    """Generate one paired dataset under the stated continuous-mixture truth."""
    n = int(n_events)
    if n < 1:
        raise ValueError("n_events must be positive")
    labels = rng.choice(4, size=n, p=[
        (1-TRUTH["fraction"])/3,
        (1-TRUTH["fraction"])/3,
        (1-TRUTH["fraction"])/3,
        TRUTH["fraction"],
    ])
    z = np.empty(n, dtype=float)
    listed = np.asarray(LISTED_Z)
    is_listed = labels < 3
    z[is_listed] = listed[labels[is_listed]]
    a, b = Z_BOUNDS
    aa = (a-TRUTH["mu_h"])/TRUTH["sigma_h"]
    bb = (b-TRUTH["mu_h"])/TRUTH["sigma_h"]
    z[~is_listed] = truncnorm.rvs(
        aa, bb, loc=TRUTH["mu_h"], scale=TRUTH["sigma_h"],
        size=np.count_nonzero(~is_listed), random_state=rng,
    )
    u = rng.normal(
        np.log(C_KM_S*z/TRUTH["H0"]), SIGMA_U, size=n
    )
    return u, z, labels


def _h0_summary(u, h0_grid, f_grid, mu_grid, sigma_grid, quadrature_order):
    post = posterior_grid(
        u, SIGMA_U, h0_grid, f_grid, mu_grid, sigma_grid,
        listed_redshifts=LISTED_Z, listed_weights=LISTED_W,
        z_bounds=Z_BOUNDS, quadrature_order=quadrature_order,
    )
    return posterior_marginals(post)["H0"]


def _interval_metrics(summary, truth):
    return {
        "median": summary["median"],
        "lower_68": summary["lower_68"],
        "upper_68": summary["upper_68"],
        "lower_95": summary["lower_95"],
        "upper_95": summary["upper_95"],
        "error": summary["median"]-truth,
        "covered_68": int(summary["lower_68"] <= truth <= summary["upper_68"]),
        "covered_95": int(summary["lower_95"] <= truth <= summary["upper_95"]),
        "width_68": summary["upper_68"]-summary["lower_68"],
        "width_95": summary["upper_95"]-summary["lower_95"],
    }


def run(replicates=2, sample_sizes=(30,), seed=20261128, h0_points=11,
        fraction_points=7, mu_points=9, quadrature_order=24,
        estimate_sigma_h=False, sigma_h_points=5):
    if replicates < 1 or any(int(n) < 1 for n in sample_sizes):
        raise ValueError("replicates and sample sizes must be positive")
    if h0_points < 5 or fraction_points < 3 or mu_points < 3:
        raise ValueError("grids are too small for a meaningful smoke test")
    rng = np.random.default_rng(seed)
    h0_grid = np.linspace(50, 90, h0_points)
    f_grid = np.linspace(0, 1, fraction_points)
    mu_grid = np.linspace(0.03, 0.095, mu_points)
    sigma_grid = (
        np.linspace(0.006, 0.020, sigma_h_points)
        if estimate_sigma_h else [TRUTH["sigma_h"]]
    )
    raw = []

    for n in sample_sizes:
        for rep in range(replicates):
            u, _, _ = generate_dataset(n, rng)
            # Fully integrated joint posterior, with sigma_h conditioned at its
            # generating value for this first recovery stage.
            joint = posterior_grid(
                u, SIGMA_U, h0_grid, f_grid, mu_grid, sigma_grid,
                listed_redshifts=LISTED_Z, listed_weights=LISTED_W,
                z_bounds=Z_BOUNDS, quadrature_order=quadrature_order,
            )
            marg = posterior_marginals(joint)
            integrated = marg["H0"]

            # Empirical-Bayes plug-in: posterior medians of f and mu_h from the
            # joint posterior, followed by an H0 posterior conditional on them.
            plugin_f = marg["fraction"]["median"]
            plugin_mu = marg["mu_h"]["median"]
            plugin_sh = marg["sigma_h"]["median"]
            plugin = _h0_summary(
                u, h0_grid, [plugin_f], [plugin_mu], [plugin_sh], quadrature_order
            )

            oracle = _h0_summary(
                u, h0_grid, [TRUTH["fraction"]], [TRUTH["mu_h"]],
                sigma_grid, quadrature_order,
            )
            fixed_wrong = _h0_summary(
                u, h0_grid, [0.10], [0.060], sigma_grid, quadrature_order
            )

            for method, summary in (
                ("oracle", oracle),
                ("fixed_wrong_f0.10_mu0.060", fixed_wrong),
                ("plug_in_joint_medians", plugin),
                ("fully_integrated_f_mu", integrated),
            ):
                row = {
                    "n_events": int(n), "replicate": int(rep), "seed": int(seed),
                    "method": method, "true_H0": TRUTH["H0"],
                    "true_fraction": TRUTH["fraction"], "true_mu_h": TRUTH["mu_h"],
                    "sigma_h_conditioned": TRUTH["sigma_h"], "sigma_u": SIGMA_U,
                    "plugin_fraction": plugin_f if method == "plug_in_joint_medians" else "",
                    "plugin_mu_h": plugin_mu if method == "plug_in_joint_medians" else "",
                    "plugin_sigma_h": plugin_sh if method == "plug_in_joint_medians" else "",
                    "sigma_h_estimated": bool(estimate_sigma_h),
                }
                row.update(_interval_metrics(summary, TRUTH["H0"]))
                raw.append(row)
    return raw


def _wilson(successes, n, z=1.959963984540054):
    p = successes/n
    den = 1 + z*z/n
    center = (p + z*z/(2*n))/den
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return center-half, center+half


def summarize(raw):
    groups = {}
    for row in raw:
        groups.setdefault((row["n_events"], row["method"]), []).append(row)
    summaries = []
    for (n_events, method), rows in sorted(groups.items()):
        errors = np.asarray([r["error"] for r in rows], dtype=float)
        c68 = int(sum(r["covered_68"] for r in rows))
        c95 = int(sum(r["covered_95"] for r in rows))
        lo68, hi68 = _wilson(c68, len(rows))
        lo95, hi95 = _wilson(c95, len(rows))
        summaries.append({
            "n_events": n_events, "method": method, "replicates": len(rows),
            "bias": float(np.mean(errors)),
            "bias_mcse": float(np.std(errors, ddof=1)/np.sqrt(len(rows))) if len(rows)>1 else np.nan,
            "rmse": float(np.sqrt(np.mean(errors**2))),
            "rmse_mcse_delta": (
                float(np.std(errors**2, ddof=1)/(2*np.sqrt(len(rows))*max(np.sqrt(np.mean(errors**2)),1e-15)))
                if len(rows)>1 else np.nan
            ),
            "coverage_68": c68/len(rows), "coverage_68_wilson_low": lo68,
            "coverage_68_wilson_high": hi68,
            "coverage_95": c95/len(rows), "coverage_95_wilson_low": lo95,
            "coverage_95_wilson_high": hi95,
            "mean_width_68": float(np.mean([r["width_68"] for r in rows])),
            "mean_width_95": float(np.mean([r["width_95"] for r in rows])),
        })
    return summaries


def paired_contrasts(raw, bootstrap_resamples=2000, seed=20261130):
    """Paired method-minus-oracle contrasts using common replicate indices."""
    if bootstrap_resamples < 100:
        raise ValueError("Use at least 100 paired bootstrap resamples")
    grouped = {}
    for row in raw:
        grouped.setdefault((row["n_events"], row["replicate"]), {})[row["method"]] = row
    by_n = {}
    for (n_events, rep), methods in grouped.items():
        if "oracle" not in methods:
            raise ValueError(f"Missing oracle result for n={n_events}, replicate={rep}")
        by_n.setdefault(n_events, []).append((rep, methods))
    rng = np.random.default_rng(seed)
    metrics = {
        "median_error": lambda r: r["error"],
        "squared_error": lambda r: r["error"]**2,
        "coverage_68": lambda r: r["covered_68"],
        "coverage_95": lambda r: r["covered_95"],
        "width_68": lambda r: r["width_68"],
        "width_95": lambda r: r["width_95"],
    }
    output = []
    for n_events, pairs in sorted(by_n.items()):
        pairs.sort(key=lambda x: x[0])
        methods = sorted(set.intersection(*[set(m) for _, m in pairs]) - {"oracle"})
        count = len(pairs)
        indices = rng.integers(0, count, size=(bootstrap_resamples, count))
        for method in methods:
            for metric, getter in metrics.items():
                differences = np.asarray([
                    getter(m[method])-getter(m["oracle"]) for _, m in pairs
                ], dtype=float)
                boot = np.mean(differences[indices], axis=1)
                output.append({
                    "n_events": n_events, "method_vs_oracle": method,
                    "metric": metric, "estimate": float(np.mean(differences)),
                    "mcse": float(np.std(differences, ddof=1)/np.sqrt(count)) if count>1 else np.nan,
                    "paired_bootstrap_95_low": float(np.quantile(boot, 0.025)),
                    "paired_bootstrap_95_high": float(np.quantile(boot, 0.975)),
                    "replicates": count, "bootstrap_resamples": bootstrap_resamples,
                    "bootstrap_seed": seed,
                })
    return output


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"Refusing to write empty output: {path}")
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=2,
                        help="Use 2 for smoke; choose >=200 for a confirmatory study.")
    parser.add_argument("--sample-sizes", default="30,50,100,200")
    parser.add_argument("--seed", type=int, default=20261128)
    parser.add_argument("--h0-points", type=int, default=41)
    parser.add_argument("--fraction-points", type=int, default=15)
    parser.add_argument("--mu-points", type=int, default=21)
    parser.add_argument("--quadrature-order", type=int, default=48)
    parser.add_argument("--estimate-sigma-h", action="store_true",
                        help="Jointly estimate hidden redshift scale instead of conditioning it at truth.")
    parser.add_argument("--sigma-h-points", type=int, default=5)
    parser.add_argument("--bootstrap-resamples", type=int, default=2000)
    parser.add_argument("--output-dir", default="results/continuous-hidden-population-recovery")
    args = parser.parse_args()
    sizes = tuple(int(x) for x in args.sample_sizes.split(","))
    raw = run(args.replicates, sizes, args.seed, args.h0_points,
              args.fraction_points, args.mu_points, args.quadrature_order,
              args.estimate_sigma_h, args.sigma_h_points)
    summary = summarize(raw)
    contrasts = paired_contrasts(raw, args.bootstrap_resamples, args.seed+1)
    root = Path(args.output_dir)
    root.mkdir(parents=True, exist_ok=True)
    write_csv(root/"per_run.csv", raw)
    write_csv(root/"summary.csv", summary)
    write_csv(root/"paired_contrasts.csv", contrasts)
    print(f"Wrote {len(raw)} paired posterior rows and {len(summary)} summaries to {root}")


if __name__ == "__main__":
    main()
