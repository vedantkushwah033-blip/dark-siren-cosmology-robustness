"""Independent adaptive-quadrature reference for grid-based two-host posterior summaries."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary

PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
GRID_SIZES = (1001, 4001, 16001)
CASES = (
    (52.0, 0.025, 0.03, "correct"),
    (52.0, 0.060, 0.10, "reversed"),
    (70.0, 0.025, 0.10, "correct"),
    (70.0, 0.040, 0.03, "reversed"),
    (70.0, 0.060, 0.10, "correct"),
    (88.0, 0.025, 0.03, "reversed"),
    (88.0, 0.060, 0.10, "correct"),
)
WEIGHTS = {"correct": (0.7, 0.3), "reversed": (0.3, 0.7)}


def reference_summary(redshifts, distance, sigma, weights):
    """Compute normalized posterior quantiles using adaptive quadrature, not a grid CDF."""
    z = np.asarray(redshifts, dtype=float)
    w = np.asarray(weights, dtype=float)

    def likelihood(h):
        means = C_KM_S * z / h
        return float(np.dot(w, norm.pdf(distance, loc=means, scale=sigma)))

    norm_const, err = quad(likelihood, PRIOR_MIN, PRIOR_MAX, epsabs=1e-12, epsrel=1e-10, limit=300)
    if not np.isfinite(norm_const) or norm_const <= 0:
        raise ValueError("Adaptive quadrature failed to normalize posterior.")

    def cdf(h):
        if h <= PRIOR_MIN:
            return 0.0
        if h >= PRIOR_MAX:
            return 1.0
        area, _ = quad(likelihood, PRIOR_MIN, h, epsabs=1e-12, epsrel=1e-10, limit=300)
        return area / norm_const

    quantiles = {}
    for name, probability in (("lower_68", 0.16), ("median", 0.5), ("upper_68", 0.84),
                              ("lower_95", 0.025), ("upper_95", 0.975)):
        quantiles[name] = float(brentq(lambda h: cdf(h) - probability, PRIOR_MIN, PRIOR_MAX,
                                       xtol=1e-10, rtol=1e-12))
    quantiles["normalization_integral"] = float(norm_const)
    quantiles["quadrature_error_estimate"] = float(err)
    return quantiles


def run_independent_reference_audit():
    rows = []
    for h0_true, z2, noise_fraction, weighting in CASES:
        redshifts = np.array([0.02, z2])
        sigma = noise_fraction * C_KM_S * 0.03 / h0_true
        # Deterministic observations exercise both host modes and avoid Monte Carlo
        # noise obscuring a purely numerical comparison.
        for host in (0, 1):
            distance = C_KM_S * redshifts[host] / h0_true + (0.35 if host == 0 else -0.35) * sigma
            weights = WEIGHTS[weighting]
            ref = reference_summary(redshifts, distance, sigma, weights)
            for size in GRID_SIZES:
                grid = np.linspace(25.0, 125.0, size)
                density = posterior_h0(grid, redshifts, distance, sigma, weights=weights,
                                       prior_min=PRIOR_MIN, prior_max=PRIOR_MAX)
                got = posterior_summary(grid, density)
                rows.append({
                    "H0_true": h0_true, "z2": z2, "noise_fraction": noise_fraction,
                    "weighting": weighting, "injected_host_index": host, "grid_size": size,
                    "distance_mpc": distance, "reference_normalization": ref["normalization_integral"],
                    "quadrature_error_estimate": ref["quadrature_error_estimate"],
                    **{f"reference_{k}": ref[k] for k in ("median", "lower_68", "upper_68", "lower_95", "upper_95")},
                    **{f"grid_{k}": got[k] for k in ("median", "lower_68", "upper_68", "lower_95", "upper_95")},
                    **{f"abs_error_{k}": abs(got[k] - ref[k]) for k in ("median", "lower_68", "upper_68", "lower_95", "upper_95")},
                })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="results/independent-reference-audit")
    args = parser.parse_args()
    rows = run_independent_reference_audit()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "grid_vs_quadrature.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"PASS: {len(rows)} grid summaries compared with adaptive-quadrature reference.")
    for size in GRID_SIZES:
        selected = [r for r in rows if r["grid_size"] == size]
        for metric in ("median", "lower_68", "upper_68", "lower_95", "upper_95"):
            errors = np.array([r[f"abs_error_{metric}"] for r in selected])
            print(f"grid={size} metric={metric} median_abs_error={np.median(errors):.8g} "
                  f"p95_abs_error={np.quantile(errors, .95):.8g} max_abs_error={np.max(errors):.8g}")


if __name__ == "__main__":
    main()
