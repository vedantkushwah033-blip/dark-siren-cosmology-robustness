"""Seeded randomized stress test of grid summaries against adaptive quadrature."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from run_independent_reference_audit import reference_summary
from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary

PRIOR_MIN, PRIOR_MAX = 50.0, 90.0


def run_stress_audit(cases=120, seed=20261024):
    if cases < 30:
        raise ValueError("Use at least 30 randomized cases.")
    rng = np.random.default_rng(seed)
    rows = []
    for case in range(cases):
        # Deliberately include truths near prior edges as well as interior values.
        h0_true = float(rng.choice([50.5, 52.0, 55.0, 70.0, 85.0, 88.0, 89.5])
                        if case % 3 == 0 else rng.uniform(50.2, 89.8))
        redshifts = np.sort(rng.uniform(0.012, 0.085, size=3))
        true_weights = rng.dirichlet([1.5, 1.5, 1.5])
        inference_weights = rng.dirichlet([1.0, 1.0, 1.0])
        host = int(rng.choice(3, p=true_weights))
        # Vary noise over a broad range relative to a reference low-z distance.
        noise_fraction = float(rng.uniform(0.015, 0.25))
        sigma = noise_fraction * C_KM_S * 0.03 / h0_true
        distance = float(rng.normal(C_KM_S * redshifts[host] / h0_true, sigma))
        while distance <= 0:
            distance = float(rng.normal(C_KM_S * redshifts[host] / h0_true, sigma))
        ref = reference_summary(redshifts, distance, sigma, inference_weights)
        for size in (1001, 4001):
            grid = np.linspace(25.0, 125.0, size)
            density = posterior_h0(grid, redshifts, distance, sigma, weights=inference_weights,
                                   prior_min=PRIOR_MIN, prior_max=PRIOR_MAX)
            got = posterior_summary(grid, density)
            rows.append({
                "case": case, "H0_true": h0_true, "redshifts": ";".join(map(str, redshifts)),
                "true_weights": ";".join(map(str, true_weights)),
                "inference_weights": ";".join(map(str, inference_weights)),
                "generated_host_index": host, "noise_fraction": noise_fraction,
                "sigma_mpc": sigma, "distance_mpc": distance, "grid_size": size,
                "reference_median": ref["median"], "grid_median": got["median"],
                "abs_error_median": abs(got["median"] - ref["median"]),
                **{f"abs_error_{k}": abs(got[k] - ref[k]) for k in
                   ("lower_68", "upper_68", "lower_95", "upper_95")},
            })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=int, default=120)
    parser.add_argument("--seed", type=int, default=20261024)
    parser.add_argument("--output-dir", default="results/randomized-numerical-stress")
    args = parser.parse_args()
    rows = run_stress_audit(args.cases, args.seed)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "randomized_grid_vs_quadrature.csv"
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"PASS: {len(rows)} grid summaries across {args.cases} seeded randomized cases.")
    for size in (1001, 4001):
        subset = [r for r in rows if r["grid_size"] == size]
        for metric in ("median", "lower_68", "upper_68", "lower_95", "upper_95"):
            errors = np.array([r[f"abs_error_{metric}"] for r in subset])
            print(f"grid={size} metric={metric} median_abs_error={np.median(errors):.8g} "
                  f"p95_abs_error={np.quantile(errors, .95):.8g} max_abs_error={np.max(errors):.8g}")


if __name__ == "__main__":
    main()
