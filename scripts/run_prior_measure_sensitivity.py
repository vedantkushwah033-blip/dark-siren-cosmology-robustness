"""Compare uniform-H0 and log-uniform-H0 measures in a two-host mixture.

This isolates prior-measure effects on integrated host-component posterior mass.
It is a synthetic low-redshift diagnostic, not a real GW analysis.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.stats import norm

C_KM_S = 299792.458
H0_TRUE = 70.0
Z1 = 0.02
Z2_VALUES = (0.025, 0.04, 0.06)
NOISE_FRACTIONS = (0.03, 0.10, 0.20)
TRUE_WEIGHTS = np.array([0.7, 0.3])
INFERENCE_WEIGHTS = {"correct": np.array([0.7, 0.3]), "reversed": np.array([0.3, 0.7])}
SUPPORTS = {"narrow_50_90": (50.0, 90.0), "broad_20_220": (20.0, 220.0)}
MEASURES = ("uniform_H0", "log_uniform_H0")
GRID = np.linspace(10.0, 250.0, 9601)


def component_masses(distance, sigma, redshifts, weights, support, measure):
    lo, hi = support
    mask = (GRID >= lo) & (GRID <= hi)
    h = GRID[mask]
    if h.size < 2:
        raise ValueError("Prior support has fewer than two grid points")
    likelihoods = np.stack([
        norm.pdf(distance, loc=C_KM_S * z / h, scale=sigma) for z in redshifts
    ])
    prior_density = np.ones_like(h) if measure == "uniform_H0" else 1.0 / h
    if measure not in MEASURES:
        raise ValueError(f"Unknown prior measure: {measure}")
    integrals = np.array([trapezoid(likelihoods[i] * prior_density, h) for i in range(2)])
    contributions = np.asarray(weights, dtype=float) * integrals
    total = float(contributions.sum())
    if not np.isfinite(total) or total <= 0:
        raise ValueError("Component evidence normalization failed")
    return contributions / total


def bootstrap_ci(values, seed, resamples):
    values = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    draws = np.empty(resamples)
    for start in range(0, resamples, 250):
        stop = min(start + 250, resamples)
        idx = rng.integers(0, len(values), size=(stop - start, len(values)))
        draws[start:stop] = values[idx].mean(axis=1)
    return float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def run(replicates=1000, seed=20261029, bootstrap_resamples=4000):
    if replicates < 200 or bootstrap_resamples < 1000:
        raise ValueError("Protocol requires at least 200 replicates and 1000 bootstrap resamples")
    rng = np.random.default_rng(seed)
    raw = []
    cell_index = 0
    summaries = []
    for z2 in Z2_VALUES:
        redshifts = np.array([Z1, z2])
        for noise in NOISE_FRACTIONS:
            sigma = noise * C_KM_S * 0.03 / H0_TRUE
            block = []
            for replicate in range(replicates):
                host = int(rng.choice(2, p=TRUE_WEIGHTS))
                distance_true = C_KM_S * redshifts[host] / H0_TRUE
                distance = float(rng.normal(distance_true, sigma))
                while distance <= 0:
                    distance = float(rng.normal(distance_true, sigma))
                row = {
                    "z1": Z1, "z2": z2, "noise_fraction": noise, "replicate": replicate,
                    "true_host_index": host, "true_host_redshift": float(redshifts[host]),
                    "true_H0": H0_TRUE, "observed_distance_mpc": distance, "sigma_mpc": sigma,
                    "base_seed": seed,
                }
                for support_name, support in SUPPORTS.items():
                    for measure in MEASURES:
                        for weight_name, weights in INFERENCE_WEIGHTS.items():
                            masses = component_masses(distance, sigma, redshifts, weights, support, measure)
                            key = f"{support_name}_{measure}_{weight_name}"
                            row[f"{key}_host1_mass"] = float(masses[0])
                            row[f"{key}_host2_mass"] = float(masses[1])
                raw.append(row)
                block.append(row)
            for support_index, (support_name, support) in enumerate(SUPPORTS.items()):
                for measure_index, measure in enumerate(MEASURES):
                    correct_key = f"{support_name}_{measure}_correct"
                    reversed_key = f"{support_name}_{measure}_reversed"
                    diffs = np.array([
                        r[f"{reversed_key}_host2_mass"] - r[f"{correct_key}_host2_mass"]
                        for r in block
                    ])
                    bootstrap_seed = seed + 104729 * cell_index + 1009 * support_index + 37 * measure_index
                    lo, hi = bootstrap_ci(diffs, bootstrap_seed, bootstrap_resamples)
                    summaries.append({
                        "z1": Z1, "z2": z2, "noise_fraction": noise,
                        "prior_support": support_name, "prior_min": support[0], "prior_max": support[1],
                        "prior_measure": measure, "replicates": replicates,
                        "correct_mean_host2_mass": float(np.mean([r[f"{correct_key}_host2_mass"] for r in block])),
                        "reversed_mean_host2_mass": float(np.mean([r[f"{reversed_key}_host2_mass"] for r in block])),
                        "mean_paired_host2_mass_difference_reversed_minus_correct": float(diffs.mean()),
                        "ci95_low": lo, "ci95_high": hi,
                        "bootstrap_resamples": bootstrap_resamples, "bootstrap_seed": bootstrap_seed,
                    })
                cell_index += 1
    return raw, summaries


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--bootstrap-resamples", type=int, default=4000)
    parser.add_argument("--output-dir", default="results/prior-measure-sensitivity")
    args = parser.parse_args()
    raw, summaries = run(args.replicates, args.seed, args.bootstrap_resamples)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_csv(out / "per_run.csv", raw)
    write_csv(out / "summary.csv", summaries)
    metadata = {
        "experiment": "uniform versus log-uniform H0 prior measure on host component masses",
        "audit_status": "PENDING", "base_seed": args.seed,
        "replicates_per_redshift_noise_cell": args.replicates,
        "bootstrap_resamples": args.bootstrap_resamples,
        "raw_records": len(raw), "summary_rows": len(summaries),
        "sha256": {name: hashlib.sha256((out / name).read_bytes()).hexdigest()
                   for name in ("per_run.csv", "summary.csv")},
    }
    (out / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Generated {len(raw)} observations and {len(summaries)} paired prior-measure summaries.")


if __name__ == "__main__":
    main()
