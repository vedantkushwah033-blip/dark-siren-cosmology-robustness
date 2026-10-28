"""Independent reconstruction audit for the paired host sensitivity surface."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.stats import norm

C = 299792.458
TRUTHS = (60.0, 70.0, 80.0)
NOISES = (0.05, 0.10, 0.20)
FRACTIONS = (0.0, 0.15, 0.30, 0.45)
HIDDEN_ZS = (0.025, 0.030, 0.035)
Z_LISTED = np.array([0.02, 0.04])
ARMS = ("base_weights", "reversed_weights", "equal_weights", "hidden_host_modelled")
COMPARISONS = (
    ("reversed_weights", "base_weights", "reversed_minus_base"),
    ("base_weights", "hidden_host_modelled", "omitted_minus_modelled_hidden"),
)
METRICS = ("mean_error", "rmse_error", "coverage_68", "coverage_95", "mean_width_68", "mean_width_95")
GRID = np.linspace(25.0, 125.0, 4001)
PRIOR = (50.0, 90.0)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def cell_specs(truths=TRUTHS, noises=NOISES, fractions=FRACTIONS, hidden_zs=HIDDEN_ZS):
    return [(float(t), float(n), float(f), float(z))
            for t in truths for n in noises for f in fractions for z in hidden_zs]


def regenerate_cell(truth, noise, fraction, hidden_z, replicates, rng):
    sigma = noise * C * 0.03 / truth
    z_values = np.array([0.02, 0.04, hidden_z])
    probabilities = np.array([(1 - fraction) * 0.7, (1 - fraction) * 0.3, fraction])
    observations = []
    for rep in range(replicates):
        host = int(rng.choice(3, p=probabilities))
        z = float(z_values[host])
        distance = float(rng.normal(C * z / truth, sigma))
        while distance <= 0:
            distance = float(rng.normal(C * z / truth, sigma))
        observations.append((rep, host, z, distance, sigma))
    return observations


def posterior_from_equation(distance, sigma, truth, fraction, hidden_z, arm):
    if arm == "base_weights":
        weights, missing = np.array([0.7, 0.3]), 0.0
    elif arm == "reversed_weights":
        weights, missing = np.array([0.3, 0.7]), 0.0
    elif arm == "equal_weights":
        weights, missing = np.array([0.5, 0.5]), 0.0
    elif arm == "hidden_host_modelled":
        weights, missing = np.array([0.7, 0.3]), fraction
    else:
        raise ValueError(f"Unexpected arm {arm}")
    mu = C * Z_LISTED[None, :] / GRID[:, None]
    density_components = norm.pdf(distance, loc=mu, scale=sigma) / norm.cdf(mu / sigma)
    likelihood = (1.0 - missing) * np.sum(weights[None, :] * density_components, axis=1)
    if missing > 0:
        mu_hidden = C * hidden_z / GRID
        hidden_density = norm.pdf(distance, loc=mu_hidden, scale=sigma) / norm.cdf(mu_hidden / sigma)
        likelihood += missing * hidden_density
    in_prior = (GRID >= PRIOR[0]) & (GRID <= PRIOR[1])
    x = GRID[in_prior]
    density = likelihood[in_prior]
    norm_const = trapezoid(density, x)
    if not np.isfinite(norm_const) or norm_const <= 0:
        raise ValueError("Invalid posterior normalization")
    density = density / norm_const
    cdf = np.concatenate(([0.0], np.cumsum((density[1:] + density[:-1]) * np.diff(x) / 2.0)))
    cdf = cdf / cdf[-1]
    quantile = {p: float(np.interp(p, cdf, x)) for p in (0.025, 0.16, 0.5, 0.84, 0.975)}
    return {
        "posterior_median": quantile[0.5], "median_error": quantile[0.5] - truth,
        "lower_68": quantile[0.16], "upper_68": quantile[0.84],
        "lower_95": quantile[0.025], "upper_95": quantile[0.975],
        "width_68": quantile[0.84] - quantile[0.16],
        "width_95": quantile[0.975] - quantile[0.025],
        "covered_68": int(quantile[0.16] <= truth <= quantile[0.84]),
        "covered_95": int(quantile[0.025] <= truth <= quantile[0.975]),
    }


def stats(records, indices=None):
    if indices is None:
        indices = np.arange(len(records))
    error = np.array([float(r["median_error"]) for r in records])[indices]
    c68 = np.array([int(r["covered_68"]) for r in records], dtype=float)[indices]
    c95 = np.array([int(r["covered_95"]) for r in records], dtype=float)[indices]
    w68 = np.array([float(r["width_68"]) for r in records])[indices]
    w95 = np.array([float(r["width_95"]) for r in records])[indices]
    return np.array([
        np.mean(error, axis=-1), np.sqrt(np.mean(error**2, axis=-1)),
        np.mean(c68, axis=-1), np.mean(c95, axis=-1),
        np.mean(w68, axis=-1), np.mean(w95, axis=-1),
    ]).T


def audit(root: Path, replicates=200, seed=20261029, resamples=1000,
          truths=TRUTHS, noises=NOISES, fractions=FRACTIONS, hidden_zs=HIDDEN_ZS, atol=2e-8):
    if replicates < 200 or resamples < 100:
        raise ValueError("Audit design below protocol minimum")
    raw = read_csv(root / "per_run.csv")
    actual_summary = read_csv(root / "summary.csv")
    actual_contrasts = read_csv(root / "paired_contrasts.csv")
    specs = cell_specs(truths, noises, fractions, hidden_zs)
    expected_n = len(specs) * replicates * len(ARMS)
    if len(raw) != expected_n or len(actual_summary) != len(specs) * len(ARMS):
        raise ValueError("Raw or summary row count does not match factorial design")
    if len(actual_contrasts) != len(specs) * len(COMPARISONS) * len(METRICS):
        raise ValueError("Paired contrast row count does not match factorial design")
    data_rng = np.random.default_rng(seed)
    boot_rng = np.random.default_rng(seed + 1)
    grouped = {}
    expected_summary = []
    expected_contrasts = []
    raw_i = 0
    summary_i = 0
    for truth, noise, fraction, hidden_z in specs:
        generated = regenerate_cell(truth, noise, fraction, hidden_z, replicates, data_rng)
        arms_data = {arm: [] for arm in ARMS}
        for rep, host, z, distance, sigma in generated:
            block = raw[raw_i:raw_i + len(ARMS)]
            raw_i += len(ARMS)
            for arm, row in zip(ARMS, block):
                if row["inference_arm"] != arm or int(row["replicate"]) != rep:
                    raise ValueError(f"Pair/row ordering mismatch in cell {(truth, noise, fraction, hidden_z)}")
                expected_observation = {
                    "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
                    "hidden_redshift": hidden_z, "true_host_index": host, "true_host_redshift": z,
                    "observed_distance_mpc": distance, "sigma_mpc": sigma,
                }
                for key, value in expected_observation.items():
                    if not math.isclose(float(row[key]), float(value), rel_tol=2e-13, abs_tol=2e-13):
                        raise ValueError(f"Seeded data provenance mismatch at row {raw_i-len(ARMS)}, field {key}")
                posterior = posterior_from_equation(distance, sigma, truth, fraction, hidden_z, arm)
                for key, value in posterior.items():
                    if key not in row or not math.isclose(float(row[key]), float(value), rel_tol=2e-10, abs_tol=atol):
                        raise ValueError(f"Posterior reconstruction mismatch at row {raw_i-len(ARMS)}, field {key}")
                arms_data[arm].append(posterior)
        if fraction == 0.0:
            for a, b in zip(arms_data["base_weights"], arms_data["hidden_host_modelled"]):
                for field in ("posterior_median", "lower_68", "upper_68", "lower_95", "upper_95"):
                    if not math.isclose(a[field], b[field], rel_tol=0, abs_tol=1e-10):
                        raise ValueError("f=0 negative control failed: hidden-host arm differs from base arm")
        for arm in ARMS:
            values = stats(arms_data[arm])
            expected_summary.append({
                "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
                "hidden_redshift": hidden_z, "inference_arm": arm, "replicates": replicates,
                "mean_error": float(values[0]), "rmse_error": float(values[1]),
                "coverage_68": float(values[2]), "coverage_95": float(values[3]),
                "mean_width_68": float(values[4]), "mean_width_95": float(values[5]),
            })
        indices = boot_rng.integers(0, replicates, size=(resamples, replicates))
        for left, right, label in COMPARISONS:
            point = stats(arms_data[left]) - stats(arms_data[right])
            boot = stats(arms_data[left], indices) - stats(arms_data[right], indices)
            for j, metric in enumerate(METRICS):
                expected_contrasts.append({
                    "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
                    "hidden_redshift": hidden_z, "contrast": label, "left_arm": left,
                    "right_arm": right, "metric": metric, "estimate": float(point[j]),
                    "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                    "ci95_high": float(np.quantile(boot[:, j], 0.975)),
                    "replicates": replicates, "bootstrap_resamples": resamples,
                    "bootstrap_seed": seed + 1,
                })
    def compare(expected_rows, actual_rows, keys, fields, label):
        actual_map = {tuple(row[k] for k in keys): row for row in actual_rows}
        if len(actual_map) != len(actual_rows):
            raise ValueError(f"Duplicate keys in {label}")
        for expected in expected_rows:
            key = tuple(str(expected[k]) for k in keys)
            if key not in actual_map:
                raise ValueError(f"Missing {label} row {key}")
            row = actual_map[key]
            for field in fields:
                if not math.isclose(float(row[field]), float(expected[field]), rel_tol=2e-10, abs_tol=atol):
                    raise ValueError(f"{label} reconstruction mismatch {key}, field {field}")
        if len(actual_map) != len(expected_rows):
            raise ValueError(f"Unexpected extra rows in {label}")
    compare(expected_summary, actual_summary,
            ("true_H0", "noise_fraction", "hidden_fraction", "hidden_redshift", "inference_arm"),
            ("replicates", "mean_error", "rmse_error", "coverage_68", "coverage_95", "mean_width_68", "mean_width_95"),
            "summary")
    compare(expected_contrasts, actual_contrasts,
            ("true_H0", "noise_fraction", "hidden_fraction", "hidden_redshift", "contrast", "metric"),
            ("estimate", "ci95_low", "ci95_high", "replicates", "bootstrap_resamples", "bootstrap_seed"),
            "paired contrast")
    return {"audit": "PASS", "factorial_cells": len(specs), "observations": len(specs) * replicates,
            "posterior_reconstructions": expected_n, "summary_rows_reconstructed": len(expected_summary),
            "paired_contrast_rows_reconstructed": len(expected_contrasts), "bootstrap_resamples": resamples}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="results/host-misspecification-sensitivity-surface")
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--bootstrap-resamples", type=int, default=1000)
    args = parser.parse_args()
    print(audit(Path(args.input_dir), args.replicates, args.seed, args.bootstrap_resamples))


if __name__ == "__main__":
    main()
