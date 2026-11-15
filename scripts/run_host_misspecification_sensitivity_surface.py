"""Paired factorial sensitivity study for host-population misspecification."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary

TRUTHS = (60.0, 70.0, 80.0)
NOISE_LEVELS = (0.05, 0.10, 0.20)
HIDDEN_FRACTIONS = (0.0, 0.15, 0.30, 0.45)
HIDDEN_REDSHIFTS = (0.025, 0.030, 0.035)
BASE_REDSHIFTS = np.array([0.02, 0.04])
ARMS = ("base_weights", "reversed_weights", "equal_weights", "hidden_host_modelled")
COMPARISONS = (
    ("reversed_weights", "base_weights", "reversed_minus_base"),
    ("base_weights", "hidden_host_modelled", "omitted_minus_modelled_hidden"),
)
METRICS = ("mean_error", "rmse_error", "coverage_68", "coverage_95", "mean_width_68", "mean_width_95")
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
GRID = np.linspace(25.0, 125.0, 4001)


def cell_specs(truths=TRUTHS, noise_levels=NOISE_LEVELS,
               hidden_fractions=HIDDEN_FRACTIONS, hidden_redshifts=HIDDEN_REDSHIFTS):
    for truth in truths:
        for noise in noise_levels:
            for fraction in hidden_fractions:
                for hidden_z in hidden_redshifts:
                    yield float(truth), float(noise), float(fraction), float(hidden_z)


def generate_cell(truth, noise_fraction, hidden_fraction, hidden_z, replicates, rng):
    sigma = noise_fraction * C_KM_S * 0.03 / truth
    population_z = np.array([0.02, 0.04, hidden_z])
    population_weights = np.array([
        (1.0 - hidden_fraction) * 0.7,
        (1.0 - hidden_fraction) * 0.3,
        hidden_fraction,
    ])
    observations = []
    for replicate in range(replicates):
        host_index = int(rng.choice(3, p=population_weights))
        z_true = float(population_z[host_index])
        distance_true = C_KM_S * z_true / truth
        observed = float(rng.normal(distance_true, sigma))
        while observed <= 0:
            observed = float(rng.normal(distance_true, sigma))
        observations.append({
            "true_H0": truth, "noise_fraction": noise_fraction,
            "hidden_fraction": hidden_fraction, "hidden_redshift": hidden_z,
            "replicate": replicate, "true_host_index": host_index,
            "true_host_redshift": z_true, "observed_distance_mpc": observed,
            "sigma_mpc": sigma,
        })
    return observations


def infer(observation, arm):
    d, sigma = observation["observed_distance_mpc"], observation["sigma_mpc"]
    hidden_fraction, hidden_z = observation["hidden_fraction"], observation["hidden_redshift"]
    common = {"prior_min": PRIOR_MIN, "prior_max": PRIOR_MAX, "positive_truncated": True}
    if arm == "base_weights":
        density = posterior_h0(GRID, BASE_REDSHIFTS, d, sigma, weights=[0.7, 0.3], **common)
    elif arm == "reversed_weights":
        density = posterior_h0(GRID, BASE_REDSHIFTS, d, sigma, weights=[0.3, 0.7], **common)
    elif arm == "equal_weights":
        density = posterior_h0(GRID, BASE_REDSHIFTS, d, sigma, weights=[0.5, 0.5], **common)
    elif arm == "hidden_host_modelled":
        density = posterior_h0(GRID, BASE_REDSHIFTS, d, sigma, weights=[0.7, 0.3],
                               missing_weight=hidden_fraction, missing_redshifts=[hidden_z], **common)
    else:
        raise ValueError(f"Unknown inference arm: {arm}")
    return posterior_summary(GRID, density)


def _metrics(rows, indices=None):
    if indices is None:
        indices = np.arange(len(rows))
    errors = np.asarray([float(r["median_error"]) for r in rows])[indices]
    covered68 = np.asarray([int(r["covered_68"]) for r in rows], dtype=float)[indices]
    covered95 = np.asarray([int(r["covered_95"]) for r in rows], dtype=float)[indices]
    width68 = np.asarray([float(r["width_68"]) for r in rows])[indices]
    width95 = np.asarray([float(r["width_95"]) for r in rows])[indices]
    return np.stack((
        np.mean(errors, axis=-1),
        np.sqrt(np.mean(errors**2, axis=-1)),
        np.mean(covered68, axis=-1),
        np.mean(covered95, axis=-1),
        np.mean(width68, axis=-1),
        np.mean(width95, axis=-1),
    ), axis=-1)


def run_experiment(replicates=200, seed=20261029, bootstrap_resamples=1000,
                   truths=TRUTHS, noise_levels=NOISE_LEVELS,
                   hidden_fractions=HIDDEN_FRACTIONS, hidden_redshifts=HIDDEN_REDSHIFTS):
    if replicates < 200:
        raise ValueError("Use at least 200 observations per factorial cell.")
    if bootstrap_resamples < 100:
        raise ValueError("Use at least 100 paired bootstrap resamples.")
    specs = list(cell_specs(truths, noise_levels, hidden_fractions, hidden_redshifts))
    if not specs:
        raise ValueError("Factorial design must contain at least one cell.")
    data_rng = np.random.default_rng(seed)
    bootstrap_rng = np.random.default_rng(seed + 1)
    raw, summary, contrasts = [], [], []
    for truth, noise, fraction, hidden_z in specs:
        observations = generate_cell(truth, noise, fraction, hidden_z, replicates, data_rng)
        by_arm = {arm: [] for arm in ARMS}
        for obs in observations:
            for arm in ARMS:
                post = infer(obs, arm)
                record = {
                    **obs, "inference_arm": arm,
                    "posterior_median": post["median"],
                    "median_error": post["median"] - truth,
                    "lower_68": post["lower_68"], "upper_68": post["upper_68"],
                    "lower_95": post["lower_95"], "upper_95": post["upper_95"],
                    "width_68": post["upper_68"] - post["lower_68"],
                    "width_95": post["upper_95"] - post["lower_95"],
                    "covered_68": int(post["lower_68"] <= truth <= post["upper_68"]),
                    "covered_95": int(post["lower_95"] <= truth <= post["upper_95"]),
                }
                by_arm[arm].append(record)
                raw.append(record)
        if fraction == 0.0:
            base = by_arm["base_weights"]
            modelled = by_arm["hidden_host_modelled"]
            for i, (a, b) in enumerate(zip(base, modelled)):
                for field in ("posterior_median", "lower_68", "upper_68", "lower_95", "upper_95"):
                    if not np.isclose(a[field], b[field], rtol=0, atol=1e-10):
                        raise AssertionError(f"f=0 negative control failed at cell {specs}: row {i}, {field}")
        for arm in ARMS:
            values = _metrics(by_arm[arm])
            summary.append({
                "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
                "hidden_redshift": hidden_z, "inference_arm": arm, "replicates": replicates,
                "mean_error": float(values[0]), "rmse_error": float(values[1]),
                "coverage_68": float(values[2]), "coverage_95": float(values[3]),
                "mean_width_68": float(values[4]), "mean_width_95": float(values[5]),
            })
        boot_indices = bootstrap_rng.integers(0, replicates, size=(bootstrap_resamples, replicates))
        for left, right, label in COMPARISONS:
            left_rows, right_rows = by_arm[left], by_arm[right]
            point = _metrics(left_rows) - _metrics(right_rows)
            boot = _metrics(left_rows, boot_indices) - _metrics(right_rows, boot_indices)
            for j, metric in enumerate(METRICS):
                contrasts.append({
                    "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
                    "hidden_redshift": hidden_z, "contrast": label,
                    "left_arm": left, "right_arm": right, "metric": metric,
                    "estimate": float(point[j]),
                    "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                    "ci95_high": float(np.quantile(boot[:, j], 0.975)),
                    "replicates": replicates, "bootstrap_resamples": bootstrap_resamples,
                    "bootstrap_seed": seed + 1,
                })
    return raw, summary, contrasts


def write_csv(path, rows):
    if not rows:
        raise ValueError("Refusing to write empty CSV.")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--bootstrap-resamples", type=int, default=1000)
    parser.add_argument("--output-dir", default="results/host-misspecification-sensitivity-surface")
    args = parser.parse_args()
    raw, summary, contrasts = run_experiment(args.replicates, args.seed, args.bootstrap_resamples)
    root = Path(args.output_dir)
    root.mkdir(parents=True, exist_ok=True)
    write_csv(root / "per_run.csv", raw)
    write_csv(root / "summary.csv", summary)
    write_csv(root / "paired_contrasts.csv", contrasts)
    print(f"Completed {len(raw)} fits from {len(raw) // 4} shared observations, "
          f"{len(summary)} arm-cell summaries, and {len(contrasts)} paired contrast rows.")
    for row in summary:
        print(row)


if __name__ == "__main__":
    main()
