"""Paired prior-support × hidden-host sensitivity factorial (synthetic model only)."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary

TRUTHS = (55.0, 70.0, 85.0)
NOISE_LEVELS = (0.05, 0.20)
HIDDEN_FRACTIONS = (0.0, 0.30, 0.45)
HIDDEN_REDSHIFTS = (0.025, 0.035)
PRIORS = ((50.0, 90.0), (40.0, 100.0), (30.0, 120.0))
BASE_REDSHIFTS = np.array([0.02, 0.04])
ARMS = ("base_weights", "reversed_weights", "equal_weights", "hidden_host_modelled")
ARM_COMPARISONS = (
    ("reversed_weights", "base_weights", "reversed_minus_base"),
    ("base_weights", "hidden_host_modelled", "omitted_minus_modelled_hidden"),
)
PRIOR_COMPARISONS = ((40.0, 100.0, "40_100_minus_50_90"), (30.0, 120.0, "30_120_minus_50_90"))
METRICS = ("mean_error", "rmse_error", "coverage_68", "coverage_95", "mean_width_68", "mean_width_95")
GRID = np.linspace(25.0, 125.0, 4001)


def cell_specs(truths=TRUTHS, noise_levels=NOISE_LEVELS,
               hidden_fractions=HIDDEN_FRACTIONS, hidden_redshifts=HIDDEN_REDSHIFTS):
    return [(float(t), float(n), float(f), float(z))
            for t in truths for n in noise_levels for f in hidden_fractions for z in hidden_redshifts]


def generate_cell(truth, noise, fraction, hidden_z, replicates, rng):
    sigma = noise * C_KM_S * 0.03 / truth
    z_values = np.array([0.02, 0.04, hidden_z])
    probabilities = np.array([(1.0 - fraction) * 0.7, (1.0 - fraction) * 0.3, fraction])
    observations = []
    for replicate in range(replicates):
        host = int(rng.choice(3, p=probabilities))
        z = float(z_values[host])
        distance = float(rng.normal(C_KM_S * z / truth, sigma))
        while distance <= 0:
            distance = float(rng.normal(C_KM_S * z / truth, sigma))
        observations.append({
            "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
            "hidden_redshift": hidden_z, "replicate": replicate, "true_host_index": host,
            "true_host_redshift": z, "observed_distance_mpc": distance, "sigma_mpc": sigma,
        })
    return observations


def infer(obs, arm, prior):
    distance, sigma = obs["observed_distance_mpc"], obs["sigma_mpc"]
    fraction, hidden_z = obs["hidden_fraction"], obs["hidden_redshift"]
    lo, hi = prior
    common = {"prior_min": lo, "prior_max": hi, "positive_truncated": True}
    if arm == "base_weights":
        density = posterior_h0(GRID, BASE_REDSHIFTS, distance, sigma, weights=[0.7, 0.3], **common)
    elif arm == "reversed_weights":
        density = posterior_h0(GRID, BASE_REDSHIFTS, distance, sigma, weights=[0.3, 0.7], **common)
    elif arm == "equal_weights":
        density = posterior_h0(GRID, BASE_REDSHIFTS, distance, sigma, weights=[0.5, 0.5], **common)
    elif arm == "hidden_host_modelled":
        density = posterior_h0(GRID, BASE_REDSHIFTS, distance, sigma, weights=[0.7, 0.3],
                               missing_weight=fraction, missing_redshifts=[hidden_z], **common)
    else:
        raise ValueError(f"Unknown inference arm: {arm}")
    return posterior_summary(GRID, density)


def metric_vector(records, indices=None):
    errors = np.asarray([float(r["median_error"]) for r in records])
    c68 = np.asarray([int(r["covered_68"]) for r in records], dtype=float)
    c95 = np.asarray([int(r["covered_95"]) for r in records], dtype=float)
    w68 = np.asarray([float(r["width_68"]) for r in records])
    w95 = np.asarray([float(r["width_95"]) for r in records])
    if indices is not None:
        errors, c68, c95, w68, w95 = (x[indices] for x in (errors, c68, c95, w68, w95))
    return np.stack((np.mean(errors, axis=-1), np.sqrt(np.mean(errors**2, axis=-1)),
                     np.mean(c68, axis=-1), np.mean(c95, axis=-1),
                     np.mean(w68, axis=-1), np.mean(w95, axis=-1)), axis=-1)


def run_experiment(replicates=200, seed=20261030, bootstrap_resamples=1000,
                   truths=TRUTHS, noise_levels=NOISE_LEVELS,
                   hidden_fractions=HIDDEN_FRACTIONS, hidden_redshifts=HIDDEN_REDSHIFTS,
                   priors=PRIORS):
    if replicates < 200:
        raise ValueError("Use at least 200 observations per base cell.")
    if bootstrap_resamples < 100:
        raise ValueError("Use at least 100 paired bootstrap resamples.")
    priors = tuple((float(lo), float(hi)) for lo, hi in priors)
    if len(priors) < 2 or len(set(priors)) != len(priors):
        raise ValueError("Need at least two unique prior intervals.")
    for lo, hi in priors:
        if not (25.0 <= lo < hi <= 125.0):
            raise ValueError("Every prior must lie inside the H0 grid.")
    specs = cell_specs(truths, noise_levels, hidden_fractions, hidden_redshifts)
    if not specs:
        raise ValueError("Factorial design must contain at least one base cell.")
    data_rng = np.random.default_rng(seed)
    bootstrap_rng = np.random.default_rng(seed + 1)
    raw, summaries, arm_contrasts, prior_contrasts = [], [], [], []

    for truth, noise, fraction, hidden_z in specs:
        observations = generate_cell(truth, noise, fraction, hidden_z, replicates, data_rng)
        grouped = {(prior, arm): [] for prior in priors for arm in ARMS}
        for obs in observations:
            for prior in priors:
                for arm in ARMS:
                    post = infer(obs, arm, prior)
                    row = {
                        **obs, "prior_min": prior[0], "prior_max": prior[1], "inference_arm": arm,
                        "posterior_median": post["median"], "median_error": post["median"] - truth,
                        "lower_68": post["lower_68"], "upper_68": post["upper_68"],
                        "lower_95": post["lower_95"], "upper_95": post["upper_95"],
                        "width_68": post["upper_68"] - post["lower_68"],
                        "width_95": post["upper_95"] - post["lower_95"],
                        "covered_68": int(post["lower_68"] <= truth <= post["upper_68"]),
                        "covered_95": int(post["lower_95"] <= truth <= post["upper_95"]),
                    }
                    grouped[(prior, arm)].append(row)
                    raw.append(row)
        if fraction == 0.0:
            for rep in range(replicates):
                base = grouped[(priors[0], "base_weights")][rep]
                modelled = grouped[(priors[0], "hidden_host_modelled")][rep]
                for field in ("posterior_median", "lower_68", "upper_68", "lower_95", "upper_95"):
                    if not np.isclose(base[field], modelled[field], rtol=0, atol=1e-10):
                        raise AssertionError(f"Zero-hidden-fraction negative control failed: {field}")
        for prior in priors:
            for arm in ARMS:
                rows = grouped[(prior, arm)]
                value = metric_vector(rows)
                summaries.append({
                    "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
                    "hidden_redshift": hidden_z, "prior_min": prior[0], "prior_max": prior[1],
                    "inference_arm": arm, "replicates": replicates,
                    **dict(zip(METRICS, [float(v) for v in value])),
                })
        indices = bootstrap_rng.integers(0, replicates, size=(bootstrap_resamples, replicates))
        for prior in priors:
            for left, right, label in ARM_COMPARISONS:
                left_rows, right_rows = grouped[(prior, left)], grouped[(prior, right)]
                point = metric_vector(left_rows) - metric_vector(right_rows)
                boot = metric_vector(left_rows, indices) - metric_vector(right_rows, indices)
                for j, metric in enumerate(METRICS):
                    arm_contrasts.append({
                        "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
                        "hidden_redshift": hidden_z, "prior_min": prior[0], "prior_max": prior[1],
                        "contrast": label, "left_arm": left, "right_arm": right, "metric": metric,
                        "estimate": float(point[j]), "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                        "ci95_high": float(np.quantile(boot[:, j], 0.975)), "replicates": replicates,
                        "bootstrap_resamples": bootstrap_resamples, "bootstrap_seed": seed + 1,
                    })
        base_prior = priors[0]
        for arm in ARMS:
            base_rows = grouped[(base_prior, arm)]
            for lo, hi, label in PRIOR_COMPARISONS:
                wide_prior = next((p for p in priors if p == (lo, hi)), None)
                if wide_prior is None:
                    continue
                wide_rows = grouped[(wide_prior, arm)]
                point = metric_vector(wide_rows) - metric_vector(base_rows)
                boot = metric_vector(wide_rows, indices) - metric_vector(base_rows, indices)
                for j, metric in enumerate(METRICS):
                    prior_contrasts.append({
                        "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
                        "hidden_redshift": hidden_z, "inference_arm": arm,
                        "comparison_prior_min": lo, "comparison_prior_max": hi,
                        "reference_prior_min": base_prior[0], "reference_prior_max": base_prior[1],
                        "contrast": label, "metric": metric, "estimate": float(point[j]),
                        "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                        "ci95_high": float(np.quantile(boot[:, j], 0.975)), "replicates": replicates,
                        "bootstrap_resamples": bootstrap_resamples, "bootstrap_seed": seed + 1,
                    })
    return raw, summaries, arm_contrasts, prior_contrasts


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--bootstrap-resamples", type=int, default=1000)
    parser.add_argument("--output-dir", default="results/host-prior-support-factorial")
    args = parser.parse_args()
    raw, summaries, arm_contrasts, prior_contrasts = run_experiment(
        args.replicates, args.seed, args.bootstrap_resamples)
    root = Path(args.output_dir)
    root.mkdir(parents=True, exist_ok=True)
    for name, rows in (("per_run.csv", raw), ("summary.csv", summaries),
                       ("arm_contrasts.csv", arm_contrasts), ("prior_contrasts.csv", prior_contrasts)):
        write_csv(root / name, rows)
    print(f"Completed {len(raw)} posterior fits from {len(raw) // (len(PRIORS) * len(ARMS))} shared observations; "
          f"{len(summaries)} summary rows, {len(arm_contrasts)} arm contrasts, "
          f"{len(prior_contrasts)} prior contrasts.")


if __name__ == "__main__":
    main()
