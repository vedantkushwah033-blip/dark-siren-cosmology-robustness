"""Independent reconstruction audit for prior-support × hidden-host factorial outputs."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.stats import norm

C = 299792.458
TRUTHS = (55.0, 70.0, 85.0)
NOISES = (0.05, 0.20)
FRACTIONS = (0.0, 0.30, 0.45)
HIDDEN_ZS = (0.025, 0.035)
PRIORS = ((50.0, 90.0), (40.0, 100.0), (30.0, 120.0))
Z_LISTED = np.array([0.02, 0.04])
ARMS = ("base_weights", "reversed_weights", "equal_weights", "hidden_host_modelled")
ARM_COMPARISONS = (
    ("reversed_weights", "base_weights", "reversed_minus_base"),
    ("base_weights", "hidden_host_modelled", "omitted_minus_modelled_hidden"),
)
PRIOR_COMPARISONS = ((40.0, 100.0, "40_100_minus_50_90"), (30.0, 120.0, "30_120_minus_50_90"))
METRICS = ("mean_error", "rmse_error", "coverage_68", "coverage_95", "mean_width_68", "mean_width_95")
GRID = np.linspace(25.0, 125.0, 4001)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Missing or empty CSV: {path}")
    return rows


def cell_specs(truths=TRUTHS, noises=NOISES, fractions=FRACTIONS, hidden_zs=HIDDEN_ZS):
    return [(float(t), float(n), float(f), float(z))
            for t in truths for n in noises for f in fractions for z in hidden_zs]


def regenerate_cell(truth, noise, fraction, hidden_z, replicates, rng):
    sigma = noise * C * 0.03 / truth
    z_values = np.array([0.02, 0.04, hidden_z])
    probabilities = np.array([(1.0 - fraction) * 0.7, (1.0 - fraction) * 0.3, fraction])
    result = []
    for rep in range(replicates):
        host = int(rng.choice(3, p=probabilities))
        z = float(z_values[host])
        distance = float(rng.normal(C * z / truth, sigma))
        while distance <= 0:
            distance = float(rng.normal(C * z / truth, sigma))
        result.append((rep, host, z, distance, sigma))
    return result


def posterior_equation(distance, sigma, truth, fraction, hidden_z, arm, prior):
    if arm == "base_weights":
        weights, missing = np.array([0.7, 0.3]), 0.0
    elif arm == "reversed_weights":
        weights, missing = np.array([0.3, 0.7]), 0.0
    elif arm == "equal_weights":
        weights, missing = np.array([0.5, 0.5]), 0.0
    elif arm == "hidden_host_modelled":
        weights, missing = np.array([0.7, 0.3]), fraction
    else:
        raise ValueError(f"Unknown inference arm: {arm}")
    mu = C * Z_LISTED[None, :] / GRID[:, None]
    listed_density = norm.pdf(distance, loc=mu, scale=sigma) / norm.cdf(mu / sigma)
    likelihood = (1.0 - missing) * np.sum(weights[None, :] * listed_density, axis=1)
    if missing > 0:
        mu_hidden = C * hidden_z / GRID
        hidden_density = norm.pdf(distance, loc=mu_hidden, scale=sigma) / norm.cdf(mu_hidden / sigma)
        likelihood += missing * hidden_density
    lo, hi = prior
    mask = (GRID >= lo) & (GRID <= hi)
    x, y = GRID[mask], likelihood[mask]
    area = float(trapezoid(y, x))
    if not math.isfinite(area) or area <= 0:
        raise ValueError("Invalid independently reconstructed posterior normalization")
    y = y / area
    cdf = np.concatenate(([0.0], np.cumsum((y[1:] + y[:-1]) * np.diff(x) / 2.0)))
    cdf /= cdf[-1]
    q = {p: float(np.interp(p, cdf, x)) for p in (0.025, 0.16, 0.5, 0.84, 0.975)}
    return {
        "posterior_median": q[0.5], "median_error": q[0.5] - truth,
        "lower_68": q[0.16], "upper_68": q[0.84], "lower_95": q[0.025], "upper_95": q[0.975],
        "width_68": q[0.84] - q[0.16], "width_95": q[0.975] - q[0.025],
        "covered_68": int(q[0.16] <= truth <= q[0.84]),
        "covered_95": int(q[0.025] <= truth <= q[0.975]),
    }


def metric_vector(records, indices=None):
    e = np.asarray([float(r["median_error"]) for r in records])
    c68 = np.asarray([int(r["covered_68"]) for r in records], dtype=float)
    c95 = np.asarray([int(r["covered_95"]) for r in records], dtype=float)
    w68 = np.asarray([float(r["width_68"]) for r in records])
    w95 = np.asarray([float(r["width_95"]) for r in records])
    if indices is not None:
        e, c68, c95, w68, w95 = (a[indices] for a in (e, c68, c95, w68, w95))
    return np.stack((np.mean(e, axis=-1), np.sqrt(np.mean(e**2, axis=-1)),
                     np.mean(c68, axis=-1), np.mean(c95, axis=-1),
                     np.mean(w68, axis=-1), np.mean(w95, axis=-1)), axis=-1)


def _compare(expected_rows, actual_rows, keys, fields, label, atol=2e-8):
    actual_map = {tuple(str(r[k]) for k in keys): r for r in actual_rows}
    if len(actual_map) != len(actual_rows):
        raise ValueError(f"Duplicate keys in {label}")
    for expected in expected_rows:
        key = tuple(str(expected[k]) for k in keys)
        if key not in actual_map:
            raise ValueError(f"Missing {label} row {key}")
        actual = actual_map[key]
        for field in fields:
            if not math.isclose(float(actual[field]), float(expected[field]), rel_tol=2e-10, abs_tol=atol):
                raise ValueError(f"{label} reconstruction mismatch {key}, field {field}")
    if len(actual_map) != len(expected_rows):
        raise ValueError(f"Unexpected extra rows in {label}")


def audit(root: Path, replicates=200, seed=20261030, resamples=1000,
          truths=TRUTHS, noises=NOISES, fractions=FRACTIONS, hidden_zs=HIDDEN_ZS,
          priors=PRIORS, atol=2e-8):
    if replicates < 200 or resamples < 100:
        raise ValueError("Audit design below protocol minimum")
    priors = tuple((float(lo), float(hi)) for lo, hi in priors)
    specs = cell_specs(truths, noises, fractions, hidden_zs)
    raw = read_csv(root / "per_run.csv")
    actual_summary = read_csv(root / "summary.csv")
    actual_arm_contrasts = read_csv(root / "arm_contrasts.csv")
    actual_prior_contrasts = read_csv(root / "prior_contrasts.csv")
    expected_raw = len(specs) * replicates * len(priors) * len(ARMS)
    expected_summary_n = len(specs) * len(priors) * len(ARMS)
    expected_arm_n = len(specs) * len(priors) * len(ARM_COMPARISONS) * len(METRICS)
    prior_pairs = [(lo, hi, label) for lo, hi, label in PRIOR_COMPARISONS if (lo, hi) in priors]
    expected_prior_n = len(specs) * len(ARMS) * len(prior_pairs) * len(METRICS)
    if len(raw) != expected_raw or len(actual_summary) != expected_summary_n:
        raise ValueError("Raw or summary row count does not match prior-support factorial design")
    if len(actual_arm_contrasts) != expected_arm_n or len(actual_prior_contrasts) != expected_prior_n:
        raise ValueError("Arm/prior contrast row count does not match factorial design")

    data_rng, boot_rng = np.random.default_rng(seed), np.random.default_rng(seed + 1)
    expected_summary, expected_arm_contrasts, expected_prior_contrasts = [], [], []
    raw_i = 0
    for truth, noise, fraction, hidden_z in specs:
        observations = regenerate_cell(truth, noise, fraction, hidden_z, replicates, data_rng)
        grouped = {(prior, arm): [] for prior in priors for arm in ARMS}
        for rep, host, z, distance, sigma in observations:
            for prior in priors:
                for arm in ARMS:
                    row = raw[raw_i]
                    raw_i += 1
                    if (row["inference_arm"] != arm or int(row["replicate"]) != rep
                            or not math.isclose(float(row["prior_min"]), prior[0], abs_tol=1e-12)
                            or not math.isclose(float(row["prior_max"]), prior[1], abs_tol=1e-12)):
                        raise ValueError(f"Raw ordering/key mismatch at row {raw_i-1}")
                    observation_fields = {
                        "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
                        "hidden_redshift": hidden_z, "true_host_index": host, "true_host_redshift": z,
                        "observed_distance_mpc": distance, "sigma_mpc": sigma,
                    }
                    for field, value in observation_fields.items():
                        if not math.isclose(float(row[field]), float(value), rel_tol=2e-13, abs_tol=2e-13):
                            raise ValueError(f"Seeded observation mismatch at row {raw_i-1}, field {field}")
                    posterior = posterior_equation(distance, sigma, truth, fraction, hidden_z, arm, prior)
                    for field, value in posterior.items():
                        if field not in row or not math.isclose(float(row[field]), float(value), rel_tol=2e-10, abs_tol=atol):
                            raise ValueError(f"Posterior reconstruction mismatch at row {raw_i-1}, field {field}")
                    grouped[(prior, arm)].append(posterior)
        if fraction == 0:
            for prior in priors:
                for a, b in zip(grouped[(prior, "base_weights")], grouped[(prior, "hidden_host_modelled")]):
                    for field in ("posterior_median", "lower_68", "upper_68", "lower_95", "upper_95"):
                        if not math.isclose(a[field], b[field], rel_tol=0, abs_tol=1e-10):
                            raise ValueError("Zero-hidden-fraction negative control failed")
        for prior in priors:
            for arm in ARMS:
                vals = metric_vector(grouped[(prior, arm)])
                expected_summary.append({
                    "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
                    "hidden_redshift": hidden_z, "prior_min": prior[0], "prior_max": prior[1],
                    "inference_arm": arm, "replicates": replicates,
                    **dict(zip(METRICS, [float(v) for v in vals])),
                })
        indices = boot_rng.integers(0, replicates, size=(resamples, replicates))
        for prior in priors:
            for left, right, label in ARM_COMPARISONS:
                left_rows, right_rows = grouped[(prior, left)], grouped[(prior, right)]
                point = metric_vector(left_rows) - metric_vector(right_rows)
                boot = metric_vector(left_rows, indices) - metric_vector(right_rows, indices)
                for j, metric in enumerate(METRICS):
                    expected_arm_contrasts.append({
                        "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
                        "hidden_redshift": hidden_z, "prior_min": prior[0], "prior_max": prior[1],
                        "contrast": label, "left_arm": left, "right_arm": right, "metric": metric,
                        "estimate": float(point[j]), "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                        "ci95_high": float(np.quantile(boot[:, j], 0.975)), "replicates": replicates,
                        "bootstrap_resamples": resamples, "bootstrap_seed": seed + 1,
                    })
        reference = priors[0]
        for arm in ARMS:
            for lo, hi, label in prior_pairs:
                wide_rows, ref_rows = grouped[((lo, hi), arm)], grouped[(reference, arm)]
                point = metric_vector(wide_rows) - metric_vector(ref_rows)
                boot = metric_vector(wide_rows, indices) - metric_vector(ref_rows, indices)
                for j, metric in enumerate(METRICS):
                    expected_prior_contrasts.append({
                        "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
                        "hidden_redshift": hidden_z, "inference_arm": arm,
                        "comparison_prior_min": lo, "comparison_prior_max": hi,
                        "reference_prior_min": reference[0], "reference_prior_max": reference[1],
                        "contrast": label, "metric": metric, "estimate": float(point[j]),
                        "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                        "ci95_high": float(np.quantile(boot[:, j], 0.975)), "replicates": replicates,
                        "bootstrap_resamples": resamples, "bootstrap_seed": seed + 1,
                    })
    _compare(expected_summary, actual_summary,
             ("true_H0", "noise_fraction", "hidden_fraction", "hidden_redshift", "prior_min", "prior_max", "inference_arm"),
             ("replicates", *METRICS), "summary", atol)
    _compare(expected_arm_contrasts, actual_arm_contrasts,
             ("true_H0", "noise_fraction", "hidden_fraction", "hidden_redshift", "prior_min", "prior_max", "contrast", "metric"),
             ("estimate", "ci95_low", "ci95_high", "replicates", "bootstrap_resamples", "bootstrap_seed"),
             "arm contrast", atol)
    _compare(expected_prior_contrasts, actual_prior_contrasts,
             ("true_H0", "noise_fraction", "hidden_fraction", "hidden_redshift", "inference_arm", "contrast", "metric"),
             ("estimate", "ci95_low", "ci95_high", "replicates", "bootstrap_resamples", "bootstrap_seed"),
             "prior contrast", atol)
    return {
        "audit": "PASS", "base_cells": len(specs), "prior_arm_cells": len(specs) * len(priors),
        "shared_observations": len(specs) * replicates, "posterior_reconstructions": expected_raw,
        "summary_rows_reconstructed": len(expected_summary),
        "arm_contrast_rows_reconstructed": len(expected_arm_contrasts),
        "prior_contrast_rows_reconstructed": len(expected_prior_contrasts), "bootstrap_resamples": resamples,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="results/host-prior-support-factorial")
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--bootstrap-resamples", type=int, default=1000)
    args = parser.parse_args()
    print(audit(Path(args.input_dir), args.replicates, args.seed, args.bootstrap_resamples))


if __name__ == "__main__":
    main()
