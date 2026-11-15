"""Independent reconstruction audit for hidden-host fraction/redshift assumption sensitivity."""
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
TRUE_FRACTIONS = (0.0, 0.15, 0.30, 0.45)
TRUE_ZS = (0.025, 0.035)
ASSUMED_FRACTIONS = (0.15, 0.30, 0.45)
ASSUMED_ZS = (0.025, 0.030, 0.035)
Z_LISTED = np.array([0.02, 0.04])
GRID = np.linspace(25.0, 125.0, 4001)
PRIOR = (50.0, 90.0)
METRICS = ("mean_error", "rmse_error", "coverage_68", "coverage_95", "mean_width_68", "mean_width_95")


def assumption_label(fraction, redshift):
    return f"assumed_f{fraction:.3f}_z{redshift:.3f}"


def arm_names(assumed_fractions=ASSUMED_FRACTIONS, assumed_redshifts=ASSUMED_ZS):
    return ("base_omitted", "oracle_hidden") + tuple(
        assumption_label(f, z) for f in assumed_fractions for z in assumed_redshifts
    )


def contrast_specs(arms):
    return (("base_omitted", "oracle_hidden", "base_omitted_minus_oracle"),) + tuple(
        (arm, "oracle_hidden", arm + "_minus_oracle") for arm in arms[2:]
    )


def cell_specs(truths=TRUTHS, noises=NOISES, fractions=TRUE_FRACTIONS, hidden_zs=TRUE_ZS):
    return [(float(t), float(n), float(f), float(z))
            for t in truths for n in noises for f in fractions for z in hidden_zs]


def regenerate_cell(truth, noise, fraction, hidden_z, replicates, rng):
    sigma = noise * C * 0.03 / truth
    z_values = np.array([0.02, 0.04, hidden_z])
    probabilities = np.array([(1.0 - fraction) * 0.7, (1.0 - fraction) * 0.3, fraction])
    output = []
    for rep in range(replicates):
        host = int(rng.choice(3, p=probabilities))
        z = float(z_values[host])
        distance = float(rng.normal(C * z / truth, sigma))
        while distance <= 0:
            distance = float(rng.normal(C * z / truth, sigma))
        output.append((rep, host, z, distance, sigma))
    return output


def posterior_equation(distance, sigma, truth, true_fraction, true_z, arm,
                       assumed_fractions, assumed_redshifts):
    if arm == "base_omitted":
        missing, hidden_z = 0.0, None
    elif arm == "oracle_hidden":
        missing, hidden_z = true_fraction, true_z
    else:
        lookup = {assumption_label(f, z): (f, z)
                  for f in assumed_fractions for z in assumed_redshifts}
        if arm not in lookup:
            raise ValueError(f"Unexpected arm {arm}")
        missing, hidden_z = lookup[arm]
    mu = C * Z_LISTED[None, :] / GRID[:, None]
    listed = norm.pdf(distance, loc=mu, scale=sigma) / norm.cdf(mu / sigma)
    likelihood = (1.0 - missing) * np.sum(np.array([0.7, 0.3])[None, :] * listed, axis=1)
    if missing > 0:
        mu_hidden = C * hidden_z / GRID
        hidden = norm.pdf(distance, loc=mu_hidden, scale=sigma) / norm.cdf(mu_hidden / sigma)
        likelihood += missing * hidden
    mask = (GRID >= PRIOR[0]) & (GRID <= PRIOR[1])
    x, y = GRID[mask], likelihood[mask]
    area = float(trapezoid(y, x))
    if not math.isfinite(area) or area <= 0:
        raise ValueError("Invalid posterior normalization")
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
        e, c68, c95, w68, w95 = (x[indices] for x in (e, c68, c95, w68, w95))
    return np.stack((np.mean(e, axis=-1), np.sqrt(np.mean(e**2, axis=-1)),
                     np.mean(c68, axis=-1), np.mean(c95, axis=-1),
                     np.mean(w68, axis=-1), np.mean(w95, axis=-1)), axis=-1)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Missing or empty CSV: {path}")
    return rows


def _compare(expected_rows, actual_rows, keys, fields, label, atol=2e-8):
    actual_map = {tuple(str(row[k]) for k in keys): row for row in actual_rows}
    if len(actual_map) != len(actual_rows):
        raise ValueError(f"Duplicate keys in {label}")
    for expected in expected_rows:
        key = tuple(str(expected[k]) for k in keys)
        if key not in actual_map:
            raise ValueError(f"Missing {label} row {key}")
        for field in fields:
            if not math.isclose(float(actual_map[key][field]), float(expected[field]), rel_tol=2e-10, abs_tol=atol):
                raise ValueError(f"{label} reconstruction mismatch {key}, field {field}")
    if len(actual_map) != len(expected_rows):
        raise ValueError(f"Unexpected extra rows in {label}")


def audit(root: Path, replicates=200, seed=20261032, resamples=1000,
          truths=TRUTHS, noises=NOISES, fractions=TRUE_FRACTIONS, true_zs=TRUE_ZS,
          assumed_fractions=ASSUMED_FRACTIONS, assumed_zs=ASSUMED_ZS, atol=2e-8):
    if replicates < 200 or resamples < 100:
        raise ValueError("Audit design below protocol minimum")
    arms = arm_names(assumed_fractions, assumed_zs)
    comparisons = contrast_specs(arms)
    specs = cell_specs(truths, noises, fractions, true_zs)
    raw, actual_summary, actual_contrasts = (
        read_csv(root / name) for name in ("per_run.csv", "summary.csv", "paired_contrasts.csv")
    )
    if len(raw) != len(specs) * replicates * len(arms):
        raise ValueError("Raw row count does not match factorial design")
    if len(actual_summary) != len(specs) * len(arms):
        raise ValueError("Summary row count does not match factorial design")
    if len(actual_contrasts) != len(specs) * len(comparisons) * len(METRICS):
        raise ValueError("Contrast row count does not match factorial design")
    data_rng, boot_rng = np.random.default_rng(seed), np.random.default_rng(seed + 1)
    expected_summary, expected_contrasts = [], []
    raw_i = 0
    for truth, noise, fraction, true_z in specs:
        generated = regenerate_cell(truth, noise, fraction, true_z, replicates, data_rng)
        by_arm = {arm: [] for arm in arms}
        for rep, host, z, distance, sigma in generated:
            for arm in arms:
                row = raw[raw_i]
                raw_i += 1
                if row["inference_arm"] != arm or int(row["replicate"]) != rep:
                    raise ValueError(f"Raw order/key mismatch at row {raw_i-1}")
                expected_obs = {
                    "true_H0": truth, "noise_fraction": noise, "true_hidden_fraction": fraction,
                    "true_hidden_redshift": true_z, "true_host_index": host, "true_host_redshift": z,
                    "observed_distance_mpc": distance, "sigma_mpc": sigma,
                }
                for field, value in expected_obs.items():
                    if not math.isclose(float(row[field]), float(value), rel_tol=2e-13, abs_tol=2e-13):
                        raise ValueError(f"Seeded observation mismatch at row {raw_i-1}, field {field}")
                post = posterior_equation(distance, sigma, truth, fraction, true_z, arm,
                                         assumed_fractions, assumed_zs)
                for field, value in post.items():
                    if field not in row or not math.isclose(float(row[field]), float(value), rel_tol=2e-10, abs_tol=atol):
                        raise ValueError(f"Posterior reconstruction mismatch at row {raw_i-1}, field {field}")
                by_arm[arm].append(post)
        if fraction == 0:
            for base, oracle in zip(by_arm["base_omitted"], by_arm["oracle_hidden"]):
                for field in ("posterior_median", "lower_68", "upper_68", "lower_95", "upper_95"):
                    if not math.isclose(base[field], oracle[field], rel_tol=0, abs_tol=1e-10):
                        raise ValueError("Zero-fraction negative control failed")
        for arm in arms:
            values = metric_vector(by_arm[arm])
            expected_summary.append({
                "true_H0": truth, "noise_fraction": noise, "true_hidden_fraction": fraction,
                "true_hidden_redshift": true_z, "inference_arm": arm, "replicates": replicates,
                **dict(zip(METRICS, [float(v) for v in values])),
            })
        indices = boot_rng.integers(0, replicates, size=(resamples, replicates))
        for left, right, label in comparisons:
            point = metric_vector(by_arm[left]) - metric_vector(by_arm[right])
            boot = metric_vector(by_arm[left], indices) - metric_vector(by_arm[right], indices)
            for j, metric in enumerate(METRICS):
                expected_contrasts.append({
                    "true_H0": truth, "noise_fraction": noise, "true_hidden_fraction": fraction,
                    "true_hidden_redshift": true_z, "contrast": label, "left_arm": left,
                    "right_arm": right, "metric": metric, "estimate": float(point[j]),
                    "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                    "ci95_high": float(np.quantile(boot[:, j], 0.975)),
                    "replicates": replicates, "bootstrap_resamples": resamples, "bootstrap_seed": seed + 1,
                })
    _compare(expected_summary, actual_summary,
             ("true_H0", "noise_fraction", "true_hidden_fraction", "true_hidden_redshift", "inference_arm"),
             ("replicates", *METRICS), "summary", atol)
    _compare(expected_contrasts, actual_contrasts,
             ("true_H0", "noise_fraction", "true_hidden_fraction", "true_hidden_redshift", "contrast", "metric"),
             ("estimate", "ci95_low", "ci95_high", "replicates", "bootstrap_resamples", "bootstrap_seed"),
             "paired contrast", atol)
    return {"audit": "PASS", "factorial_cells": len(specs), "shared_observations": len(specs) * replicates,
            "posterior_reconstructions": len(raw), "summary_rows_reconstructed": len(expected_summary),
            "paired_contrast_rows_reconstructed": len(expected_contrasts), "bootstrap_resamples": resamples}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="results/hidden-host-assumption-sensitivity")
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20261032)
    parser.add_argument("--bootstrap-resamples", type=int, default=1000)
    args = parser.parse_args()
    print(audit(Path(args.input_dir), args.replicates, args.seed, args.bootstrap_resamples))


if __name__ == "__main__":
    main()
