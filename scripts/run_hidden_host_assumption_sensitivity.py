"""Paired sensitivity to misspecified hidden-host fraction and redshift (synthetic only)."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary

TRUTHS = (55.0, 70.0, 85.0)
NOISE_LEVELS = (0.05, 0.20)
TRUE_FRACTIONS = (0.0, 0.15, 0.30, 0.45)
TRUE_HIDDEN_ZS = (0.025, 0.035)
ASSUMED_FRACTIONS = (0.15, 0.30, 0.45)
ASSUMED_HIDDEN_ZS = (0.025, 0.030, 0.035)
BASE_ZS = np.array([0.02, 0.04])
GRID = np.linspace(25.0, 125.0, 4001)
PRIOR = (50.0, 90.0)
METRICS = ("mean_error", "rmse_error", "coverage_68", "coverage_95", "mean_width_68", "mean_width_95")


def assumption_label(fraction, redshift):
    return f"assumed_f{fraction:.3f}_z{redshift:.3f}"


def arm_names(assumed_fractions=ASSUMED_FRACTIONS, assumed_redshifts=ASSUMED_HIDDEN_ZS):
    return ("base_omitted", "oracle_hidden") + tuple(
        assumption_label(f, z) for f in assumed_fractions for z in assumed_redshifts
    )


def cell_specs(truths=TRUTHS, noises=NOISE_LEVELS, fractions=TRUE_FRACTIONS, hidden_zs=TRUE_HIDDEN_ZS):
    return [(float(t), float(n), float(f), float(z))
            for t in truths for n in noises for f in fractions for z in hidden_zs]


def generate_cell(truth, noise, fraction, hidden_z, replicates, rng):
    sigma = noise * C_KM_S * 0.03 / truth
    z_values = np.array([0.02, 0.04, hidden_z])
    probabilities = np.array([(1.0 - fraction) * 0.7, (1.0 - fraction) * 0.3, fraction])
    records = []
    for rep in range(replicates):
        host = int(rng.choice(3, p=probabilities))
        z = float(z_values[host])
        distance = float(rng.normal(C_KM_S * z / truth, sigma))
        while distance <= 0:
            distance = float(rng.normal(C_KM_S * z / truth, sigma))
        records.append({
            "true_H0": truth, "noise_fraction": noise, "true_hidden_fraction": fraction,
            "true_hidden_redshift": hidden_z, "replicate": rep, "true_host_index": host,
            "true_host_redshift": z, "observed_distance_mpc": distance, "sigma_mpc": sigma,
        })
    return records


def infer(obs, arm, assumed_fractions=ASSUMED_FRACTIONS, assumed_redshifts=ASSUMED_HIDDEN_ZS):
    truth_fraction, truth_z = obs["true_hidden_fraction"], obs["true_hidden_redshift"]
    if arm == "base_omitted":
        missing_fraction, missing_z = 0.0, None
    elif arm == "oracle_hidden":
        missing_fraction, missing_z = truth_fraction, truth_z
    else:
        lookup = {assumption_label(f, z): (f, z)
                  for f in assumed_fractions for z in assumed_redshifts}
        if arm not in lookup:
            raise ValueError(f"Unknown inference arm {arm}")
        missing_fraction, missing_z = lookup[arm]
    common = {"prior_min": PRIOR[0], "prior_max": PRIOR[1], "positive_truncated": True}
    kwargs = {"weights": [0.7, 0.3], **common}
    if missing_fraction > 0:
        kwargs.update(missing_weight=missing_fraction, missing_redshifts=[missing_z])
    density = posterior_h0(GRID, BASE_ZS, obs["observed_distance_mpc"], obs["sigma_mpc"], **kwargs)
    return posterior_summary(GRID, density)


def metric_vector(records, indices=None):
    error = np.asarray([float(r["median_error"]) for r in records])
    c68 = np.asarray([int(r["covered_68"]) for r in records], dtype=float)
    c95 = np.asarray([int(r["covered_95"]) for r in records], dtype=float)
    w68 = np.asarray([float(r["width_68"]) for r in records])
    w95 = np.asarray([float(r["width_95"]) for r in records])
    if indices is not None:
        error, c68, c95, w68, w95 = (a[indices] for a in (error, c68, c95, w68, w95))
    return np.stack((np.mean(error, axis=-1), np.sqrt(np.mean(error**2, axis=-1)),
                     np.mean(c68, axis=-1), np.mean(c95, axis=-1),
                     np.mean(w68, axis=-1), np.mean(w95, axis=-1)), axis=-1)


def run_experiment(replicates=200, seed=20261032, bootstrap_resamples=1000,
                   truths=TRUTHS, noises=NOISE_LEVELS, fractions=TRUE_FRACTIONS,
                   hidden_zs=TRUE_HIDDEN_ZS, assumed_fractions=ASSUMED_FRACTIONS,
                   assumed_redshifts=ASSUMED_HIDDEN_ZS):
    if replicates < 200:
        raise ValueError("Use at least 200 observations per true cell.")
    if bootstrap_resamples < 100:
        raise ValueError("Use at least 100 paired bootstrap resamples.")
    assumed_fractions = tuple(float(f) for f in assumed_fractions)
    assumed_redshifts = tuple(float(z) for z in assumed_redshifts)
    if not assumed_fractions or any(not 0.0 < f < 1.0 for f in assumed_fractions):
        raise ValueError("Assumed hidden fractions must lie strictly between zero and one.")
    if not assumed_redshifts or any(z <= 0 or not np.isfinite(z) for z in assumed_redshifts):
        raise ValueError("Assumed hidden redshifts must be finite and positive.")
    arms = arm_names(assumed_fractions, assumed_redshifts)
    contrasts = (("base_omitted", "oracle_hidden", "base_omitted_minus_oracle"),) + tuple(
        (arm, "oracle_hidden", arm + "_minus_oracle") for arm in arms[2:]
    )
    specs = cell_specs(truths, noises, fractions, hidden_zs)
    if not specs:
        raise ValueError("Factorial design must contain at least one cell.")
    data_rng, boot_rng = np.random.default_rng(seed), np.random.default_rng(seed + 1)
    raw, summaries, contrast_rows = [], [], []
    for truth, noise, fraction, hidden_z in specs:
        observations = generate_cell(truth, noise, fraction, hidden_z, replicates, data_rng)
        by_arm = {arm: [] for arm in arms}
        for obs in observations:
            for arm in arms:
                post = infer(obs, arm, assumed_fractions, assumed_redshifts)
                row = {
                    **obs, "inference_arm": arm,
                    "posterior_median": post["median"], "median_error": post["median"] - truth,
                    "lower_68": post["lower_68"], "upper_68": post["upper_68"],
                    "lower_95": post["lower_95"], "upper_95": post["upper_95"],
                    "width_68": post["upper_68"] - post["lower_68"],
                    "width_95": post["upper_95"] - post["lower_95"],
                    "covered_68": int(post["lower_68"] <= truth <= post["upper_68"]),
                    "covered_95": int(post["lower_95"] <= truth <= post["upper_95"]),
                }
                by_arm[arm].append(row)
                raw.append(row)
        if fraction == 0.0:
            for base, oracle in zip(by_arm["base_omitted"], by_arm["oracle_hidden"]):
                for field in ("posterior_median", "lower_68", "upper_68", "lower_95", "upper_95"):
                    if not np.isclose(base[field], oracle[field], rtol=0, atol=1e-10):
                        raise AssertionError(f"Zero-fraction oracle negative control failed: {field}")
        for arm in arms:
            values = metric_vector(by_arm[arm])
            summaries.append({
                "true_H0": truth, "noise_fraction": noise, "true_hidden_fraction": fraction,
                "true_hidden_redshift": hidden_z, "inference_arm": arm, "replicates": replicates,
                **dict(zip(METRICS, [float(v) for v in values])),
            })
        indices = boot_rng.integers(0, replicates, size=(bootstrap_resamples, replicates))
        for left, right, label in contrasts:
            point = metric_vector(by_arm[left]) - metric_vector(by_arm[right])
            boot = metric_vector(by_arm[left], indices) - metric_vector(by_arm[right], indices)
            for j, metric in enumerate(METRICS):
                contrast_rows.append({
                    "true_H0": truth, "noise_fraction": noise, "true_hidden_fraction": fraction,
                    "true_hidden_redshift": hidden_z, "contrast": label, "left_arm": left,
                    "right_arm": right, "metric": metric, "estimate": float(point[j]),
                    "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                    "ci95_high": float(np.quantile(boot[:, j], 0.975)),
                    "replicates": replicates, "bootstrap_resamples": bootstrap_resamples,
                    "bootstrap_seed": seed + 1,
                })
    return raw, summaries, contrast_rows


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
    parser.add_argument("--seed", type=int, default=20261032)
    parser.add_argument("--bootstrap-resamples", type=int, default=1000)
    parser.add_argument("--output-dir", default="results/hidden-host-assumption-sensitivity")
    args = parser.parse_args()
    raw, summaries, contrasts = run_experiment(args.replicates, args.seed, args.bootstrap_resamples)
    root = Path(args.output_dir)
    root.mkdir(parents=True, exist_ok=True)
    for name, rows in (("per_run.csv", raw), ("summary.csv", summaries), ("paired_contrasts.csv", contrasts)):
        write_csv(root / name, rows)
    print(f"Completed {len(raw)} fits from {len(raw)//len(arm_names())} shared observations; "
          f"{len(summaries)} summaries and {len(contrasts)} paired contrast rows.")


if __name__ == "__main__":
    main()
