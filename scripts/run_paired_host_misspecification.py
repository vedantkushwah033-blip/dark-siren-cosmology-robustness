"""Paired host-misspecification experiment on identical synthetic observations."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary

H0_TRUE = 70.0
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
Z_OBS = np.array([0.02, 0.04])
Z_HIDDEN = 0.03
POPULATION_Z = np.array([0.02, 0.04, 0.03])
POPULATION_WEIGHTS = np.array([0.35, 0.35, 0.30])
SIGMA_MPC = 0.10 * C_KM_S * 0.03 / H0_TRUE
GRID = np.linspace(25.0, 125.0, 4001)
ARMS = ("correct_weights", "reversed_weights", "hidden_host_ignored", "hidden_host_modelled")
CONTRASTS = (
    ("reversed_weights", "correct_weights", "reversed_minus_correct"),
    ("hidden_host_ignored", "hidden_host_modelled", "omitted_minus_modelled_hidden"),
)
METRICS = ("mean_error", "rmse_error", "coverage_68", "coverage_95", "mean_width_68", "mean_width_95")


def generate_observations(replicates: int, seed: int):
    if replicates < 200:
        raise ValueError("Use at least 200 paired observations.")
    rng = np.random.default_rng(seed)
    out = []
    for replicate in range(replicates):
        host_idx = int(rng.choice(len(POPULATION_Z), p=POPULATION_WEIGHTS))
        z_true = float(POPULATION_Z[host_idx])
        distance_true = C_KM_S * z_true / H0_TRUE
        observed = float(rng.normal(distance_true, SIGMA_MPC))
        while observed <= 0:
            observed = float(rng.normal(distance_true, SIGMA_MPC))
        out.append({
            "replicate": replicate,
            "true_H0": H0_TRUE,
            "true_host_index": host_idx,
            "true_host_redshift": z_true,
            "observed_distance_mpc": observed,
            "sigma_mpc": SIGMA_MPC,
        })
    return out


def infer_arm(observation, arm):
    d, sigma = observation["observed_distance_mpc"], observation["sigma_mpc"]
    kwargs = {"prior_min": PRIOR_MIN, "prior_max": PRIOR_MAX, "positive_truncated": True}
    if arm == "correct_weights":
        density = posterior_h0(GRID, Z_OBS, d, sigma, weights=[0.7, 0.3], **kwargs)
    elif arm == "reversed_weights":
        density = posterior_h0(GRID, Z_OBS, d, sigma, weights=[0.3, 0.7], **kwargs)
    elif arm == "hidden_host_ignored":
        density = posterior_h0(GRID, Z_OBS, d, sigma, weights=[0.5, 0.5], **kwargs)
    elif arm == "hidden_host_modelled":
        density = posterior_h0(GRID, Z_OBS, d, sigma, weights=[0.5, 0.5],
                               missing_weight=0.30, missing_redshifts=[Z_HIDDEN], **kwargs)
    else:
        raise ValueError(f"Unknown inference arm: {arm}")
    return posterior_summary(GRID, density)


def run_experiment(replicates=1000, seed=20261027, bootstrap_resamples=4000):
    if bootstrap_resamples < 100:
        raise ValueError("Use at least 100 paired bootstrap resamples.")
    observations = generate_observations(replicates, seed)
    raw = []
    for observation in observations:
        for arm in ARMS:
            posterior = infer_arm(observation, arm)
            raw.append({
                **observation,
                "inference_arm": arm,
                "posterior_median": posterior["median"],
                "median_error": posterior["median"] - H0_TRUE,
                "lower_68": posterior["lower_68"],
                "upper_68": posterior["upper_68"],
                "lower_95": posterior["lower_95"],
                "upper_95": posterior["upper_95"],
                "width_68": posterior["upper_68"] - posterior["lower_68"],
                "width_95": posterior["upper_95"] - posterior["lower_95"],
                "covered_68": int(posterior["lower_68"] <= H0_TRUE <= posterior["upper_68"]),
                "covered_95": int(posterior["lower_95"] <= H0_TRUE <= posterior["upper_95"]),
            })
    summary = []
    by_arm = {arm: [r for r in raw if r["inference_arm"] == arm] for arm in ARMS}
    for arm in ARMS:
        rows = by_arm[arm]
        errors = np.asarray([r["median_error"] for r in rows])
        summary.append({
            "inference_arm": arm, "replicates": replicates, "true_H0": H0_TRUE,
            "mean_error": float(np.mean(errors)),
            "rmse_error": float(np.sqrt(np.mean(errors**2))),
            "coverage_68": float(np.mean([r["covered_68"] for r in rows])),
            "coverage_95": float(np.mean([r["covered_95"] for r in rows])),
            "mean_width_68": float(np.mean([r["width_68"] for r in rows])),
            "mean_width_95": float(np.mean([r["width_95"] for r in rows])),
        })
    rng = np.random.default_rng(seed + 1)
    indices = rng.integers(0, replicates, size=(bootstrap_resamples, replicates))
    contrasts = []
    for left, right, label in CONTRASTS:
        a, b = by_arm[left], by_arm[right]
        a_errors = np.asarray([r["median_error"] for r in a])
        b_errors = np.asarray([r["median_error"] for r in b])
        a68 = np.asarray([r["covered_68"] for r in a], dtype=float)
        b68 = np.asarray([r["covered_68"] for r in b], dtype=float)
        a95 = np.asarray([r["covered_95"] for r in a], dtype=float)
        b95 = np.asarray([r["covered_95"] for r in b], dtype=float)
        a_w68 = np.asarray([r["width_68"] for r in a])
        b_w68 = np.asarray([r["width_68"] for r in b])
        a_w95 = np.asarray([r["width_95"] for r in a])
        b_w95 = np.asarray([r["width_95"] for r in b])
        point = np.asarray([
            np.mean(a_errors - b_errors),
            np.sqrt(np.mean(a_errors**2)) - np.sqrt(np.mean(b_errors**2)),
            np.mean(a68 - b68), np.mean(a95 - b95),
            np.mean(a_w68 - b_w68), np.mean(a_w95 - b_w95),
        ])
        boot = np.stack([
            np.mean(a_errors[indices] - b_errors[indices], axis=1),
            np.sqrt(np.mean(a_errors[indices]**2, axis=1)) - np.sqrt(np.mean(b_errors[indices]**2, axis=1)),
            np.mean(a68[indices] - b68[indices], axis=1),
            np.mean(a95[indices] - b95[indices], axis=1),
            np.mean(a_w68[indices] - b_w68[indices], axis=1),
            np.mean(a_w95[indices] - b_w95[indices], axis=1),
        ], axis=1)
        for j, metric in enumerate(METRICS):
            contrasts.append({
                "contrast": label, "left_arm": left, "right_arm": right, "metric": metric,
                "estimate": float(point[j]),
                "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                "ci95_high": float(np.quantile(boot[:, j], 0.975)),
                "replicates": replicates, "bootstrap_resamples": bootstrap_resamples,
                "bootstrap_seed": seed + 1,
            })
    return raw, summary, contrasts


def write_csv(path: Path, rows):
    if not rows:
        raise ValueError("Refusing to write empty CSV.")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261027)
    parser.add_argument("--bootstrap-resamples", type=int, default=4000)
    parser.add_argument("--output-dir", default="results/paired-host-misspecification")
    args = parser.parse_args()
    raw, summary, contrasts = run_experiment(args.replicates, args.seed, args.bootstrap_resamples)
    root = Path(args.output_dir)
    root.mkdir(parents=True, exist_ok=True)
    write_csv(root / "per_run.csv", raw)
    write_csv(root / "summary.csv", summary)
    write_csv(root / "paired_contrasts.csv", contrasts)
    print(f"PASS: {len(raw)} posterior fits on {args.replicates} shared observations; "
          f"{len(contrasts)} paired contrast rows.")
    for row in summary:
        print(row)
    for row in contrasts:
        print(row)


if __name__ == "__main__":
    main()
