"""Independent raw-data, posterior, and paired-bootstrap audit."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.stats import norm

C = 299792.458
H0 = 70.0
Z_OBS = np.array([0.02, 0.04])
Z_HIDDEN = 0.03
Z_POP = np.array([0.02, 0.04, 0.03])
P_POP = np.array([0.35, 0.35, 0.30])
SIGMA = 0.10 * C * 0.03 / H0
GRID = np.linspace(25.0, 125.0, 4001)
PRIOR = (50.0, 90.0)
ARMS = ("correct_weights", "reversed_weights", "hidden_host_ignored", "hidden_host_modelled")
COMPARISONS = (
    ("reversed_weights", "correct_weights", "reversed_minus_correct"),
    ("hidden_host_ignored", "hidden_host_modelled", "omitted_minus_modelled_hidden"),
)
METRICS = ("mean_error", "rmse_error", "coverage_68", "coverage_95", "mean_width_68", "mean_width_95")


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def regenerate(replicates, seed):
    rng = np.random.default_rng(seed)
    out = []
    for replicate in range(replicates):
        index = int(rng.choice(3, p=P_POP))
        z = float(Z_POP[index])
        observed = float(rng.normal(C * z / H0, SIGMA))
        while observed <= 0:
            observed = float(rng.normal(C * z / H0, SIGMA))
        out.append((replicate, index, z, observed))
    return out


def posterior_reconstruction(observed, arm):
    if arm == "correct_weights":
        weights, missing = np.array([0.7, 0.3]), 0.0
    elif arm == "reversed_weights":
        weights, missing = np.array([0.3, 0.7]), 0.0
    elif arm == "hidden_host_ignored":
        weights, missing = np.array([0.5, 0.5]), 0.0
    elif arm == "hidden_host_modelled":
        weights, missing = np.array([0.5, 0.5]), 0.30
    else:
        raise ValueError(f"Unexpected inference arm: {arm}")
    means = C * Z_OBS[None, :] / GRID[:, None]
    components = norm.pdf(observed, loc=means, scale=SIGMA) / norm.cdf(means / SIGMA)
    likelihood = (1.0 - missing) * np.sum(weights[None, :] * components, axis=1)
    if missing:
        hidden_mean = C * Z_HIDDEN / GRID
        hidden_component = norm.pdf(observed, loc=hidden_mean, scale=SIGMA) / norm.cdf(hidden_mean / SIGMA)
        likelihood = likelihood + missing * hidden_component
    keep = (GRID >= PRIOR[0]) & (GRID <= PRIOR[1])
    x = GRID[keep]
    density = likelihood[keep]
    area = trapezoid(density, x)
    if not np.isfinite(area) or area <= 0:
        raise ValueError("Invalid reconstructed posterior normalization")
    density = density / area
    cdf = np.concatenate(([0.0], np.cumsum((density[1:] + density[:-1]) * np.diff(x) / 2.0)))
    cdf = cdf / cdf[-1]
    q = {p: float(np.interp(p, cdf, x)) for p in (0.025, 0.16, 0.5, 0.84, 0.975)}
    return {
        "posterior_median": q[0.5], "median_error": q[0.5] - H0,
        "lower_68": q[0.16], "upper_68": q[0.84],
        "lower_95": q[0.025], "upper_95": q[0.975],
        "width_68": q[0.84] - q[0.16], "width_95": q[0.975] - q[0.025],
        "covered_68": int(q[0.16] <= H0 <= q[0.84]),
        "covered_95": int(q[0.025] <= H0 <= q[0.975]),
    }


def reconstructed_contrasts(rows_by_arm, replicates, seed, resamples):
    rng = np.random.default_rng(seed + 1)
    indices = rng.integers(0, replicates, size=(resamples, replicates))
    results = []
    for left, right, label in COMPARISONS:
        a, b = rows_by_arm[left], rows_by_arm[right]
        ae = np.array([float(r["median_error"]) for r in a])
        be = np.array([float(r["median_error"]) for r in b])
        ac68 = np.array([int(r["covered_68"]) for r in a], dtype=float)
        bc68 = np.array([int(r["covered_68"]) for r in b], dtype=float)
        ac95 = np.array([int(r["covered_95"]) for r in a], dtype=float)
        bc95 = np.array([int(r["covered_95"]) for r in b], dtype=float)
        aw68 = np.array([float(r["width_68"]) for r in a])
        bw68 = np.array([float(r["width_68"]) for r in b])
        aw95 = np.array([float(r["width_95"]) for r in a])
        bw95 = np.array([float(r["width_95"]) for r in b])
        point = np.array([
            np.mean(ae - be), np.sqrt(np.mean(ae**2)) - np.sqrt(np.mean(be**2)),
            np.mean(ac68 - bc68), np.mean(ac95 - bc95),
            np.mean(aw68 - bw68), np.mean(aw95 - bw95),
        ])
        boot = np.column_stack((
            np.mean(ae[indices] - be[indices], axis=1),
            np.sqrt(np.mean(ae[indices]**2, axis=1)) - np.sqrt(np.mean(be[indices]**2, axis=1)),
            np.mean(ac68[indices] - bc68[indices], axis=1),
            np.mean(ac95[indices] - bc95[indices], axis=1),
            np.mean(aw68[indices] - bw68[indices], axis=1),
            np.mean(aw95[indices] - bw95[indices], axis=1),
        ))
        for j, metric in enumerate(METRICS):
            results.append({
                "contrast": label, "left_arm": left, "right_arm": right, "metric": metric,
                "estimate": float(point[j]), "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                "ci95_high": float(np.quantile(boot[:, j], 0.975)),
                "replicates": replicates, "bootstrap_resamples": resamples, "bootstrap_seed": seed + 1,
            })
    return results


def audit(root: Path, replicates=1000, seed=20261027, resamples=4000, atol=2e-8):
    if replicates < 200 or resamples < 100:
        raise ValueError("Audit design below protocol minimum")
    raw = read_csv(root / "per_run.csv")
    summaries = read_csv(root / "summary.csv")
    actual_contrasts = read_csv(root / "paired_contrasts.csv")
    expected_obs = regenerate(replicates, seed)
    if len(raw) != replicates * len(ARMS) or len(summaries) != len(ARMS) or len(actual_contrasts) != len(COMPARISONS) * len(METRICS):
        raise ValueError("Unexpected raw or summary row counts")
    grouped = {arm: [] for arm in ARMS}
    for i, (row, obs) in enumerate(zip(raw, [(o, arm) for o in expected_obs for arm in ARMS])):
        (replicate, host_idx, z, observed), arm = obs
        if row["inference_arm"] != arm or int(row["replicate"]) != replicate:
            raise ValueError(f"Pair order/replicate mismatch at raw row {i}")
        expected = {"true_H0": H0, "true_host_index": host_idx, "true_host_redshift": z,
                    "observed_distance_mpc": observed, "sigma_mpc": SIGMA}
        for key, value in expected.items():
            if not math.isclose(float(row[key]), float(value), rel_tol=2e-13, abs_tol=2e-13):
                raise ValueError(f"Observation provenance mismatch at row {i}, {key}")
        reconstructed = posterior_reconstruction(observed, arm)
        for key, value in reconstructed.items():
            if key not in row or not math.isclose(float(row[key]), float(value), rel_tol=2e-10, abs_tol=atol):
                raise ValueError(f"Posterior/derived-field mismatch at row {i}, {key}")
        grouped[arm].append(reconstructed)

    for i, arm in enumerate(ARMS):
        rows = grouped[arm]
        errors = np.array([r["median_error"] for r in rows])
        expected = {
            "inference_arm": arm, "replicates": replicates, "true_H0": H0,
            "mean_error": float(np.mean(errors)), "rmse_error": float(np.sqrt(np.mean(errors**2))),
            "coverage_68": float(np.mean([r["covered_68"] for r in rows])),
            "coverage_95": float(np.mean([r["covered_95"] for r in rows])),
            "mean_width_68": float(np.mean([r["width_68"] for r in rows])),
            "mean_width_95": float(np.mean([r["width_95"] for r in rows])),
        }
        for key, value in expected.items():
            if key == "inference_arm":
                if summaries[i][key] != value:
                    raise ValueError(f"Summary arm/order mismatch for {arm}")
            elif not math.isclose(float(summaries[i][key]), float(value), rel_tol=2e-10, abs_tol=atol):
                raise ValueError(f"Summary reconstruction mismatch {arm}.{key}")

    expected_contrasts = reconstructed_contrasts(grouped, replicates, seed, resamples)
    actual_index = {(r["contrast"], r["metric"]): r for r in actual_contrasts}
    if len(actual_index) != len(actual_contrasts):
        raise ValueError("Duplicate paired-contrast keys")
    for expected in expected_contrasts:
        key = (expected["contrast"], expected["metric"])
        if key not in actual_index:
            raise ValueError(f"Missing paired contrast {key}")
        actual = actual_index[key]
        for field in ("left_arm", "right_arm"):
            if actual[field] != expected[field]:
                raise ValueError(f"Contrast arm mismatch for {key}")
        for field in ("estimate", "ci95_low", "ci95_high", "replicates", "bootstrap_resamples", "bootstrap_seed"):
            if not math.isclose(float(actual[field]), float(expected[field]), rel_tol=2e-10, abs_tol=atol):
                raise ValueError(f"Paired-bootstrap reconstruction mismatch {key}.{field}")
        if float(actual["ci95_low"]) > float(actual["ci95_high"]):
            raise ValueError(f"Reversed confidence interval for {key}")
    return {"audit": "PASS", "observations": replicates, "posterior_reconstructions": len(raw),
            "paired_contrasts_reconstructed": len(expected_contrasts), "bootstrap_resamples": resamples}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="results/paired-host-misspecification")
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261027)
    parser.add_argument("--bootstrap-resamples", type=int, default=4000)
    args = parser.parse_args()
    print(audit(Path(args.input_dir), args.replicates, args.seed, args.bootstrap_resamples))


if __name__ == "__main__":
    main()
