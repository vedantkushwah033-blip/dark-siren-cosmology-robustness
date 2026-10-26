"""Independent reconstruction audit for host-misspecification stress-test CSVs.

This auditor intentionally does not import the production runner or inference
implementation. It regenerates the declared RNG stream and independently
reconstructs posterior summaries and aggregate metrics for the toy model.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.stats import norm

C = 299792.458
H0_TRUE = 70.0
PRIOR = (50.0, 90.0)
Z_OBS = np.array([0.02, 0.04])
Z_HIDDEN = 0.03
SIGMA = 0.10 * C * 0.03 / H0_TRUE
GRID = np.linspace(25.0, 125.0, 4001)
SCENARIOS = ("correct_weights", "reversed_weights", "hidden_host_ignored", "hidden_host_modelled")
SUMMARY_FIELDS = ("mean_median_error", "median_median_error", "rmse_median",
                  "coverage_68", "coverage_95", "mean_width_68", "mean_width_95")


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def regenerate(replicates: int, seed: int):
    """Reproduce the published data-generating process without importing it."""
    rng = np.random.default_rng(seed)
    designs = {
        "correct_weights": (np.array([0.02, 0.04]), np.array([0.7, 0.3])),
        "reversed_weights": (np.array([0.02, 0.04]), np.array([0.7, 0.3])),
        "hidden_host_ignored": (np.array([0.02, 0.04, 0.03]), np.array([0.35, 0.35, 0.30])),
        "hidden_host_modelled": (np.array([0.02, 0.04, 0.03]), np.array([0.35, 0.35, 0.30])),
    }
    out = []
    for scenario, (z_values, weights) in designs.items():
        for replicate in range(replicates):
            idx = int(rng.choice(len(z_values), p=weights))
            z = float(z_values[idx])
            observed = float(rng.normal(C * z / H0_TRUE, SIGMA))
            while observed <= 0:
                observed = float(rng.normal(C * z / H0_TRUE, SIGMA))
            out.append((scenario, replicate, idx, z, observed))
    return out


def independent_summary(scenario: str, observed: float):
    """Rebuild the posterior from the model equation, not production helpers."""
    z = Z_OBS
    if scenario == "correct_weights":
        weights, missing_weight = np.array([0.7, 0.3]), 0.0
    elif scenario == "reversed_weights":
        weights, missing_weight = np.array([0.3, 0.7]), 0.0
    elif scenario == "hidden_host_ignored":
        weights, missing_weight = np.array([0.5, 0.5]), 0.0
    elif scenario == "hidden_host_modelled":
        weights, missing_weight = np.array([0.5, 0.5]), 0.30
    else:
        raise ValueError(f"Unknown scenario {scenario}")
    means = C * z[None, :] / GRID[:, None]
    likelihood = (1.0 - missing_weight) * np.sum(
        weights[None, :] * norm.pdf(observed, loc=means, scale=SIGMA), axis=1
    )
    if missing_weight:
        hidden_mean = C * Z_HIDDEN / GRID
        likelihood += missing_weight * norm.pdf(observed, loc=hidden_mean, scale=SIGMA)
    mask = (GRID >= PRIOR[0]) & (GRID <= PRIOR[1])
    x, density = GRID[mask], likelihood[mask]
    density = density / trapezoid(density, x)
    cdf = np.concatenate(([0.0], np.cumsum((density[1:] + density[:-1]) * np.diff(x) / 2)))
    cdf /= cdf[-1]
    quantiles = {q: float(np.interp(q, cdf, x)) for q in (0.025, 0.16, 0.5, 0.84, 0.975)}
    return {
        "posterior_median": quantiles[0.5],
        "median_error": quantiles[0.5] - H0_TRUE,
        "lower_68": quantiles[0.16], "upper_68": quantiles[0.84],
        "lower_95": quantiles[0.025], "upper_95": quantiles[0.975],
    }


def audit_directory(root: Path, replicates: int = 1000, seed: int = 20261027):
    raw = read_csv(root / "per_run.csv")
    actual = read_csv(root / "summary.csv")
    if replicates < 200:
        raise ValueError("At least 200 replicates are required by the protocol.")
    if len(raw) != replicates * len(SCENARIOS) or len(actual) != len(SCENARIOS):
        raise ValueError("Raw/summary row count does not match the declared design.")
    expected_data = regenerate(replicates, seed)
    reconstructed = {name: [] for name in SCENARIOS}
    for i, (row, expected) in enumerate(zip(raw, expected_data)):
        scenario, rep, host_idx, z, observed = expected
        if row["scenario"] != scenario or int(row["replicate"]) != rep:
            raise ValueError(f"Scenario/order/replicate mismatch at raw row {i}")
        checks = {
            "true_H0": H0_TRUE, "true_host_index": host_idx,
            "true_host_redshift": z, "observed_distance_mpc": observed, "sigma_mpc": SIGMA,
        }
        for key, value in checks.items():
            got = float(row[key])
            if not math.isclose(got, float(value), rel_tol=2e-13, abs_tol=2e-13):
                raise ValueError(f"RNG/data provenance mismatch at row {i}, field {key}")
        post = independent_summary(scenario, observed)
        for key, value in post.items():
            if not math.isclose(float(row[key]), value, rel_tol=2e-10, abs_tol=2e-8):
                raise ValueError(f"Posterior reconstruction mismatch at row {i}, field {key}")
        derived = {
            **post,
            "covered_68": int(float(row["lower_68"]) <= H0_TRUE <= float(row["upper_68"])),
            "covered_95": int(float(row["lower_95"]) <= H0_TRUE <= float(row["upper_95"])),
            "width_68": post["upper_68"] - post["lower_68"],
            "width_95": post["upper_95"] - post["lower_95"],
        }
        for key in ("covered_68", "covered_95", "width_68", "width_95"):
            if not math.isclose(float(row[key]), float(derived[key]), rel_tol=2e-10, abs_tol=2e-8):
                raise ValueError(f"Per-run derived-field mismatch at row {i}, field {key}")
        reconstructed[scenario].append(derived)

    for summary_row, scenario in zip(actual, SCENARIOS):
        if summary_row["scenario"] != scenario or int(summary_row["replicates"]) != replicates:
            raise ValueError(f"Summary design mismatch for {scenario}")
        rows = reconstructed[scenario]
        errors = np.array([r["median_error"] for r in rows])
        expected_metrics = {
            "mean_median_error": float(errors.mean()),
            "median_median_error": float(np.median(errors)),
            "rmse_median": float(np.sqrt(np.mean(errors**2))),
            "coverage_68": float(np.mean([r["covered_68"] for r in rows])),
            "coverage_95": float(np.mean([r["covered_95"] for r in rows])),
            "mean_width_68": float(np.mean([r["width_68"] for r in rows])),
            "mean_width_95": float(np.mean([r["width_95"] for r in rows])),
        }
        for field, expected_value in expected_metrics.items():
            if not math.isclose(float(summary_row[field]), expected_value, rel_tol=2e-10, abs_tol=2e-8):
                raise ValueError(f"Summary reconstruction mismatch for {scenario}.{field}")
    return {"status": "PASS", "raw_rows": len(raw), "summary_rows": len(actual),
            "posterior_reconstructions": len(raw), "seed": seed, "replicates_per_scenario": replicates}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="results/host-misspecification-audit")
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261027)
    args = parser.parse_args()
    result = audit_directory(Path(args.input_dir), args.replicates, args.seed)
    print("PASS: independent RNG provenance, per-fit posterior reconstruction, and summary reconstruction")
    print(result)


if __name__ == "__main__":
    main()
