"""Matched host-structure misspecification control; synthetic diagnostic only.

Generate detected distances from a two-host redshift mixture, then analyse each
same observation with the generating two-host model and a collapsed single-host
model at the mixture-weighted mean redshift. Compare with and without the toy
selection normalization. This isolates analysis host-structure mismatch better
than comparing separately generated one-host and two-host datasets.
"""
from __future__ import annotations
import argparse, csv, json, platform
from pathlib import Path
import numpy as np
import scipy
from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary
from dark_siren_robustness.selection import posterior_h0_selected

HOST_REDSHIFTS = np.array([0.02, 0.04], dtype=float)
HOST_WEIGHTS = np.array([0.70, 0.30], dtype=float)
COLLAPSED_REDSHIFT = float(np.dot(HOST_REDSHIFTS, HOST_WEIGHTS))
SIGMA_MPC = 20.0
DETECTION_MIN_MPC, DETECTION_MAX_MPC = 0.0, 120.0
ASSUMED_MAX_MPC = 120.0
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
TRUTHS = (52.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 88.0)
METHODS = ("two_host_uncorrected", "two_host_corrected",
           "collapsed_uncorrected", "collapsed_corrected")


def draw_detected(truth, rng, max_attempts=100000):
    if not np.isfinite(truth) or truth <= 0:
        raise ValueError("truth must be finite and positive")
    for _ in range(max_attempts):
        host = int(rng.choice(len(HOST_REDSHIFTS), p=HOST_WEIGHTS))
        distance = float(rng.normal(C_KM_S * HOST_REDSHIFTS[host] / truth, SIGMA_MPC))
        if DETECTION_MIN_MPC < distance < DETECTION_MAX_MPC:
            return host, distance
    raise RuntimeError("Detection rejection sampler exceeded limit")


def wilson(successes, n, z=1.959963984540054):
    if n <= 0 or not 0 <= successes <= n:
        raise ValueError("invalid binomial counts")
    p = successes / n
    den = 1 + z*z/n
    center = (p + z*z/(2*n)) / den
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return (0.0 if successes == 0 else max(0.0, float(center-half)),
            1.0 if successes == n else min(1.0, float(center+half)))


def run(replicates=300, seed=20261031, grid_size=4001, bootstrap_draws=5000):
    if replicates < 100:
        raise ValueError("Use at least 100 replicates per truth")
    if grid_size < 1001:
        raise ValueError("grid_size must be at least 1001")
    if bootstrap_draws < 1000:
        raise ValueError("Use at least 1000 bootstrap draws")
    rng = np.random.default_rng(seed)
    grid = np.linspace(25.0, 125.0, grid_size)
    rows = []
    for truth in TRUTHS:
        for replicate in range(replicates):
            host, distance = draw_detected(truth, rng)
            posteriors = {
                "two_host_uncorrected": posterior_h0(
                    grid, HOST_REDSHIFTS, distance, SIGMA_MPC, weights=HOST_WEIGHTS,
                    prior_min=PRIOR_MIN, prior_max=PRIOR_MAX),
                "two_host_corrected": posterior_h0_selected(
                    grid, HOST_REDSHIFTS, distance, SIGMA_MPC, ASSUMED_MAX_MPC,
                    weights=HOST_WEIGHTS, prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
                    detection_distance_min=DETECTION_MIN_MPC),
                "collapsed_uncorrected": posterior_h0(
                    grid, np.array([COLLAPSED_REDSHIFT]), distance, SIGMA_MPC,
                    weights=np.array([1.0]), prior_min=PRIOR_MIN, prior_max=PRIOR_MAX),
                "collapsed_corrected": posterior_h0_selected(
                    grid, np.array([COLLAPSED_REDSHIFT]), distance, SIGMA_MPC,
                    ASSUMED_MAX_MPC, weights=np.array([1.0]), prior_min=PRIOR_MIN,
                    prior_max=PRIOR_MAX, detection_distance_min=DETECTION_MIN_MPC),
            }
            row = {"true_H0": truth, "replicate": replicate,
                   "generated_host_index": host, "observed_distance_mpc": distance}
            for method, posterior in posteriors.items():
                summary = posterior_summary(grid, posterior)
                for key in ("median", "lower_68", "upper_68", "lower_95", "upper_95"):
                    row[f"{method}_{key}"] = float(summary[key])
                error = float(summary["median"] - truth)
                row[f"{method}_error"] = error
                row[f"{method}_absolute_error"] = abs(error)
            rows.append(row)

    summaries = []
    bootstrap_rng = np.random.default_rng(seed + 1)
    for truth in TRUTHS:
        cell = [r for r in rows if r["true_H0"] == truth]
        for method in METHODS:
            errors = np.array([r[f"{method}_error"] for r in cell])
            hit68 = np.array([r[f"{method}_lower_68"] <= truth <= r[f"{method}_upper_68"] for r in cell])
            hit95 = np.array([r[f"{method}_lower_95"] <= truth <= r[f"{method}_upper_95"] for r in cell])
            lo68, hi68 = wilson(int(hit68.sum()), len(cell))
            lo95, hi95 = wilson(int(hit95.sum()), len(cell))
            summaries.append({
                "true_H0": truth, "method": method, "n": len(cell),
                "bias": float(errors.mean()), "mae": float(np.abs(errors).mean()),
                "rmse": float(np.sqrt(np.mean(errors**2))),
                "coverage_68": float(hit68.mean()), "coverage_68_wilson_low": lo68,
                "coverage_68_wilson_high": hi68, "coverage_95": float(hit95.mean()),
                "coverage_95_wilson_low": lo95, "coverage_95_wilson_high": hi95,
                "mean_width_68": float(np.mean([r[f"{method}_upper_68"]-r[f"{method}_lower_68"] for r in cell])),
                "mean_width_95": float(np.mean([r[f"{method}_upper_95"]-r[f"{method}_lower_95"] for r in cell])),
            })
        for correction in ("uncorrected", "corrected"):
            collapsed = f"collapsed_{correction}"
            mixture = f"two_host_{correction}"
            delta = np.array([r[f"{collapsed}_absolute_error"] - r[f"{mixture}_absolute_error"] for r in cell])
            boot = np.empty(bootstrap_draws)
            for start in range(0, bootstrap_draws, 256):
                stop = min(start + 256, bootstrap_draws)
                indices = bootstrap_rng.integers(0, len(delta), size=(stop-start, len(delta)))
                boot[start:stop] = delta[indices].mean(axis=1)
            summaries.append({
                "true_H0": truth,
                "method": f"paired_collapsed_minus_two_host_absolute_error_{correction}",
                "n": len(cell), "mean_paired_change": float(delta.mean()),
                "bootstrap_95_low": float(np.quantile(boot, 0.025)),
                "bootstrap_95_high": float(np.quantile(boot, 0.975)),
                "bootstrap_draws": bootstrap_draws, "bootstrap_seed": seed + 1,
            })
    manifest = {
        "experiment": "matched_host_structure_misspecification_control",
        "scope": "synthetic; two-host mixture generates data; not an observational measurement",
        "seed": seed, "replicates_per_truth": replicates, "grid_size": grid_size,
        "grid_bounds": [25.0, 125.0], "prior": [PRIOR_MIN, PRIOR_MAX],
        "truths": list(TRUTHS), "host_redshifts": HOST_REDSHIFTS.tolist(),
        "host_weights": HOST_WEIGHTS.tolist(), "collapsed_redshift_weighted_mean": COLLAPSED_REDSHIFT,
        "sigma_mpc": SIGMA_MPC, "true_detection_interval_mpc": [DETECTION_MIN_MPC, DETECTION_MAX_MPC],
        "assumed_detection_interval_mpc": [DETECTION_MIN_MPC, ASSUMED_MAX_MPC],
        "methods": list(METHODS), "bootstrap_draws": bootstrap_draws,
        "bootstrap_seed": seed + 1, "raw_rows": len(rows), "summary_rows": len(summaries),
        "python_version": platform.python_version(), "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
    }
    return rows, summaries, manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--grid-size", type=int, default=4001)
    parser.add_argument("--bootstrap-draws", type=int, default=5000)
    parser.add_argument("--output-dir", default="results/matched-host-structure-control")
    args = parser.parse_args()
    rows, summaries, manifest = run(args.replicates, args.seed, args.grid_size, args.bootstrap_draws)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for filename, data in (("per_run.csv", rows), ("summary.csv", summaries)):
        with (out / filename).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for row in data for k in row)))
            writer.writeheader()
            writer.writerows(data)
    (out / "experiment_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: {len(rows)} paired observations; {len(summaries)} summary rows.")
    for item in summaries:
        print(json.dumps(item, sort_keys=True))


if __name__ == "__main__":
    main()
