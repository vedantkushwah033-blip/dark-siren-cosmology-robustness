"""Run a seeded SBC and fixed-truth calibration experiment for the H0-dependent toy.

This is a synthetic distance-only model, not a detector-level gravitational-wave
analysis. The raw per-replicate truth and detected distances are preserved.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.stats import chi2

from dark_siren_robustness.selection_conditioned_h0_dependent_noise import (
    dataset_posterior_h0_selection_conditioned,
    simulate_detected_distance_events,
)

PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
FRACTIONAL_SIGMA, SIGMA_FLOOR_MPC = 0.12, 8.0
DETECTION_DISTANCE_MIN_MPC, N_EVENTS = 300.0, 3
N_POPULATION_NODES = 12
H0_POINTS = 161
MODE_VALUES = {"prior": None, "fixed_60": 60.0, "fixed_70": 70.0, "fixed_80": 80.0}


def wilson_interval(successes: int, trials: int) -> list[float]:
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("invalid Wilson interval counts")
    z = 1.959963984540054
    p = successes / trials
    z2 = z * z
    denominator = 1.0 + z2 / trials
    center = (p + z2 / (2.0 * trials)) / denominator
    half = z * np.sqrt(p * (1.0 - p) / trials + z2 / (4.0 * trials * trials)) / denominator
    return [float(max(0.0, center - half)), float(min(1.0, center + half))]


def _posterior_metrics(h0_grid, density, h0_true):
    area = float(trapezoid(density, h0_grid))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("posterior density has invalid normalization")
    density = np.asarray(density, dtype=float) / area
    cdf = np.r_[0.0, cumulative_trapezoid(density, h0_grid)]
    cdf /= cdf[-1]
    q05, q16, median, q84, q95 = np.interp(
        [0.05, 0.16, 0.50, 0.84, 0.95], cdf, h0_grid
    )
    return {
        "posterior_mean": float(trapezoid(h0_grid * density, h0_grid)),
        "posterior_median": float(median),
        "ci68_low": float(q16), "ci68_high": float(q84),
        "ci90_low": float(q05), "ci90_high": float(q95),
        "rank": float(np.interp(h0_true, h0_grid, cdf)),
        "posterior_integral": float(trapezoid(density, h0_grid)),
    }


def _summarise_mode(rows, mode):
    selected = [r for r in rows if r["mode"] == mode]
    n = len(selected)
    truths = np.asarray([r["h0_true"] for r in selected], dtype=float)
    means = np.asarray([r["posterior_mean"] for r in selected], dtype=float)
    errors = means - truths
    hit68 = int(sum(r["covered_68"] for r in selected))
    hit90 = int(sum(r["covered_90"] for r in selected))
    ranks = np.asarray([r["rank"] for r in selected], dtype=float)
    counts, _ = np.histogram(ranks, bins=10, range=(0.0, 1.0))
    expected = n / 10.0
    chi2_stat = float(np.sum((counts - expected) ** 2 / expected))
    result = {
        "mode": mode, "n_replicates": n,
        "bias": float(np.mean(errors)),
        "bias_standard_error": float(np.std(errors, ddof=1) / np.sqrt(n)),
        "rmse": float(np.sqrt(np.mean(errors ** 2))),
        "coverage68": hit68 / n, "coverage68_count": hit68,
        "coverage68_wilson95": wilson_interval(hit68, n),
        "coverage90": hit90 / n, "coverage90_count": hit90,
        "coverage90_wilson95": wilson_interval(hit90, n),
        "mean_rank": float(np.mean(ranks)),
        "rank_standard_error": float(np.std(ranks, ddof=1) / np.sqrt(n)),
        "rank_decile_counts": counts.tolist(),
        "rank_chi2_9df": chi2_stat,
        "rank_chi2_p_descriptive": float(chi2.sf(chi2_stat, 9)),
        "mean_interval68_width": float(np.mean([r["ci68_high"] - r["ci68_low"] for r in selected])),
        "mean_interval90_width": float(np.mean([r["ci90_high"] - r["ci90_low"] for r in selected])),
    }
    return result


def run_experiment(*, replicates=300, seed=20261030):
    if isinstance(replicates, bool) or not isinstance(replicates, int) or replicates < 30:
        raise ValueError("replicates must be an integer of at least 30")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("seed must be a non-negative integer")
    redshifts = np.linspace(0.025, 0.16, N_POPULATION_NODES)
    weights_raw = np.linspace(1.0, 2.5, N_POPULATION_NODES)
    h0_grid = np.linspace(PRIOR_MIN, PRIOR_MAX, H0_POINTS)
    rows = []
    for mode_index, (mode, fixed_h0) in enumerate(MODE_VALUES.items()):
        for replicate in range(replicates):
            replicate_seed = seed + mode_index * 100_000 + replicate
            rng = np.random.default_rng(replicate_seed)
            h0_true = float(rng.uniform(PRIOR_MIN, PRIOR_MAX)) if fixed_h0 is None else float(fixed_h0)
            events = simulate_detected_distance_events(
                rng=rng, n_events=N_EVENTS, population_redshifts=redshifts,
                population_weights=weights_raw, H0_true=h0_true,
                detection_distance_min=DETECTION_DISTANCE_MIN_MPC,
                fractional_sigma=FRACTIONAL_SIGMA, sigma_floor_mpc=SIGMA_FLOOR_MPC,
            )
            posterior = dataset_posterior_h0_selection_conditioned(
                h0_grid, events, prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
                population_redshifts=redshifts, population_weights=weights_raw,
                detection_distance_min=DETECTION_DISTANCE_MIN_MPC,
                fractional_sigma=FRACTIONAL_SIGMA, sigma_floor_mpc=SIGMA_FLOOR_MPC,
            )
            metrics = _posterior_metrics(h0_grid, posterior["posterior_density"], h0_true)
            rows.append({
                "mode": mode, "replicate": replicate, "seed": replicate_seed,
                "h0_true": h0_true,
                "distances_mpc_json": json.dumps([float(e["gw_distance"]) for e in events], separators=(",", ":")),
                **metrics,
                "covered_68": bool(metrics["ci68_low"] <= h0_true <= metrics["ci68_high"]),
                "covered_90": bool(metrics["ci90_low"] <= h0_true <= metrics["ci90_high"]),
            })
    summaries = [_summarise_mode(rows, mode) for mode in MODE_VALUES]
    configuration = {
        "prior": [PRIOR_MIN, PRIOR_MAX], "population_redshifts": redshifts.tolist(),
        "population_weights_raw": weights_raw.tolist(), "h0_points": H0_POINTS,
        "fractional_sigma": FRACTIONAL_SIGMA, "sigma_floor_mpc": SIGMA_FLOOR_MPC,
        "detection_distance_min_mpc": DETECTION_DISTANCE_MIN_MPC,
        "n_detected_events": N_EVENTS, "replicates_per_mode": replicates,
        "seed_start": seed, "seed_end_by_mode": {
            mode: seed + i * 100_000 + replicates - 1 for i, mode in enumerate(MODE_VALUES)
        },
        "modes": {mode: fixed for mode, fixed in MODE_VALUES.items()},
        "model_scope": "synthetic_discrete_host_distance_only_model_not_real_GW_inference",
    }
    summary = {
        "status": "synthetic_h0_dependent_noise_selection_sbc_not_real_h0_measurement",
        "git_commit": os.environ.get("GITHUB_SHA"),
        "configuration": configuration,
        "summaries": summaries,
        "interpretation": (
            "Prior-drawn ranks and coverage assess internal Bayesian calibration under the stipulated toy model. "
            "Fixed-H0 modes assess conditional bias and coverage. Neither validates realistic detector response."
        ),
    }
    return rows, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    rows, summary = run_experiment(replicates=args.replicates, seed=args.seed)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "replicates.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary_path = args.output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    file_records = {}
    for path in (csv_path, summary_path):
        file_records[path.name] = {"sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                                   "bytes": path.stat().st_size}
    manifest = {
        "schema_version": 1, "git_commit": summary["git_commit"],
        "replicate_count": len(rows), "files": file_records,
        "configuration": summary["configuration"],
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
