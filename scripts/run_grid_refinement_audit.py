"""Grid-refinement audit for posterior summaries in representative two-host cases.

The same synthetic observations are analyzed on 1,001-, 2,001-, and 4,001-point
H0 grids. This is a numerical sensitivity check, not astrophysical validation.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary

H0_VALUES = (52.0, 70.0, 88.0)
Z1 = 0.02
Z2_VALUES = (0.025, 0.04, 0.06)
NOISE_FRACTIONS = (0.03, 0.10)
WEIGHTS = {"correct": np.array([0.7, 0.3]), "reversed": np.array([0.3, 0.7])}
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
GRID_SIZES = (1001, 2001, 4001)


def run_grid_audit(replicates=200, seed=20261030):
    if replicates < 100:
        raise ValueError("Use at least 100 replicates per design cell.")
    rng = np.random.default_rng(seed)
    per_run, summaries = [], []
    for h0 in H0_VALUES:
        for z2 in Z2_VALUES:
            redshifts = np.array([Z1, z2])
            for noise in NOISE_FRACTIONS:
                sigma = noise * C_KM_S * 0.03 / h0
                for replicate in range(replicates):
                    host = int(rng.choice(2, p=[0.7, 0.3]))
                    d_true = C_KM_S * redshifts[host] / h0
                    d_obs = float(rng.normal(d_true, sigma))
                    while d_obs <= 0:
                        d_obs = float(rng.normal(d_true, sigma))
                    for weighting, weights in WEIGHTS.items():
                        summaries_by_grid = {}
                        for size in GRID_SIZES:
                            grid = np.linspace(25.0, 125.0, size)
                            density = posterior_h0(
                                grid, redshifts, d_obs, sigma, weights=weights,
                                prior_min=PRIOR_MIN, prior_max=PRIOR_MAX,
                            )
                            post = posterior_summary(grid, density)
                            summaries_by_grid[size] = post
                            per_run.append({
                                "H0_true": h0, "z2": z2, "noise_fraction": noise,
                                "replicate": replicate, "inference_weighting": weighting,
                                "grid_size": size, "observed_distance_mpc": d_obs,
                                "posterior_median": post["median"],
                                "lower_68": post["lower_68"], "upper_68": post["upper_68"],
                                "lower_95": post["lower_95"], "upper_95": post["upper_95"],
                            })
                        coarse, medium, fine = (summaries_by_grid[x] for x in GRID_SIZES)
                        summaries.append({
                            "H0_true": h0, "z2": z2, "noise_fraction": noise,
                            "inference_weighting": weighting, "replicates": replicates,
                            "median_abs_median_diff_1001_vs_4001": abs(coarse["median"] - fine["median"]),
                            "abs_lower68_diff_1001_vs_4001": abs(coarse["lower_68"] - fine["lower_68"]),
                            "abs_upper68_diff_1001_vs_4001": abs(coarse["upper_68"] - fine["upper_68"]),
                            "abs_lower95_diff_1001_vs_4001": abs(coarse["lower_95"] - fine["lower_95"]),
                            "abs_upper95_diff_1001_vs_4001": abs(coarse["upper_95"] - fine["upper_95"]),
                            "median_abs_median_diff_2001_vs_4001": abs(medium["median"] - fine["median"]),
                        })
    cells = []
    keys = ("H0_true", "z2", "noise_fraction", "inference_weighting")
    for h0 in H0_VALUES:
        for z2 in Z2_VALUES:
            for noise in NOISE_FRACTIONS:
                for weighting in WEIGHTS:
                    group = [r for r in summaries if r["H0_true"] == h0 and r["z2"] == z2
                             and r["noise_fraction"] == noise and r["inference_weighting"] == weighting]
                    cell = {k: group[0][k] for k in keys}
                    cell["replicates"] = len(group)
                    for metric in (
                        "median_abs_median_diff_1001_vs_4001",
                        "abs_lower68_diff_1001_vs_4001", "abs_upper68_diff_1001_vs_4001",
                        "abs_lower95_diff_1001_vs_4001", "abs_upper95_diff_1001_vs_4001",
                        "median_abs_median_diff_2001_vs_4001",
                    ):
                        vals = np.array([r[metric] for r in group])
                        cell[f"median_{metric}"] = float(np.median(vals))
                        cell[f"p95_{metric}"] = float(np.quantile(vals, 0.95))
                        cell[f"max_{metric}"] = float(np.max(vals))
                    cells.append(cell)
    return per_run, cells


def largest_upper95_discrepancy(per_run):
    """Return the observation with the largest 1,001-vs-4,001 upper-95% shift."""
    grouped = {}
    identity = ("H0_true", "z2", "noise_fraction", "replicate", "inference_weighting")
    for row in per_run:
        key = tuple(row[name] for name in identity)
        grouped.setdefault(key, {})[int(row["grid_size"])] = row
    candidates = []
    for key, by_grid in grouped.items():
        if not all(size in by_grid for size in GRID_SIZES):
            raise ValueError(f"Incomplete grid resolutions for case {key}.")
        coarse, medium, fine = (by_grid[size] for size in GRID_SIZES)
        candidates.append({
            **dict(zip(identity, key)),
            "observed_distance_mpc": coarse["observed_distance_mpc"],
            "upper95_1001": coarse["upper_95"],
            "upper95_2001": medium["upper_95"],
            "upper95_4001": fine["upper_95"],
            "abs_diff_1001_vs_4001": abs(coarse["upper_95"] - fine["upper_95"]),
            "abs_diff_2001_vs_4001": abs(medium["upper_95"] - fine["upper_95"]),
        })
    if not candidates:
        raise ValueError("No complete cases found in per_run results.")
    return max(candidates, key=lambda row: row["abs_diff_1001_vs_4001"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--output-dir", default="results/grid-refinement-audit")
    args = parser.parse_args()
    runs, cells = run_grid_audit(args.replicates, args.seed)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, data in (("per_run.csv", runs), ("summary.csv", cells)):
        with (out / name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    print(f"PASS: {len(runs)} posterior fits across {len(cells)} design cells and 3 grid resolutions.")
    print({"largest_upper95_grid_discrepancy": largest_upper95_discrepancy(runs)})
    for metric in ("median_abs_median_diff_1001_vs_4001", "abs_lower68_diff_1001_vs_4001",
                   "abs_upper68_diff_1001_vs_4001", "abs_lower95_diff_1001_vs_4001",
                   "abs_upper95_diff_1001_vs_4001", "median_abs_median_diff_2001_vs_4001"):
        print({"metric": metric,
               "median_of_cell_medians": float(np.median([r["median_" + metric] for r in cells])),
               "max_cell_p95": float(max(r["p95_" + metric] for r in cells)),
               "max_observed_difference": float(max(r["max_" + metric] for r in cells))})


if __name__ == "__main__":
    main()
