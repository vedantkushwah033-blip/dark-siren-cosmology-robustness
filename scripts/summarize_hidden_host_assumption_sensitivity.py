"""Aggregate hidden-host assumption effects over cells with paired bootstrap uncertainty."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

METRICS = ("mean_error", "rmse_error", "coverage_68", "coverage_95", "mean_width_68", "mean_width_95")
CELL_KEYS = ("true_H0", "noise_fraction", "true_hidden_fraction", "true_hidden_redshift")


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Missing or empty CSV: {path}")
    return rows


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def metric_arrays(rows):
    return {
        "error": np.asarray([float(r["median_error"]) for r in rows]),
        "covered_68": np.asarray([int(r["covered_68"]) for r in rows], dtype=float),
        "covered_95": np.asarray([int(r["covered_95"]) for r in rows], dtype=float),
        "width_68": np.asarray([float(r["width_68"]) for r in rows]),
        "width_95": np.asarray([float(r["width_95"]) for r in rows]),
    }


def stats(arr, indices=None):
    if indices is None:
        e, c68, c95, w68, w95 = (arr[k] for k in ("error", "covered_68", "covered_95", "width_68", "width_95"))
    else:
        e, c68, c95, w68, w95 = (
            np.take_along_axis(arr[k][None, :], indices, axis=1)
            for k in ("error", "covered_68", "covered_95", "width_68", "width_95")
        )
    return np.stack((np.mean(e, axis=-1), np.sqrt(np.mean(e**2, axis=-1)),
                     np.mean(c68, axis=-1), np.mean(c95, axis=-1),
                     np.mean(w68, axis=-1), np.mean(w95, axis=-1)), axis=-1)


def summarize(input_dir: Path, output_dir: Path, resamples=1000, seed=20261033):
    if resamples < 100:
        raise ValueError("Use at least 100 bootstrap resamples.")
    raw = read_csv(input_dir / "per_run.csv")
    summaries = read_csv(input_dir / "summary.csv")
    contrasts = read_csv(input_dir / "paired_contrasts.csv")
    cells = sorted({tuple(float(r[k]) for k in CELL_KEYS) for r in raw})
    arms = tuple(dict.fromkeys(r["inference_arm"] for r in raw))
    if "base_omitted" not in arms or "oracle_hidden" not in arms or len(arms) < 3:
        raise ValueError("Raw output is missing base/oracle arms or assumed-host arms")
    if len(raw) % (len(cells) * len(arms)):
        raise ValueError("Raw row count is not divisible by base-cell × arm count")
    reps = len(raw) // (len(cells) * len(arms))
    if reps < 200:
        raise ValueError("Need at least 200 replicates per base cell")
    comparisons = [("base_omitted", "oracle_hidden", "base_omitted_minus_oracle")]
    comparisons.extend((arm, "oracle_hidden", arm + "_minus_oracle") for arm in arms[2:])
    if len(summaries) != len(cells) * len(arms):
        raise ValueError("Summary row count does not match raw factorial design")
    if len(contrasts) != len(cells) * len(comparisons) * len(METRICS):
        raise ValueError("Paired contrast row count does not match raw factorial design")

    groups, observations = {}, {}
    for row in raw:
        cell = tuple(float(row[k]) for k in CELL_KEYS)
        arm = row["inference_arm"]
        if arm not in arms:
            raise ValueError(f"Unexpected arm {arm}")
        groups.setdefault((cell, arm), []).append(row)
        obs_key = (cell, int(row["replicate"]))
        obs = (float(row["true_host_index"]), float(row["true_host_redshift"]),
               float(row["observed_distance_mpc"]), float(row["sigma_mpc"]))
        if obs_key in observations and observations[obs_key] != obs:
            raise ValueError(f"Observation pairing mismatch at {obs_key}")
        observations[obs_key] = obs
    for key, rows in groups.items():
        rows.sort(key=lambda r: int(r["replicate"]))
        if [int(r["replicate"]) for r in rows] != list(range(reps)):
            raise ValueError(f"Missing/duplicate replicate IDs in {key}")

    actual_map = {}
    for row in summaries:
        key = (tuple(float(row[k]) for k in CELL_KEYS), row["inference_arm"])
        if key in actual_map:
            raise ValueError(f"Duplicate summary key {key}")
        actual_map[key] = row
    cell_stats = {}
    for cell in cells:
        for arm in arms:
            values = stats(metric_arrays(groups[(cell, arm)]))
            cell_stats[(cell, arm)] = values
            actual = actual_map.get((cell, arm))
            if actual is None:
                raise ValueError(f"Missing summary row {(cell, arm)}")
            for metric, value in zip(METRICS, values):
                if not np.isclose(float(actual[metric]), float(value), rtol=2e-10, atol=2e-8):
                    raise ValueError(f"Raw-to-summary mismatch at {(cell, arm)}: {metric}")
    if len(actual_map) != len(summaries):
        raise ValueError("Unexpected extra summary rows")

    fractions = sorted({cell[2] for cell in cells})
    rng = np.random.default_rng(seed)
    marginal_summary, marginal_contrasts = [], []
    for fraction in fractions:
        fraction_cells = [cell for cell in cells if cell[2] == fraction]
        indices = {cell: rng.integers(0, reps, size=(resamples, reps)) for cell in fraction_cells}
        boot_stats = {}
        for cell in fraction_cells:
            for arm in arms:
                boot_stats[(cell, arm)] = stats(metric_arrays(groups[(cell, arm)]), indices[cell])
        for arm in arms:
            values = np.stack([cell_stats[(cell, arm)] for cell in fraction_cells])
            row = {"true_hidden_fraction": fraction, "inference_arm": arm,
                   "cells": len(fraction_cells), "replicates_per_cell": reps}
            for j, metric in enumerate(METRICS):
                row["equal_cell_mean_" + metric] = float(np.mean(values[:, j]))
                row["cell_min_" + metric] = float(np.min(values[:, j]))
                row["cell_max_" + metric] = float(np.max(values[:, j]))
            marginal_summary.append(row)
        for left, right, label in comparisons:
            point = np.mean(np.stack([
                cell_stats[(cell, left)] - cell_stats[(cell, right)] for cell in fraction_cells
            ]), axis=0)
            boot = np.mean(np.stack([
                boot_stats[(cell, left)] - boot_stats[(cell, right)] for cell in fraction_cells
            ]), axis=0)
            for j, metric in enumerate(METRICS):
                marginal_contrasts.append({
                    "true_hidden_fraction": fraction, "contrast": label, "left_arm": left,
                    "right_arm": right, "metric": metric, "equal_cell_mean_difference": float(point[j]),
                    "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                    "ci95_high": float(np.quantile(boot[:, j], 0.975)),
                    "cells": len(fraction_cells), "replicates_per_cell": reps,
                    "bootstrap_resamples": resamples, "bootstrap_seed": seed,
                })
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "marginal_summary.csv", marginal_summary)
    write_csv(output_dir / "marginal_contrasts.csv", marginal_contrasts)
    return {"base_cells": len(cells), "true_fraction_strata": len(fractions),
            "marginal_summary_rows": len(marginal_summary),
            "marginal_contrast_rows": len(marginal_contrasts),
            "bootstrap_resamples": resamples, "bootstrap_seed": seed}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="results/hidden-host-assumption-sensitivity")
    parser.add_argument("--output-dir", default="results/hidden-host-assumption-sensitivity")
    parser.add_argument("--bootstrap-resamples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261033)
    args = parser.parse_args()
    print(summarize(Path(args.input_dir), Path(args.output_dir), args.bootstrap_resamples, args.seed))


if __name__ == "__main__":
    main()
