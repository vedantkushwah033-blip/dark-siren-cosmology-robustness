"""Produce equal-cell marginal results and stratified paired-bootstrap intervals."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

PRIORS = ((50.0, 90.0), (40.0, 100.0), (30.0, 120.0))
ARMS = ("base_weights", "reversed_weights", "equal_weights", "hidden_host_modelled")
ARM_COMPARISONS = (
    ("reversed_weights", "base_weights", "reversed_minus_base"),
    ("base_weights", "hidden_host_modelled", "omitted_minus_modelled_hidden"),
)
PRIOR_COMPARISONS = ((40.0, 100.0, "40_100_minus_50_90"), (30.0, 120.0, "30_120_minus_50_90"))
METRICS = ("mean_error", "rmse_error", "coverage_68", "coverage_95", "mean_width_68", "mean_width_95")
CELL_KEYS = ("true_H0", "noise_fraction", "hidden_fraction", "hidden_redshift")


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


def _metric_arrays(records):
    return {
        "error": np.asarray([float(r["median_error"]) for r in records]),
        "covered_68": np.asarray([int(r["covered_68"]) for r in records], dtype=float),
        "covered_95": np.asarray([int(r["covered_95"]) for r in records], dtype=float),
        "width_68": np.asarray([float(r["width_68"]) for r in records]),
        "width_95": np.asarray([float(r["width_95"]) for r in records]),
    }


def _stats(arr, indices=None):
    if indices is None:
        e, c68, c95, w68, w95 = (arr[k] for k in ("error", "covered_68", "covered_95", "width_68", "width_95"))
    else:
        e, c68, c95, w68, w95 = (
            np.take_along_axis(arr[k][None, :], indices, axis=1)
            for k in ("error", "covered_68", "covered_95", "width_68", "width_95")
        )
    return np.stack((
        np.mean(e, axis=-1), np.sqrt(np.mean(e**2, axis=-1)),
        np.mean(c68, axis=-1), np.mean(c95, axis=-1),
        np.mean(w68, axis=-1), np.mean(w95, axis=-1),
    ), axis=-1)


def summarize(input_dir: Path, output_dir: Path, resamples=1000, seed=20261031):
    if resamples < 100:
        raise ValueError("Use at least 100 bootstrap resamples.")
    raw = read_csv(input_dir / "per_run.csv")
    summary = read_csv(input_dir / "summary.csv")
    arm_contrasts = read_csv(input_dir / "arm_contrasts.csv")
    prior_contrasts = read_csv(input_dir / "prior_contrasts.csv")
    if len(raw) % (len(PRIORS) * len(ARMS)):
        raise ValueError("Raw row count is not divisible by prior × arm count")
    reps = len(raw) // (len({tuple(float(r[k]) for k in CELL_KEYS) for r in raw}) * len(PRIORS) * len(ARMS))
    if reps < 200:
        raise ValueError("Raw data has fewer than 200 replicates per cell")
    cells = sorted({tuple(float(r[k]) for k in CELL_KEYS) for r in raw})
    expected_summary_n = len(cells) * len(PRIORS) * len(ARMS)
    expected_arm_n = len(cells) * len(PRIORS) * len(ARM_COMPARISONS) * len(METRICS)
    fractions = sorted({c[2] for c in cells})
    expected_prior_n = len(cells) * len(ARMS) * len(PRIOR_COMPARISONS) * len(METRICS)
    if len(summary) != expected_summary_n or len(arm_contrasts) != expected_arm_n or len(prior_contrasts) != expected_prior_n:
        raise ValueError("Summary/contrast row count does not match raw factorial design")

    groups = {}
    observations = {}
    for row in raw:
        cell = tuple(float(row[k]) for k in CELL_KEYS)
        prior = (float(row["prior_min"]), float(row["prior_max"]))
        arm = row["inference_arm"]
        if prior not in PRIORS or arm not in ARMS:
            raise ValueError(f"Unexpected prior or inference arm: {prior}/{arm}")
        groups.setdefault((cell, prior, arm), []).append(row)
        obs_key = (cell, int(row["replicate"]))
        obs_value = (float(row["true_host_index"]), float(row["true_host_redshift"]),
                     float(row["observed_distance_mpc"]), float(row["sigma_mpc"]))
        if obs_key in observations and observations[obs_key] != obs_value:
            raise ValueError(f"Prior/arm pairing mismatch for {obs_key}")
        observations[obs_key] = obs_value
    for key, rows in groups.items():
        rows.sort(key=lambda r: int(r["replicate"]))
        if [int(r["replicate"]) for r in rows] != list(range(reps)):
            raise ValueError(f"Missing/duplicate replicates in group {key}")

    # Reconstruct cell-level summaries from raw per-run metrics before aggregation.
    actual_map = {}
    for row in summary:
        key = (tuple(float(row[k]) for k in CELL_KEYS), (float(row["prior_min"]), float(row["prior_max"])), row["inference_arm"])
        if key in actual_map:
            raise ValueError(f"Duplicate summary key {key}")
        actual_map[key] = row
    cell_stats = {}
    for cell in cells:
        for prior in PRIORS:
            for arm in ARMS:
                rows = groups[(cell, prior, arm)]
                values = _stats(_metric_arrays(rows))
                key = (cell, prior, arm)
                cell_stats[key] = values
                actual = actual_map.get(key)
                if actual is None:
                    raise ValueError(f"Missing summary row {key}")
                for metric, value in zip(METRICS, values):
                    if not np.isclose(float(actual[metric]), float(value), rtol=2e-10, atol=2e-8):
                        raise ValueError(f"Raw-to-summary mismatch at {key}: {metric}")
    if len(actual_map) != len(summary):
        raise ValueError("Unexpected extra summary rows")

    marginal_summary = []
    marginal_arm_contrasts = []
    marginal_prior_contrasts = []
    rng = np.random.default_rng(seed)
    for fraction in fractions:
        fraction_cells = [c for c in cells if c[2] == fraction]
        # One replicate-index array per base cell; reuse it across all priors and arms.
        boot_indices = {cell: rng.integers(0, reps, size=(resamples, reps)) for cell in fraction_cells}
        for prior in PRIORS:
            for arm in ARMS:
                values = np.stack([cell_stats[(cell, prior, arm)] for cell in fraction_cells])
                row = {"hidden_fraction": fraction, "prior_min": prior[0], "prior_max": prior[1],
                       "inference_arm": arm, "cells": len(fraction_cells), "replicates_per_cell": reps}
                for j, metric in enumerate(METRICS):
                    row["equal_cell_mean_" + metric] = float(np.mean(values[:, j]))
                    row["cell_min_" + metric] = float(np.min(values[:, j]))
                    row["cell_max_" + metric] = float(np.max(values[:, j]))
                marginal_summary.append(row)

        # Calculate bootstrap metric arrays per cell/prior/arm, then contrast paired arms/priors.
        boot_stats = {}
        for cell in fraction_cells:
            indices = boot_indices[cell]
            for prior in PRIORS:
                for arm in ARMS:
                    boot_stats[(cell, prior, arm)] = _stats(_metric_arrays(groups[(cell, prior, arm)]), indices)

        for prior in PRIORS:
            for left, right, label in ARM_COMPARISONS:
                point = np.mean(np.stack([
                    cell_stats[(cell, prior, left)] - cell_stats[(cell, prior, right)]
                    for cell in fraction_cells
                ]), axis=0)
                boot = np.mean(np.stack([
                    boot_stats[(cell, prior, left)] - boot_stats[(cell, prior, right)]
                    for cell in fraction_cells
                ]), axis=0)
                for j, metric in enumerate(METRICS):
                    marginal_arm_contrasts.append({
                        "hidden_fraction": fraction, "prior_min": prior[0], "prior_max": prior[1],
                        "contrast": label, "metric": metric, "equal_cell_mean_difference": float(point[j]),
                        "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                        "ci95_high": float(np.quantile(boot[:, j], 0.975)),
                        "cells": len(fraction_cells), "replicates_per_cell": reps,
                        "bootstrap_resamples": resamples, "bootstrap_seed": seed,
                    })
        reference = PRIORS[0]
        for arm in ARMS:
            for lo, hi, label in PRIOR_COMPARISONS:
                point = np.mean(np.stack([
                    cell_stats[(cell, (lo, hi), arm)] - cell_stats[(cell, reference, arm)]
                    for cell in fraction_cells
                ]), axis=0)
                boot = np.mean(np.stack([
                    boot_stats[(cell, (lo, hi), arm)] - boot_stats[(cell, reference, arm)]
                    for cell in fraction_cells
                ]), axis=0)
                for j, metric in enumerate(METRICS):
                    marginal_prior_contrasts.append({
                        "hidden_fraction": fraction, "inference_arm": arm, "contrast": label,
                        "metric": metric, "equal_cell_mean_difference": float(point[j]),
                        "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                        "ci95_high": float(np.quantile(boot[:, j], 0.975)),
                        "cells": len(fraction_cells), "replicates_per_cell": reps,
                        "bootstrap_resamples": resamples, "bootstrap_seed": seed,
                    })
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "marginal_summary.csv", marginal_summary)
    write_csv(output_dir / "marginal_arm_contrasts.csv", marginal_arm_contrasts)
    write_csv(output_dir / "marginal_prior_contrasts.csv", marginal_prior_contrasts)
    return {"base_cells": len(cells), "hidden_fraction_strata": len(fractions),
            "marginal_summary_rows": len(marginal_summary),
            "marginal_arm_contrast_rows": len(marginal_arm_contrasts),
            "marginal_prior_contrast_rows": len(marginal_prior_contrasts),
            "bootstrap_resamples": resamples, "bootstrap_seed": seed}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="results/host-prior-support-factorial")
    parser.add_argument("--output-dir", default="results/host-prior-support-factorial")
    parser.add_argument("--bootstrap-resamples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261031)
    args = parser.parse_args()
    print(summarize(Path(args.input_dir), Path(args.output_dir), args.bootstrap_resamples, args.seed))


if __name__ == "__main__":
    main()
