"""Summarize audited factorial results with equal-cell means and paired bootstrap CIs."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

ARMS = ("base_weights", "reversed_weights", "equal_weights", "hidden_host_modelled")
COMPARISONS = (
    ("reversed_weights", "base_weights", "reversed_minus_base"),
    ("base_weights", "hidden_host_modelled", "omitted_minus_modelled_hidden"),
)
METRICS = ("mean_error", "rmse_error", "coverage_68", "coverage_95", "mean_width_68", "mean_width_95")
FACTOR_KEYS = ("true_H0", "noise_fraction", "hidden_fraction", "hidden_redshift")


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"Refusing to write empty file: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _arrays(rows):
    return {
        "error": np.asarray([float(r["median_error"]) for r in rows]),
        "covered_68": np.asarray([float(r["covered_68"]) for r in rows]),
        "covered_95": np.asarray([float(r["covered_95"]) for r in rows]),
        "width_68": np.asarray([float(r["width_68"]) for r in rows]),
        "width_95": np.asarray([float(r["width_95"]) for r in rows]),
    }


def _stats(arr, indices=None):
    if indices is None:
        e, c68, c95, w68, w95 = (arr[k] for k in ("error", "covered_68", "covered_95", "width_68", "width_95"))
    else:
        e, c68, c95, w68, w95 = (
            np.take_along_axis(arr[k][:, None, :], indices, axis=2)
            for k in ("error", "covered_68", "covered_95", "width_68", "width_95")
        )
    return np.stack((
        np.mean(e, axis=-1), np.sqrt(np.mean(e**2, axis=-1)),
        np.mean(c68, axis=-1), np.mean(c95, axis=-1),
        np.mean(w68, axis=-1), np.mean(w95, axis=-1),
    ), axis=-1)


def summarize(root: Path, output_dir: Path, resamples=1000, seed=20261031):
    if resamples < 100:
        raise ValueError("Use at least 100 bootstrap resamples.")
    raw = read_csv(root / "per_run.csv")
    summaries = read_csv(root / "summary.csv")
    contrasts = read_csv(root / "paired_contrasts.csv")
    if len(raw) % len(ARMS):
        raise ValueError("Raw row count is not divisible by inference-arm count")
    cell_map = {}
    for row in raw:
        cell = tuple(float(row[k]) for k in FACTOR_KEYS)
        arm = row["inference_arm"]
        if arm not in ARMS:
            raise ValueError(f"Unexpected arm {arm}")
        cell_map.setdefault(cell, {name: [] for name in ARMS})[arm].append(row)
    if not cell_map:
        raise ValueError("No factorial cells found")
    reps_set = set()
    for cell, arm_map in cell_map.items():
        for arm in ARMS:
            block = sorted(arm_map[arm], key=lambda r: int(r["replicate"]))
            arm_map[arm] = block
            ids = [int(r["replicate"]) for r in block]
            if ids != list(range(len(ids))):
                raise ValueError(f"Missing, duplicated, or unordered replicates in {cell}/{arm}")
            reps_set.add(len(ids))
        if len(reps_set) != 1:
            raise ValueError("Unequal replicate counts across factorial arms/cells")
        for i in range(len(arm_map[ARMS[0]])):
            values = {arm_map[arm][i]["observed_distance_mpc"] for arm in ARMS}
            if len(values) != 1:
                raise ValueError(f"Observation pairing mismatch in {cell}, replicate {i}")
    replicates = next(iter(reps_set))
    expected_summary_count = len(cell_map) * len(ARMS)
    if len(summaries) != expected_summary_count:
        raise ValueError("Summary row count does not match raw factorial design")
    if len(contrasts) != len(cell_map) * len(COMPARISONS) * len(METRICS):
        raise ValueError("Cell contrast row count does not match raw factorial design")

    # Reconstruct cell metrics from raw records and verify the supplied summary.
    expected_summary = []
    for cell in sorted(cell_map):
        truth, noise, fraction, hidden_z = cell
        for arm in ARMS:
            stats = _stats(_arrays(cell_map[cell][arm]))
            expected_summary.append({
                "true_H0": truth, "noise_fraction": noise, "hidden_fraction": fraction,
                "hidden_redshift": hidden_z, "inference_arm": arm, "replicates": replicates,
                **dict(zip(METRICS, [float(v) for v in stats])),
            })
    summary_index = {(tuple(float(r[k]) for k in FACTOR_KEYS), r["inference_arm"]): r for r in summaries}
    if len(summary_index) != len(summaries):
        raise ValueError("Duplicate cell-summary keys")
    for expected in expected_summary:
        key = (tuple(float(expected[k]) for k in FACTOR_KEYS), expected["inference_arm"])
        actual = summary_index[key]
        for metric in METRICS:
            if not np.isclose(float(actual[metric]), float(expected[metric]), rtol=2e-10, atol=2e-8):
                raise ValueError(f"Raw-to-summary mismatch at {key}, {metric}")

    marginal_arms = []
    fractions = sorted({cell[2] for cell in cell_map})
    for fraction in fractions:
        cells = [cell for cell in sorted(cell_map) if cell[2] == fraction]
        for arm in ARMS:
            rows = [r for r in expected_summary if float(r["hidden_fraction"]) == fraction
                    and r["inference_arm"] == arm]
            record = {"hidden_fraction": fraction, "inference_arm": arm, "cells": len(cells),
                      "replicates_per_cell": replicates}
            for metric in METRICS:
                vals = np.asarray([float(r[metric]) for r in rows])
                record["equal_cell_mean_" + metric] = float(np.mean(vals))
                record["cell_min_" + metric] = float(np.min(vals))
                record["cell_max_" + metric] = float(np.max(vals))
            marginal_arms.append(record)

    rng = np.random.default_rng(seed)
    marginal_contrasts = []
    for fraction in fractions:
        cells = [cell for cell in sorted(cell_map) if cell[2] == fraction]
        arrays_by_arm = {
            arm: {k: np.stack([_arrays(cell_map[cell][arm])[k] for cell in cells], axis=0)
                  for k in ("error", "covered_68", "covered_95", "width_68", "width_95")}
            for arm in ARMS
        }
        indices = rng.integers(0, replicates, size=(len(cells), resamples, replicates))
        for left, right, label in COMPARISONS:
            point = np.mean(_stats(arrays_by_arm[left]) - _stats(arrays_by_arm[right]), axis=0)
            boot = _stats(arrays_by_arm[left], indices) - _stats(arrays_by_arm[right], indices)
            boot = np.mean(boot, axis=0) if boot.ndim == 2 else boot
            # For bootstrap arrays, _stats returns (cell, resample, metric).
            if boot.shape == (len(cells), resamples, len(METRICS)):
                boot = np.mean(boot, axis=0)
            for j, metric in enumerate(METRICS):
                marginal_contrasts.append({
                    "hidden_fraction": fraction, "contrast": label, "metric": metric,
                    "equal_cell_mean_difference": float(point[j]),
                    "ci95_low": float(np.quantile(boot[:, j], 0.025)),
                    "ci95_high": float(np.quantile(boot[:, j], 0.975)),
                    "cells": len(cells), "replicates_per_cell": replicates,
                    "bootstrap_resamples": resamples, "bootstrap_seed": seed,
                })
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "marginal_arm_summary.csv", marginal_arms)
    write_csv(output_dir / "marginal_paired_contrasts.csv", marginal_contrasts)
    return {"cells": len(cell_map), "replicates_per_cell": replicates,
            "marginal_arm_rows": len(marginal_arms), "marginal_contrast_rows": len(marginal_contrasts),
            "bootstrap_resamples": resamples, "bootstrap_seed": seed}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", default="results/host-misspecification-sensitivity-surface")
    parser.add_argument("--output-dir", default="results/host-misspecification-sensitivity-surface")
    parser.add_argument("--bootstrap-resamples", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261031)
    args = parser.parse_args()
    print(summarize(Path(args.input_dir), Path(args.output_dir), args.bootstrap_resamples, args.seed))


if __name__ == "__main__":
    main()
