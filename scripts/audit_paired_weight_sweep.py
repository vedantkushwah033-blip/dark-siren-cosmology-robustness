"""Independent integrity and summary audit for paired host-weight sweep outputs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

C_KM_S = 299792.458
H0_TRUE = 70.0
Z1 = 0.02
Z2_VALUES = (0.025, 0.04, 0.06)
NOISE_FRACTIONS = (0.03, 0.10, 0.20)
PRIOR_SUPPORTS = {
    "narrow_50_90": (50.0, 90.0),
    "broad_20_220": (20.0, 220.0),
}
NUMERIC_FIELDS = (
    "prior_min", "prior_max", "z1", "z2", "redshift_separation", "noise_fraction",
    "correct_mean_median_error", "reversed_mean_median_error", "correct_rmse", "reversed_rmse",
    "rmse_diff_reversed_minus_correct", "rmse_diff_ci95_low", "rmse_diff_ci95_high",
    "correct_coverage_68", "reversed_coverage_68", "coverage_68_diff_reversed_minus_correct",
    "coverage_68_diff_ci95_low", "coverage_68_diff_ci95_high",
    "correct_mean_width_68", "reversed_mean_width_68",
    "mean_paired_abs_error_diff_reversed_minus_correct",
    "median_paired_abs_error_diff_reversed_minus_correct",
    "fraction_reversed_has_larger_abs_error",
    "paired_abs_error_diff_ci95_low", "paired_abs_error_diff_ci95_high",
)
INTEGER_FIELDS = (
    "replicates", "bootstrap_resamples", "bootstrap_seed_abs_error",
    "bootstrap_seed_coverage", "bootstrap_seed_rmse",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def percentile_ci(values, rng, resamples, statistic="mean"):
    values = np.asarray(values, dtype=float)
    draws = np.empty(resamples)
    for start in range(0, resamples, 250):
        stop = min(start + 250, resamples)
        idx = rng.integers(0, len(values), size=(stop - start, len(values)))
        if statistic == "mean":
            draws[start:stop] = values[idx].mean(axis=1)
        elif statistic == "rmse_delta":
            draws[start:stop] = (
                np.sqrt(np.mean(values[idx, 1] ** 2, axis=1))
                - np.sqrt(np.mean(values[idx, 0] ** 2, axis=1))
            )
        else:
            raise ValueError(f"Unknown statistic: {statistic}")
    lo, hi = np.quantile(draws, [0.025, 0.975])
    return float(lo), float(hi)


def reconstruct_summary(rows, support_name, bounds, replicates, bootstrap_resamples,
                        base_seed, cell_index, support_index):
    correct_error = np.array([float(r[f"{support_name}_correct_median_error"]) for r in rows])
    reversed_error = np.array([float(r[f"{support_name}_reversed_median_error"]) for r in rows])
    abs_diffs = np.array([float(r[f"paired_abs_error_difference_{support_name}"]) for r in rows])
    coverage_diffs = np.array([
        float(r[f"{support_name}_reversed_covered_68"]) - float(r[f"{support_name}_correct_covered_68"])
        for r in rows
    ])
    bootstrap_seed = base_seed + 104729 * cell_index + 1009 * support_index
    abs_ci = percentile_ci(abs_diffs, np.random.default_rng(bootstrap_seed), bootstrap_resamples)
    coverage_ci = percentile_ci(coverage_diffs, np.random.default_rng(bootstrap_seed + 1), bootstrap_resamples)
    rmse_ci = percentile_ci(
        np.column_stack([correct_error, reversed_error]),
        np.random.default_rng(bootstrap_seed + 2), bootstrap_resamples, statistic="rmse_delta",
    )
    correct_coverage = np.array([float(r[f"{support_name}_correct_covered_68"]) for r in rows])
    reversed_coverage = np.array([float(r[f"{support_name}_reversed_covered_68"]) for r in rows])
    correct_width = np.array([float(r[f"{support_name}_correct_width_68"]) for r in rows])
    reversed_width = np.array([float(r[f"{support_name}_reversed_width_68"]) for r in rows])
    z2 = float(rows[0]["z2"])
    noise = float(rows[0]["noise_fraction"])
    rmse_correct = float(np.sqrt(np.mean(correct_error ** 2)))
    rmse_reversed = float(np.sqrt(np.mean(reversed_error ** 2)))
    return {
        "prior_support": support_name, "prior_min": bounds[0], "prior_max": bounds[1],
        "z1": Z1, "z2": z2, "redshift_separation": z2 - Z1,
        "noise_fraction": noise, "replicates": replicates,
        "correct_mean_median_error": float(correct_error.mean()),
        "reversed_mean_median_error": float(reversed_error.mean()),
        "correct_rmse": rmse_correct, "reversed_rmse": rmse_reversed,
        "rmse_diff_reversed_minus_correct": rmse_reversed - rmse_correct,
        "rmse_diff_ci95_low": rmse_ci[0], "rmse_diff_ci95_high": rmse_ci[1],
        "correct_coverage_68": float(correct_coverage.mean()),
        "reversed_coverage_68": float(reversed_coverage.mean()),
        "coverage_68_diff_reversed_minus_correct": float(coverage_diffs.mean()),
        "coverage_68_diff_ci95_low": coverage_ci[0], "coverage_68_diff_ci95_high": coverage_ci[1],
        "correct_mean_width_68": float(correct_width.mean()),
        "reversed_mean_width_68": float(reversed_width.mean()),
        "mean_paired_abs_error_diff_reversed_minus_correct": float(abs_diffs.mean()),
        "median_paired_abs_error_diff_reversed_minus_correct": float(np.median(abs_diffs)),
        "fraction_reversed_has_larger_abs_error": float(np.mean(abs_diffs > 0)),
        "paired_abs_error_diff_ci95_low": abs_ci[0],
        "paired_abs_error_diff_ci95_high": abs_ci[1],
        "bootstrap_resamples": bootstrap_resamples,
        "bootstrap_seed_abs_error": bootstrap_seed,
        "bootstrap_seed_coverage": bootstrap_seed + 1,
        "bootstrap_seed_rmse": bootstrap_seed + 2,
    }


def audit_directory(root: Path) -> dict[str, object]:
    raw_path, summary_path, metadata_path = (
        root / "per_run.csv", root / "summary.csv", root / "run_metadata.json"
    )
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("audit_status") != "PENDING":
        raise ValueError("run metadata must start in PENDING audit state")
    for name in ("per_run.csv", "summary.csv"):
        if metadata.get("sha256", {}).get(name) != sha256(root / name):
            raise ValueError(f"metadata hash mismatch for {name}")
    raw_rows, actual_summaries = read_csv(raw_path), read_csv(summary_path)
    replicates = int(metadata["replicates_per_redshift_noise_cell"])
    base_seed = int(metadata["base_seed"])
    bootstrap_resamples = int(metadata["bootstrap_resamples"])
    if replicates < 200 or bootstrap_resamples < 1000:
        raise ValueError("replicate/bootstrap counts are below protocol minimums")
    if len(raw_rows) != len(Z2_VALUES) * len(NOISE_FRACTIONS) * replicates:
        raise ValueError(f"unexpected raw row count: {len(raw_rows)}")
    if len(actual_summaries) != len(Z2_VALUES) * len(NOISE_FRACTIONS) * len(PRIOR_SUPPORTS):
        raise ValueError(f"unexpected summary row count: {len(actual_summaries)}")

    groups = {}
    for index, row in enumerate(raw_rows):
        z2, noise, rep = float(row["z2"]), float(row["noise_fraction"]), int(row["replicate"])
        if z2 not in Z2_VALUES or noise not in NOISE_FRACTIONS or not 0 <= rep < replicates:
            raise ValueError(f"invalid design key in raw row {index}")
        if int(row["base_seed"]) != base_seed:
            raise ValueError(f"base seed mismatch in raw row {index}")
        if float(row["true_H0"]) != H0_TRUE:
            raise ValueError(f"unexpected injected H0 in raw row {index}")
        host_index = int(row["true_host_index"])
        expected_z = Z1 if host_index == 0 else z2 if host_index == 1 else None
        if expected_z is None or not math.isclose(float(row["true_host_redshift"]), expected_z, abs_tol=1e-12):
            raise ValueError(f"host index/redshift mismatch in raw row {index}")
        expected_sigma = noise * C_KM_S * 0.03 / H0_TRUE
        if not math.isclose(float(row["sigma_mpc"]), expected_sigma, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError(f"sigma mismatch in raw row {index}")
        key = (z2, noise)
        group = groups.setdefault(key, {})
        if rep in group:
            raise ValueError(f"duplicate replicate {rep} in cell {key}")
        group[rep] = row
        for support_name in PRIOR_SUPPORTS:
            for label in ("correct", "reversed"):
                for field in ("median", "median_error", "lower_68", "upper_68", "covered_68", "width_68"):
                    value = float(row[f"{support_name}_{label}_{field}"])
                    if not math.isfinite(value):
                        raise ValueError(f"non-finite {support_name}_{label}_{field} in raw row {index}")
            for label in ("correct", "reversed"):
                median = float(row[f"{support_name}_{label}_median"])
                error = float(row[f"{support_name}_{label}_median_error"])
                lower = float(row[f"{support_name}_{label}_lower_68"])
                upper = float(row[f"{support_name}_{label}_upper_68"])
                covered = int(row[f"{support_name}_{label}_covered_68"])
                width = float(row[f"{support_name}_{label}_width_68"])
                if not math.isclose(error, median - H0_TRUE, rel_tol=1e-10, abs_tol=1e-10):
                    raise ValueError(f"median error inconsistent with median in raw row {index}")
                if not math.isclose(width, upper - lower, rel_tol=1e-10, abs_tol=1e-10):
                    raise ValueError(f"interval width inconsistent with endpoints in raw row {index}")
                if covered != int(lower <= H0_TRUE <= upper):
                    raise ValueError(f"coverage indicator inconsistent with interval in raw row {index}")
            for support_name in PRIOR_SUPPORTS:
                expected_abs_diff = (
                    abs(float(row[f"{support_name}_reversed_median"]) - H0_TRUE)
                    - abs(float(row[f"{support_name}_correct_median"]) - H0_TRUE)
                )
                if not math.isclose(float(row[f"paired_abs_error_difference_{support_name}"]),
                                    expected_abs_diff, rel_tol=1e-10, abs_tol=1e-10):
                    raise ValueError(f"paired absolute-error difference inconsistent in raw row {index}")

    expected_keys = {(z2, noise) for z2 in Z2_VALUES for noise in NOISE_FRACTIONS}
    if set(groups) != expected_keys or any(set(g) != set(range(replicates)) for g in groups.values()):
        raise ValueError("raw data do not contain complete replicate blocks for every design cell")

    expected_summaries = []
    cell_index = 0
    for z2 in Z2_VALUES:
        for noise in NOISE_FRACTIONS:
            cell_rows = [groups[(z2, noise)][rep] for rep in range(replicates)]
            for support_index, (support_name, bounds) in enumerate(PRIOR_SUPPORTS.items()):
                expected_summaries.append(reconstruct_summary(
                    cell_rows, support_name, bounds, replicates, bootstrap_resamples,
                    base_seed, cell_index, support_index,
                ))
            cell_index += 1

    if len(actual_summaries) != len(expected_summaries):
        raise ValueError("summary count mismatch")
    for i, (actual, expected) in enumerate(zip(actual_summaries, expected_summaries)):
        if actual["prior_support"] != expected["prior_support"]:
            raise ValueError(f"prior-support/order mismatch at summary row {i}")
        for field in NUMERIC_FIELDS:
            observed, wanted = float(actual[field]), float(expected[field])
            if not math.isfinite(observed) or not math.isclose(observed, wanted, rel_tol=1e-10, abs_tol=1e-10):
                raise ValueError(f"{field} mismatch at summary row {i}: {observed} != {wanted}")
        for field in INTEGER_FIELDS:
            if int(actual[field]) != int(expected[field]):
                raise ValueError(f"{field} mismatch at summary row {i}")
    metadata["audit_status"] = "PASS"
    metadata["independent_audit"] = {
        "raw_records": len(raw_rows),
        "summary_rows": len(actual_summaries),
        "design_cells": len(expected_keys),
        "prior_support_cells": len(actual_summaries),
        "posterior_fits": 4 * len(raw_rows),
        "hashes_verified": True,
        "all_summary_metrics_and_bootstrap_intervals_reconstructed": True,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    report = {
        "audit": "PASS",
        "raw_records": len(raw_rows),
        "summary_rows": len(actual_summaries),
        "design_cells": len(expected_keys),
        "prior_support_cells": len(actual_summaries),
        "posterior_fits": 4 * len(raw_rows),
        "bootstrap_resamples_per_metric_cell": bootstrap_resamples,
        "per_run_sha256": sha256(raw_path),
        "summary_sha256": sha256(summary_path),
        "independent_reconstruction": "All summary point estimates and deterministic paired-bootstrap intervals reconstructed from per-run CSV.",
    }
    (root / "audit_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", nargs="?", default="results/paired-weight-sweep")
    args = parser.parse_args()
    print(json.dumps(audit_directory(Path(args.output_dir)), indent=2))


if __name__ == "__main__":
    main()
