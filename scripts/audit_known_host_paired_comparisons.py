"""Independent verification of paired known-host calibration contrasts.

This module intentionally does not import the primary paired-analysis script.
It reconstructs replicate-level contrasts from the audited raw CSV, independently
recomputes point estimates and deterministic percentile-bootstrap intervals, and
checks the output's provenance hashes before accepting the paired table.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

TRUTHS = (50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0)
NOISES = (0.05, 0.10, 0.15, 0.20)
METRICS = (
    "absolute_median_error",
    "width_68",
    "width_95",
    "coverage_68_indicator",
    "coverage_95_indicator",
)
NUMERIC_FIELDS = (
    "mean_difference",
    "median_difference",
    "bootstrap_mean_difference_ci95_low",
    "bootstrap_mean_difference_ci95_high",
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def metric_value(row: dict[str, str], truth: float, metric: str) -> float:
    median = float(row["median"])
    if metric == "absolute_median_error":
        return abs(median - truth)
    if metric == "width_68":
        return float(row["upper_68"]) - float(row["lower_68"])
    if metric == "width_95":
        return float(row["upper_95"]) - float(row["lower_95"])
    if metric == "coverage_68_indicator":
        return float(float(row["lower_68"]) <= truth <= float(row["upper_68"]))
    if metric == "coverage_95_indicator":
        return float(float(row["lower_95"]) <= truth <= float(row["upper_95"]))
    raise ValueError(f"unknown metric: {metric}")


def reconstruct_differences(
    raw_rows: list[dict[str, str]], replicates: int, seed: int
) -> dict[tuple[float, float, str], np.ndarray]:
    """Rebuild paired replicate differences without using analysis output."""
    cells: dict[tuple[float, int], dict[float, dict[str, str]]] = {}
    for index, row in enumerate(raw_rows):
        try:
            truth = float(row["true_h0"])
            noise = float(row["fractional_distance_error"])
            run = int(row["run"])
            row_seed = int(row["seed"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid raw record at row {index}") from exc
        if truth not in TRUTHS or noise not in NOISES:
            raise ValueError(f"unexpected truth/noise combination at row {index}")
        truth_index = TRUTHS.index(truth)
        if not 0 <= run < replicates or row_seed != seed + truth_index * replicates + run:
            raise ValueError(f"incorrect run index or seed at row {index}")
        group = cells.setdefault((truth, run), {})
        if noise in group:
            raise ValueError(f"duplicate paired member at row {index}")
        group[noise] = row

    expected_groups = {(truth, run) for truth in TRUTHS for run in range(replicates)}
    if set(cells) != expected_groups or any(set(group) != set(NOISES) for group in cells.values()):
        raise ValueError("raw data are not a complete truth × replicate × noise design")

    differences: dict[tuple[float, float, str], np.ndarray] = {}
    for truth in TRUTHS:
        for comparison in NOISES[1:]:
            for metric in METRICS:
                values = []
                for run in range(replicates):
                    group = cells[(truth, run)]
                    reference_row, candidate_row = group[NOISES[0]], group[comparison]
                    for field in ("host_redshift", "true_distance_mpc"):
                        if not math.isclose(
                            float(reference_row[field]), float(candidate_row[field]),
                            rel_tol=1e-12, abs_tol=1e-12
                        ):
                            raise ValueError(f"pairing mismatch in {field} at truth={truth}, run={run}")
                    reference_value = metric_value(reference_row, truth, metric)
                    candidate_value = metric_value(candidate_row, truth, metric)
                    values.append(candidate_value - reference_value)
                differences[(truth, comparison, metric)] = np.asarray(values, dtype=float)
    return differences


def reconstruct_expected_table(
    differences: dict[tuple[float, float, str], np.ndarray],
    replicates: int,
    bootstrap_resamples: int,
    base_seed: int,
) -> list[dict[str, object]]:
    rng = np.random.default_rng(base_seed + 7919)
    expected = []
    for truth in TRUTHS:
        for comparison in NOISES[1:]:
            for metric in METRICS:
                values = differences[(truth, comparison, metric)]
                draws = rng.integers(0, replicates, size=(bootstrap_resamples, replicates))
                bootstrap_means = values[draws].mean(axis=1)
                low, high = np.quantile(bootstrap_means, [0.025, 0.975])
                expected.append({
                    "true_h0": truth,
                    "reference_fractional_error": NOISES[0],
                    "comparison_fractional_error": comparison,
                    "metric": metric,
                    "n_paired_replicates": replicates,
                    "mean_difference": float(values.mean()),
                    "median_difference": float(np.median(values)),
                    "bootstrap_mean_difference_ci95_low": float(low),
                    "bootstrap_mean_difference_ci95_high": float(high),
                    "bootstrap_resamples": bootstrap_resamples,
                    "bootstrap_seed": base_seed + 7919,
                    "difference_definition": "comparison minus 5% reference",
                })
    return expected


def validate_pair_table(
    raw_rows: list[dict[str, str]],
    actual_rows: list[dict[str, str]],
    replicates: int,
    bootstrap_resamples: int,
    base_seed: int,
    tolerance: float = 1e-10,
) -> dict[str, int]:
    expected_count = len(TRUTHS) * (len(NOISES) - 1) * len(METRICS)
    if len(actual_rows) != expected_count:
        raise ValueError(f"expected {expected_count} paired rows, found {len(actual_rows)}")
    differences = reconstruct_differences(raw_rows, replicates, base_seed)
    expected = reconstruct_expected_table(differences, replicates, bootstrap_resamples, base_seed)
    seen = set()
    for index, (actual, target) in enumerate(zip(actual_rows, expected)):
        key = (float(actual["true_h0"]), float(actual["comparison_fractional_error"]), actual["metric"])
        if key in seen:
            raise ValueError(f"duplicate paired comparison key: {key}")
        seen.add(key)
        expected_key = (target["true_h0"], target["comparison_fractional_error"], target["metric"])
        if key != expected_key:
            raise ValueError(f"unexpected row order/key at paired row {index}: {key} != {expected_key}")
        for field in NUMERIC_FIELDS:
            observed = float(actual[field])
            wanted = float(target[field])
            if not math.isfinite(observed) or not math.isclose(
                observed, wanted, rel_tol=tolerance, abs_tol=tolerance
            ):
                raise ValueError(
                    f"{field} mismatch at paired row {index}: observed={observed}, expected={wanted}"
                )
        for field in ("reference_fractional_error", "n_paired_replicates",
                      "bootstrap_resamples", "bootstrap_seed"):
            if float(actual[field]) != float(target[field]):
                raise ValueError(f"{field} mismatch at paired row {index}")
        if actual["difference_definition"] != target["difference_definition"]:
            raise ValueError(f"difference definition mismatch at paired row {index}")
    return {"raw_records": len(raw_rows), "paired_comparisons": len(actual_rows),
            "paired_groups": len(TRUTHS) * replicates}


def audit_directory(directory: Path) -> dict[str, object]:
    audit_path = directory / "audit_report.json"
    metadata_path = directory / "cell_summary.json"
    paired_metadata_path = directory / "paired_comparisons.json"
    raw_path = directory / "per_run_results.csv"
    paired_path = directory / "paired_comparisons.csv"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    paired_meta = json.loads(paired_metadata_path.read_text(encoding="utf-8"))
    if audit.get("audit") != "PASS":
        raise ValueError("raw calibration audit is not PASS")
    if audit.get("posterior_summaries_independently_reconstructed_from_raw_observations") is not True:
        raise ValueError("raw audit does not confirm independent posterior reconstruction")
    for filename in ("per_run_results.csv", "cell_summary.csv", "cell_summary.json"):
        expected_hash = audit.get("input_sha256", {}).get(filename)
        path = directory / filename
        if not isinstance(expected_hash, str) or sha256(path) != expected_hash:
            raise ValueError(f"raw audit hash mismatch for {filename}")
    raw_hash = sha256(raw_path)
    if paired_meta.get("raw_input_sha256") != raw_hash:
        raise ValueError("paired metadata raw-input hash does not match raw CSV")
    if paired_meta.get("paired_comparisons_sha256") != sha256(paired_path):
        raise ValueError("paired-comparison CSV hash does not match paired metadata")
    if paired_meta.get("audit") != "PASS" or paired_meta.get("source_audit_run_id") != audit.get("workflow_run_id"):
        raise ValueError("paired metadata does not point to the passing source audit")
    replicates = int(metadata["replicates_per_cell"])
    seed = int(metadata["seed"])
    bootstrap_resamples = int(paired_meta["bootstrap_resamples"])
    if int(paired_meta["replicates_per_truth"]) != replicates or int(paired_meta["base_seed"]) != seed:
        raise ValueError("paired metadata does not match raw experiment design")
    raw_rows = read_csv(raw_path)
    paired_rows = read_csv(paired_path)
    result = validate_pair_table(raw_rows, paired_rows, replicates, bootstrap_resamples, seed)
    return {
        "audit": "PASS",
        "raw_calibration_audit": "PASS",
        "raw_records": result["raw_records"],
        "paired_comparisons": result["paired_comparisons"],
        "paired_truth_replicate_groups": result["paired_groups"],
        "raw_input_sha256": raw_hash,
        "paired_comparisons_sha256": sha256(paired_path),
        "bootstrap_resamples": bootstrap_resamples,
        "bootstrap_seed": seed + 7919,
        "verification": "paired point estimates and deterministic percentile-bootstrap intervals independently reconstructed",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", nargs="?", default="results/known-host-calibration-surface")
    args = parser.parse_args()
    result = audit_directory(Path(args.output_dir))
    out_path = Path(args.output_dir) / "paired_comparison_audit.json"
    out_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
