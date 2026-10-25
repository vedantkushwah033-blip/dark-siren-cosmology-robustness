"""Paired replicate-level comparisons across distance-noise levels.

Run only after the independent raw-record audit. Differences are each noise
level minus the 5% reference within the same truth/replicate block.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import statistics

TRUTHS = (50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0)
NOISES = (0.05, 0.10, 0.15, 0.20)
REFERENCE = 0.05
METRICS = {
    "absolute_median_error": lambda r, t: abs(float(r["median"]) - t),
    "width_68": lambda r, t: float(r["upper_68"]) - float(r["lower_68"]),
    "width_95": lambda r, t: float(r["upper_95"]) - float(r["lower_95"]),
    "covers_68": lambda r, t: float(r["lower_68"]) <= t <= float(r["upper_68"]),
    "covers_95": lambda r, t: float(r["lower_95"]) <= t <= float(r["upper_95"]),
}


def percentile(values: list[float], q: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("cannot take percentile of empty values")
    pos = (len(ordered) - 1) * q
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (pos - lo)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_and_validate(out: Path, expected_replicates: int | None) -> tuple[list[dict[str, str]], dict]:
    audit_path = out / "audit_report.json"
    metadata_path = out / "cell_summary.json"
    raw_path = out / "per_run_results.csv"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if audit.get("audit") != "PASS":
        raise ValueError("paired analysis requires an independent audit PASS")
    if audit.get("posterior_summaries_independently_reconstructed_from_raw_observations") is not True:
        raise ValueError("paired analysis requires independently reconstructed posterior summaries")
    actual_hash = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    expected_hash = audit.get("input_sha256", {}).get("per_run_results.csv")
    if not isinstance(expected_hash, str) or actual_hash != expected_hash:
        raise ValueError("raw records changed since independent audit")
    if metadata.get("true_h0_values") != list(TRUTHS) or metadata.get("fractional_distance_errors") != list(NOISES):
        raise ValueError("metadata design does not match the prespecified truth/noise grid")
    if metadata.get("original_commit") != "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f":
        raise ValueError("unexpected source simulator revision")
    n = metadata.get("replicates_per_cell")
    if isinstance(n, bool) or not isinstance(n, int) or n < 1:
        raise ValueError("invalid replicate count in metadata")
    if expected_replicates is not None and n != expected_replicates:
        raise ValueError(f"metadata replicates {n} do not match requested {expected_replicates}")
    rows = read_csv(raw_path)
    if len(rows) != 9 * 4 * n or len(rows) != audit.get("raw_records"):
        raise ValueError("raw row count disagrees with the audited design")
    seen: dict[tuple[float, int, float], dict[str, str]] = {}
    for row in rows:
        truth = float(row["true_h0"])
        noise = float(row["fractional_distance_error"])
        run = int(row["run"])
        if truth not in TRUTHS or noise not in NOISES or not 0 <= run < n:
            raise ValueError("raw record has an unexpected truth, noise level, or replicate ID")
        key = (truth, run, noise)
        if key in seen:
            raise ValueError(f"duplicate truth/replicate/noise record: {key}")
        seen[key] = row
        for field in ("median", "lower_68", "upper_68", "lower_95", "upper_95"):
            if not math.isfinite(float(row[field])):
                raise ValueError(f"non-finite {field} in raw record {key}")
    expected_keys = {(t, run, noise) for t in TRUTHS for run in range(n) for noise in NOISES}
    if set(seen) != expected_keys:
        raise ValueError("raw records do not form complete four-level replicate blocks")
    # The common-random-number design must preserve host data and standardized
    # noise across uncertainty levels within each truth/replicate.
    for truth in TRUTHS:
        for run in range(n):
            block = [seen[(truth, run, noise)] for noise in NOISES]
            for field in ("seed", "host_redshift", "true_distance_mpc"):
                if len({r[field] for r in block}) != 1:
                    raise ValueError(f"pairing mismatch in {field} for truth={truth}, replicate={run}")
            standardized = [
                (float(r["observed_distance_mpc"]) - float(r["true_distance_mpc"])) / float(r["distance_sigma_mpc"])
                for r in block
            ]
            if max(standardized) - min(standardized) > 1e-8:
                raise ValueError(f"standardized event noise mismatch for truth={truth}, replicate={run}")
    return rows, metadata


def analyze(rows: list[dict[str, str]], replicates: int, bootstrap_replicates: int, seed: int) -> list[dict]:
    by_key = {(float(r["true_h0"]), int(r["run"]), float(r["fractional_distance_error"])): r for r in rows}
    rng = random.Random(seed)
    output = []
    # One block bootstrap index vector per truth/resample is reused across all
    # comparisons and metrics, retaining the complete four-level pairing.
    bootstrap_indices = {
        truth: [[rng.randrange(replicates) for _ in range(replicates)] for _ in range(bootstrap_replicates)]
        for truth in TRUTHS
    }
    for truth in TRUTHS:
        for noise in NOISES[1:]:
            for metric_name, metric_fn in METRICS.items():
                diffs = [
                    float(metric_fn(by_key[(truth, run, noise)], truth))
                    - float(metric_fn(by_key[(truth, run, REFERENCE)], truth))
                    for run in range(replicates)
                ]
                boot_means = [
                    statistics.fmean(diffs[i] for i in indices)
                    for indices in bootstrap_indices[truth]
                ]
                boot_medians = [
                    statistics.median(diffs[i] for i in indices)
                    for indices in bootstrap_indices[truth]
                ]
                output.append({
                    "true_h0": truth,
                    "reference_fractional_error": REFERENCE,
                    "comparison_fractional_error": noise,
                    "metric": metric_name,
                    "difference_definition": "comparison minus 5% reference; positive means larger metric at comparison level",
                    "replicates": replicates,
                    "bootstrap_replicates": bootstrap_replicates,
                    "bootstrap_seed": seed,
                    "mean_paired_difference": statistics.fmean(diffs),
                    "median_paired_difference": statistics.median(diffs),
                    "mean_bootstrap_ci95_low": percentile(boot_means, 0.025),
                    "mean_bootstrap_ci95_high": percentile(boot_means, 0.975),
                    "median_bootstrap_ci95_low": percentile(boot_medians, 0.025),
                    "median_bootstrap_ci95_high": percentile(boot_medians, 0.975),
                })
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", nargs="?", default="results/known-host-calibration-surface")
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261026)
    args = parser.parse_args()
    if args.bootstrap_replicates < 200:
        parser.error("--bootstrap-replicates must be at least 200")
    if args.seed < 0:
        parser.error("--seed must be non-negative")
    out = Path(args.output_dir)
    rows, metadata = load_and_validate(out, args.replicates)
    result = analyze(rows, int(metadata["replicates_per_cell"]), args.bootstrap_replicates, args.seed)
    path = out / "paired_comparison_summary.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result[0]))
        writer.writeheader()
        writer.writerows(result)
    report = {
        "analysis": "paired replicate-level differences across distance uncertainty",
        "status": "PASS",
        "truth_count": len(TRUTHS),
        "comparison_levels": list(NOISES[1:]),
        "metrics": list(METRICS),
        "rows": len(result),
        "replicates_per_truth": int(metadata["replicates_per_cell"]),
        "bootstrap_replicates": args.bootstrap_replicates,
        "bootstrap_seed": args.seed,
        "reference_fractional_error": REFERENCE,
        "raw_records_sha256": hashlib.sha256((out / "per_run_results.csv").read_bytes()).hexdigest(),
        "interpretation": "Descriptive paired differences under this synthetic matched model; not causal or observational inference. Coverage differences are paired binary-outcome differences.",
    }
    (out / "paired_comparison_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
