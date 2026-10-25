"""Paired replicate-level comparisons across distance-noise levels.

Requires a passing raw-record audit. Differences are defined as (higher-noise
metric minus the 5% reference) within the same truth and replicate ID.
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
REFERENCE = 0.05
METRICS = {
    "absolute_median_error": lambda r, t: abs(float(r["median"]) - t),
    "width_68": lambda r, t: float(r["upper_68"]) - float(r["lower_68"]),
    "width_95": lambda r, t: float(r["upper_95"]) - float(r["lower_95"]),
    "coverage_68_indicator": lambda r, t: float(float(r["lower_68"]) <= t <= float(r["upper_68"])),
    "coverage_95_indicator": lambda r, t: float(float(r["lower_95"]) <= t <= float(r["upper_95"])),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def verify_audit(out: Path, audit: dict) -> None:
    if audit.get("audit") != "PASS":
        raise ValueError("paired analysis requires audit_report.json with audit=PASS")
    if audit.get("posterior_summaries_independently_reconstructed_from_raw_observations") is not True:
        raise ValueError("audit did not independently reconstruct posterior summaries")
    hashes = audit.get("input_sha256", {})
    for name in ("per_run_results.csv", "cell_summary.csv", "cell_summary.json"):
        expected = hashes.get(name)
        if not isinstance(expected, str) or len(expected) != 64:
            raise ValueError(f"audit is missing SHA-256 for {name}")
        actual = hashlib.sha256((out / name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f"{name} changed after audit")


def percentile_ci(values: np.ndarray, rng: np.random.Generator, n_boot: int) -> tuple[float, float]:
    n = len(values)
    draws = rng.integers(0, n, size=(n_boot, n))
    means = values[draws].mean(axis=1)
    low, high = np.quantile(means, [0.025, 0.975])
    return float(low), float(high)


def analyze(rows: list[dict[str, str]], replicates: int, n_boot: int, seed: int) -> list[dict]:
    keyed: dict[tuple[float, int], dict[float, dict[str, str]]] = {}
    for i, row in enumerate(rows):
        truth, noise = float(row["true_h0"]), float(row["fractional_distance_error"])
        run, row_seed = int(row["run"]), int(row["seed"])
        if truth not in TRUTHS or noise not in NOISES:
            raise ValueError(f"unexpected truth/noise at raw row {i}")
        if not 0 <= run < replicates or row_seed != seed + run:
            raise ValueError(f"unexpected run or seed at raw row {i}")
        cell = keyed.setdefault((truth, run), {})
        if noise in cell:
            raise ValueError(f"duplicate truth/run/noise record: {(truth, run, noise)}")
        cell[noise] = row

    expected_keys = {(t, r) for t in TRUTHS for r in range(replicates)}
    if set(keyed) != expected_keys or any(set(v) != set(NOISES) for v in keyed.values()):
        raise ValueError("raw records do not form a complete paired 9-truth × replicate × 4-noise design")

    rng = np.random.default_rng(seed + 7919)
    output = []
    for truth in TRUTHS:
        for noise in NOISES[1:]:
            for metric, get_value in METRICS.items():
                diffs = []
                for run in range(replicates):
                    block = keyed[(truth, run)]
                    ref_value = get_value(block[REFERENCE], truth)
                    candidate_value = get_value(block[noise], truth)
                    diffs.append(candidate_value - ref_value)
                    # Catalogue/host quantities should remain fixed within paired blocks.
                    for field in ("host_redshift", "true_distance_mpc"):
                        if not math.isclose(float(block[REFERENCE][field]), float(block[noise][field]),
                                            rel_tol=1e-12, abs_tol=1e-12):
                            raise ValueError(f"pairing mismatch in {field}, truth={truth}, run={run}")
                arr = np.asarray(diffs, dtype=float)
                lo, hi = percentile_ci(arr, rng, n_boot)
                output.append({
                    "true_h0": truth,
                    "reference_fractional_error": REFERENCE,
                    "comparison_fractional_error": noise,
                    "metric": metric,
                    "n_paired_replicates": replicates,
                    "mean_difference": float(arr.mean()),
                    "median_difference": float(np.median(arr)),
                    "bootstrap_mean_difference_ci95_low": lo,
                    "bootstrap_mean_difference_ci95_high": hi,
                    "bootstrap_resamples": n_boot,
                    "bootstrap_seed": seed + 7919,
                    "difference_definition": "comparison minus 5% reference",
                })
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", nargs="?", default="results/known-host-calibration-surface")
    parser.add_argument("--replicates", type=int, default=300)
    parser.add_argument("--bootstrap-resamples", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=20261025)
    args = parser.parse_args()
    if args.replicates < 100 or args.bootstrap_resamples < 1000 or args.seed < 0:
        parser.error("replicates must be >=100, bootstrap-resamples >=1000, and seed >=0")
    out = Path(args.output_dir)
    audit = json.loads((out / "audit_report.json").read_text(encoding="utf-8"))
    verify_audit(out, audit)
    metadata = json.loads((out / "cell_summary.json").read_text(encoding="utf-8"))
    if metadata.get("replicates_per_cell") != args.replicates or metadata.get("seed") != args.seed:
        raise ValueError("CLI replicate count/seed must match audited experiment metadata")
    rows = read_csv(out / "per_run_results.csv")
    expected = len(TRUTHS) * len(NOISES) * args.replicates
    if len(rows) != expected or audit.get("raw_records") != expected:
        raise ValueError(f"expected {expected} audited raw rows")
    results = analyze(rows, args.replicates, args.bootstrap_resamples, args.seed)
    fields = list(results[0])
    with (out / "paired_comparisons.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)
    report = {
        "analysis": "paired replicate-level noise comparisons",
        "audit": "PASS",
        "source_audit_run_id": audit.get("workflow_run_id"),
        "source_simulator_commit": metadata["original_commit"],
        "base_seed": args.seed,
        "replicates_per_truth": args.replicates,
        "bootstrap_resamples": args.bootstrap_resamples,
        "bootstrap_seed": args.seed + 7919,
        "pairing_unit": "replicate ID within fixed true H0; all four noise-level records retained as a block",
        "reference_fractional_error": REFERENCE,
        "difference_definition": "comparison minus reference; negative absolute-error difference favours the comparison level, positive width difference means wider intervals, positive coverage-indicator difference means inclusion changed from 0 to 1",
        "metrics": list(METRICS),
        "comparison_count": len(results),
        "raw_input_sha256": audit["input_sha256"]["per_run_results.csv"],
        "paired_comparisons_sha256": hashlib.sha256((out / "paired_comparisons.csv").read_bytes()).hexdigest(),
        "interpretation": "Descriptive paired comparisons; percentile bootstrap intervals quantify Monte Carlo uncertainty in mean paired differences. They are not a model-validity test and no multiplicity-adjusted confirmatory claims are made.",
    }
    (out / "paired_comparisons.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    lines = [
        "# Paired noise-level comparisons",
        "",
        "All differences are **comparison level minus the 5% reference**, paired by replicate ID within each fixed true-H₀ stratum. Negative absolute-error differences favour the comparison level; positive width differences mean wider intervals. Coverage-indicator differences are paired changes in truth inclusion.",
        "",
        f"- Replicates per truth: {args.replicates}",
        f"- Percentile bootstrap resamples: {args.bootstrap_resamples}; seed: {args.seed + 7919}",
        "- Bootstrap resampling unit: replicate ID; resampling retains all four noise-level records as a block.",
        "- These are descriptive comparisons, without multiplicity-adjusted confirmatory claims.",
        "",
        "| True H₀ | Comparison error | Metric | Mean paired difference | Median paired difference | Bootstrap 95% CI for mean |",
        "|---:|---:|---|---:|---:|---:|",
    ]
    for r in results:
        lines.append(f"| {r['true_h0']:.0f} | {100*r['comparison_fractional_error']:.0f}% | {r['metric']} | {r['mean_difference']:+.6g} | {r['median_difference']:+.6g} | [{r['bootstrap_mean_difference_ci95_low']:+.6g}, {r['bootstrap_mean_difference_ci95_high']:+.6g}] |")
    lines.extend(["", "The two prior-boundary truths (50 and 90) are retained but must be interpreted separately because equal-tailed intervals on a continuous posterior supported on [50, 90] do not include the exact boundary. These comparisons do not validate realistic dark-siren inference.", ""])
    (out / "paired_comparisons_report.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"comparison_rows": len(results), "paired_replicates_per_truth": args.replicates,
                      "bootstrap_resamples": args.bootstrap_resamples}, indent=2))


if __name__ == "__main__":
    main()
