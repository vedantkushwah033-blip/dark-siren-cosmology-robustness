"""Independent integrity checks for prior × selection-threshold factorial outputs.

Checks the raw simulation table, derived summaries, and paired contrasts as
separate artifacts. This is a synthetic-model validation utility, not an
astrophysical validity check.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import math
from collections import defaultdict
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve().with_name("run_prior_threshold_factorial.py")
SPEC = importlib.util.spec_from_file_location("factorial_module", SCRIPT)
factorial = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(factorial)

SUMMARY_KEYS = ("posterior_median", "lower_68", "upper_68", "lower_95", "upper_95")
RAW_KEY = ("true_H0", "replicate", "prior_min", "prior_max", "assumed_detection_max_mpc", "method")
SUMMARY_KEY = ("true_H0", "prior_min", "prior_max", "assumed_detection_max_mpc", "method")
CONTRAST_KEY = ("true_H0", "method", "metric", "comparison_prior", "comparison_threshold_mpc")


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"Missing CSV header: {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty CSV: {path}")
    for line, row in enumerate(rows, 2):
        for key, value in row.items():
            if value is None:
                raise ValueError(f"Missing field {key} at {path}:{line}")
            try:
                number = float(value)
            except (ValueError, TypeError):
                continue
            if not math.isfinite(number):
                raise ValueError(f"Non-finite numeric field {key} at {path}:{line}")
    return rows


def as_float(row, key):
    return float(row[key])


def assert_close(actual, expected, label, atol=1e-9):
    if not math.isclose(float(actual), float(expected), rel_tol=1e-10, abs_tol=atol):
        raise ValueError(f"{label}: expected {expected:.12g}, got {actual:.12g}")


def validate(per_run_path: Path, summary_path: Path, contrasts_path: Path):
    raw = read_csv(per_run_path)
    summaries = read_csv(summary_path)
    contrasts = read_csv(contrasts_path)
    reps_by_truth = {}
    for truth in factorial.TRUTHS:
        values = {int(r["replicate"]) for r in raw if float(r["true_H0"]) == truth}
        if not values or values != set(range(max(values) + 1)):
            raise ValueError(f"Replicate indices are incomplete for H0={truth}")
        reps_by_truth[truth] = len(values)
    if len(set(reps_by_truth.values())) != 1:
        raise ValueError(f"Different replicate counts by truth: {reps_by_truth}")
    reps = next(iter(reps_by_truth.values()))
    expected_raw = len(factorial.TRUTHS) * reps * len(factorial.PRIORS) * len(factorial.ASSUMED_MAXIMA) * len(factorial.METHODS)
    expected_summaries = len(factorial.TRUTHS) * len(factorial.PRIORS) * len(factorial.ASSUMED_MAXIMA) * len(factorial.METHODS)
    expected_contrasts = len(factorial.TRUTHS) * (len(factorial.PRIORS) * len(factorial.ASSUMED_MAXIMA) - 1) * len(factorial.METHODS) * 6
    if len(raw) != expected_raw:
        raise ValueError(f"Raw row count: expected {expected_raw}, got {len(raw)}")
    if len(summaries) != expected_summaries:
        raise ValueError(f"Summary row count: expected {expected_summaries}, got {len(summaries)}")
    if len(contrasts) != expected_contrasts:
        raise ValueError(f"Contrast row count: expected {expected_contrasts}, got {len(contrasts)}")

    raw_lookup = {}
    observations = defaultdict(set)
    for row in raw:
        key = (float(row["true_H0"]), int(row["replicate"]), float(row["prior_min"]),
               float(row["prior_max"]), float(row["assumed_detection_max_mpc"]), row["method"])
        if key in raw_lookup:
            raise ValueError(f"Duplicate raw configuration: {key}")
        raw_lookup[key] = row
        sim_key = (float(row["true_H0"]), int(row["replicate"]))
        observations[sim_key].add((int(row["host_index"]), float(row["observed_distance_mpc"])))
        if float(row["true_detection_max_mpc"]) != factorial.TRUE_MAX:
            raise ValueError("Unexpected data-generating detection threshold")
        interval = [float(row[k]) for k in ("lower_95", "lower_68", "posterior_median", "upper_68", "upper_95")]
        if interval != sorted(interval):
            raise ValueError(f"Credible intervals are not ordered: {key}")
    if any(len(v) != 1 for v in observations.values()):
        raise ValueError("The same replicate does not retain a single host and observation across configurations")

    summary_lookup = {}
    grouped = defaultdict(list)
    for row in raw:
        key = (float(row["true_H0"]), float(row["prior_min"]), float(row["prior_max"]),
               float(row["assumed_detection_max_mpc"]), row["method"])
        grouped[key].append(row)
    for row in summaries:
        key = (float(row["true_H0"]), float(row["prior_min"]), float(row["prior_max"]),
               float(row["assumed_detection_max_mpc"]), row["method"])
        if key in summary_lookup:
            raise ValueError(f"Duplicate summary configuration: {key}")
        summary_lookup[key] = row
        cell = grouped.get(key)
        if cell is None or len(cell) != reps:
            raise ValueError(f"Summary lacks complete raw replicate group: {key}")
        truth = key[0]
        errors = np.array([float(x["posterior_median"]) - truth for x in cell])
        coverage68 = np.array([float(x["lower_68"]) <= truth <= float(x["upper_68"]) for x in cell])
        coverage95 = np.array([float(x["lower_95"]) <= truth <= float(x["upper_95"]) for x in cell])
        checks = {
            "bias": float(errors.mean()),
            "mae": float(np.abs(errors).mean()),
            "rmse": float(np.sqrt(np.mean(errors**2))),
            "coverage_68": float(coverage68.mean()),
            "coverage_95": float(coverage95.mean()),
            "mean_width_68": float(np.mean([float(x["upper_68"])-float(x["lower_68"]) for x in cell])),
            "mean_width_95": float(np.mean([float(x["upper_95"])-float(x["lower_95"]) for x in cell])),
        }
        for metric, expected in checks.items():
            assert_close(row[metric], expected, f"{key} {metric}")
        for level, success in ((68, int(coverage68.sum())), (95, int(coverage95.sum()))):
            low, high = factorial.wilson(success, reps)
            assert_close(row[f"coverage_{level}_wilson_low"], low, f"{key} Wilson {level}% low")
            assert_close(row[f"coverage_{level}_wilson_high"], high, f"{key} Wilson {level}% high")

    if len(summary_lookup) != expected_summaries:
        raise ValueError("Summary configuration set is incomplete")

    seen_contrasts = set()
    for row in contrasts:
        key = (float(row["true_H0"]), row["method"], row["metric"],
               row["comparison_prior"], float(row["comparison_threshold_mpc"]))
        if key in seen_contrasts:
            raise ValueError(f"Duplicate contrast row: {key}")
        seen_contrasts.add(key)
        truth, method, metric, prior_text, threshold = key
        try:
            pmin, pmax = (float(x) for x in prior_text.strip("[]").split(","))
        except Exception as exc:
            raise ValueError(f"Malformed comparison_prior {prior_text}") from exc
        if metric not in {"median", "absolute_error", "width_68", "width_95", "coverage_68", "coverage_95"}:
            raise ValueError(f"Unknown contrast metric: {metric}")
        baseline = [raw_lookup[(truth, rep, 50.0, 90.0, 120.0, method)] for rep in range(reps)]
        comparison = [raw_lookup[(truth, rep, pmin, pmax, threshold, method)] for rep in range(reps)]
        def value(item):
            if metric == "median": return float(item["posterior_median"])
            if metric == "absolute_error": return abs(float(item["posterior_median"]) - truth)
            if metric == "width_68": return float(item["upper_68"]) - float(item["lower_68"])
            if metric == "width_95": return float(item["upper_95"]) - float(item["lower_95"])
            if metric == "coverage_68": return float(float(item["lower_68"]) <= truth <= float(item["upper_68"]))
            return float(float(item["lower_95"]) <= truth <= float(item["upper_95"]))
        differences = np.array([value(b) - value(a) for a, b in zip(baseline, comparison)])
        if int(row["n_pairs"]) != reps:
            raise ValueError(f"Wrong n_pairs for contrast {key}")
        assert_close(row["mean_difference_vs_baseline"], differences.mean(), f"{key} paired mean")
        if float(row["bootstrap_95_low"]) > float(row["bootstrap_95_high"]):
            raise ValueError(f"Reversed bootstrap interval for {key}")
        if int(row["bootstrap_draws"]) < 100:
            raise ValueError(f"Insufficient bootstrap draws for {key}")
    if len(seen_contrasts) != expected_contrasts:
        raise ValueError("Paired contrast configuration set is incomplete")

    return {
        "status": "PASS_FACTORIAL_ARTIFACT_INTEGRITY",
        "replicates_per_truth": reps,
        "raw_rows": len(raw),
        "summary_rows": len(summaries),
        "paired_contrast_rows": len(contrasts),
        "unique_observations": len(observations),
        "summary_metrics_recomputed": 7 * len(summaries),
        "paired_mean_contrasts_recomputed": len(contrasts),
        "limitations": [
            "Checks numerical consistency of exported artifacts; does not validate scientific model realism.",
            "Bootstrap interval endpoints are checked for ordering, while the bootstrap itself is reproducible from the run seed recorded in each row.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-run", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--contrasts", required=True)
    args = parser.parse_args()
    result = validate(Path(args.per_run), Path(args.summary), Path(args.contrasts))
    import json
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
