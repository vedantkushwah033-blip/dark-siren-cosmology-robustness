"""Bootstrap uncertainty for factorial dark-siren scenario contrasts.

Input is per_run_results.csv from run_factorial_sensitivity.py. Each contrast
is scenario B minus scenario A. Intervals quantify Monte Carlo variability
conditional on this synthetic design, not uncertainty in real H0.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from pathlib import Path

TRUE_H0 = 70.0
SCENARIO_CONTRASTS = (
    ("complete_uniform_matched", "complete_luminosity_matched", "Complete catalogue: matched luminosity minus uniform"),
    ("complete_luminosity_matched", "complete_luminosity_host_uniform_inference", "Complete catalogue: mismatched minus matched luminosity"),
    ("random70_uniform_matched", "random70_luminosity_matched", "Random 70% catalogue: matched luminosity minus uniform"),
    ("random70_luminosity_matched", "random70_luminosity_host_uniform_inference", "Random 70% catalogue: mismatched minus matched luminosity"),
    ("redshift70_uniform_matched", "redshift70_luminosity_matched", "Redshift-dependent 70% catalogue: matched luminosity minus uniform"),
    ("redshift70_luminosity_matched", "redshift70_luminosity_host_uniform_inference", "Redshift-dependent 70% catalogue: mismatched minus matched luminosity"),
    ("random70_uniform_matched", "redshift70_uniform_matched", "Uniform weights: redshift-dependent minus random missingness"),
    ("random70_luminosity_matched", "redshift70_luminosity_matched", "Matched luminosity weights: redshift-dependent minus random missingness"),
)


def value(row: dict, metric: str) -> float:
    median = float(row["median"])
    mean = float(row["mean"])
    lower68, upper68 = float(row["lower_68"]), float(row["upper_68"])
    lower95, upper95 = float(row["lower_95"]), float(row["upper_95"])
    if metric == "coverage_68":
        return float(lower68 <= TRUE_H0 <= upper68)
    if metric == "coverage_95":
        return float(lower95 <= TRUE_H0 <= upper95)
    if metric == "width_68":
        return upper68 - lower68
    if metric == "width_95":
        return upper95 - lower95
    if metric == "median_error":
        return median - TRUE_H0
    if metric == "absolute_median_error":
        return abs(median - TRUE_H0)
    if metric == "squared_median_error":
        return (median - TRUE_H0) ** 2
    if metric == "mean_error":
        return mean - TRUE_H0
    if metric == "pull_median":
        pull_raw = row.get("pull_median")
        if pull_raw is not None and str(pull_raw).strip() != "":
            return float(pull_raw)
        return (median - TRUE_H0) / max(float(row["posterior_sd"]), 1e-12)
    raise ValueError(f"Unknown metric: {metric}")


METRICS = ("coverage_68", "coverage_95", "width_68", "width_95",
           "median_error", "squared_median_error", "pull_median")


def coverage_mc_se(coverage: float, n: int) -> float:
    """Binomial Monte Carlo standard error for empirical coverage."""
    if n <= 0:
        raise ValueError("n must be positive")
    if not 0.0 <= coverage <= 1.0:
        raise ValueError("coverage must be between 0 and 1")
    return math.sqrt(coverage * (1.0 - coverage) / n)


def wilson_interval(successes: int, n: int, z: float = 1.959963984540054) -> list[float]:
    """Wilson score interval for a binomial coverage proportion."""
    if n <= 0:
        raise ValueError("n must be positive")
    if successes < 0 or successes > n:
        raise ValueError("successes must be between 0 and n")
    if z <= 0.0:
        raise ValueError("z must be positive")
    p = successes / n
    z2 = z * z
    denom = 1.0 + z2 / n
    center = (p + z2 / (2.0 * n)) / denom
    half = z * math.sqrt(p * (1.0 - p) / n + z2 / (4.0 * n * n)) / denom
    lower = max(0.0, center - half)
    upper = min(1.0, center + half)
    # Avoid tiny floating-point residuals at mathematically exact boundaries.
    if successes == 0:
        lower = 0.0
    if successes == n:
        upper = 1.0
    return [lower, upper]


def percentile(xs: list[float], q: float) -> float:
    xs = sorted(xs)
    p = (len(xs) - 1) * q
    lo = int(math.floor(p))
    hi = int(math.ceil(p))
    return xs[lo] + (xs[hi] - xs[lo]) * (p - lo)


def bootstrap_mean_difference(a: list[float], b: list[float], rng: random.Random,
                              repetitions: int) -> dict:
    observed = sum(b) / len(b) - sum(a) / len(a)
    diffs = []
    for _ in range(repetitions):
        aa = [a[rng.randrange(len(a))] for _ in a]
        bb = [b[rng.randrange(len(b))] for _ in b]
        diffs.append(sum(bb) / len(bb) - sum(aa) / len(aa))
    return {
        "difference_B_minus_A": observed,
        "bootstrap_95_interval": [percentile(diffs, 0.025), percentile(diffs, 0.975)],
        "n_A": len(a), "n_B": len(b),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="results/factorial-sensitivity/per_run_results.csv")
    parser.add_argument("--output", default="results/factorial-sensitivity/contrast_analysis.json")
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261023)
    args = parser.parse_args()
    if args.bootstrap < 100:
        parser.error("--bootstrap must be at least 100")

    with Path(args.input).open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit("Input CSV has no rows")
    grouped: dict[str, list[dict]] = {}
    for row in rows:
        grouped.setdefault(row["scenario"], []).append(row)

    expected_scenarios = {
        scenario
        for scenario_a, scenario_b, _ in SCENARIO_CONTRASTS
        for scenario in (scenario_a, scenario_b)
    }
    actual_scenarios = set(grouped)
    if actual_scenarios != expected_scenarios:
        missing = sorted(expected_scenarios - actual_scenarios)
        unexpected = sorted(actual_scenarios - expected_scenarios)
        raise SystemExit(
            f"Scenario set mismatch; missing={missing}, unexpected={unexpected}"
        )
    replicate_counts = {name: len(group) for name, group in grouped.items()}
    if len(set(replicate_counts.values())) != 1:
        raise SystemExit(f"Inconsistent replicate counts across scenarios: {replicate_counts}")
    for name, group in grouped.items():
        for index, row in enumerate(group):
            try:
                lo68, hi68 = float(row["lower_68"]), float(row["upper_68"])
                lo95, hi95 = float(row["lower_95"]), float(row["upper_95"])
            except (KeyError, TypeError, ValueError) as exc:
                raise SystemExit(
                    f"Invalid credible-interval fields for scenario {name}, row {index}"
                ) from exc
            if not all(math.isfinite(x) for x in (lo68, hi68, lo95, hi95)):
                raise SystemExit(
                    f"Non-finite credible interval for scenario {name}, row {index}"
                )
            if lo68 > hi68 or lo95 > hi95:
                raise SystemExit(
                    f"Reversed credible interval for scenario {name}, row {index}"
                )
            try:
                median = float(row["median"])
                mean = float(row["mean"])
                posterior_sd = float(row["posterior_sd"])
            except (KeyError, TypeError, ValueError) as exc:
                raise SystemExit(
                    f"Invalid posterior summary fields for scenario {name}, row {index}"
                ) from exc
            if not all(math.isfinite(x) for x in (median, mean, posterior_sd)):
                raise SystemExit(
                    f"Non-finite posterior summary for scenario {name}, row {index}"
                )
            if posterior_sd <= 0.0:
                raise SystemExit(
                    f"Non-positive posterior_sd for scenario {name}, row {index}"
                )
            if "pull_median" in row and str(row["pull_median"]).strip() != "":
                try:
                    supplied_pull = float(row["pull_median"])
                except (TypeError, ValueError) as exc:
                    raise SystemExit(
                        f"Invalid pull_median for scenario {name}, row {index}"
                    ) from exc
                if not math.isfinite(supplied_pull):
                    raise SystemExit(
                        f"Non-finite pull_median for scenario {name}, row {index}"
                    )

    rng = random.Random(args.seed)
    result = {
        "input": args.input,
        "true_h0_km_s_mpc": TRUE_H0,
        "bootstrap_replicates": args.bootstrap,
        "bootstrap_seed": args.seed,
        "definition": "Each difference is scenario B minus scenario A. Bootstrap intervals quantify Monte Carlo uncertainty under this synthetic design, not real-world H0 uncertainty.",
        "contrasts": [],
        "scenario_metrics": [],
    }

    for name, group in sorted(grouped.items()):
        if len(group) < 2:
            raise SystemExit(f"Scenario {name} has fewer than two rows")
        summary = {"scenario": name, "n": len(group)}
        for metric in METRICS:
            vals = [value(row, metric) for row in group]
            summary[metric] = sum(vals) / len(vals)
        # Binomial Monte Carlo standard errors quantify finite-replicate
        # noise in empirical credible-interval coverage (not real-H0 uncertainty).
        for metric in ("coverage_68", "coverage_95"):
            p = summary[metric]
            successes = round(p * len(group))
            summary[f"{metric}_mc_se"] = coverage_mc_se(p, len(group))
            summary[f"{metric}_wilson_95_interval"] = wilson_interval(successes, len(group))
        summary["rmse_median"] = math.sqrt(summary["squared_median_error"])
        result["scenario_metrics"].append(summary)

    for scenario_a, scenario_b, label in SCENARIO_CONTRASTS:
        if scenario_a not in grouped or scenario_b not in grouped:
            raise SystemExit(f"Missing scenario: {scenario_a} or {scenario_b}")
        for metric in METRICS:
            a = [value(row, metric) for row in grouped[scenario_a]]
            b = [value(row, metric) for row in grouped[scenario_b]]
            contrast = bootstrap_mean_difference(a, b, rng, args.bootstrap)
            result["contrasts"].append({
                "label": label, "scenario_A": scenario_a,
                "scenario_B": scenario_b, "metric": metric, **contrast,
            })

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"Wrote {output} with {len(result['contrasts'])} contrasts.")
    print("Interpretation: synthetic-model sensitivity only; inspect coverage, bias, RMSE and widths together.")


if __name__ == "__main__":
    main()
