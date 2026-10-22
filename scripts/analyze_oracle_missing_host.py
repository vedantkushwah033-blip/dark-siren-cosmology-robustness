"""Bootstrap contrasts for the missing-host oracle-control experiment."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from pathlib import Path
import random

ANALYSIS_PATH = Path(__file__).resolve().parent / "analyze_factorial_contrasts.py"
SPEC = importlib.util.spec_from_file_location("factorial_contrasts", ANALYSIS_PATH)
analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(analysis)

PAIRS = (
    ("random70_observed_only", "random70_missing_term_inference_weights", "Random missingness: inference-weight missing term vs observed-only"),
    ("random70_missing_term_inference_weights", "random70_missing_term_oracle_weights", "Random missingness: oracle vs inference missing-host weights"),
    ("redshift70_observed_only", "redshift70_missing_term_inference_weights", "Redshift-dependent missingness: inference-weight missing term vs observed-only"),
    ("redshift70_missing_term_inference_weights", "redshift70_missing_term_oracle_weights", "Redshift-dependent missingness: oracle vs inference missing-host weights"),
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="results/oracle-missing-host/per_run_results.csv")
    parser.add_argument("--output", default="results/oracle-missing-host/contrast_analysis.json")
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261024)
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

    rng = random.Random(args.seed)
    metrics = ("coverage_68", "coverage_95", "width_68", "width_95",
               "median_error", "squared_median_error", "pull_median")
    result = {
        "input": args.input,
        "true_h0_km_s_mpc": 70.0,
        "bootstrap_replicates": args.bootstrap,
        "bootstrap_seed": args.seed,
        "warning": "Synthetic-model control only. Oracle weighting uses hidden simulated population information unavailable to real surveys. Contrasts are independent-scenario bootstrap comparisons, not paired effects.",
        "contrasts": [],
    }
    for scenario_a, scenario_b, label in PAIRS:
        if scenario_a not in grouped or scenario_b not in grouped:
            raise SystemExit(f"Missing scenario: {scenario_a} or {scenario_b}")
        for metric in metrics:
            a = [analysis.value(row, metric) for row in grouped[scenario_a]]
            b = [analysis.value(row, metric) for row in grouped[scenario_b]]
            result["contrasts"].append({
                "label": label,
                "scenario_A": scenario_a,
                "scenario_B": scenario_b,
                "metric": metric,
                **analysis.bootstrap_mean_difference(a, b, rng, args.bootstrap),
            })
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"Wrote {output} with {len(result['contrasts'])} contrasts.")


if __name__ == "__main__":
    main()
