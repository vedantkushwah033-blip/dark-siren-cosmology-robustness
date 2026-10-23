"""Paired-seed sensitivity experiment using shared mock universes.

This is a separate extension; it does not replace the independent-scenario baseline.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

SCENARIOS = (
    ("complete_uniform", 1.0, "complete", "uniform", "uniform"),
    ("random70_uniform", 0.7, "random", "uniform", "uniform"),
    ("redshift70_uniform", 0.7, "redshift", "uniform", "uniform"),
    ("complete_luminosity_matched", 1.0, "complete", "luminosity", "luminosity"),
    ("complete_luminosity_host_uniform_inference", 1.0, "complete", "luminosity", "uniform"),
)
CONTRASTS = (
    ("random70_vs_complete_uniform", "random70_uniform", "complete_uniform"),
    ("redshift70_vs_random70_uniform", "redshift70_uniform", "random70_uniform"),
    ("uniform_vs_matched_luminosity_inference", "complete_luminosity_host_uniform_inference", "complete_luminosity_matched"),
)
METRICS = ("median", "mean", "posterior_sd", "width_68", "width_95", "pull_median", "pull_mean", "covered_68", "covered_95")


def mean_se(values: list[float]) -> tuple[float, float]:
    if not values:
        raise ValueError("Cannot summarize an empty paired contrast.")
    if not all(math.isfinite(v) for v in values):
        raise ValueError("Contrast values must all be finite.")
    mean = sum(values) / len(values)
    if len(values) < 2:
        return mean, float("nan")
    variance = sum((v - mean) ** 2 for v in values) / (len(values) - 1)
    return mean, math.sqrt(variance / len(values))


def paired_differences(a_rows: list[dict], b_rows: list[dict], metric: str) -> list[float]:
    """Return B - A contrasts only when replicate identities are explicitly aligned."""
    if len(a_rows) != len(b_rows):
        raise ValueError(
            f"Cannot pair scenarios with different row counts: {len(a_rows)} != {len(b_rows)}"
        )
    differences: list[float] = []
    for index, (a, b) in enumerate(zip(a_rows, b_rows)):
        if a.get("run") != b.get("run"):
            raise ValueError(
                f"Pair alignment mismatch at row {index}: "
                f"{a.get('run')!r} != {b.get('run')!r}"
            )
        try:
            value_a = float(a[metric])
            value_b = float(b[metric])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Invalid or missing metric {metric!r} at pair {index}") from exc
        if not math.isfinite(value_a) or not math.isfinite(value_b):
            raise ValueError(f"Non-finite {metric!r} value at pair {index}")
        differences.append(value_b - value_a)
    return differences


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", default="original-project")
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--output-dir", default="results/paired-seed")
    args = parser.parse_args()
    if args.replicates < 2:
        parser.error("--replicates must be at least 2")

    original = Path(args.original_dir).resolve()
    if not (original / "src" / "experiments.py").is_file():
        raise FileNotFoundError(f"Original source tree not found under {original}")
    sys.path.insert(0, str(original))
    from src.experiments import run_repeated

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    by_scenario: dict[str, list[dict]] = {}
    run_rows: list[dict] = []

    for name, completeness, missingness, host_generation, inference_weighting in SCENARIOS:
        # Crucial: every scenario receives the same base seed. In run_repeated,
        # each replicate derives catalogue, host, GW and mask seeds from seed + r.
        # This shares the parent catalogue/host/GW draw when host_generation matches.
        rows = run_repeated(
            n_runs=args.replicates,
            n_galaxies=500,
            H0_true=70.0,
            completeness=completeness,
            missingness=missingness,
            host_generation=host_generation,
            inference_weighting=inference_weighting,
            seed=args.seed,
        )
        if len(rows) != args.replicates:
            raise RuntimeError(f"{name}: expected {args.replicates} rows, got {len(rows)}")
        # Coverage is a per-universe indicator; paired contrasts then estimate
        # the difference in empirical coverage for the same simulated universes.
        for row in rows:
            row["covered_68"] = int(row["lower_68"] <= 70.0 <= row["upper_68"])
            row["covered_95"] = int(row["lower_95"] <= 70.0 <= row["upper_95"])
        by_scenario[name] = rows
        for row in rows:
            run_rows.append({
                "scenario": name,
                "pair_id": row["run"],
                "base_seed": args.seed,
                "target_completeness": completeness,
                "missingness": missingness,
                "host_generation": host_generation,
                "inference_weighting": inference_weighting,
                **row,
            })

    contrast_rows: list[dict] = []
    contrast_summary: list[dict] = []
    for contrast_name, scenario_b, scenario_a in CONTRASTS:
        a_rows, b_rows = by_scenario[scenario_a], by_scenario[scenario_b]
        for metric in METRICS:
            differences = paired_differences(a_rows, b_rows, metric)
            avg, se = mean_se(differences)
            contrast_summary.append({
                "contrast": contrast_name,
                "scenario_b_minus_a": f"{scenario_b} - {scenario_a}",
                "metric": metric,
                "n_pairs": len(differences),
                "mean_difference": avg,
                "paired_difference_se": se,
            })
            for a, b, diff in zip(a_rows, b_rows, differences):
                contrast_rows.append({
                    "contrast": contrast_name,
                    "scenario_a": scenario_a,
                    "scenario_b": scenario_b,
                    "pair_id": a["run"],
                    "metric": metric,
                    "value_a": float(a[metric]),
                    "value_b": float(b[metric]),
                    "difference_b_minus_a": diff,
                })

    def write_csv(path: Path, rows: list[dict]) -> None:
        if not rows:
            raise ValueError(f"No rows to write for {path}")
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    write_csv(out / "paired_per_run_results.csv", run_rows)
    write_csv(out / "paired_contrasts.csv", contrast_rows)
    write_csv(out / "paired_contrast_summary.csv", contrast_summary)
    metadata = {
        "experiment": "paired-seed sensitivity",
        "status": "implementation; not yet executed or scientifically validated",
        "original_commit": "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f",
        "base_seed": args.seed,
        "replicates_per_scenario": args.replicates,
        "n_galaxies": 500,
        "true_h0": 70.0,
        "scenario_count": len(SCENARIOS),
        "contrast_count": len(CONTRASTS),
        "pairing_note": (
            "Scenarios share the same seed and replicate index. The pinned simulator derives "
            "catalogue, host, GW and mask seeds from these values. Host/GW draws are shared only "
            "when host_generation is identical; masks are scenario-specific transformations. "
            "Verify simulator implementation and outputs before interpreting contrasts."
        ),
        "scenarios": [
            {"name": n, "completeness": c, "missingness": m, "host_generation": h, "inference_weighting": w}
            for n, c, m, h, w in SCENARIOS
        ],
        "contrasts": [
            {"name": n, "scenario_b": b, "scenario_a": a}
            for n, b, a in CONTRASTS
        ],
    }
    (out / "paired_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
