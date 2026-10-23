"""Paired mock-universe sensitivity experiment.

Each replicate reuses one parent galaxy catalogue. Within a host-generation
setting, host selection and GW noise are shared across catalogue-selection
scenarios. Inference-weight comparisons also share the same event and mask.
This is synthetic-model analysis, not LVK observational inference.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

SCENARIOS = (
    ("complete_uniform_uniform", "complete", 1.0, "uniform", "uniform"),
    ("random70_uniform_uniform", "random", 0.7, "uniform", "uniform"),
    ("redshift70_uniform_uniform", "redshift", 0.7, "uniform", "uniform"),
    ("complete_luminosity_luminosity", "complete", 1.0, "luminosity", "luminosity"),
    ("complete_luminosity_uniform", "complete", 1.0, "luminosity", "uniform"),
    ("random70_luminosity_luminosity", "random", 0.7, "luminosity", "luminosity"),
    ("random70_luminosity_uniform", "random", 0.7, "luminosity", "uniform"),
    ("redshift70_luminosity_luminosity", "redshift", 0.7, "luminosity", "luminosity"),
    ("redshift70_luminosity_uniform", "redshift", 0.7, "luminosity", "uniform"),
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", default="original-project")
    parser.add_argument("--replicates", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261023)
    parser.add_argument("--output-dir", default="results/paired-factorial-sensitivity")
    args = parser.parse_args()
    if args.replicates < 2:
        parser.error("--replicates must be at least 2")

    original = Path(args.original_dir).resolve()
    if not (original / "src" / "experiments.py").is_file():
        raise FileNotFoundError(f"Original source tree not found under {original}")
    sys.path.insert(0, str(original))
    import numpy as np
    from src.simulations_v2 import simulate_galaxy_catalogue, choose_host, simulate_gw_event
    from src.experiments import _mask, _weights
    from src.bayesian_inference import posterior_h0, posterior_summary
    from src.metrics_v2 import summarize_runs

    grid = np.linspace(50.0, 90.0, 1601)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    per_run: list[dict] = []

    for r in range(args.replicates):
        universe_seed = args.seed + r * 1_000_000
        cat = simulate_galaxy_catalogue(n_galaxies=500, H0_true=70.0, seed=universe_seed)
        for host_generation in ("uniform", "luminosity"):
            host = choose_host(cat, weighting=host_generation, seed=universe_seed + 100_000)
            gw = simulate_gw_event(
                cat["true_distance_mpc"][host],
                seed=universe_seed + 200_000,
            )
            for scenario_name, missingness, completeness, gen_weight, infer_weight in SCENARIOS:
                if gen_weight != host_generation:
                    continue
                mask = _mask(cat, completeness, missingness, universe_seed + 300_000 + {
                    "complete": 0, "random": 1, "redshift": 2
                }[missingness])
                if not mask.any():
                    raise RuntimeError(f"Empty observed catalogue in {scenario_name}, replicate {r}")
                weights = _weights(cat, infer_weight)[mask]
                posterior = posterior_h0(
                    grid,
                    cat["redshift"][mask],
                    gw["observed_distance_mpc"],
                    gw["sigma_mpc"],
                    weights=weights,
                )
                summary = posterior_summary(grid, posterior)
                row = {
                    "universe_id": r,
                    "universe_seed": universe_seed,
                    "scenario": scenario_name,
                    "run": r,
                    "H0_true": 70.0,
                    "host_generation": gen_weight,
                    "inference_weighting": infer_weight,
                    "missingness": missingness,
                    "target_completeness": completeness,
                    "realized_completeness": float(mask.mean()),
                    "host_observed": bool(mask[host]),
                    "host_redshift": float(cat["redshift"][host]),
                    "observed_distance_mpc": float(gw["observed_distance_mpc"]),
                    "sigma_mpc": float(gw["sigma_mpc"]),
                    "median": float(summary["median"]),
                    "mean": float(summary["mean"]),
                    "posterior_sd": float(summary["posterior_sd"]),
                    "lower_68": float(summary["lower_68"]),
                    "upper_68": float(summary["upper_68"]),
                    "lower_95": float(summary["lower_95"]),
                    "upper_95": float(summary["upper_95"]),
                    "width_68": float(summary["upper_68"] - summary["lower_68"]),
                    "width_95": float(summary["upper_95"] - summary["lower_95"]),
                    "median_error": float(summary["median"] - 70.0),
                    "mean_error": float(summary["mean"] - 70.0),
                }
                if not all(math.isfinite(float(row[k])) for k in (
                    "median", "mean", "posterior_sd", "lower_68", "upper_68",
                    "lower_95", "upper_95", "realized_completeness"
                )):
                    raise RuntimeError(f"Non-finite result in {scenario_name}, replicate {r}")
                per_run.append(row)

    scenario_rows = []
    for scenario_name, missingness, completeness, gen_weight, infer_weight in SCENARIOS:
        group = [row for row in per_run if row["scenario"] == scenario_name]
        if len(group) != args.replicates:
            raise RuntimeError(f"Wrong replicate count for {scenario_name}: {len(group)}")
        # Summarizer uses the same public row schema as the original simulator.
        metric_rows = [{
            **row,
            "pull_median": row["median_error"] / max(row["posterior_sd"], 1e-12),
            "pull_mean": row["mean_error"] / max(row["posterior_sd"], 1e-12),
        } for row in group]
        scenario_rows.append({
            "scenario": scenario_name,
            "missingness": missingness,
            "target_completeness": completeness,
            "host_generation": gen_weight,
            "inference_weighting": infer_weight,
            **summarize_runs(metric_rows, true_h0=70.0),
            "median_bias": float(np.median([row["median_error"] for row in group])),
            "mean_median_error": float(np.mean([row["median_error"] for row in group])),
        })

    # Contrasts are paired by universe_id. Only compare scenarios that share
    # the same host-generation model; weighting contrasts additionally share
    # the same mask. B-minus-A signs are stored explicitly.
    contrast_specs = (
        ("complete_uniform_uniform", "random70_uniform_uniform", "Uniform hosts: random 70% minus complete"),
        ("complete_uniform_uniform", "redshift70_uniform_uniform", "Uniform hosts: redshift-dependent 70% minus complete"),
        ("random70_uniform_uniform", "redshift70_uniform_uniform", "Uniform hosts: redshift-dependent minus random 70%"),
        ("complete_luminosity_luminosity", "complete_luminosity_uniform", "Complete: uniform inference minus matched luminosity inference"),
        ("random70_luminosity_luminosity", "random70_luminosity_uniform", "Random 70%: uniform inference minus matched luminosity inference"),
        ("redshift70_luminosity_luminosity", "redshift70_luminosity_uniform", "Redshift 70%: uniform inference minus matched luminosity inference"),
        ("random70_luminosity_uniform", "redshift70_luminosity_uniform", "Luminosity-generated hosts, uniform inference: redshift minus random"),
    )
    contrasts = []
    metrics = ("median_error", "mean_error", "width_68", "width_95", "coverage_68", "coverage_95", "host_observed")
    by_scenario = {name: {int(row["universe_id"]): row for row in per_run if row["scenario"] == name}
                   for name, *_ in SCENARIOS}
    rng = np.random.default_rng(args.seed + 909_090)
    for a_name, b_name, label in contrast_specs:
        a_map, b_map = by_scenario[a_name], by_scenario[b_name]
        if set(a_map) != set(b_map):
            raise RuntimeError(f"Unpaired universe IDs for {a_name} and {b_name}")
        for metric in metrics:
            differences = []
            for universe_id in sorted(a_map):
                a, b = a_map[universe_id], b_map[universe_id]
                if metric in ("coverage_68", "coverage_95"):
                    low, high = ("lower_68", "upper_68") if metric == "coverage_68" else ("lower_95", "upper_95")
                    av = float(a[low] <= 70.0 <= a[high])
                    bv = float(b[low] <= 70.0 <= b[high])
                elif metric == "host_observed":
                    av, bv = float(a[metric]), float(b[metric])
                else:
                    av, bv = float(a[metric]), float(b[metric])
                differences.append(bv - av)
            array = np.asarray(differences, dtype=float)
            boot = np.empty(5000)
            for j in range(len(boot)):
                boot[j] = rng.choice(array, size=len(array), replace=True).mean()
            contrasts.append({
                "label": label, "scenario_A": a_name, "scenario_B": b_name,
                "metric": metric, "definition": "B minus A, paired by universe_id",
                "n_pairs": len(array), "mean_difference_B_minus_A": float(array.mean()),
                "median_difference_B_minus_A": float(np.median(array)),
                "bootstrap_95_interval": [float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))],
            })

    def write_csv(path: Path, rows: list[dict]) -> None:
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    write_csv(out / "scenario_summary.csv", scenario_rows)
    write_csv(out / "per_run_results.csv", per_run)
    write_csv(out / "paired_contrasts.csv", contrasts)
    metadata = {
        "original_repository": "vedantkushwah033-blip/bayesian-dark-siren-cosmology",
        "original_commit": "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f",
        "true_h0": 70.0, "n_galaxies": 500,
        "requested_replicates_per_scenario": args.replicates,
        "base_seed": args.seed, "scenario_count": len(SCENARIOS),
        "contrast_count": len(contrasts),
        "pairing": "One parent catalogue per universe; shared host/GW realization within host-generation setting; contrasts are paired by universe_id. Mask seeds are scenario-specific.",
        "interpretation_warning": "Synthetic model only; not LVK data or real-H0 measurement. Host weights are synthetic proxies.",
        "scenarios": scenario_rows,
    }
    (out / "scenario_summary.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps({"scenario_count": len(scenario_rows), "contrast_count": len(contrasts),
                      "replicates_per_scenario": args.replicates,
                      "per_run_rows": len(per_run), "paired_contrast_rows": len(contrasts)}, indent=2))


if __name__ == "__main__":
    main()
