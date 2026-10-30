from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from dark_siren_robustness.finite_population_mean_bootstrap import (
    studentized_pseudopopulation_srswor_bootstrap,
)


def _population(rng: np.random.Generator, kind: str, size: int) -> np.ndarray:
    if kind == "normal":
        return rng.normal(0.0, 1.0, size)
    if kind == "gamma":
        return rng.gamma(shape=0.7, scale=1.0 / 0.7, size=size)
    if kind == "lognormal":
        return rng.lognormal(mean=0.0, sigma=1.0, size=size)
    raise ValueError(f"unknown population kind: {kind}")


def _sha256_array(x: np.ndarray) -> str:
    arr = np.ascontiguousarray(np.asarray(x, dtype="<f8"))
    h = hashlib.sha256()
    h.update(str(arr.shape).encode("ascii"))
    h.update(arr.tobytes())
    return h.hexdigest()


def run_benchmark(
    *,
    scenarios: int = 3,
    repetitions: int = 100,
    population_size: int = 1000,
    sample_sizes: tuple[int, ...] = (100, 400),
    bootstrap_replicates: int = 199,
    seed: int = 20261031,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    for name, value, minimum in (
        ("scenarios", scenarios, 1),
        ("repetitions", repetitions, 2),
        ("population_size", population_size, 20),
        ("bootstrap_replicates", bootstrap_replicates, 20),
        ("seed", seed, 0),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")
    if not sample_sizes or len(set(sample_sizes)) != len(sample_sizes):
        raise ValueError("sample_sizes must be nonempty and unique")
    if any(isinstance(n, bool) or not isinstance(n, int) or n < 4 or n >= population_size
           for n in sample_sizes):
        raise ValueError("each sample size must be an integer in [4, population_size)")
    kinds = ("normal", "gamma", "lognormal")
    root = np.random.SeedSequence(seed)
    scenario_seeds = root.spawn(scenarios)
    rows: list[dict[str, object]] = []
    replicate_rows: list[dict[str, object]] = []
    frame_rows: list[dict[str, object]] = []

    for s, child in enumerate(scenario_seeds):
        kind = kinds[s % len(kinds)]
        rng = np.random.default_rng(child)
        population = _population(rng, kind, population_size)
        target = float(np.mean(population))
        frame_hash = _sha256_array(population)
        frame_rows.extend({"scenario": s, "unit": i, "value": float(value), "frame_sha256": frame_hash} for i, value in enumerate(population))
        sample_seeds = child.spawn(repetitions * len(sample_sizes))
        seed_index = 0
        for n in sample_sizes:
            method_stats = {
                "studentized_pseudopop_srswor": [],
                "wald_srswor": [],
                "ordinary_percentile_bootstrap": [],
            }
            method_widths = {key: [] for key in method_stats}
            method_failures = {key: 0 for key in method_stats}
            for r in range(repetitions):
                sample_seed = int(sample_seeds[seed_index].generate_state(1)[0])
                seed_index += 1
                srng = np.random.default_rng(sample_seed)
                indices = srng.choice(population_size, size=n, replace=False)
                sample = population[indices]
                estimate = float(np.mean(sample))
                sample_var = float(np.var(sample, ddof=1))
                design_se = float(np.sqrt((1.0 - n / population_size) * sample_var / n))

                # Conventional normal/Wald interval using the SRSWOR variance estimate.
                z = 1.959963984540054
                wald = (estimate - z * design_se, estimate + z * design_se)
                method_stats["wald_srswor"].append(wald[0] <= target <= wald[1])
                method_widths["wald_srswor"].append(wald[1] - wald[0])
                replicate_rows.append({"scenario": s, "population_kind": kind, "population_size": population_size, "sample_size": n, "repetition": r, "method": "wald_srswor", "target": target, "estimate": estimate, "interval_lower": wald[0], "interval_upper": wald[1], "covered": int(wald[0] <= target <= wald[1]), "failure_code": "", "sample_indices": json.dumps(indices.tolist()), "sample_seed": sample_seed, "ordinary_bootstrap_seed": ordinary_bootstrap_seed, "studentized_bootstrap_seed": studentized_bootstrap_seed, "frame_sha256": frame_hash})

                # Ordinary with-replacement percentile bootstrap: diagnostic comparator.
                ordinary_bootstrap_seed = int(np.random.SeedSequence([seed, s, n, r, 41031]).generate_state(1)[0])
                brng = np.random.default_rng(ordinary_bootstrap_seed)
                draws = brng.integers(0, n, size=(bootstrap_replicates, n))
                means = np.mean(sample[draws], axis=1)
                lo, hi = np.quantile(means, [0.025, 0.975])
                method_stats["ordinary_percentile_bootstrap"].append(
                    bool(lo <= target <= hi)
                )
                method_widths["ordinary_percentile_bootstrap"].append(float(hi - lo))
                replicate_rows.append({"scenario": s, "population_kind": kind, "population_size": population_size, "sample_size": n, "repetition": r, "method": "ordinary_percentile_bootstrap", "target": target, "estimate": estimate, "interval_lower": float(lo), "interval_upper": float(hi), "covered": int(lo <= target <= hi), "failure_code": "", "sample_indices": json.dumps(indices.tolist()), "frame_sha256": frame_hash})

                # Studentized pseudo-population bootstrap, with deterministic per-draw seed.
                studentized_bootstrap_seed = int(np.random.SeedSequence([seed, s, n, r, 41032]).generate_state(1)[0])
                try:
                    result = studentized_pseudopopulation_srswor_bootstrap(
                        sample, population_size,
                        n_bootstrap=bootstrap_replicates,
                        seed=studentized_bootstrap_seed,
                    )
                    interval = result["interval"]
                    if result["status"] != "ok" or not np.all(np.isfinite(interval)):
                        raise ArithmeticError("studentized bootstrap did not yield a finite interval")
                    method_stats["studentized_pseudopop_srswor"].append(
                        bool(interval[0] <= target <= interval[1])
                    )
                    method_widths["studentized_pseudopop_srswor"].append(
                        float(interval[1] - interval[0])
                    )
                    replicate_rows.append({"scenario": s, "population_kind": kind, "population_size": population_size, "sample_size": n, "repetition": r, "method": "studentized_pseudopop_srswor", "target": target, "estimate": estimate, "interval_lower": float(interval[0]), "interval_upper": float(interval[1]), "covered": int(interval[0] <= target <= interval[1]), "failure_code": "", "sample_indices": json.dumps(indices.tolist()), "frame_sha256": frame_hash})
                except (ValueError, ArithmeticError, FloatingPointError) as exc:
                    method_failures["studentized_pseudopop_srswor"] += 1
                    replicate_rows.append({"scenario": s, "population_kind": kind, "population_size": population_size, "sample_size": n, "repetition": r, "method": "studentized_pseudopop_srswor", "target": target, "estimate": estimate, "interval_lower": float("nan"), "interval_upper": float("nan"), "covered": "", "failure_code": type(exc).__name__ + ":" + str(exc), "sample_indices": json.dumps(indices.tolist()), "frame_sha256": frame_hash})

            for method, covered in method_stats.items():
                successful = len(covered)
                coverage = float(np.mean(covered)) if successful else float("nan")
                # Binomial Monte Carlo SE conditional on successful interval construction.
                mcse = float(np.sqrt(coverage * (1.0 - coverage) / successful)) if successful else float("nan")
                rows.append({
                    "scenario": s,
                    "population_kind": kind,
                    "population_size": population_size,
                    "sample_size": n,
                    "sampling_fraction": n / population_size,
                    "repetitions": repetitions,
                    "method": method,
                    "target_finite_population_mean": target,
                    "frame_sha256": frame_hash,
                    "successful_intervals": successful,
                    "failed_intervals": method_failures[method],
                    "coverage": coverage,
                    "coverage_mcse_conditional_on_success": mcse,
                    "mean_interval_width": float(np.mean(method_widths[method])) if method_widths[method] else float("nan"),
                })
    return rows, replicate_rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Repeated-sampling linear finite-population bootstrap benchmark")
    parser.add_argument("--scenarios", type=int, default=3)
    parser.add_argument("--repetitions", type=int, default=100)
    parser.add_argument("--population-size", type=int, default=1000)
    parser.add_argument("--sample-sizes", type=int, nargs="+", default=[100, 400])
    parser.add_argument("--bootstrap-replicates", type=int, default=199)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    rows = run_benchmark(
        scenarios=args.scenarios,
        repetitions=args.repetitions,
        population_size=args.population_size,
        sample_sizes=tuple(args.sample_sizes),
        bootstrap_replicates=args.bootstrap_replicates,
        seed=args.seed,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "linear_finite_population_bootstrap.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    frame_path = args.output_dir / "finite_population_frames.csv"
    with frame_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(frame_rows[0]))
        writer.writeheader()
        writer.writerows(frame_rows)
    replicate_path = args.output_dir / "linear_finite_population_bootstrap_replicates.csv"
    with replicate_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(replicate_rows[0]))
        writer.writeheader()
        writer.writerows(replicate_rows)
    manifest = {
        "schema": "linear_finite_population_bootstrap_benchmark_v1",
        "status": "PILOT_NOT_CONFIRMATORY",
        "seed": args.seed,
        "scenarios": args.scenarios,
        "repetitions_per_scenario_and_sample_size": args.repetitions,
        "population_size": args.population_size,
        "sample_sizes": args.sample_sizes,
        "bootstrap_replicates": args.bootstrap_replicates,
        "methods": sorted({str(row["method"]) for row in rows}),
        "rows": len(rows),
        "interpretation": "Coverage is conditional on successful interval construction; failures are reported separately. Scenarios are distinct fixed finite frames. No result validates a nonlinear posterior functional.",
        "summary_csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        "replicate_csv_sha256": hashlib.sha256(replicate_path.read_bytes()).hexdigest(),
        "frame_csv_sha256": hashlib.sha256(frame_path.read_bytes()).hexdigest(),
        "frame_rows": len(frame_rows),
        "replicate_rows": len(replicate_rows),
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
