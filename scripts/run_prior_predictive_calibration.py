"""Prior-predictive calibration diagnostic for the pinned dark-siren simulator.

Unlike fixed-H0 repeated-sampling coverage, this diagnostic draws H0 from the
same uniform prior used by posterior_h0. It reports interval coverage and
posterior-CDF ranks separately for the known-host and unknown-host cases.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np


LEVELS = (0.50, 0.68, 0.80, 0.90, 0.95)


def interval(grid: np.ndarray, posterior: np.ndarray, level: float) -> tuple[float, float, np.ndarray]:
    cdf = np.concatenate(
        [[0.0], np.cumsum((posterior[1:] + posterior[:-1]) * np.diff(grid) / 2.0)]
    )
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ValueError("Posterior CDF normalization failed.")
    cdf /= cdf[-1]
    alpha = (1.0 - level) / 2.0
    low = float(np.interp(alpha, cdf, grid))
    high = float(np.interp(1.0 - alpha, cdf, grid))
    return low, high, cdf


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", default="original-project")
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261023)
    parser.add_argument("--output-dir", default="results/prior-predictive-calibration")
    args = parser.parse_args()
    if args.replicates < 2:
        parser.error("--replicates must be at least 2")

    original = Path(args.original_dir).resolve()
    sys.path.insert(0, str(original))
    from src.simulations_v2 import simulate_galaxy_catalogue, choose_host, simulate_gw_event
    from src.bayesian_inference import posterior_h0

    h0_min, h0_max = 50.0, 90.0
    grid = np.linspace(h0_min, h0_max, 1601)
    rng = np.random.default_rng(args.seed)
    per_run: list[dict] = []
    aggregate = {
        scenario: {
            level: {"covered": 0, "widths": []}
            for level in LEVELS
        }
        for scenario in ("unknown_host", "known_host")
    }

    for r in range(args.replicates):
        h0_true = float(rng.uniform(h0_min, h0_max))
        sim_seed = args.seed + 1_000_000 + r
        cat = simulate_galaxy_catalogue(n_galaxies=500, H0_true=h0_true, seed=sim_seed)
        host = choose_host(cat, weighting="uniform", seed=sim_seed + 100000)
        gw = simulate_gw_event(
            cat["true_distance_mpc"][host],
            seed=sim_seed + 200000,
        )
        candidates = {
            "unknown_host": (cat["redshift"], np.ones(len(cat["redshift"]))),
            "known_host": (np.array([cat["redshift"][host]]), np.array([1.0])),
        }

        for scenario, (redshifts, weights) in candidates.items():
            posterior = posterior_h0(
                grid,
                redshifts,
                gw["observed_distance_mpc"],
                gw["sigma_mpc"],
                weights=weights,
                prior_min=h0_min,
                prior_max=h0_max,
            )
            cdf = np.concatenate(
                [[0.0], np.cumsum((posterior[1:] + posterior[:-1]) * np.diff(grid) / 2.0)]
            )
            cdf /= cdf[-1]
            rank = float(np.interp(h0_true, grid, cdf))
            row = {
                "replicate": r,
                "H0_true": h0_true,
                "scenario": scenario,
                "posterior_cdf_rank": rank,
            }
            for level in LEVELS:
                low, high, _ = interval(grid, posterior, level)
                covered = int(low <= h0_true <= high)
                aggregate[scenario][level]["covered"] += covered
                aggregate[scenario][level]["widths"].append(high - low)
                suffix = int(round(level * 100))
                row[f"lower_{suffix}"] = low
                row[f"upper_{suffix}"] = high
                row[f"covered_{suffix}"] = covered
            per_run.append(row)

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "prior_predictive_per_run.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(per_run[0]))
        writer.writeheader()
        writer.writerows(per_run)

    summary_rows = []
    for scenario, by_level in aggregate.items():
        for level, stats in by_level.items():
            coverage = stats["covered"] / args.replicates
            summary_rows.append({
                "scenario": scenario,
                "nominal_level": level,
                "replicates": args.replicates,
                "coverage": coverage,
                "coverage_mc_se": float(np.sqrt(coverage * (1.0 - coverage) / args.replicates)),
                "mean_interval_width": float(np.mean(stats["widths"])),
            })
    with (output / "prior_predictive_coverage.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)

    ranks = {
        scenario: [row["posterior_cdf_rank"] for row in per_run if row["scenario"] == scenario]
        for scenario in ("unknown_host", "known_host")
    }
    rank_summary = {
        scenario: {
            "mean_rank": float(np.mean(values)),
            "rank_sd": float(np.std(values, ddof=1)),
            "n": len(values),
        }
        for scenario, values in ranks.items()
    }
    metadata = {
        "replicates": args.replicates,
        "seed": args.seed,
        "H0_prior": f"Uniform({h0_min}, {h0_max})",
        "H0_grid_points": len(grid),
        "true_H0_draw_method": "independent uniform draws from the same bounds as the posterior prior",
        "rank_summary": rank_summary,
        "interpretation": (
            "Prior-predictive diagnostic, distinct from fixed-H0 frequentist coverage. "
            "Finite-replicate deviations are expected; do not infer calibration from a single metric."
        ),
    }
    (output / "prior_predictive_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))
    print(f"Wrote {len(per_run)} per-run rows and {len(summary_rows)} coverage rows to {output}")


if __name__ == "__main__":
    main()
