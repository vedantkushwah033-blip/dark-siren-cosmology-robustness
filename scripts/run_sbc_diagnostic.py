"""Prior-predictive simulation-based calibration for the synthetic dark-siren model.

The measurement uncertainty is drawn independently of H0 and host identity so
that the data-generating model matches the likelihood's fixed-sigma assumption.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

import numpy as np
from scipy.stats import kstest

C_KM_S = 299792.458
H0_MIN, H0_MAX = 50.0, 90.0
N_BINS = 10


def posterior_cdf(grid: np.ndarray, posterior: np.ndarray) -> np.ndarray:
    cdf = np.concatenate(
        [[0.0], np.cumsum((posterior[1:] + posterior[:-1]) * np.diff(grid) / 2)]
    )
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ValueError("Posterior CDF normalization failed")
    return cdf / cdf[-1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original-dir", default="original-project")
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=99173)
    parser.add_argument("--output-dir", default="results/sbc")
    args = parser.parse_args()
    if args.replicates < 100:
        parser.error("--replicates must be at least 100 for a useful rank histogram")

    original = Path(args.original_dir).resolve()
    sys.path.insert(0, str(original))
    from src.simulations_v2 import simulate_galaxy_catalogue, choose_host
    from src.bayesian_inference import posterior_h0

    grid = np.linspace(H0_MIN, H0_MAX, 1601)
    ranks: dict[str, list[dict]] = {"known_host": [], "unknown_host": []}

    for r in range(args.replicates):
        rng = np.random.default_rng(args.seed + r)
        h0_true = float(rng.uniform(H0_MIN, H0_MAX))
        catalogue_seed = args.seed + 10_000_000 + r
        host_seed = args.seed + 20_000_000 + r
        cat = simulate_galaxy_catalogue(
            n_galaxies=500, H0_true=h0_true, seed=catalogue_seed
        )
        host = choose_host(cat, weighting="uniform", seed=host_seed)

        # Sigma is a data-design quantity drawn independently of H0 and host,
        # matching the likelihood's assumption that the reported sigma is fixed.
        sigma_mpc = float(rng.uniform(15.0, 35.0))
        true_distance = C_KM_S * float(cat["redshift"][host]) / h0_true
        observed_distance = float(rng.normal(true_distance, sigma_mpc))

        scenarios = {
            "known_host": (
                np.array([cat["redshift"][host]], dtype=float),
                np.array([1.0], dtype=float),
            ),
            "unknown_host": (
                np.asarray(cat["redshift"], dtype=float),
                np.ones(len(cat["redshift"]), dtype=float),
            ),
        }
        for scenario, (redshifts, weights) in scenarios.items():
            posterior = posterior_h0(
                grid,
                redshifts,
                observed_distance,
                sigma_mpc,
                weights=weights,
                prior_min=H0_MIN,
                prior_max=H0_MAX,
            )
            cdf = posterior_cdf(grid, posterior)
            rank = float(np.interp(h0_true, grid, cdf))
            ranks[scenario].append({
                "replicate": r,
                "H0_true": h0_true,
                "observed_distance_mpc": observed_distance,
                "sigma_mpc": sigma_mpc,
                "posterior_cdf_at_truth": rank,
            })

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    for scenario, rows in ranks.items():
        with (output_dir / f"{scenario}_sbc_ranks.csv").open(
            "w", newline="", encoding="utf-8"
        ) as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        values = np.asarray([row["posterior_cdf_at_truth"] for row in rows])
        counts, _ = np.histogram(values, bins=np.linspace(0.0, 1.0, N_BINS + 1))
        ks = kstest(values, "uniform")
        summary_rows.append({
            "scenario": scenario,
            "replicates": len(values),
            "ks_statistic_vs_uniform": float(ks.statistic),
            "ks_p_value_descriptive": float(ks.pvalue),
            "mean_rank": float(values.mean()),
            "rank_histogram_counts": ";".join(str(int(x)) for x in counts),
            "minimum_rank": float(values.min()),
            "maximum_rank": float(values.max()),
        })
    with (output_dir / "sbc_summary.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Wrote SBC ranks and summaries to {output_dir}")


if __name__ == "__main__":
    main()
