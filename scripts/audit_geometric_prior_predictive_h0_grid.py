#!/usr/bin/env python3
"""Exact geometric H0-grid convergence audit for the 20-replicate prior-predictive pilot.

The archived KDE densities are held fixed. This script reconstructs the saved
1,121-point selected likelihoods, then exactly re-evaluates event numerators and
selection efficiencies on a 2,241-point geometric H0 grid. Arm E is re-integrated
over its 25 shared coefficient nodes. Posterior integration remains with respect
to dH0 under a uniform-H0 prior.
"""
from __future__ import annotations

import argparse
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import log_ndtr, logsumexp
from scipy.stats import norm

ARMS = (
    "A_catalogue_only",
    "B_parent_oracle",
    "C_known_probability_ipw",
    "D_map_plugin_ipw",
    "E_uncertainty_integrated_ipw",
)
GAMMAS = ("gm1", "gp0", "gp1")
SIZES = (100, 300, 500)
QUANTILES = np.array([0.025, 0.05, 0.10, 0.16, 0.50, 0.84, 0.90, 0.95, 0.975])


def _read_csv(archive: zipfile.ZipFile, path: str) -> pd.DataFrame:
    return pd.read_csv(io.BytesIO(archive.read(path)))


def _posterior(grid: np.ndarray, log_likelihood: np.ndarray, truth: float):
    shifted = log_likelihood - np.max(log_likelihood)
    density = np.exp(shifted)
    density /= trapezoid(density, grid)
    cdf = np.r_[0.0, cumulative_trapezoid(density, grid)]
    cdf /= cdf[-1]
    return density, np.interp(QUANTILES, cdf, grid), float(np.interp(truth, grid, cdf))


def _evaluate(grid, densities, redshift_grid, distances, sigma, threshold):
    delta = np.diff(redshift_grid)
    weights = np.empty_like(redshift_grid)
    weights[0] = delta[0] / 2.0
    weights[-1] = delta[-1] / 2.0
    weights[1:-1] = (delta[:-1] + delta[1:]) / 2.0
    result = np.zeros((densities.shape[0], grid.size), dtype=float)
    for start in range(0, grid.size, 128):
        stop = min(start + 128, grid.size)
        means = 299792.458 * redshift_grid[None, :] / grid[start:stop, None]
        efficiency_kernel = np.exp(log_ndtr((means - threshold) / sigma)) * weights
        efficiency = efficiency_kernel @ densities.T
        for distance in distances:
            numerator_kernel = norm.pdf(distance, loc=means, scale=sigma) * weights
            numerator = numerator_kernel @ densities.T
            if np.any(numerator <= 0) or np.any(efficiency <= 0):
                raise ArithmeticError("non-positive selected-likelihood integral")
            result[:, start:stop] += np.log(numerator.T) - np.log(efficiency.T)
    if not np.all(np.isfinite(result)):
        raise ArithmeticError("non-finite selected log likelihood")
    return result


def _find_root(archive: zipfile.ZipFile) -> str:
    candidates = [name for name in archive.namelist()
                  if name.endswith("replicates/rep000_gm1/config.json")]
    if len(candidates) != 1:
        raise ValueError(f"expected one rep000_gm1 config, found {len(candidates)}")
    return candidates[0].removesuffix("replicates/rep000_gm1/config.json")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-zip", required=True)
    parser.add_argument("--output-dir", default="results/geometric-h0-grid-audit")
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--source-points", type=int, default=1121)
    parser.add_argument("--fine-points", type=int, default=2241)
    args = parser.parse_args()
    if args.replicates < 2 or args.source_points < 2 or args.fine_points <= args.source_points:
        raise ValueError("invalid replicate or grid-point counts")
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    source_grid = np.geomspace(30.0, 110.0, args.source_points)
    fine_grid = np.geomspace(30.0, 110.0, args.fine_points)
    rows = []
    reconstruction_max = 0.0
    with zipfile.ZipFile(args.artifact_zip) as archive:
        root_prefix = _find_root(archive)
        for replicate in range(args.replicates):
            for gamma in GAMMAS:
                root = f"{root_prefix}replicates/rep{replicate:03d}_{gamma}/"
                cfg = json.loads(archive.read(root + "config.json"))
                if (cfg.get("n_detected") != 30
                    or cfg.get("h0_prior_min") != 30.0
                    or cfg.get("h0_prior_max") != 110.0
                    or cfg.get("h0_grid_type") != "geometric"
                    or cfg.get("h0_grid_points") != args.source_points):
                    raise ValueError(f"unexpected frozen design metadata in {root}")
                events = _read_csv(archive, root + "observable/gw_events_observed.csv")
                if len(events) != 30 or not np.all(events["detected"].to_numpy() == 1):
                    raise ValueError(f"event table integrity failure in {root}")
                distances = events["distance_mpc"].to_numpy(dtype=float)
                sigma = float(cfg["sigma_mpc"])
                threshold = float(cfg["detection_threshold_mpc"])
                truth = float(cfg["h0_true"])
                data = np.load(io.BytesIO(archive.read(root + "results/log_likelihoods.npz")))
                zgrid = data["redshift_grid"]
                if not np.allclose(data["H0_grid"], source_grid, rtol=0, atol=1e-12):
                    raise ValueError(f"archived geometric H0 grid mismatch in {root}")
                for n in SIZES:
                    density_stack = np.vstack([
                        data[f"density_{name}_n{n}"] for name in ARMS[:4]
                    ] + [data[f"density_E_nodes_n{n}"]])
                    if density_stack.shape != (29, zgrid.size):
                        raise ValueError(f"unexpected density stack shape in {root}, n={n}")
                    qweights = data[f"posterior_weights_n{n}"]
                    if qweights.shape != (25,) or not np.isclose(qweights.sum(), 1.0):
                        raise ValueError(f"invalid E quadrature weights in {root}, n={n}")
                    likelihoods = {}
                    for label, grid in (("source", source_grid), ("fine", fine_grid)):
                        node_ll = _evaluate(grid, density_stack, zgrid, distances, sigma, threshold)
                        marginal_e = logsumexp(node_ll[4:] + np.log(qweights[:, None]), axis=0)
                        likelihoods[label] = list(node_ll[:4]) + [marginal_e]
                        if label == "source":
                            for arm_index, arm in enumerate(ARMS[:4]):
                                reconstruction_max = max(
                                    reconstruction_max,
                                    float(np.max(np.abs(likelihoods[label][arm_index] - data[f"{arm}_n{n}"]))),
                                )
                            reconstruction_max = max(
                                reconstruction_max,
                                float(np.max(np.abs(marginal_e - data[f"{ARMS[4]}_n{n}"]))),
                            )
                    for arm_index, arm in enumerate(ARMS):
                        old_density, old_q, old_rank = _posterior(source_grid, likelihoods["source"][arm_index], truth)
                        fine_density, fine_q, fine_rank = _posterior(fine_grid, likelihoods["fine"][arm_index], truth)
                        l1 = float(trapezoid(
                            np.abs(np.interp(fine_grid, source_grid, old_density) - fine_density),
                            fine_grid,
                        ))
                        rows.append({
                            "replicate": replicate,
                            "gamma": gamma,
                            "validation_size": n,
                            "arm": arm,
                            "rank_source": old_rank,
                            "rank_fine": fine_rank,
                            "rank_shift": abs(old_rank - fine_rank),
                            "median_shift": abs(old_q[4] - fine_q[4]),
                            "max_nine_quantile_shift": float(np.max(np.abs(old_q - fine_q))),
                            "posterior_density_l1": l1,
                        })
    result = pd.DataFrame(rows)
    expected = args.replicates * len(GAMMAS) * len(SIZES) * len(ARMS)
    keys = ["replicate", "gamma", "validation_size", "arm"]
    if len(result) != expected or result.duplicated(keys).any():
        raise ValueError(f"expected {expected} unique rows, got {len(result)}")
    result.to_csv(output / "geometric_h0_grid_refinement_rows.csv", index=False)
    refined_pass = (
        reconstruction_max <= 1e-8
        and result["rank_shift"].max() <= 0.002
        and result["max_nine_quantile_shift"].max() <= 0.05
        and result["posterior_density_l1"].max() <= 0.001
    )
    summary = {
        "status": "PASS" if refined_pass else "FAIL",
        "rows": len(result),
        "replicates": args.replicates,
        "source_grid_points": args.source_points,
        "fine_grid_points": args.fine_points,
        "grid_type": "geometric",
        "likelihood_reconstruction_max_abs_error": reconstruction_max,
        "max_rank_shift": float(result["rank_shift"].max()),
        "max_median_shift": float(result["median_shift"].max()),
        "max_nine_quantile_shift": float(result["max_nine_quantile_shift"].max()),
        "max_posterior_density_l1": float(result["posterior_density_l1"].max()),
        "integration_measure": "dH0",
        "method": "exact selected-likelihood recomputation conditional on archived KDE densities; 25 shared E coefficient nodes",
        "limitations": "This audit checks H0-grid integration error conditional on archived KDE densities. It does not regenerate KDEs from raw galaxy weights or establish calibration.",
    }
    (output / "geometric_h0_grid_refinement_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    if not refined_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
