#!/usr/bin/env python3
"""Exact selected-likelihood H0-grid refinement audit for the five-arm SBC artifact.

The archived KDE host-density arrays are held fixed. For every archived replicate,
curvature and validation size, this script re-evaluates event numerators and
selection efficiencies at 561, 1121 and 2241 H0 points. Arm E is integrated over
its archived shared coefficient-node quadrature weights after evaluating each
node's complete multi-event likelihood. Posterior integration is with respect to
dH0 (not d log H0).
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
    quantiles = np.interp(QUANTILES, cdf, grid)
    rank = float(np.interp(truth, grid, cdf))
    return density, cdf, quantiles, rank


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
        for distance in distances:
            numerator_kernel = norm.pdf(distance, loc=means, scale=sigma) * weights
            numerator = numerator_kernel @ densities.T
            efficiency = efficiency_kernel @ densities.T
            if np.any(numerator <= 0) or np.any(efficiency <= 0):
                raise ArithmeticError("non-positive selected-likelihood integral")
            result[:, start:stop] += np.log(numerator.T) - np.log(efficiency.T)
    if not np.all(np.isfinite(result)):
        raise ArithmeticError("non-finite selected log likelihood")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-zip", required=True)
    parser.add_argument("--output-dir", default="results/exact-h0-grid-audit")
    args = parser.parse_args()
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    grids = {
        561: np.linspace(30.0, 110.0, 561),
        1121: np.geomspace(30.0, 110.0, 1121),
        2241: np.geomspace(30.0, 110.0, 2241),
    }
    rows = []
    reconstruction_max = 0.0
    with zipfile.ZipFile(args.artifact_zip) as archive:
        for replicate in range(50):
            for gamma in GAMMAS:
                root = f"five_arm_prior_predictive_50rep/replicates/rep{replicate:03d}_{gamma}/"
                cfg = json.loads(archive.read(root + "config.json"))
                if cfg.get("n_detected") != 30 or cfg.get("h0_prior_min") != 30.0 or cfg.get("h0_prior_max") != 110.0:
                    raise ValueError(f"unexpected design metadata in {root}")
                events = _read_csv(archive, root + "observable/gw_events_observed.csv")
                if len(events) != 30 or not np.all(events["detected"].to_numpy() == 1):
                    raise ValueError(f"event table integrity failure in {root}")
                distances = events["distance_mpc"].to_numpy(dtype=float)
                sigma = float(cfg["sigma_mpc"])
                threshold = float(cfg["detection_threshold_mpc"])
                truth = float(cfg["h0_true"])
                data = np.load(io.BytesIO(archive.read(root + "results/log_likelihoods.npz")))
                zgrid = data["redshift_grid"]
                if not np.allclose(data["H0_grid"], grids[561], rtol=0, atol=1e-12):
                    raise ValueError(f"archived 561-point grid mismatch in {root}")
                for n in SIZES:
                    density_stack = np.vstack([
                        data[f"density_{name}_n{n}"] for name in ARMS[:4]
                    ] + [data[f"density_E_nodes_n{n}"]])
                    if density_stack.shape != (29, zgrid.size):
                        raise ValueError(f"unexpected density stack shape in {root}, n={n}")
                    qweights = data[f"posterior_weights_n{n}"]
                    if qweights.shape != (25,) or not np.isclose(qweights.sum(), 1.0):
                        raise ValueError(f"invalid E quadrature weights in {root}, n={n}")
                    posteriors = {}
                    for points, grid in grids.items():
                        node_ll = _evaluate(grid, density_stack, zgrid, distances, sigma, threshold)
                        marginal_e = logsumexp(node_ll[4:] + np.log(qweights[:, None]), axis=0)
                        likelihoods = list(node_ll[:4]) + [marginal_e]
                        posteriors[points] = []
                        for arm_index, arm in enumerate(ARMS):
                            density, cdf, qs, rank = _posterior(grid, likelihoods[arm_index], truth)
                            posteriors[points].append((density, cdf, qs, rank))
                        if points == 561:
                            for arm_index, arm in enumerate(ARMS[:4]):
                                archived = data[f"{arm}_n{n}"]
                                reconstruction_max = max(
                                    reconstruction_max,
                                    float(np.max(np.abs(likelihoods[arm_index] - archived))),
                                )
                            archived_e = data[f"{ARMS[4]}_n{n}"]
                            reconstruction_max = max(
                                reconstruction_max, float(np.max(np.abs(marginal_e - archived_e)))
                            )
                    for arm_index, arm in enumerate(ARMS):
                        old = posteriors[561][arm_index]
                        mid = posteriors[1121][arm_index]
                        fine = posteriors[2241][arm_index]
                        # Compare normalized posterior densities on a common H0
                        # coordinate using interpolation only for the L1 diagnostic.
                        common = grids[2241]
                        l1_561_1121 = float(trapezoid(
                            np.abs(np.interp(common, grids[561], old[0]) - np.interp(common, grids[1121], mid[0])),
                            common,
                        ))
                        l1_1121_2241 = float(trapezoid(np.abs(np.interp(common, grids[1121], mid[0]) - fine[0]), common))
                        rows.append({
                            "replicate": replicate, "gamma": gamma, "validation_size": n, "arm": arm,
                            "rank_561": old[3], "rank_1121": mid[3], "rank_2241": fine[3],
                            "rank_shift_561_1121": abs(old[3] - mid[3]),
                            "rank_shift_1121_2241": abs(mid[3] - fine[3]),
                            "max_quantile_shift_561_1121": float(np.max(np.abs(old[2] - mid[2]))),
                            "max_quantile_shift_1121_2241": float(np.max(np.abs(mid[2] - fine[2]))),
                            "median_shift_561_1121": abs(old[2][4] - mid[2][4]),
                            "median_shift_1121_2241": abs(mid[2][4] - fine[2][4]),
                            "posterior_density_l1_561_1121": l1_561_1121,
                            "posterior_density_l1_1121_2241": l1_1121_2241,
                        })
    result = pd.DataFrame(rows)
    expected = 50 * 3 * 3 * 5
    if len(result) != expected or result.duplicated(["replicate", "gamma", "validation_size", "arm"]).any():
        raise ValueError(f"expected {expected} unique rows, got {len(result)}")
    result.to_csv(output / "exact_h0_grid_refinement_rows.csv", index=False)
    metrics = [c for c in result.columns if "shift_" in c or "_l1_" in c]
    refined_pass = (
        reconstruction_max <= 1e-8
        and result["rank_shift_1121_2241"].max() <= 0.002
        and result["max_quantile_shift_1121_2241"].max() <= 0.05
        and result["posterior_density_l1_1121_2241"].max() <= 0.001
    )
    legacy_pass = (
        result["rank_shift_561_1121"].max() <= 0.002
        and result["max_quantile_shift_561_1121"].max() <= 0.05
        and result["posterior_density_l1_561_1121"].max() <= 0.001
    )
    summary = {
        "status": "PASS_REFINED_REFERENCE_WITH_LEGACY_GRID_WARNING" if refined_pass and not legacy_pass
                  else ("PASS" if refined_pass and legacy_pass else "FAIL"),
        "refined_1121_to_2241_status": "PASS" if refined_pass else "FAIL",
        "legacy_561_to_1121_status": "PASS" if legacy_pass else "FAIL",
        "rows": len(result),
        "replicates": 50,
        "gamma_strata": list(GAMMAS),
        "validation_sizes": list(SIZES),
        "arms": list(ARMS),
        "likelihood_reconstruction_max_abs_error": reconstruction_max,
        "maxima": {column: float(result[column].max()) for column in metrics},
        "integration_measure": "dH0",
        "method": "exact selected-likelihood recomputation from frozen archived KDE densities; 25 shared E coefficient nodes",
        "limitations": "This isolates H0-grid integration error conditional on the archived KDE density estimates. It does not independently regenerate KDE densities from raw galaxy weights.",
    }
    (output / "exact_h0_grid_refinement_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if not refined_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
