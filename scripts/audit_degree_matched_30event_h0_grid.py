#!/usr/bin/env python3
"""Exact H0-grid refinement audit for the paired 30-event degree pilot.

Re-evaluates every event numerator and selection-efficiency denominator from
archived KDE densities at 281 uniform, 561 geometric, and 1121 geometric H0
points. The 281-point curves are reconstructed first and must match the
archived selected log likelihoods before any refinement comparison is accepted.
Posterior normalization remains with respect to dH0 (uniform H0 prior).
"""
from __future__ import annotations

import argparse
import hashlib
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
GAMMAS = ("gamma_m1", "gamma_0", "gamma_p1")
SIZES = (100, 500)
PROBS = np.array([0.025, 0.16, 0.50, 0.84, 0.975])


def _read_csv(zf: zipfile.ZipFile, path: str) -> pd.DataFrame:
    return pd.read_csv(io.BytesIO(zf.read(path)))


def _bundle_prefix(zf: zipfile.ZipFile, gamma_dir: str) -> str:
    """Resolve the uploaded artifact's root layout without assuming a parent folder."""
    names = set(zf.namelist())
    for candidate in (f"{gamma_dir}/", f"degree_quadrature/{gamma_dir}/"):
        if candidate + "manifest.json" in names:
            return candidate
    raise FileNotFoundError(f"no manifest found for {gamma_dir} at expected artifact roots")


def _verify_manifest(zf: zipfile.ZipFile, prefix: str) -> int:
    manifest = json.loads(zf.read(prefix + "manifest.json"))
    if manifest.get("status") != "PASS":
        raise ValueError(f"{prefix}: bundle manifest is not PASS")
    files = manifest.get("files", {})
    if manifest.get("file_count") != len(files):
        raise ValueError(f"{prefix}: manifest file_count mismatch")
    checked = 0
    for rel, meta in files.items():
        payload = zf.read(prefix + rel)
        if len(payload) != meta["bytes"]:
            raise ValueError(f"{prefix}{rel}: byte count mismatch")
        if hashlib.sha256(payload).hexdigest() != meta["sha256"]:
            raise ValueError(f"{prefix}{rel}: SHA-256 mismatch")
        checked += 1
    return checked


def _selected_loglike(grid: np.ndarray, densities: np.ndarray,
                      zgrid: np.ndarray, distances: np.ndarray,
                      sigma: float, threshold: float) -> np.ndarray:
    """Independent trapezoid integration of selected event likelihood terms."""
    delta = np.diff(zgrid)
    dz = np.empty_like(zgrid)
    dz[0], dz[-1] = delta[0] / 2.0, delta[-1] / 2.0
    dz[1:-1] = (delta[:-1] + delta[1:]) / 2.0
    out = np.zeros((densities.shape[0], grid.size), dtype=float)
    for start in range(0, grid.size, 128):
        stop = min(start + 128, grid.size)
        means = 299792.458 * zgrid[None, :] / grid[start:stop, None]
        efficiency_kernel = np.exp(log_ndtr((means - threshold) / sigma)) * dz
        # The selection denominator depends on H0 and the host density, not on
        # the observed event distance; compute it once per H0 block and reuse it.
        efficiency = efficiency_kernel @ densities.T
        if np.any(efficiency <= 0):
            raise ArithmeticError("non-positive selection-efficiency integral")
        log_efficiency = np.log(efficiency.T)
        for distance in distances:
            numerator_kernel = norm.pdf(distance, loc=means, scale=sigma) * dz
            numerator = numerator_kernel @ densities.T
            if np.any(numerator <= 0):
                raise ArithmeticError("non-positive event-numerator integral")
            out[:, start:stop] += np.log(numerator.T) - log_efficiency
    if not np.all(np.isfinite(out)):
        raise ArithmeticError("non-finite selected log likelihood")
    return out


def _posterior(grid: np.ndarray, loglike: np.ndarray, truth: float):
    density = np.exp(loglike - np.max(loglike))
    area = float(trapezoid(density, grid))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("posterior normalization failed")
    density /= area
    cdf = np.r_[0.0, cumulative_trapezoid(density, grid)]
    cdf /= cdf[-1]
    qs = np.interp(PROBS, cdf, grid)
    return density, qs, float(np.interp(truth, grid, cdf))


def run(artifact_zip: Path, output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise ValueError("output directory must be empty")
    grids = {
        281: np.linspace(30.0, 110.0, 281),
        561: np.geomspace(30.0, 110.0, 561),
        1121: np.geomspace(30.0, 110.0, 1121),
    }
    rows = []
    manifest_files = 0
    max_reconstruction = 0.0
    with zipfile.ZipFile(artifact_zip) as zf:
        for gamma_dir in GAMMAS:
            prefix = _bundle_prefix(zf, gamma_dir)
            manifest_files += _verify_manifest(zf, prefix)
            base = prefix + "degree1_order9/"
            cfg = json.loads(zf.read(base + "config.json"))
            if (cfg.get("n_parent") != 1200 or cfg.get("n_detected") != 30
                or cfg.get("h0_true") != 70.0
                or cfg.get("h0_prior_support") != [30.0, 110.0]
                or cfg.get("h0_grid_points") != 281
                or cfg.get("quadrature_order") != 9):
                raise ValueError(f"{gamma_dir}: unexpected frozen configuration")
            truth = float(cfg["h0_true"])
            sigma = float(cfg["sigma_mpc"])
            threshold = float(cfg["detection_threshold_mpc"])
            events = _read_csv(zf, base + "observable/gw_events_observed.csv")
            if len(events) != 30 or not np.all(events["detected"].to_numpy() == 1):
                raise ValueError(f"{gamma_dir}: expected 30 detected events")
            distances = events["distance_mpc"].to_numpy(float)

            for degree in (1, 2):
                child_prefix = prefix + f"degree{degree}_order9/"
                child_cfg = json.loads(zf.read(child_prefix + "config.json"))
                if child_cfg.get("inclusion_model_degree") != degree:
                    raise ValueError(f"{gamma_dir}: degree metadata mismatch")
                npz = np.load(io.BytesIO(zf.read(child_prefix + "results/log_likelihoods.npz")))
                zgrid = np.asarray(npz["redshift_grid"], float)
                if not np.array_equal(np.asarray(npz["H0_grid"], float), grids[281]):
                    raise ValueError(f"{gamma_dir}: archived 281-point H0 grid mismatch")

                for n in SIZES:
                    qweights = np.asarray(npz[f"posterior_weights_n{n}"], float)
                    nodes = np.asarray(npz[f"posterior_nodes_n{n}"])
                    node_density = np.asarray(npz[f"density_E_nodes_n{n}"], float)
                    if node_density.ndim != 2 or node_density.shape[0] != len(qweights):
                        raise ValueError(f"{gamma_dir}: E node density/weight shape mismatch")
                    if nodes.shape[0] != len(qweights) or not np.isclose(qweights.sum(), 1.0):
                        raise ValueError(f"{gamma_dir}: invalid shared E coefficient quadrature")
                    if np.any(qweights <= 0) or np.any(node_density < 0):
                        raise ValueError(f"{gamma_dir}: invalid E quadrature weights/densities")

                    density_by_arm = {}
                    for arm in ARMS[:4]:
                        density = np.asarray(npz[f"density_{arm}_n{n}"], float)
                        if density.shape != zgrid.shape:
                            raise ValueError(f"{gamma_dir}: {arm} density shape mismatch")
                        density_by_arm[arm] = density[None, :]
                    density_by_arm[ARMS[4]] = node_density

                    loglikes_by_grid = {}
                    for points, grid in grids.items():
                        arm_loglikes = {}
                        for arm, density_stack in density_by_arm.items():
                            node_ll = _selected_loglike(
                                grid, density_stack, zgrid, distances, sigma, threshold
                            )
                            if arm == ARMS[4]:
                                marginal = logsumexp(node_ll + np.log(qweights[:, None]), axis=0)
                                arm_loglikes[arm] = marginal
                            else:
                                arm_loglikes[arm] = node_ll[0]
                        loglikes_by_grid[points] = arm_loglikes

                    for arm in ARMS:
                        archived = np.asarray(npz[f"{arm}_n{n}"], float)
                        diff = float(np.max(np.abs(loglikes_by_grid[281][arm] - archived)))
                        max_reconstruction = max(max_reconstruction, diff)
                        if diff > 1e-8:
                            raise ValueError(
                                f"{gamma_dir} degree{degree} n={n} {arm}: "
                                f"281-grid likelihood reconstruction error {diff}"
                            )
                        posterior = {
                            points: _posterior(grid, loglikes_by_grid[points][arm], truth)
                            for points, grid in grids.items()
                        }
                        common = grids[1121]
                        dens281 = np.interp(common, grids[281], posterior[281][0])
                        dens561 = np.interp(common, grids[561], posterior[561][0])
                        dens1121 = posterior[1121][0]
                        q281, q561, q1121 = (
                            posterior[p][1] for p in (281, 561, 1121)
                        )
                        rank281, rank561, rank1121 = (
                            posterior[p][2] for p in (281, 561, 1121)
                        )
                        rows.append({
                            "gamma": gamma_dir, "degree": degree,
                            "validation_size": n, "arm": arm,
                            "rank_281": rank281, "rank_561": rank561, "rank_1121": rank1121,
                            "rank_shift_281_561": abs(rank281-rank561),
                            "rank_shift_561_1121": abs(rank561-rank1121),
                            "median_shift_281_561": abs(q281[2]-q561[2]),
                            "median_shift_561_1121": abs(q561[2]-q1121[2]),
                            "max_quantile_shift_281_561": float(np.max(np.abs(q281-q561))),
                            "max_quantile_shift_561_1121": float(np.max(np.abs(q561-q1121))),
                            "density_l1_281_561": float(trapezoid(np.abs(dens281-dens561), common)),
                            "density_l1_561_1121": float(trapezoid(np.abs(dens561-dens1121), common)),
                        })
                npz.close()

    result = pd.DataFrame(rows)
    expected = len(GAMMAS) * 2 * len(SIZES) * len(ARMS)
    if len(result) != expected or result.duplicated(
        ["gamma", "degree", "validation_size", "arm"]
    ).any():
        raise ValueError(f"expected {expected} unique curves, got {len(result)}")
    result.to_csv(output_dir / "h0_grid_refinement_rows.csv", index=False)
    fine = result[["rank_shift_561_1121", "median_shift_561_1121",
                   "max_quantile_shift_561_1121", "density_l1_561_1121"]]
    legacy = result[["rank_shift_281_561", "median_shift_281_561",
                     "max_quantile_shift_281_561", "density_l1_281_561"]]
    fine_max = fine.max().to_dict()
    legacy_max = legacy.max().to_dict()
    fine_pass = (
        fine_max["rank_shift_561_1121"] <= 0.002
        and fine_max["median_shift_561_1121"] <= 0.01
        and fine_max["max_quantile_shift_561_1121"] <= 0.05
        and fine_max["density_l1_561_1121"] <= 0.001
    )
    legacy_pass = (
        legacy_max["rank_shift_281_561"] <= 0.002
        and legacy_max["median_shift_281_561"] <= 0.01
        and legacy_max["max_quantile_shift_281_561"] <= 0.05
        and legacy_max["density_l1_281_561"] <= 0.001
    )
    summary = {
        "status": "PASS_REFINED_REFERENCE_WITH_LEGACY_GRID_WARNING"
                  if fine_pass and not legacy_pass else ("PASS" if fine_pass else "FAIL"),
        "source_artifact": "degree-matched-30event-quadrature-sensitivity",
        "source_manifest_payloads_verified": manifest_files,
        "curvature_bundles": len(GAMMAS),
        "unique_posterior_curves": len(result),
        "likelihood_reconstruction_max_abs_error": max_reconstruction,
        "legacy_281_to_561_status": "PASS" if legacy_pass else "FAIL",
        "refined_561_to_1121_status": "PASS" if fine_pass else "FAIL",
        "legacy_maxima": legacy_max,
        "refined_maxima": fine_max,
        "integration_measure": "dH0",
        "method": "exact selected-likelihood recomputation from archived KDE density arrays and shared E coefficient nodes",
        "limitations": [
            "Conditional on the archived KDE densities and coefficient posterior quadrature nodes.",
            "Does not regenerate the KDEs from raw galaxy weights or assess bandwidth uncertainty.",
            "Single fixed-truth realization; not calibration, SBC, or estimator superiority evidence.",
        ],
    }
    (output_dir / "h0_grid_refinement_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    if not fine_pass:
        raise SystemExit(1)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-zip", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    run(args.artifact_zip, args.output_dir)


if __name__ == "__main__":
    main()
