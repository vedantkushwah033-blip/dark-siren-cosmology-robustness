"""Finite-parent SRSWOR benchmark for validation-based host weighting.

Synthetic conditional design: each scenario fixes its parent frame, thinned
catalogue, and GW distances; repeated validation samples are drawn exactly
without replacement. This is not a real LVK or survey analysis.
"""
from __future__ import annotations

from collections.abc import Sequence
import hashlib

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.special import expit

from dark_siren_robustness.catalogue_density import estimate_host_density
from dark_siren_robustness.estimated_inclusion import estimate_host_density_from_validation
from dark_siren_robustness.generative_selection import (
    generate_synthetic_population,
    latent_galaxy_intensity,
)
from dark_siren_robustness.host_density_benchmark import posterior_summary
from dark_siren_robustness.selection_toy import selection_aware_toy_likelihood


def run_finite_parent_benchmark(
    *,
    n_scenarios: int = 10,
    n_validation_draws: int = 20,
    n_parent: int = 4000,
    n_detected: int = 30,
    validation_sizes: Sequence[int] = (100, 300),
    seed: int = 20261032,
    h0_true: float = 70.0,
    sigma_mpc: float = 140.0,
    detection_threshold_mpc: float = 250.0,
    h0_points: int = 81,
    z_points: int = 401,
    true_intercept: float = 0.25,
    true_slope: float = -1.2,
    true_quadratic: float = 0.0,
) -> pd.DataFrame:
    """Run a nested exact-SRSWOR benchmark and preserve every failed attempt."""
    sizes = tuple(sorted(set(int(x) for x in validation_sizes)))
    if n_scenarios < 1 or n_validation_draws < 2 or n_parent < 100:
        raise ValueError("invalid scenario, draw, or parent counts")
    if not sizes or min(sizes) < 10 or max(sizes) > n_parent:
        raise ValueError("validation sizes must lie in [10, n_parent]")
    if h0_points < 31 or z_points < 101:
        raise ValueError("grids are too small for the declared benchmark")

    z_grid = np.linspace(0.01, 0.40, z_points)
    h0_grid = np.linspace(50.0, 90.0, h0_points)
    zmin, zmax = float(z_grid[0]), float(z_grid[-1])
    zmid, zhalf = (zmin + zmax) / 2, (zmax - zmin) / 2
    rate = 1.0 + 1.5 * (z_grid - zmin) / (zmax - zmin)
    target = latent_galaxy_intensity(z_grid) * rate
    target = target / trapezoid(target, z_grid)
    rows: list[dict] = []

    for scenario in range(n_scenarios):
        root = np.random.SeedSequence([int(seed), scenario])
        sim_ss, thin_ss, val_ss = root.spawn(3)
        sim_seed = int(sim_ss.generate_state(1, dtype=np.uint32)[0])
        thin_seed = int(thin_ss.generate_state(1, dtype=np.uint32)[0])
        validation_seed = int(val_ss.generate_state(1, dtype=np.uint32)[0])
        sim = generate_synthetic_population(
            seed=sim_seed, n_parent=n_parent, n_detected=n_detected,
            h0_true=h0_true, sigma_mpc=sigma_mpc,
            detection_threshold_mpc=detection_threshold_mpc,
            z_min=zmin, z_max=zmax,
        )
        parent_z = sim["parent_truth"]["redshift"].to_numpy(dtype=float)
        distances = sim["gw_events_observed"]["distance_mpc"].to_numpy(dtype=float)
        x = (parent_z - zmid) / zhalf
        p_true = expit(true_intercept + true_slope*x + true_quadratic*(x*x - 1/3))
        if np.any((p_true <= 0) | (p_true >= 1)) or not np.all(np.isfinite(p_true)):
            raise ValueError("true inclusion probabilities saturated")
        thin_rng = np.random.default_rng(thin_seed)
        retained = thin_rng.random(n_parent) < p_true
        observed_z = parent_z[retained]
        parent_hash = hashlib.sha256(np.ascontiguousarray(parent_z).view(np.uint8)).hexdigest()
        catalogue_hash = hashlib.sha256(np.ascontiguousarray(observed_z).view(np.uint8)).hexdigest()
        gw_hash = hashlib.sha256(np.ascontiguousarray(distances).view(np.uint8)).hexdigest()
        if observed_z.size < 5 or np.unique(observed_z).size < 3:
            for size in sizes:
                for draw in range(n_validation_draws):
                    rows.append({"scenario": scenario, "draw": draw,
                                 "validation_size": size, "status": "failed",
                                 "failure_reason": "catalogue_minimum_size_gate",
                                 "parent_frame_sha256": parent_hash,
                                 "catalogue_sha256": catalogue_hash,
                                 "gw_distance_sha256": gw_hash,
                                 "covered_68": None})
            continue

        # Oracle and catalogue-only are fixed within a scenario; the primary
        # repeated-sampling estimand is the variation across validation draws.
        oracle = estimate_host_density(
            z_grid, observed_z, inclusion_probability=p_true[retained],
            merger_weight=rate,
        )
        oracle_post = posterior_summary(
            h0_grid, selection_aware_toy_likelihood(
                z_grid, oracle, h0_grid, distances, sigma=sigma_mpc,
                detection_threshold=detection_threshold_mpc,
            )["selected_log_likelihood"],
        )
        catalogue = estimate_host_density(z_grid, observed_z, merger_weight=rate)
        cat_post = posterior_summary(
            h0_grid, selection_aware_toy_likelihood(
                z_grid, catalogue, h0_grid, distances, sigma=sigma_mpc,
                detection_threshold=detection_threshold_mpc,
            )["selected_log_likelihood"],
        )
        rng = np.random.default_rng(validation_seed)
        for size in sizes:
            for draw in range(n_validation_draws):
                indices = rng.choice(n_parent, size=size, replace=False)
                # SRSWOR is checked at generation time, not inferred afterward.
                if np.unique(indices).size != size:
                    raise RuntimeError("SRSWOR generator returned duplicate indices")
                try:
                    result = estimate_host_density_from_validation(
                        z_grid, observed_z, parent_z[indices], retained[indices].astype(int),
                        merger_weight=rate,
                    )
                    density = np.asarray(result["density"], dtype=float)
                    likelihood = selection_aware_toy_likelihood(
                        z_grid, density, h0_grid, distances, sigma=sigma_mpc,
                        detection_threshold=detection_threshold_mpc,
                    )["selected_log_likelihood"]
                    post = posterior_summary(h0_grid, likelihood)
                    rows.append({
                        "scenario": scenario, "draw": draw,
                        "validation_size": size, "status": "success",
                        "failure_reason": "",
                        "validation_seed": validation_seed,
                        "validation_indices": ",".join(map(str, indices.tolist())),
                        "n_validation_included": int(np.sum(retained[indices])),
                        "fit_status": result["model"]["fit_status"],
                        "fit_boundary_limited": bool(result["model"]["coefficient_bound_hit"]),
                        "weight_ess": float(result["inverse_weight_effective_sample_size"]),
                        "max_inverse_weight": float(result["inverse_weight_maximum"]),
                        "p99_inverse_weight": float(result["inverse_weight_p99"]),
                        "inclusion_probability_mae": float(np.mean(np.abs(
                            expit(float(result["model"]["intercept"]) +
                                  float(result["model"]["slope"])*x) - p_true
                        ))),
                        "density_ise": float(trapezoid((density-target)**2, z_grid)),
                        "posterior_mean": float(post["posterior_mean"]),
                        "posterior_median": float(post["posterior_median"]),
                        "ci68_low": float(post["ci68_low"]),
                        "ci68_high": float(post["ci68_high"]),
                        "ci68_width": float(post["ci68_width"]),
                        "covered_68": bool(post["ci68_low"] <= h0_true <= post["ci68_high"]),
                        "absolute_error": float(abs(post["posterior_mean"]-h0_true)),
                        "oracle_posterior_mean": float(oracle_post["posterior_mean"]),
                        "catalogue_only_posterior_mean": float(cat_post["posterior_mean"]),
                        "h0_true": float(h0_true),
                        "n_catalogue": int(observed_z.size),
                        "n_detected": int(distances.size),
                        "parent_size": int(n_parent),
                        "parent_frame_sha256": parent_hash,
                        "catalogue_sha256": catalogue_hash,
                        "gw_distance_sha256": gw_hash,
                    })
                except Exception as exc:
                    rows.append({
                        "scenario": scenario, "draw": draw,
                        "validation_size": size, "status": "failed",
                        "failure_reason": f"{type(exc).__name__}: {exc}",
                        "validation_seed": validation_seed,
                        "validation_indices": ",".join(map(str, indices.tolist())),
                        "n_validation_included": int(np.sum(retained[indices])),
                        "covered_68": None, "h0_true": float(h0_true),
                        "n_catalogue": int(observed_z.size),
                        "n_detected": int(distances.size), "parent_size": int(n_parent),
                        "parent_frame_sha256": parent_hash,
                        "catalogue_sha256": catalogue_hash,
                        "gw_distance_sha256": gw_hash,
                    })
    return pd.DataFrame(rows)
