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


def _require_integer(name: str, value: object, minimum: int) -> int:
    """Validate an experiment-design integer without silently truncating inputs."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer")
    result = int(value)
    if result < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return result


def nested_srswor_sample_indices(
    rng: np.random.Generator,
    parent_size: int,
    sample_sizes: Sequence[int],
) -> dict[int, np.ndarray]:
    """Draw nested exact-SRSWOR samples from one random permutation.

    Each requested sample is a prefix of the same uniformly random
    permutation. Consequently every individual sample is an SRSWOR sample,
    and smaller samples are nested in larger ones for paired comparisons.
    """
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng must be a numpy.random.Generator")
    parent_size = _require_integer("parent_size", parent_size, 1)
    sizes_raw = tuple(sample_sizes)
    if not sizes_raw:
        raise ValueError("sample_sizes must not be empty")
    sizes = tuple(
        _require_integer(f"sample_sizes[{i}]", value, 1)
        for i, value in enumerate(sizes_raw)
    )
    if len(set(sizes)) != len(sizes):
        raise ValueError("sample_sizes must not contain duplicates")
    if max(sizes) > parent_size:
        raise ValueError("sample sizes cannot exceed parent_size")
    permutation = rng.permutation(parent_size)
    return {size: permutation[:size].copy() for size in sorted(sizes)}


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
    n_scenarios = _require_integer("n_scenarios", n_scenarios, 1)
    n_validation_draws = _require_integer("n_validation_draws", n_validation_draws, 2)
    n_parent = _require_integer("n_parent", n_parent, 100)
    n_detected = _require_integer("n_detected", n_detected, 1)
    seed = _require_integer("seed", seed, 0)
    h0_points = _require_integer("h0_points", h0_points, 31)
    z_points = _require_integer("z_points", z_points, 101)
    raw_sizes = tuple(validation_sizes)
    if not raw_sizes:
        raise ValueError("validation_sizes must contain at least one size")
    checked_sizes = tuple(
        _require_integer(f"validation_sizes[{index}]", value, 10)
        for index, value in enumerate(raw_sizes)
    )
    if len(set(checked_sizes)) != len(checked_sizes):
        raise ValueError("validation_sizes must not contain duplicates")
    sizes = tuple(sorted(checked_sizes))
    if max(sizes) > n_parent:
        raise ValueError("validation sizes must lie in [10, n_parent]")
    for name, value, minimum in (
        ("h0_true", h0_true, 0.0),
        ("sigma_mpc", sigma_mpc, 0.0),
        ("detection_threshold_mpc", detection_threshold_mpc, 0.0),
    ):
        if not np.isfinite(value) or value <= minimum:
            raise ValueError(f"{name} must be finite and positive")
    for name, value in (
        ("true_intercept", true_intercept),
        ("true_slope", true_slope),
        ("true_quadratic", true_quadratic),
    ):
        if not np.isfinite(value):
            raise ValueError(f"{name} must be finite")

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
                                 "validation_seed": validation_seed,
                                 "simulation_seed": sim_seed, "thinning_seed": thin_seed,
                                 "h0_true": float(h0_true), "sigma_mpc": float(sigma_mpc),
                                 "detection_threshold_mpc": float(detection_threshold_mpc),
                                 "true_intercept": float(true_intercept),
                                 "true_slope": float(true_slope),
                                 "true_quadratic": float(true_quadratic),
                                 "n_detected": int(distances.size),
                                 "n_catalogue": int(observed_z.size),
                                 "parent_size": int(n_parent),
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

        # The finite-parent reference is the same estimated-inclusion pipeline
        # applied to every paired record in this scenario's actual N-record
        # validation frame. It is distinct from the generating H0 truth and
        # from any pseudo-population reference used by a bootstrap.
        finite_parent_result = estimate_host_density_from_validation(
            z_grid, observed_z, parent_z, retained.astype(int),
            merger_weight=rate,
        )
        finite_parent_post = posterior_summary(
            h0_grid, selection_aware_toy_likelihood(
                z_grid, finite_parent_result["density"], h0_grid, distances,
                sigma=sigma_mpc, detection_threshold=detection_threshold_mpc,
            )["selected_log_likelihood"],
        )
        rng = np.random.default_rng(validation_seed)
        # One random permutation per draw creates nested validation samples:
        # every smaller-size sample is a prefix/subset of the larger-size sample.
        # This enables paired, higher-precision comparisons across validation sizes.
        for draw in range(n_validation_draws):
            nested_samples = nested_srswor_sample_indices(rng, n_parent, sizes)
            for size in sizes:
                indices = nested_samples[size]
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
                        "simulation_seed": sim_seed, "thinning_seed": thin_seed,
                        "true_intercept": float(true_intercept),
                        "true_slope": float(true_slope),
                        "true_quadratic": float(true_quadratic),
                        "sigma_mpc": float(sigma_mpc),
                        "detection_threshold_mpc": float(detection_threshold_mpc),
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
                        "finite_parent_target_posterior_mean": float(finite_parent_post["posterior_mean"]),
                        "finite_parent_target_posterior_median": float(finite_parent_post["posterior_median"]),
                        "finite_parent_target_ci68_low": float(finite_parent_post["ci68_low"]),
                        "finite_parent_target_ci68_high": float(finite_parent_post["ci68_high"]),
                        "error_vs_finite_parent_target": float(
                            post["posterior_mean"] - finite_parent_post["posterior_mean"]
                        ),
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
                        "simulation_seed": sim_seed, "thinning_seed": thin_seed,
                        "true_intercept": float(true_intercept),
                        "true_slope": float(true_slope),
                        "true_quadratic": float(true_quadratic),
                        "sigma_mpc": float(sigma_mpc),
                        "detection_threshold_mpc": float(detection_threshold_mpc),
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
