"""End-to-end synthetic population generator with separate catalogue/GW selection.

Controlled simulation primitive, not a real survey or detector model. Truth
tables are deliberately separated from observable tables.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

C_KM_S = 299_792.458

def latent_galaxy_intensity(redshift: np.ndarray) -> np.ndarray:
    """Positive baseline plus two smooth overdensity features (arbitrary units)."""
    z = np.asarray(redshift, dtype=float)
    if z.ndim != 1 or not np.all(np.isfinite(z)) or np.any(z <= 0):
        raise ValueError("redshift must be a finite, positive one-dimensional array")
    return 1.0 + 2.2 * np.exp(-0.5 * ((z - 0.105) / 0.025) ** 2) + 1.6 * np.exp(-0.5 * ((z - 0.275) / 0.040) ** 2)

def catalogue_inclusion_probability(redshift: np.ndarray, *, floor: float = 0.10, midpoint: float = 0.22, width: float = 0.05) -> np.ndarray:
    """Smooth redshift-dependent catalogue inclusion probability."""
    z = np.asarray(redshift, dtype=float)
    if z.ndim != 1 or not np.all(np.isfinite(z)) or np.any(z <= 0):
        raise ValueError("redshift must be a finite, positive one-dimensional array")
    if not np.isfinite(floor) or not 0 <= floor <= 1:
        raise ValueError("floor must be in [0, 1]")
    if not np.isfinite(midpoint) or midpoint <= 0:
        raise ValueError("midpoint must be finite and positive")
    if not np.isfinite(width) or width <= 0:
        raise ValueError("width must be finite and positive")
    return floor + (1.0 - floor) / (1.0 + np.exp((z - midpoint) / width))

def _sample_from_density(n: int, rng: np.random.Generator, *, z_min: float, z_max: float, grid_size: int = 20001) -> np.ndarray:
    if n < 1:
        raise ValueError("n must be positive")
    grid = np.linspace(z_min, z_max, grid_size)
    intensity = latent_galaxy_intensity(grid)
    increments = 0.5 * (intensity[:-1] + intensity[1:]) * np.diff(grid)
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    cdf /= cdf[-1]
    return np.interp(rng.random(n), cdf, grid)

def generate_synthetic_population(*, seed: int = 20261028, n_parent: int = 12000, n_detected: int = 250, h0_true: float = 70.0, sigma_mpc: float = 140.0, detection_threshold_mpc: float = 250.0, z_min: float = 0.01, z_max: float = 0.40, catalogue_floor: float = 0.10, catalogue_midpoint: float = 0.22, catalogue_width: float = 0.05) -> dict[str, object]:
    """Generate parent galaxies, an incomplete catalogue, and detected GW distances.

    catalogue_observed and gw_events_observed are observable tables. parent_truth
    and gw_event_truth are scoring tables and must not enter non-oracle inference.
    Hosts are sampled from the complete parent population, not the thinned catalogue.
    """
    if not isinstance(seed, (int, np.integer)) or isinstance(seed, (bool, np.bool_)) or seed < 0:
        raise ValueError("seed must be a non-negative integer")
    if n_parent < 100 or n_detected < 1:
        raise ValueError("n_parent must be >= 100 and n_detected must be positive")
    if not np.isfinite(h0_true) or h0_true <= 0:
        raise ValueError("h0_true must be finite and positive")
    if not np.isfinite(sigma_mpc) or sigma_mpc <= 0:
        raise ValueError("sigma_mpc must be finite and positive")
    if not np.isfinite(detection_threshold_mpc):
        raise ValueError("detection_threshold_mpc must be finite")
    if not (0 < z_min < z_max):
        raise ValueError("redshift bounds must satisfy 0 < z_min < z_max")
    rng = np.random.default_rng(int(seed))
    z_parent = _sample_from_density(n_parent, rng, z_min=z_min, z_max=z_max)
    p_catalogue = catalogue_inclusion_probability(
        z_parent, floor=catalogue_floor,
        midpoint=catalogue_midpoint, width=catalogue_width,
    )
    retained = rng.random(n_parent) < p_catalogue
    parent_truth = pd.DataFrame({
        "object_id": np.arange(n_parent, dtype=np.int64),
        "redshift": z_parent,
        "inclusion_probability": p_catalogue,
        "retained": retained.astype(np.int8),
    })
    catalogue_observed = parent_truth.loc[retained, ["object_id", "redshift"]].reset_index(drop=True)
    if len(catalogue_observed) < 2:
        raise RuntimeError("Catalogue thinning retained too few galaxies for this realization")
    host_weight = 1.0 + 1.5 * (z_parent - z_min) / (z_max - z_min)
    host_probability = host_weight / host_weight.sum()
    accepted_ids: list[np.ndarray] = []
    accepted_distances: list[np.ndarray] = []
    n_accepted = 0
    attempts = 0
    n_detected_candidates = 0
    n_rejected_candidates = 0
    n_excess_detected_candidates_discarded = 0
    generation_trace_parts: list[pd.DataFrame] = []
    max_attempts = max(100_000, n_detected * 2_000)
    while n_accepted < n_detected and attempts < max_attempts:
        batch = min(max(1000, 2 * (n_detected - n_accepted)), max_attempts - attempts)
        host_ids = rng.choice(n_parent, size=batch, replace=True, p=host_probability)
        host_z = z_parent[host_ids]
        measured_distance = rng.normal(C_KM_S * host_z / h0_true, sigma_mpc)
        detected = measured_distance >= detection_threshold_mpc
        detected_count = int(np.count_nonzero(detected))
        n_detected_candidates += detected_count
        n_rejected_candidates += batch - detected_count
        selected_as_event = np.zeros(batch, dtype=np.int8)
        if detected_count:
            detected_positions = np.flatnonzero(detected)
            ids, distances = host_ids[detected], measured_distance[detected]
            take = min(n_detected - n_accepted, len(ids))
            accepted_ids.append(ids[:take])
            accepted_distances.append(distances[:take])
            selected_as_event[detected_positions[:take]] = 1
            n_accepted += take
            n_excess_detected_candidates_discarded += len(ids) - take
        generation_trace_parts.append(pd.DataFrame({
            "attempt_id": np.arange(attempts, attempts + batch, dtype=np.int64),
            "host_object_id": host_ids.astype(np.int64),
            "distance_mpc": measured_distance.astype(float),
            "detected": detected.astype(np.int8),
            "selected_as_event": selected_as_event,
        }))
        attempts += batch
    if n_accepted != n_detected:
        raise RuntimeError(f"Only generated {n_accepted}/{n_detected} detected events within {attempts} attempts")
    host_ids = np.concatenate(accepted_ids)
    distances = np.concatenate(accepted_distances)
    host_z = z_parent[host_ids]
    gw_events_observed = pd.DataFrame({
        "event_id": np.arange(n_detected, dtype=np.int64),
        "distance_mpc": distances,
        "detected": np.ones(n_detected, dtype=np.int8),
    })
    gw_event_truth = pd.DataFrame({
        "event_id": np.arange(n_detected, dtype=np.int64),
        "host_object_id": host_ids.astype(np.int64),
        "host_redshift": host_z,
        "latent_mean_distance_mpc": C_KM_S * host_z / h0_true,
    })
    gw_generation_trace_truth = pd.concat(generation_trace_parts, ignore_index=True)
    config = {
        "experiment": "synthetic_catalogue_and_gw_selection_separated",
        "status": "synthetic_design_not_real_data_or_cosmological_measurement",
        "seed": int(seed), "n_parent": int(n_parent),
        "n_catalogue_observed": int(len(catalogue_observed)),
        "n_detected_events": int(n_detected),
        "n_generation_attempts": int(attempts),
        "n_detected_candidates": int(n_detected_candidates),
        "n_rejected_candidates": int(n_rejected_candidates),
        "n_excess_detected_candidates_discarded": int(n_excess_detected_candidates_discarded),
        "h0_true_km_s_mpc": float(h0_true), "sigma_mpc": float(sigma_mpc),
        "detection_threshold_mpc": float(detection_threshold_mpc),
        "redshift_domain": [float(z_min), float(z_max)],
        "catalogue_selection": {
            "formula": "floor + (1-floor)/(1+exp((z-midpoint)/width))",
            "floor": float(catalogue_floor),
            "midpoint": float(catalogue_midpoint),
            "width": float(catalogue_width),
        },
        "merger_weight": "1 + 1.5*(z-z_min)/(z_max-z_min)",
        "distance_model": "d_obs ~ Normal(c*z/H0_true, sigma_mpc^2); detected iff d_obs >= threshold",
        "observable_files": ["observable/catalogue_observed.csv", "observable/gw_events_observed.csv"],
        "truth_only_files": [
            "truth/parent_truth.csv", "truth/gw_event_truth.csv",
            "truth/gw_generation_trace_truth.csv",
        ],
        "scientific_limit": "Low-redshift Euclidean distance toy with fixed H0 and fixed Gaussian noise; not realistic GW inference or real survey data.",
    }
    return {
        "config": config,
        "parent_truth": parent_truth,
        "catalogue_observed": catalogue_observed,
        "gw_events_observed": gw_events_observed,
        "gw_event_truth": gw_event_truth,
        "gw_generation_trace_truth": gw_generation_trace_truth,
    }


def generate_independent_validation_sample(
    *,
    seed: int,
    n_validation_frame: int = 4000,
    n_validation: int = 300,
    z_min: float = 0.01,
    z_max: float = 0.40,
    catalogue_floor: float = 0.10,
    catalogue_midpoint: float = 0.22,
    catalogue_width: float = 0.05,
) -> dict[str, object]:
    """Draw an SRS validation sample from a separately seeded parent frame.

    The returned observable validation table contains only redshift and the
    Bernoulli catalogue-retention outcome. Inclusion probabilities and the
    unobserved frame remain in a separate truth table for simulation auditing.
    This frame is independent of the analysis-population generator when a
    distinct seed is supplied.
    """
    if not isinstance(seed, (int, np.integer)) or isinstance(seed, (bool, np.bool_)) or seed < 0:
        raise ValueError("seed must be a non-negative integer")
    if n_validation_frame < 5 or n_validation < 5 or n_validation > n_validation_frame:
        raise ValueError("validation sizes must satisfy 5 <= n_validation <= n_validation_frame")
    if not (0 < z_min < z_max):
        raise ValueError("redshift bounds must satisfy 0 < z_min < z_max")
    rng = np.random.default_rng(int(seed))
    z_frame = _sample_from_density(
        n_validation_frame, rng, z_min=z_min, z_max=z_max
    )
    probability = catalogue_inclusion_probability(
        z_frame, floor=catalogue_floor,
        midpoint=catalogue_midpoint, width=catalogue_width,
    )
    retained = (rng.random(n_validation_frame) < probability).astype(np.int8)
    # A single random permutation makes samples nested across requested sizes
    # when the same seed and frame size are reused for paired validation studies.
    selected = rng.permutation(n_validation_frame)[:n_validation]
    observable = pd.DataFrame({
        "validation_id": selected.astype(np.int64),
        "redshift": z_frame[selected],
        "retained": retained[selected],
    })
    truth = pd.DataFrame({
        "frame_id": np.arange(n_validation_frame, dtype=np.int64),
        "redshift": z_frame,
        "inclusion_probability": probability,
        "retained": retained,
        "selected_for_validation": np.isin(
            np.arange(n_validation_frame), selected
        ).astype(np.int8),
    })
    config = {
        "seed": int(seed),
        "n_validation_frame": int(n_validation_frame),
        "n_validation": int(n_validation),
        "sampling_design": "simple_random_sample_without_replacement_from_independent_parent_frame",
        "validation_observable_columns": ["validation_id", "redshift", "retained"],
        "truth_only_columns": [
            "frame_id", "redshift", "inclusion_probability",
            "retained", "selected_for_validation",
        ],
        "independent_frame": True,
        "catalogue_selection": {
            "floor": float(catalogue_floor),
            "midpoint": float(catalogue_midpoint),
            "width": float(catalogue_width),
        },
    }
    return {
        "config": config,
        "validation_observed": observable,
        "validation_frame_truth": truth,
    }

def quadratic_logit_parent_measure(*, z_min: float = 0.01, z_max: float = 0.40,
                                   coordinate_center: float = 0.205,
                                   coordinate_scale: float = 0.195,
                                   grid_size: int = 20001) -> dict[str, float]:
    """Return the parent-intensity centering and normalization constants.

    The parent sampling density is proportional to latent_galaxy_intensity;
    merger/host weights are deliberately not included in this measure.
    """
    from scipy.integrate import trapezoid

    if not (0 < z_min < z_max) or coordinate_scale <= 0 or grid_size < 1001:
        raise ValueError("invalid redshift bounds, coordinate scale, or grid size")
    grid = np.linspace(z_min, z_max, int(grid_size))
    intensity = latent_galaxy_intensity(grid)
    normalizer = float(trapezoid(intensity, grid))
    x = (grid - coordinate_center) / coordinate_scale
    x2_center = float(trapezoid(intensity * x**2, grid) / normalizer)
    return {"x2_center": x2_center, "intensity_normalizer": normalizer}


def calibrate_quadratic_logit_intercept(*, gamma: float, beta: float = -3.9,
                                        target_rate: float | None = None,
                                        alpha_reference: float = 0.0,
                                        z_min: float = 0.01, z_max: float = 0.40,
                                        coordinate_center: float = 0.205,
                                        coordinate_scale: float = 0.195,
                                        grid_size: int = 20001) -> dict[str, float]:
    """Calibrate alpha_gamma to preserve the gamma=0 parent-averaged rate.

    If target_rate is omitted, it is defined by alpha_reference and gamma=0
    under the actual parent intensity measure. This avoids importing a rate
    from a different parent distribution or an undocumented numerical grid.
    """
    from scipy.integrate import trapezoid
    from scipy.optimize import brentq
    from scipy.special import expit

    if not np.isfinite(gamma) or not np.isfinite(beta) or not np.isfinite(alpha_reference):
        raise ValueError("gamma, beta, and reference intercept must be finite")
    constants = quadratic_logit_parent_measure(
        z_min=z_min, z_max=z_max, coordinate_center=coordinate_center,
        coordinate_scale=coordinate_scale, grid_size=grid_size,
    )
    grid = np.linspace(z_min, z_max, int(grid_size))
    intensity = latent_galaxy_intensity(grid)
    x = (grid - coordinate_center) / coordinate_scale
    centered_x2 = x**2 - constants["x2_center"]
    normalizer = constants["intensity_normalizer"]
    if target_rate is None:
        target_rate = float(trapezoid(
            intensity * expit(alpha_reference + beta * x), grid
        ) / normalizer)
    if not np.isfinite(target_rate) or not 0 < target_rate < 1:
        raise ValueError("target_rate must be finite and strictly between zero and one")
    if gamma == 0.0:
        alpha = float(alpha_reference)
    else:
        def rate_difference(candidate: float) -> float:
            probability = expit(candidate + beta * x + gamma * centered_x2)
            return float(trapezoid(intensity * probability, grid) / normalizer - target_rate)
        alpha = float(brentq(rate_difference, -40.0, 40.0, xtol=1e-13, rtol=1e-13))
    achieved = float(trapezoid(
        intensity * expit(alpha + beta * x + gamma * centered_x2), grid
    ) / normalizer)
    return {
        "alpha": alpha, "beta": float(beta), "gamma": float(gamma),
        "target_rate": float(target_rate), "achieved_rate": achieved,
        "x2_center": constants["x2_center"],
        "coordinate_center": float(coordinate_center),
        "coordinate_scale": float(coordinate_scale),
    }


def quadratic_logit_inclusion_probability(
    redshift: np.ndarray, *, gamma: float, beta: float = -3.9,
    alpha: float | None = None, target_rate: float | None = None,
    z_min: float = 0.01, z_max: float = 0.40,
    coordinate_center: float = 0.205, coordinate_scale: float = 0.195,
) -> np.ndarray:
    """Quadratic-logit catalogue inclusion probability for curvature stress tests."""
    from scipy.special import expit

    z = np.asarray(redshift, dtype=float)
    if z.ndim != 1 or not np.all(np.isfinite(z)) or np.any(z <= 0):
        raise ValueError("redshift must be a finite, positive one-dimensional array")
    if np.any((z < z_min) | (z > z_max)):
        raise ValueError("redshift values must lie within the declared parent domain")
    calibration = calibrate_quadratic_logit_intercept(
        gamma=gamma, beta=beta, target_rate=target_rate,
        z_min=z_min, z_max=z_max, coordinate_center=coordinate_center,
        coordinate_scale=coordinate_scale,
    )
    alpha_use = calibration["alpha"] if alpha is None else float(alpha)
    if not np.isfinite(alpha_use):
        raise ValueError("alpha must be finite")
    x = (z - coordinate_center) / coordinate_scale
    return expit(alpha_use + beta * x + gamma * (x**2 - calibration["x2_center"]))


def generate_quadratic_logit_inclusion_samples(
    *, population_seed: int, validation_seed: int, n_parent: int = 4000,
    n_validation_frame: int = 1500, n_validation: int = 1000,
    gamma: float = 0.0, beta: float = -3.9,
    z_min: float = 0.01, z_max: float = 0.40,
    coordinate_center: float = 0.205, coordinate_scale: float = 0.195,
) -> dict[str, object]:
    """Generate paired parent-catalogue and independent validation samples.

    The observable validation table contains only redshift and retention.
    Inclusion probabilities and the unobserved validation frame are truth-only.
    """
    for name, seed in (("population_seed", population_seed), ("validation_seed", validation_seed)):
        if not isinstance(seed, (int, np.integer)) or isinstance(seed, (bool, np.bool_)) or seed < 0:
            raise ValueError(f"{name} must be a non-negative integer")
    if n_parent < 100 or n_validation_frame < 5 or not 5 <= n_validation <= n_validation_frame:
        raise ValueError("invalid parent or validation sample sizes")
    if not (0 < z_min < z_max):
        raise ValueError("invalid redshift bounds")

    calibration = calibrate_quadratic_logit_intercept(
        gamma=gamma, beta=beta, z_min=z_min, z_max=z_max,
        coordinate_center=coordinate_center, coordinate_scale=coordinate_scale,
    )
    rng_parent = np.random.default_rng(int(population_seed))
    z_parent = _sample_from_density(n_parent, rng_parent, z_min=z_min, z_max=z_max)
    p_parent = quadratic_logit_inclusion_probability(
        z_parent, gamma=gamma, beta=beta, alpha=calibration["alpha"],
        target_rate=calibration["target_rate"], z_min=z_min, z_max=z_max,
        coordinate_center=coordinate_center, coordinate_scale=coordinate_scale,
    )
    retained_parent = (rng_parent.random(n_parent) < p_parent).astype(np.int8)
    parent_truth = pd.DataFrame({
        "object_id": np.arange(n_parent, dtype=np.int64),
        "redshift": z_parent, "inclusion_probability": p_parent,
        "retained": retained_parent,
    })
    catalogue_observed = parent_truth.loc[
        retained_parent.astype(bool), ["object_id", "redshift"]
    ].reset_index(drop=True)

    rng_validation = np.random.default_rng(int(validation_seed))
    z_frame = _sample_from_density(
        n_validation_frame, rng_validation, z_min=z_min, z_max=z_max
    )
    p_frame = quadratic_logit_inclusion_probability(
        z_frame, gamma=gamma, beta=beta, alpha=calibration["alpha"],
        target_rate=calibration["target_rate"], z_min=z_min, z_max=z_max,
        coordinate_center=coordinate_center, coordinate_scale=coordinate_scale,
    )
    retained_frame = (rng_validation.random(n_validation_frame) < p_frame).astype(np.int8)
    selected = rng_validation.permutation(n_validation_frame)[:n_validation]
    validation_observed = pd.DataFrame({
        "validation_id": selected.astype(np.int64),
        "redshift": z_frame[selected], "retained": retained_frame[selected],
    })
    validation_frame_truth = pd.DataFrame({
        "frame_id": np.arange(n_validation_frame, dtype=np.int64),
        "redshift": z_frame, "inclusion_probability": p_frame,
        "retained": retained_frame,
        "selected_for_validation": np.isin(
            np.arange(n_validation_frame), selected
        ).astype(np.int8),
    })
    config = {
        "population_seed": int(population_seed), "validation_seed": int(validation_seed),
        "n_parent": int(n_parent), "n_validation_frame": int(n_validation_frame),
        "n_validation": int(n_validation), "sampling_design": "independent_parent_frames_with_nested_srs_validation_prefix",
        "catalogue_selection_model": "quadratic_logit",
        "validation_observable_columns": ["validation_id", "redshift", "retained"],
        "truth_only_columns": ["frame_id", "redshift", "inclusion_probability", "retained", "selected_for_validation"],
        **calibration,
    }
    return {
        "config": config, "parent_truth": parent_truth,
        "catalogue_observed": catalogue_observed,
        "validation_observed": validation_observed,
        "validation_frame_truth": validation_frame_truth,
    }

