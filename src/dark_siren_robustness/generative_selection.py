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
    if not np.isfinite(floor) or not 0 < floor <= 1:
        raise ValueError("floor must be in (0, 1]")
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

def generate_synthetic_population(*, seed: int = 20261028, n_parent: int = 12000, n_detected: int = 250, h0_true: float = 70.0, sigma_mpc: float = 140.0, detection_threshold_mpc: float = 250.0, z_min: float = 0.01, z_max: float = 0.40) -> dict[str, object]:
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
    p_catalogue = catalogue_inclusion_probability(z_parent)
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
    max_attempts = max(100_000, n_detected * 2_000)
    while n_accepted < n_detected and attempts < max_attempts:
        batch = min(max(1000, 2 * (n_detected - n_accepted)), max_attempts - attempts)
        host_ids = rng.choice(n_parent, size=batch, replace=True, p=host_probability)
        host_z = z_parent[host_ids]
        measured_distance = rng.normal(C_KM_S * host_z / h0_true, sigma_mpc)
        detected = measured_distance >= detection_threshold_mpc
        if np.any(detected):
            ids, distances = host_ids[detected], measured_distance[detected]
            take = min(n_detected - n_accepted, len(ids))
            accepted_ids.append(ids[:take])
            accepted_distances.append(distances[:take])
            n_accepted += take
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
    config = {
        "experiment": "synthetic_catalogue_and_gw_selection_separated",
        "status": "synthetic_design_not_real_data_or_cosmological_measurement",
        "seed": int(seed), "n_parent": int(n_parent),
        "n_catalogue_observed": int(len(catalogue_observed)),
        "n_detected_events": int(n_detected),
        "h0_true_km_s_mpc": float(h0_true), "sigma_mpc": float(sigma_mpc),
        "detection_threshold_mpc": float(detection_threshold_mpc),
        "redshift_domain": [float(z_min), float(z_max)],
        "catalogue_selection": {"formula": "floor + (1-floor)/(1+exp((z-midpoint)/width))", "floor": 0.10, "midpoint": 0.22, "width": 0.05},
        "merger_weight": "1 + 1.5*(z-z_min)/(z_max-z_min)",
        "distance_model": "d_obs ~ Normal(c*z/H0_true, sigma_mpc^2); detected iff d_obs >= threshold",
        "observable_files": ["catalogue_observed.csv", "gw_events_observed.csv"],
        "truth_only_files": ["parent_truth.csv", "gw_event_truth.csv"],
        "scientific_limit": "Low-redshift Euclidean distance toy with fixed H0 and fixed Gaussian noise; not realistic GW inference or real survey data.",
    }
    return {"config": config, "parent_truth": parent_truth, "catalogue_observed": catalogue_observed, "gw_events_observed": gw_events_observed, "gw_event_truth": gw_event_truth}
