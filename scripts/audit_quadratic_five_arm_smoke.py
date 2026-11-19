"""Independent raw-input and likelihood audit for quadratic-logit five-arm smoke artifacts.

This auditor does not import the five-arm runner, weight builder, likelihood
implementation, or H0 posterior normalizer. It reconstructs weights and the
selected likelihood directly from saved observable/truth inputs and persisted
validation-posterior nodes. It audits likelihood propagation, not the posterior
quadrature construction itself.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from itertools import product
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import root
import pandas as pd
from scipy.special import expit, log_ndtr, logsumexp, roots_hermitenorm
from scipy.stats import gaussian_kde, norm

C_KM_S = 299_792.458
ARM_KEYS = {
    "A_catalogue_only": "A_catalogue_only",
    "B_parent_oracle": "B_parent_oracle",
    "C_known_probability_ipw": "C_known_probability_ipw",
    "D_map_plugin_ipw": "D_map_plugin_ipw",
    "E_uncertainty_integrated_ipw": "E_uncertainty_integrated_ipw",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _assert_close(actual, expected, label: str, *, atol: float = 2e-10) -> None:
    if not np.allclose(actual, expected, rtol=2e-10, atol=atol, equal_nan=False):
        difference = float(np.max(np.abs(np.asarray(actual) - np.asarray(expected))))
        raise AssertionError(f"{label} differs; maximum absolute difference={difference}")


def _independent_density(
    z_grid: np.ndarray, sample: np.ndarray, raw_weights: np.ndarray,
    rate: np.ndarray, scott_multiplier: float,
) -> np.ndarray:
    """Independent reflected weighted KDE reconstruction; no production imports."""
    raw = np.asarray(raw_weights, dtype=float)
    if raw.shape != sample.shape or not np.all(np.isfinite(raw)) or np.any(raw <= 0):
        raise AssertionError("independent raw host weights are invalid")
    scaled = raw / np.max(raw)
    reflected = np.concatenate((
        sample, 2.0 * z_grid[0] - sample, 2.0 * z_grid[-1] - sample,
    ))
    reflected_weights = np.tile(scaled / 3.0, 3)
    original_probability = scaled / np.sum(scaled)
    original_denominator = 1.0 - float(original_probability @ original_probability)
    original_mean = float(original_probability @ sample)
    original_variance = float(
        np.sum(original_probability * (sample - original_mean) ** 2)
        / original_denominator
    )
    augmented_probability = reflected_weights / np.sum(reflected_weights)
    augmented_denominator = 1.0 - float(augmented_probability @ augmented_probability)
    augmented_mean = float(augmented_probability @ reflected)
    augmented_variance = float(
        np.sum(augmented_probability * (reflected - augmented_mean) ** 2)
        / augmented_denominator
    )
    if original_variance <= 0 or augmented_variance <= 0:
        raise AssertionError("independent weighted covariance is invalid")
    neff = float(np.sum(scaled) ** 2 / np.sum(scaled ** 2))
    covariance_adjustment = np.sqrt(original_variance / augmented_variance)
    factor = float(scott_multiplier) * neff ** (-1.0 / 5.0) * covariance_adjustment
    estimator = gaussian_kde(
        reflected, weights=reflected_weights, bw_method=lambda _kde: factor,
    )
    unnormalized = estimator(z_grid) * rate
    area = float(trapezoid(unnormalized, z_grid))
    if not np.isfinite(area) or area <= 0:
        raise AssertionError("independent KDE has invalid normalization")
    density = unnormalized / area
    if not np.isclose(trapezoid(density, z_grid), 1.0, rtol=1e-10, atol=1e-10):
        raise AssertionError("independent KDE failed unit-integral check")
    return density



def _independent_bandwidth(
    sample: np.ndarray, raw_weights: np.ndarray,
    lower: float, upper: float, multiplier: float,
) -> float:
    raw = np.asarray(raw_weights, dtype=float)
    if raw.shape != sample.shape or np.any(raw <= 0) or not np.all(np.isfinite(raw)):
        raise AssertionError("invalid weights for independent bandwidth reconstruction")
    scaled = raw / np.max(raw)
    probability = scaled / np.sum(scaled)
    denominator = 1.0 - float(probability @ probability)
    mean = float(probability @ sample)
    variance = float(np.sum(probability * (sample - mean) ** 2) / denominator)
    neff = float(np.sum(scaled) ** 2 / np.sum(scaled ** 2))
    width = float(multiplier) * neff ** (-1.0 / 5.0) * np.sqrt(variance)
    if not np.isfinite(width) or width <= 0:
        raise AssertionError("independent physical bandwidth is invalid")
    return width


def _selected_kde_log_likelihood(
    h0_grid: np.ndarray, z_grid: np.ndarray, density: np.ndarray,
    distances: np.ndarray, sigma: float, threshold: float,
) -> tuple[np.ndarray, np.ndarray]:
    means = C_KM_S * z_grid[None, :] / h0_grid[:, None]
    event_log_mix, event_log_beta = [], []
    for distance in distances:
        numerator = density[None, :] * norm.pdf(distance, loc=means, scale=sigma)
        efficiency = density[None, :] * np.exp(log_ndtr((means - threshold) / sigma))
        event_log_mix.append(np.log(trapezoid(numerator, z_grid, axis=1)))
        event_log_beta.append(np.log(trapezoid(efficiency, z_grid, axis=1)))
    event_log_mix = np.asarray(event_log_mix)
    event_log_beta = np.asarray(event_log_beta)
    return np.sum(event_log_mix - event_log_beta, axis=0), event_log_beta

def _posterior_quantile(grid: np.ndarray, log_likelihood: np.ndarray, q: float) -> float:
    shifted = np.exp(log_likelihood - np.max(log_likelihood))
    area = float(trapezoid(shifted, grid))
    density = shifted / area
    increments = 0.5 * (density[:-1] + density[1:]) * np.diff(grid)
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    cdf /= cdf[-1]
    return float(np.interp(q, cdf, grid))



def _independent_inclusion_quadrature(
    validation_z: np.ndarray, retained: np.ndarray, *, order: int,
    center: float, scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Rebuild the proper logistic posterior and corrected GH rule independently.

    Uses a damped Newton optimizer rather than the production trust-exact
    optimizer, then reconstructs the Laplace proposal and target/proposal
    importance correction from the Bernoulli likelihood and Gaussian prior.
    """
    z = np.asarray(validation_z, dtype=float)
    y = np.asarray(retained, dtype=float)
    x = (z - center) / scale
    design = np.column_stack((np.ones(z.size), x))
    prior_mean = np.zeros(2)
    prior_sd = np.array([5.0, 2.5])
    precision = 1.0 / prior_sd**2

    def objective(theta):
        eta = design @ theta
        return float(
            np.sum(np.logaddexp(0.0, eta) - y * eta)
            + 0.5 * np.sum(((theta - prior_mean) / prior_sd) ** 2)
        )

    def gradient(theta):
        p = expit(design @ theta)
        return design.T @ (p - y) + (theta - prior_mean) * precision

    def information(theta):
        p = expit(design @ theta)
        return design.T @ ((p * (1.0 - p))[:, None] * design) + np.diag(precision)

    mode = np.zeros(2)
    for _ in range(100):
        score = gradient(mode)
        if np.linalg.norm(score, ord=np.inf) < 2e-11:
            break
        hessian = information(mode)
        direction = np.linalg.solve(hessian, score)
        decrement = float(score @ direction)
        old_value = objective(mode)
        step = 1.0
        while step > 2.0**-35:
            candidate = mode - step * direction
            if objective(candidate) <= old_value - 1e-4 * step * decrement:
                mode = candidate
                break
            step *= 0.5
        else:
            raise AssertionError("independent logistic Newton solver could not descend")
    residual = float(np.linalg.norm(gradient(mode), ord=np.inf))
    if residual > 1e-8:
        # Rare ill-conditioned validation frames can stall the damped Newton
        # line search despite a well-defined, strictly convex posterior. Fall
        # back to an independent score-equation root solve; this is still
        # separate from the production trust-exact optimizer.
        solved = root(
            gradient, mode, jac=information, method="hybr",
            options={"xtol": 1e-12},
        )
        candidate_residual = float(np.linalg.norm(gradient(solved.x), ord=np.inf))
        if np.all(np.isfinite(solved.x)) and candidate_residual <= 1e-8:
            mode = np.asarray(solved.x, dtype=float)
            residual = candidate_residual
        else:
            raise AssertionError(
                "independent logistic Newton/root solver did not converge; "
                f"Newton residual={residual:.6g}, root residual={candidate_residual:.6g}, "
                f"root status={solved.message}"
            )

    covariance = np.linalg.inv(information(mode))
    chol = np.linalg.cholesky(covariance)
    gh_nodes, gh_weights = roots_hermitenorm(order)
    indices = list(product(range(order), repeat=2))
    standard = np.asarray([[gh_nodes[i] for i in pair] for pair in indices], dtype=float)
    base_log_weights = np.asarray([
        sum(math.log(gh_weights[i]) for i in pair) - math.log(2.0 * math.pi)
        for pair in indices
    ])
    nodes = mode[None, :] + standard @ chol.T
    log_proposal = (
        -math.log(2.0 * math.pi)
        -0.5 * np.sum(standard**2, axis=1)
        -float(np.sum(np.log(np.diag(chol))))
    )
    eta = design @ nodes.T
    log_likelihood = np.sum(y[:, None] * eta - np.logaddexp(0.0, eta), axis=0)
    log_prior = (
        -0.5 * np.sum(((nodes - prior_mean) / prior_sd) ** 2, axis=1)
        -float(np.sum(np.log(prior_sd)))
        -math.log(2.0 * math.pi)
    )
    log_importance = base_log_weights + log_likelihood + log_prior - log_proposal
    weights = np.exp(log_importance - logsumexp(log_importance))
    if not np.all(np.isfinite(weights)) or not np.isclose(weights.sum(), 1.0, atol=1e-12):
        raise AssertionError("independent logistic quadrature weights are invalid")
    return mode, covariance, nodes, weights


def audit(root: Path) -> dict:
    root = root.resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    checked_files = 0
    for relative, metadata in manifest["files"].items():
        path = root / relative
        if not path.is_file():
            raise AssertionError(f"manifest file is missing: {relative}")
        if path.stat().st_size != int(metadata["bytes"]):
            raise AssertionError(f"file size mismatch: {relative}")
        if _sha256(path) != metadata["sha256"]:
            raise AssertionError(f"SHA-256 mismatch: {relative}")
        checked_files += 1
    if checked_files != int(manifest["file_count"]):
        raise AssertionError("manifest file_count does not match entries")

    config = json.loads((root / "config.json").read_text(encoding="utf-8"))
    if config["status"] != "engineering_smoke_not_confirmatory_not_observational":
        raise AssertionError("artifact claim-status label is missing or unexpected")
    floor = float(config["catalogue_floor"])
    midpoint, width = 0.22, 0.05
    center, scale = 0.205, 0.195
    catalogue_model = config.get("catalogue_model", "floor_logistic")
    quadratic_calibration = config.get("quadratic_calibration", {})
    threshold = float(config["detection_threshold_mpc"])
    sigma = float(config["sigma_mpc"])
    h0_true = float(config["h0_true"])
    prior_min, prior_max = map(float, config.get("h0_prior_support", [50.0, 90.0]))
    if not np.isfinite(prior_min) or not np.isfinite(prior_max) or prior_min >= prior_max:
        raise AssertionError("H0 prior support must be finite and increasing")
    if not prior_min <= h0_true <= prior_max:
        raise AssertionError("H0 truth lies outside the declared prior support")

    parent = pd.read_csv(root / "truth" / "parent_truth.csv")
    catalogue = pd.read_csv(root / "observable" / "catalogue_observed.csv")
    events = pd.read_csv(root / "observable" / "gw_events_observed.csv")
    validation_truth = pd.read_csv(root / "truth" / "validation_frame_truth.csv")
    event_truth = pd.read_csv(root / "truth" / "gw_event_truth.csv")
    generation_trace = pd.read_csv(root / "truth" / "gw_generation_trace_truth.csv")
    selected_trace = generation_trace.loc[
        generation_trace["selected_as_event"] == 1
    ].sort_values("attempt_id")
    if len(events) != int(config["n_detected"]) or len(event_truth) != len(events):
        raise AssertionError("detected-event observable/truth row counts differ")
    if np.any(events["distance_mpc"].to_numpy(float) < threshold):
        raise AssertionError("an observed event fails the declared GW detection threshold")
    if len(catalogue) != int(config["n_catalogue_observed"]):
        raise AssertionError("catalogue row count differs from config")
    expected_catalogue_ids = set(parent.loc[parent["retained"] == 1, "object_id"].astype(int))
    if set(catalogue["object_id"].astype(int)) != expected_catalogue_ids:
        raise AssertionError("observable catalogue does not match retained parent IDs")
    if not np.array_equal(events["event_id"].to_numpy(int), event_truth["event_id"].to_numpy(int)):
        raise AssertionError("GW event observable/truth IDs differ")
    if len(selected_trace) != len(events):
        raise AssertionError("selected generation-trace rows do not match detected-event count")
    _assert_close(
        events["distance_mpc"].to_numpy(float),
        selected_trace["distance_mpc"].to_numpy(float),
        "observed GW distances versus selected generation trace",
    )
    if not np.array_equal(
        event_truth["host_object_id"].to_numpy(int),
        selected_trace["host_object_id"].to_numpy(int),
    ):
        raise AssertionError("GW event host IDs disagree with selected generation trace")
    event_join = event_truth.merge(
        parent[["object_id", "redshift"]],
        left_on="host_object_id", right_on="object_id",
        suffixes=("_event", "_parent"), validate="many_to_one",
    )
    _assert_close(
        event_join["host_redshift"].to_numpy(float),
        event_join["redshift"].to_numpy(float),
        "GW host-redshift truth join",
    )
    _assert_close(
        event_join["latent_mean_distance_mpc"].to_numpy(float),
        C_KM_S * event_join["redshift"].to_numpy(float) / h0_true,
        "GW latent distance truth",
    )
    joined = catalogue.merge(
        parent[["object_id", "redshift", "inclusion_probability"]],
        on="object_id", suffixes=("_observed", "_truth"), validate="one_to_one",
    )
    _assert_close(
        joined["redshift_observed"].to_numpy(float),
        joined["redshift_truth"].to_numpy(float),
        "catalogue redshift join",
    )
    catalogue_z = joined["redshift_observed"].to_numpy(float)
    if catalogue_model == "quadratic_logit":
        alpha = float(quadratic_calibration["alpha"])
        beta = float(quadratic_calibration["beta"])
        gamma = float(quadratic_calibration["gamma"])
        x2_center = float(quadratic_calibration["x2_center"])
        # Independently reconstruct the parent-intensity centering and the
        # calibrated parent-average retention rate used by the quadratic DGP.
        z_cal = np.linspace(0.01, 0.40, 20001)
        intensity = (1.0 + 2.2*np.exp(-0.5*((z_cal-0.105)/0.025)**2)
                     + 1.6*np.exp(-0.5*((z_cal-0.275)/0.040)**2))
        x_cal = (z_cal-center)/scale
        intensity_norm = float(trapezoid(intensity,z_cal))
        rebuilt_x2 = float(trapezoid(intensity*x_cal**2,z_cal)/intensity_norm)
        _assert_close(x2_center, rebuilt_x2, "quadratic parent-law x2 centering")
        rebuilt_rate = float(trapezoid(
            intensity*expit(alpha+beta*x_cal+gamma*(x_cal**2-x2_center)),z_cal
        )/intensity_norm)
        _assert_close(rebuilt_rate,float(quadratic_calibration["target_rate"]),
                      "quadratic parent-average retention calibration")
        _assert_close(rebuilt_rate,float(quadratic_calibration["achieved_rate"]),
                      "quadratic achieved retention rate")
        x_catalogue = (catalogue_z-center)/scale
        true_probability = expit(alpha+beta*x_catalogue+gamma*(x_catalogue**2-x2_center))
        x_parent = (parent["redshift"].to_numpy(float)-center)/scale
        expected_parent_probability = expit(alpha+beta*x_parent+gamma*(x_parent**2-x2_center))
        _assert_close(parent["inclusion_probability"].to_numpy(float),
                      expected_parent_probability,"quadratic parent inclusion probabilities")
        x_frame=(validation_truth["redshift"].to_numpy(float)-center)/scale
        expected_frame_probability=expit(alpha+beta*x_frame+gamma*(x_frame**2-x2_center))
        _assert_close(validation_truth["inclusion_probability"].to_numpy(float),
                      expected_frame_probability,"quadratic validation inclusion probabilities")
    elif catalogue_model == "floor_logistic":
        true_probability = floor + (1.0 - floor) / (
            1.0 + np.exp((catalogue_z - midpoint) / width)
        )
    else:
        raise AssertionError(f"unknown catalogue model: {catalogue_model}")
    _assert_close(
        joined["inclusion_probability"].to_numpy(float),
        true_probability,
        "catalogue truth inclusion probabilities",
    )
    parent_z = parent["redshift"].to_numpy(float)
    parent_base = 1.0 + 1.5 * (parent_z - 0.01) / 0.39
    catalogue_base = 1.0 + 1.5 * (catalogue_z - 0.01) / 0.39
    distances = events["distance_mpc"].to_numpy(float)
    grid = np.linspace(prior_min, prior_max, int(config["h0_grid_points"]))
    z_grid = np.linspace(0.01, 0.40, int(config["redshift_grid_points"]))
    rate_grid = np.ones_like(z_grid)
    scott_multiplier = float(config["scott_multiplier"])
    archive = np.load(root / "results" / "log_likelihoods.npz")
    _assert_close(archive["H0_grid"], grid, "H0 grid")
    _assert_close(archive["redshift_grid"], z_grid, "redshift grid")
    _assert_close(archive["merger_weight_grid"], rate_grid, "merger weight grid")
    stored_summary = pd.read_csv(root / "results" / "posterior_summaries.csv")
    stored_bandwidths = pd.read_csv(root / "results" / "realized_bandwidths.csv")

    validation_sizes = [int(value) for value in config["validation_sizes"]]
    expected_bandwidth_rows = 0
    for n_value in validation_sizes:
        expected_bandwidth_rows += 4 + len(archive[f"posterior_nodes_n{n_value}"])
    if len(stored_bandwidths) != expected_bandwidth_rows:
        raise AssertionError(
            f"realized bandwidth table has {len(stored_bandwidths)} rows; "
            f"expected {expected_bandwidth_rows}"
        )
    if stored_bandwidths.duplicated(["validation_size", "arm", "node"]).any():
        raise AssertionError("realized bandwidth table contains duplicate arm/node keys")
    validation_ids = {}
    posterior_rows = 0
    max_likelihood_error = 0.0
    for n in validation_sizes:
        validation = pd.read_csv(root / "observable" / f"validation_n{n}.csv")
        validation_ids[n] = validation["validation_id"].to_numpy(int)
        if len(validation) != n:
            raise AssertionError(f"validation size mismatch at n={n}")
        validation_join = validation.merge(
            validation_truth[["frame_id", "redshift", "retained"]],
            left_on="validation_id", right_on="frame_id",
            suffixes=("_observable", "_truth"), validate="one_to_one",
        )
        _assert_close(
            validation_join["redshift_observable"].to_numpy(float),
            validation_join["redshift_truth"].to_numpy(float),
            f"validation redshift join n={n}",
        )
        _assert_close(
            validation_join["retained_observable"].to_numpy(int),
            validation_join["retained_truth"].to_numpy(int),
            f"validation retention join n={n}",
            atol=0.0,
        )
        y = validation["retained"].to_numpy(int)
        z_validation = validation["redshift"].to_numpy(float)
        x_validation = (z_validation - center) / scale
        design_validation = np.column_stack([np.ones(n), x_validation])
        nodes = archive[f"posterior_nodes_n{n}"]
        mode = archive[f"posterior_mode_n{n}"]
        quadrature = archive[f"posterior_weights_n{n}"]
        if nodes.shape[1] != 2 or mode.shape != (2,):
            raise AssertionError(f"unexpected inclusion model dimension at n={n}")
        if not np.isclose(np.sum(quadrature), 1.0) or np.any(quadrature < 0):
            raise AssertionError(f"invalid coefficient quadrature weights at n={n}")
        if not np.all(np.isfinite(design_validation @ mode)):
            raise AssertionError(f"non-finite validation MAP predictor at n={n}")

        rebuilt_mode, rebuilt_covariance, rebuilt_nodes, rebuilt_weights = (
            _independent_inclusion_quadrature(
                validation["redshift"].to_numpy(float),
                validation["retained"].to_numpy(int),
                order=int(config["quadrature_order"]),
                center=center,
                scale=scale,
            )
        )
        if not np.allclose(mode, rebuilt_mode, rtol=2e-8, atol=2e-9):
            raise AssertionError(f"independent inclusion posterior mode mismatch n={n}")
        if not np.allclose(nodes, rebuilt_nodes, rtol=2e-7, atol=2e-8):
            raise AssertionError(f"independent inclusion quadrature nodes mismatch n={n}")
        if not np.allclose(quadrature, rebuilt_weights, rtol=2e-7, atol=2e-9):
            raise AssertionError(f"independent inclusion quadrature weights mismatch n={n}")

        x_catalogue = (catalogue_z - center) / scale
        design_catalogue = np.column_stack([np.ones(len(catalogue_z)), x_catalogue])
        p_map = expit(design_catalogue @ mode)
        p_nodes = expit(design_catalogue @ nodes.T).T
        if np.any(p_map <= 0) or np.any(p_nodes <= 0):
            raise AssertionError(f"non-positive estimated inclusion probability at n={n}")

        arm_inputs = {
            "A_catalogue_only": (catalogue_z, catalogue_base),
            "B_parent_oracle": (parent_z, parent_base),
            "C_known_probability_ipw": (catalogue_z, catalogue_base / true_probability),
            "D_map_plugin_ipw": (catalogue_z, catalogue_base / p_map),
        }
        reconstructed = {}
        expected_bandwidths = {}
        for arm, (z_host, raw_weights) in arm_inputs.items():
            expected_bandwidths[(arm, -1)] = _independent_bandwidth(
                z_host, raw_weights, 0.01, 0.40, scott_multiplier
            )
            p_host = _independent_density(
                z_grid, z_host, raw_weights, rate_grid, scott_multiplier
            )
            _assert_close(
                archive[f"density_{arm}_n{n}"], p_host,
                f"{arm} independently reconstructed KDE n={n}",
            )
            reconstructed[arm], _ = _selected_kde_log_likelihood(
                grid, z_grid, p_host, distances, sigma, threshold
            )
        node_likelihoods = []
        node_efficiencies = []
        node_densities = []
        for node_index, node_probability in enumerate(p_nodes):
            node_raw_weights = catalogue_base / node_probability
            expected_bandwidths[("E_uncertainty_integrated_ipw", node_index)] = _independent_bandwidth(
                catalogue_z, node_raw_weights, 0.01, 0.40, scott_multiplier
            )
            p_host = _independent_density(
                z_grid, catalogue_z, node_raw_weights, rate_grid, scott_multiplier
            )
            node_densities.append(p_host)
            ll, efficiency = _selected_kde_log_likelihood(
                grid, z_grid, p_host, distances, sigma, threshold,
            )
            node_likelihoods.append(ll)
            node_efficiencies.append(efficiency)
        _assert_close(
            archive[f"density_E_nodes_n{n}"], np.asarray(node_densities),
            f"E node KDEs n={n}",
        )
        node_likelihoods = np.asarray(node_likelihoods)
        node_efficiencies = np.asarray(node_efficiencies)
        positive = quadrature > 0
        e_log_likelihood = logsumexp(
            node_likelihoods[positive] + np.log(quadrature[positive, None]),
            axis=0,
        )
        reconstructed["E_uncertainty_integrated_ipw"] = e_log_likelihood

        for (arm_name, node_index), expected_width in expected_bandwidths.items():
            row = stored_bandwidths.loc[
                (stored_bandwidths["validation_size"] == n)
                & (stored_bandwidths["arm"] == arm_name)
                & (stored_bandwidths["node"] == node_index)
            ]
            if len(row) != 1 or not np.isclose(
                float(row.iloc[0]["realized_bandwidth_redshift"]),
                expected_width, rtol=2e-10, atol=2e-12,
            ):
                raise AssertionError(
                    f"realized KDE bandwidth mismatch for {arm_name}, n={n}, node={node_index}"
                )

        for arm, expected in reconstructed.items():
            key = f"{arm}_n{n}"
            _assert_close(archive[key], expected, f"{arm} log likelihood n={n}")
            max_likelihood_error = max(
                max_likelihood_error,
                float(np.max(np.abs(archive[key] - expected))),
            )
            row = stored_summary.loc[
                (stored_summary["validation_size"] == n)
                & (stored_summary["arm"] == arm)
            ]
            if len(row) != 1 or row.iloc[0]["status"] != "SUCCESS":
                raise AssertionError(f"posterior summary row missing/invalid: {arm}, n={n}")
            expected_quantiles = {
                "posterior_median": _posterior_quantile(grid, expected, 0.50),
                "q025": _posterior_quantile(grid, expected, 0.025),
                "q05": _posterior_quantile(grid, expected, 0.05),
                "q16": _posterior_quantile(grid, expected, 0.16),
                "q84": _posterior_quantile(grid, expected, 0.84),
                "q95": _posterior_quantile(grid, expected, 0.95),
                "q975": _posterior_quantile(grid, expected, 0.975),
            }
            for field, value in expected_quantiles.items():
                if not np.isclose(float(row.iloc[0][field]), value, rtol=2e-9, atol=2e-9):
                    raise AssertionError(f"{arm} n={n}: posterior summary {field} mismatch")
            if not np.isclose(float(row.iloc[0]["median_error"]), expected_quantiles["posterior_median"] - h0_true):
                raise AssertionError(f"{arm} n={n}: median error mismatch")
            posterior_rows += 1

        _assert_close(
            archive[f"E_node_log_likelihood_n{n}"],
            node_likelihoods,
            f"E node likelihoods n={n}",
        )
        _assert_close(
            archive[f"E_node_log_efficiency_n{n}"],
            node_efficiencies,
            f"E node detection efficiencies n={n}",
        )

    for smaller, larger in zip(validation_sizes, validation_sizes[1:]):
        if not np.array_equal(
            validation_ids[smaller], validation_ids[larger][:len(validation_ids[smaller])]
        ):
            raise AssertionError("validation samples are not nested prefixes")

    return {
        "status": "PASS",
        "manifest_files_verified": checked_files,
        "posterior_summary_rows_reconstructed": posterior_rows,
        "maximum_log_likelihood_absolute_difference": max_likelihood_error,
        "validation_sizes_verified": validation_sizes,
        "scope": (
            "Independently reconstructs input joins, inclusion probabilities, "
            "A-D/E weights, selected likelihoods, posterior summaries and file hashes. "
            "Independently reconstructs the logistic posterior mode, Laplace "
            "covariance, Gauss-Hermite nodes, and target-corrected quadrature "
            "weights using a separate damped-Newton implementation."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = audit(args.artifact_dir)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
