"""Independent raw-input and likelihood audit for quadratic-logit five-arm smoke artifacts.

This auditor does not import the five-arm runner, weight builder, likelihood
implementation, or H0 posterior normalizer. It independently reconstructs the
logistic coefficient posterior mode, Laplace covariance, Gauss-Hermite nodes,
quadrature weights, host weights, and selected likelihood from saved inputs.
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
    center: float, scale: float, degree: int = 1, quadratic_center: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Rebuild the proper logistic posterior and corrected GH rule independently.

    Uses a damped Newton optimizer rather than the production trust-exact
    optimizer, then reconstructs the Laplace proposal and target/proposal
    importance correction from the Bernoulli likelihood and Gaussian prior.
    """
    z = np.asarray(validation_z, dtype=float)
    y = np.asarray(retained, dtype=float)
    x = (z - center) / scale
    if degree not in (1, 2):
        raise ValueError("independent audit degree must be 1 or 2")
    dimension = degree + 1
    features = [np.ones_like(x), x]
    if degree == 2:
        features.append(x * x - quadratic_center)
    design = np.column_stack(features)
    prior_mean = np.zeros(dimension)
    prior_sd = np.asarray([5.0] + [2.5] * (dimension - 1), dtype=float)
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

    mode = np.zeros(dimension)
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
    indices = list(product(range(order), repeat=dimension))
    standard = np.asarray([[gh_nodes[i] for i in index] for index in indices], dtype=float)
    base_log_weights = np.asarray([
        sum(math.log(gh_weights[i]) for i in index) - 0.5 * dimension * math.log(2.0 * math.pi)
        for index in indices
    ])
    nodes = mode[None, :] + standard @ chol.T
    log_proposal = (
        -0.5 * dimension * math.log(2.0 * math.pi)
        -0.5 * np.sum(standard**2, axis=1)
        -float(np.sum(np.log(np.diag(chol))))
    )
    eta = design @ nodes.T
    log_likelihood = np.sum(y[:, None] * eta - np.logaddexp(0.0, eta), axis=0)
    log_prior = (
        -0.5 * np.sum(((nodes - prior_mean) / prior_sd) ** 2, axis=1)
        -float(np.sum(np.log(prior_sd)))
        -0.5 * dimension * math.log(2.0 * math.pi)
    )
    log_importance = base_log_weights + log_likelihood + log_prior - log_proposal
    weights = np.exp(log_importance - logsumexp(log_importance))
    if not np.all(np.isfinite(weights)) or not np.isclose(weights.sum(), 1.0, atol=1e-12):
        raise AssertionError("independent logistic quadrature weights are invalid")
    return mode, covariance, nodes, weights


def audit_cell(input_dir: Path, output_dir: Path) -> dict:
    root, cell = input_dir.resolve(), output_dir.resolve()
    config_path, source_manifest_path = root / "config.json", root / "manifest.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    provenance = json.loads((cell / "provenance.json").read_text(encoding="utf-8"))
    summary = json.loads((cell / "summary.json").read_text(encoding="utf-8"))
    output_manifest = json.loads((cell / "manifest.json").read_text(encoding="utf-8"))
    if source_manifest.get("status") != "PASS" or output_manifest.get("status") != "PASS":
        raise AssertionError("source and cell manifests must both be PASS")
    if provenance["source_config_sha256"] != _sha256(config_path):
        raise AssertionError("source config hash mismatch")
    if provenance["source_manifest_sha256"] != _sha256(source_manifest_path):
        raise AssertionError("source manifest hash mismatch")
    for name, entry in output_manifest["files"].items():
        path = cell / name
        if not path.is_file() or path.stat().st_size != int(entry["bytes"]) or _sha256(path) != entry["sha256"]:
            raise AssertionError(f"cell output manifest mismatch: {name}")
    arm, nval = provenance["arm"], int(provenance["validation_size"])
    if arm not in ARM_KEYS or arm != summary["arm"] or nval != int(summary["validation_size"]):
        raise AssertionError("cell identity differs between provenance and summary")
    if int(provenance["replicate"]) != int(summary["replicate"]):
        raise AssertionError("replicate identity mismatch")
    verified_inputs = provenance["verified_inputs"]
    for relative, expected_hash in verified_inputs.items():
        path = root / relative
        entry = source_manifest.get("files", {}).get(relative)
        if entry is None or not path.is_file():
            raise AssertionError(f"source manifest is missing required input {relative}")
        if path.stat().st_size != int(entry["bytes"]) or _sha256(path) != expected_hash or expected_hash != entry["sha256"]:
            raise AssertionError(f"source input provenance mismatch: {relative}")

    events_frame = pd.read_csv(root / "observable" / "gw_events_observed.csv")
    distances = events_frame["distance_mpc"].to_numpy(float)
    sigma = float(config["sigma_mpc"])
    threshold = float(config["detection_threshold_mpc"])
    events_sigma = np.full(distances.shape, sigma)
    prior_min, prior_max = map(float, config["h0_prior_support"])
    h0_grid = np.linspace(prior_min, prior_max, int(config["h0_grid_points"]))
    z_grid = np.linspace(0.01, 0.40, int(config["redshift_grid_points"]))
    rate = np.ones_like(z_grid)
    scott = float(config["scott_multiplier"])
    parent = catalogue = None
    parent_base = catalogue_base = None
    if "truth/parent_truth.csv" in verified_inputs:
        parent = pd.read_csv(root / "truth/parent_truth.csv")
        parent_z = parent["redshift"].to_numpy(float)
        parent_base = 1.0 + 1.5 * (parent_z - 0.01) / 0.39
    if "observable/catalogue_observed.csv" in verified_inputs:
        catalogue_frame = pd.read_csv(root / "observable/catalogue_observed.csv")
        catalogue_z = catalogue_frame["redshift"].to_numpy(float)
        catalogue_base = 1.0 + 1.5 * (catalogue_z - 0.01) / 0.39

    mode = nodes = quadrature = None
    if arm == "B_parent_oracle":
        sample, raw_components = parent_z, parent_base[None, :]
    else:
        sample = catalogue_z
        if arm == "A_catalogue_only":
            raw_components = catalogue_base[None, :]
        elif arm == "C_known_probability_ipw":
            probability = catalogue_frame["object_id"].map(
                parent.set_index("object_id")["inclusion_probability"]
            ).to_numpy(float)
            raw_components = (catalogue_base / probability)[None, :]
        elif arm in {"D_map_plugin_ipw", "E_uncertainty_integrated_ipw"}:
            validation = pd.read_csv(root / f"observable/validation_n{nval}.csv")
            degree = int(config["inclusion_model_degree"])
            center, scale = 0.205, 0.195
            quadratic_center = (
                float(config["quadratic_calibration"]["x2_center"]) if degree == 2 else 0.0
            )
            mode, _cov, nodes, quadrature = _independent_inclusion_quadrature(
                validation["redshift"].to_numpy(float),
                validation["retained"].to_numpy(int),
                order=int(config["quadrature_order"]), center=center, scale=scale,
                degree=degree, quadratic_center=quadratic_center,
            )
            x = (catalogue_z - center) / scale
            features = [np.ones_like(x), x]
            if degree == 2:
                features.append(x * x - quadratic_center)
            design = np.column_stack(features)
            if arm == "D_map_plugin_ipw":
                raw_components = (catalogue_base / expit(design @ mode))[None, :]
            else:
                probabilities = expit(design @ nodes.T).T
                raw_components = catalogue_base[None, :] / probabilities

    if not np.all(np.isfinite(raw_components)) or np.any(raw_components <= 0):
        raise AssertionError("independently reconstructed raw weights are invalid")
    if arm == "E_uncertainty_integrated_ipw":
        mixture = quadrature
    else:
        mixture = np.ones(1)
    archive = np.load(cell / "arm_cell_arrays.npz")
    _assert_close(archive["H0_grid"], h0_grid, "H0 grid")
    _assert_close(archive["redshift_grid"], z_grid, "redshift grid")
    _assert_close(archive["component_weights"], raw_components, "component raw weights")
    _assert_close(archive["raw_weights"], raw_components[0], "primary raw weights")
    _assert_close(archive["mixture_weights"], mixture, "mixture weights")
    if mode is None:
        if archive["inclusion_mode"].size or archive["inclusion_nodes"].size or archive["inclusion_quadrature_weights"].size:
            raise AssertionError("A/B/C must not persist inclusion-posterior fit outputs")
    else:
        _assert_close(archive["inclusion_mode"], mode, "independent inclusion mode")
        _assert_close(archive["inclusion_nodes"], nodes, "independent inclusion nodes")
        _assert_close(archive["inclusion_quadrature_weights"], quadrature, "independent inclusion weights")

    densities, likelihoods, efficiencies, bandwidths = [], [], [], []
    for weights in raw_components:
        density = _independent_density(z_grid, sample, weights, rate, scott)
        ll, log_efficiency = _selected_kde_log_likelihood(
            h0_grid, z_grid, density, distances, sigma, threshold
        )
        densities.append(density)
        likelihoods.append(ll)
        efficiencies.append(log_efficiency)
        bandwidths.append(_independent_bandwidth(sample, weights, 0.01, 0.40, scott))
    densities = np.asarray(densities)
    likelihoods = np.asarray(likelihoods)
    efficiencies = np.asarray(efficiencies)
    bandwidths = np.asarray(bandwidths)
    _assert_close(archive["component_density"], densities, "independent component KDE")
    _assert_close(archive["component_log_likelihood"], likelihoods, "independent component likelihood")
    _assert_close(archive["component_log_detection_efficiency"], efficiencies, "independent detection efficiency")
    _assert_close(archive["realized_bandwidths"], bandwidths, "realized KDE bandwidth")
    if len(mixture) > 1:
        positive = mixture > 0
        marginal = logsumexp(likelihoods[positive] + np.log(mixture[positive, None]), axis=0)
    else:
        marginal = likelihoods[0]
    _assert_close(archive["marginal_log_likelihood"], marginal, "marginal selected likelihood")
    quantiles = {
        "posterior_median": _posterior_quantile(h0_grid, marginal, 0.50),
        "q025": _posterior_quantile(h0_grid, marginal, 0.025),
        "q05": _posterior_quantile(h0_grid, marginal, 0.05),
        "q16": _posterior_quantile(h0_grid, marginal, 0.16),
        "q84": _posterior_quantile(h0_grid, marginal, 0.84),
        "q95": _posterior_quantile(h0_grid, marginal, 0.95),
        "q975": _posterior_quantile(h0_grid, marginal, 0.975),
    }
    for key, value in quantiles.items():
        if not np.isclose(float(summary[key]), value, rtol=2e-8, atol=2e-8):
            raise AssertionError(f"posterior summary mismatch: {key}")
    if not np.isclose(float(summary["median_error"]), quantiles["posterior_median"] - float(config["h0_true"]), rtol=2e-8, atol=2e-8):
        raise AssertionError("posterior median error mismatch")
    return {
        "status": "PASS",
        "arm": arm,
        "validation_size": nval,
        "replicate": int(provenance["replicate"]),
        "component_count": int(raw_components.shape[0]),
        "maximum_log_likelihood_absolute_difference": float(np.max(np.abs(archive["marginal_log_likelihood"] - marginal))),
        "output_files_verified": int(output_manifest["file_count"]),
        "source_inputs_verified": len(verified_inputs),
        "claim_limit": "Independent single-cell reconstruction; not repeated-sampling calibration.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--cell-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = audit_cell(args.input_dir, args.cell_dir)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
