from pathlib import Path
import sys
import zipfile

import numpy as np
from scipy.integrate import trapezoid
from scipy.special import log_ndtr
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_geometric_prior_predictive_h0_grid import (
    _evaluate,
    _find_root,
    _posterior,
)


def test_posterior_rank_and_quantiles_are_normalized_on_dh0_grid():
    grid = np.geomspace(30.0, 110.0, 2241)
    log_likelihood = -0.5 * ((grid - 70.0) / 8.0) ** 2
    density, quantiles, rank = _posterior(grid, log_likelihood, 70.0)
    assert np.isclose(trapezoid(density, grid), 1.0, atol=1e-10)
    assert quantiles.shape == (9,)
    assert np.all(np.diff(quantiles) > 0)
    assert 0.45 < rank < 0.55


def test_selected_likelihood_recomputation_is_finite_and_shaped():
    z = np.linspace(0.01, 0.04, 101)
    densities = np.vstack([
        np.ones_like(z) / (z[-1] - z[0]),
        (z - z[0]) / trapezoid(z - z[0], z),
    ])
    h0 = np.geomspace(30.0, 110.0, 51)
    result = _evaluate(
        h0,
        densities,
        z,
        np.array([350.0, 500.0]),
        sigma=140.0,
        threshold=250.0,
    )
    assert result.shape == (2, h0.size)
    assert np.all(np.isfinite(result))


def test_batched_selected_likelihood_matches_eventwise_reference():
    z = np.linspace(0.01, 0.04, 83)
    densities = np.vstack([
        np.ones_like(z) / (z[-1] - z[0]),
        (z - z[0]) / trapezoid(z - z[0], z),
        np.exp(-0.5 * ((z - 0.025) / 0.004) ** 2),
    ])
    densities[2] /= trapezoid(densities[2], z)
    h0 = np.geomspace(30.0, 110.0, 37)
    distances = np.array([320.0, 410.0, 560.0])
    sigma = 140.0
    threshold = 250.0

    actual = _evaluate(h0, densities, z, distances, sigma, threshold)

    delta = np.diff(z)
    weights = np.empty_like(z)
    weights[0] = delta[0] / 2.0
    weights[-1] = delta[-1] / 2.0
    weights[1:-1] = (delta[:-1] + delta[1:]) / 2.0
    means = 299792.458 * z[None, :] / h0[:, None]
    efficiency = (np.exp(log_ndtr((means - threshold) / sigma)) * weights) @ densities.T
    expected = np.zeros_like(actual)
    for distance in distances:
        numerator = (norm.pdf(distance, loc=means, scale=sigma) * weights) @ densities.T
        expected += np.log(numerator.T) - np.log(efficiency.T)

    np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-12)


def test_source_root_detection_requires_expected_replicate_path(tmp_path):
    archive_path = tmp_path / "source.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr(
            "experiment/replicates/rep000_gm1/config.json",
            "{}",
        )
    with zipfile.ZipFile(archive_path) as archive:
        assert _find_root(archive) == "experiment/"
