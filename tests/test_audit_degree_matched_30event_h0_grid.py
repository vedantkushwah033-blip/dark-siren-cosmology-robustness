import numpy as np
import pytest
from scipy.special import logsumexp

from scripts.audit_degree_matched_30event_h0_grid import (
    _posterior,
    _selected_loglike,
)


def test_posterior_normalizes_and_returns_truth_rank():
    grid = np.linspace(30.0, 110.0, 281)
    loglike = -0.5 * ((grid - 70.0) / 5.0) ** 2
    density, quantiles, rank = _posterior(grid, loglike, 70.0)
    assert np.isclose(np.trapz(density, grid), 1.0)
    assert quantiles[2] == pytest.approx(70.0, abs=0.01)
    assert rank == pytest.approx(0.5, abs=0.01)


def test_selected_loglike_is_finite_and_shape_stable():
    z = np.linspace(0.01, 0.40, 201)
    density = np.exp(-0.5 * ((z - 0.2) / 0.06) ** 2)
    density /= np.trapz(density, z)
    grid = np.linspace(30.0, 110.0, 41)
    result = _selected_loglike(
        grid, density[None, :], z, np.array([300.0, 500.0]),
        sigma=140.0, threshold=250.0,
    )
    assert result.shape == (1, 41)
    assert np.isfinite(result).all()


def test_degree_two_e_nodes_can_be_marginalized():
    grid = np.linspace(30.0, 110.0, 41)
    z = np.linspace(0.01, 0.40, 201)
    base = np.exp(-0.5 * ((z - 0.2) / 0.06) ** 2)
    base /= np.trapz(base, z)
    densities = np.vstack([base, base * (1 + 0.1 * (z - 0.2))])
    densities /= np.trapz(densities, z, axis=1)[:, None]
    node_loglikes = _selected_loglike(
        grid, densities, z, np.array([400.0]), 140.0, 250.0
    )
    weights = np.array([0.25, 0.75])
    marginal = logsumexp(node_loglikes + np.log(weights[:, None]), axis=0)
    assert marginal.shape == (41,)
    assert np.isfinite(marginal).all()
