import numpy as np
import pytest
from scipy.special import logsumexp

from dark_siren_robustness.bayesian_inclusion import InclusionPosteriorApproximation
from dark_siren_robustness.finite_validation_h0_marginalization import (
    marginalize_selected_draws,
    propagate_inclusion_posterior_to_h0,
    selection_correct_each_draw,
)


def test_selection_correction_is_per_draw_before_weighted_marginalization():
    unconditional = np.array([[-2.0, -1.0], [-0.5, -3.0]])
    beta = np.array([[0.2, 0.4], [0.8, 0.9]])
    weights = np.array([0.25, 0.75])
    log_weights = np.log(weights)
    selected = selection_correct_each_draw(
        unconditional, beta, n_detected=3
    )
    expected = logsumexp(selected + log_weights[:, None], axis=0) - logsumexp(log_weights)
    wrong = (
        logsumexp(unconditional + log_weights[:, None], axis=0)
        - logsumexp(log_weights)
        - 3 * np.log(np.average(beta, axis=0, weights=weights))
    )
    result = marginalize_selected_draws(selected, log_weights)
    np.testing.assert_allclose(selected, unconditional - 3 * np.log(beta))
    np.testing.assert_allclose(result, expected)
    assert not np.allclose(result, wrong)


def test_equal_weights_recover_unweighted_logmeanexp():
    selected = np.array([[-1.0, -2.0], [-3.0, -0.5], [-0.2, -1.2]])
    result = marginalize_selected_draws(selected, np.zeros(3))
    expected = logsumexp(selected, axis=0) - np.log(3)
    np.testing.assert_allclose(result, expected)


def test_invalid_draw_or_efficiency_is_not_silently_removed():
    log_likelihood = np.zeros((2, 3))
    with pytest.raises(ValueError):
        selection_correct_each_draw(
            log_likelihood, np.array([[0.5, 0.5, 0.5], [0.5, 0.0, 0.5]]),
            n_detected=2,
        )
    with pytest.raises(ValueError):
        selection_correct_each_draw(
            log_likelihood, np.array([[0.5, 0.5, 0.5], [0.5, np.nan, 0.5]]),
            n_detected=2,
        )
    with pytest.raises(ValueError):
        marginalize_selected_draws(log_likelihood, np.array([0.0, np.nan]))


def test_pipeline_uses_each_proposal_and_its_importance_weight():
    proposals = np.array([[-0.5, -0.8], [0.2, -1.1]])
    weights = np.array([0.3, 0.7])
    posterior = InclusionPosteriorApproximation(
        status="SUCCESS",
        message="test fixture",
        degree=1,
        z_min=0.01,
        z_max=0.40,
        prior_mean=np.zeros(2),
        prior_sd=np.ones(2),
        map_coefficients=np.zeros(2),
        laplace_covariance=np.eye(2),
        coefficient_proposals=proposals,
        normalized_importance_weights=weights,
        log_importance_weights=np.log(weights),
        posterior_mean=np.zeros(2),
        posterior_covariance=np.eye(2),
        importance_ess=float(1 / np.sum(weights**2)),
        max_normalized_weight=float(np.max(weights)),
        log_evidence_estimate=0.0,
        log_evidence_relative_mcse=0.0,
        n_draws=2,
        newton_iterations=1,
    )
    catalogue = np.linspace(0.02, 0.39, 80)
    z_grid = np.linspace(0.01, 0.40, 101)
    h0_grid = np.linspace(55, 85, 31)
    merger = 1.0 + 1.5 * (z_grid - 0.01) / 0.39
    distances = np.array([330.0, 410.0, 520.0])
    result = propagate_inclusion_posterior_to_h0(
        posterior,
        catalogue,
        z_grid,
        h0_grid,
        distances,
        merger_weight=merger,
        scott_multiplier=0.3,
        sigma=120.0,
        detection_threshold=250.0,
    )
    assert result.host_density_draws.shape == (2, 101)
    assert result.selected_log_likelihood_by_draw.shape == (2, 31)
    assert not np.allclose(
        result.detection_efficiency_by_draw[0],
        result.detection_efficiency_by_draw[1],
    )
    expected = logsumexp(
        result.selected_log_likelihood_by_draw + np.log(weights)[:, None],
        axis=0,
    ) - logsumexp(np.log(weights))
    np.testing.assert_allclose(result.marginalized_selected_log_likelihood, expected)
    assert np.isclose(result.importance_ess, 1 / np.sum(weights**2))


def test_pipeline_rejects_truncated_redshift_domain_and_oversized_output():
    proposals = np.array([[-0.5, -0.8], [0.2, -1.1]])
    weights = np.array([0.3, 0.7])
    posterior = InclusionPosteriorApproximation(
        status="SUCCESS",
        message="test fixture",
        degree=1,
        z_min=0.01,
        z_max=0.40,
        prior_mean=np.zeros(2),
        prior_sd=np.ones(2),
        map_coefficients=np.zeros(2),
        laplace_covariance=np.eye(2),
        coefficient_proposals=proposals,
        normalized_importance_weights=weights,
        log_importance_weights=np.log(weights),
        posterior_mean=np.zeros(2),
        posterior_covariance=np.eye(2),
        importance_ess=float(1 / np.sum(weights**2)),
        max_normalized_weight=float(np.max(weights)),
        log_evidence_estimate=0.0,
        log_evidence_relative_mcse=0.0,
        n_draws=2,
        newton_iterations=1,
    )
    catalogue = np.linspace(0.02, 0.39, 80)
    z_grid = np.linspace(0.01, 0.40, 101)
    h0_grid = np.linspace(55, 85, 31)
    merger = 1.0 + 1.5 * (z_grid - 0.01) / 0.39
    arguments = dict(
        posterior=posterior,
        analysis_redshift=catalogue,
        redshift_grid=z_grid,
        H0_grid=h0_grid,
        distances=np.array([330.0, 410.0, 520.0]),
        merger_weight=merger,
        scott_multiplier=0.3,
        sigma=120.0,
        detection_threshold=250.0,
    )
    with pytest.raises(ValueError, match="full fixed inclusion-prior domain"):
        propagate_inclusion_posterior_to_h0(
            **{**arguments, "redshift_grid": z_grid[1:], "merger_weight": merger[1:]}
        )
    with pytest.raises(MemoryError, match="reduce proposal draws or grid size"):
        propagate_inclusion_posterior_to_h0(**arguments, max_output_bytes=100)
