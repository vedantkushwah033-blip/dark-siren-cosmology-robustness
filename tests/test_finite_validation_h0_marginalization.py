import numpy as np
import pandas as pd
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


def test_real_importance_posterior_flows_through_density_and_selected_likelihood():
    from scipy.special import expit
    from dark_siren_robustness.bayesian_inclusion import (
        importance_sample_inclusion_posterior,
    )

    rng = np.random.default_rng(910)
    validation_z = rng.uniform(0.01, 0.40, 160)
    validation_u = (validation_z - 0.01) / 0.39
    validation_probability = expit(-0.8 - 0.7 * validation_u)
    validation_labels = rng.binomial(1, validation_probability)
    posterior = importance_sample_inclusion_posterior(
        validation_z,
        validation_labels,
        degree=1,
        prior_sd=1.0,
        n_draws=1000,
        seed=12,
        min_ess=30,
        min_ess_fraction=0.01,
    )
    assert posterior.status == "SUCCESS", posterior.message

    analysis_parent_z = rng.uniform(0.01, 0.40, 400)
    analysis_u = (analysis_parent_z - 0.01) / 0.39
    analysis_probability = expit(-0.8 - 0.7 * analysis_u)
    analysis_retained = rng.binomial(1, analysis_probability)
    analysis_catalogue = analysis_parent_z[analysis_retained == 1]
    assert analysis_catalogue.size >= 20

    redshift_grid = np.linspace(0.01, 0.40, 101)
    h0_grid = np.linspace(55.0, 85.0, 41)
    merger = 1.0 + 1.5 * (redshift_grid - 0.01) / 0.39
    distances = np.array([330.0, 410.0, 520.0])
    result = propagate_inclusion_posterior_to_h0(
        posterior,
        analysis_catalogue,
        redshift_grid,
        h0_grid,
        distances,
        merger_weight=merger,
        scott_multiplier=0.3,
        sigma=120.0,
        detection_threshold=250.0,
    )
    assert result.host_density_draws.shape == (1000, 101)
    assert result.selected_log_likelihood_by_draw.shape == (1000, 41)
    assert np.all(np.isfinite(result.marginalized_selected_log_likelihood))
    assert result.importance_ess >= 30
    assert result.max_normalized_importance_weight < 1.0


def test_loop_reference_pipeline_matches_independent_vectorized_h0_likelihood():
    from dark_siren_robustness.finite_validation_h0 import (
        selected_log_likelihood_by_draw as vectorized_selected_likelihood,
    )

    proposals = np.array([[-0.5, -0.8], [0.2, -1.1]])
    weights = np.array([0.3, 0.7])
    posterior = InclusionPosteriorApproximation(
        status="SUCCESS",
        message="cross-implementation fixture",
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
    h0_grid = np.linspace(55.0, 85.0, 31)
    merger = 1.0 + 1.5 * (z_grid - 0.01) / 0.39
    distances = np.array([330.0, 410.0, 520.0])
    loop_result = propagate_inclusion_posterior_to_h0(
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
    vectorized = vectorized_selected_likelihood(
        z_grid,
        loop_result.host_density_draws,
        h0_grid,
        distances,
        sigma=120.0,
        detection_threshold=250.0,
        h0_batch_size=7,
    )
    np.testing.assert_allclose(
        loop_result.unconditional_log_likelihood_by_draw,
        vectorized["unconditional_log_likelihoods"],
        rtol=1e-10,
        atol=1e-10,
    )
    np.testing.assert_allclose(
        loop_result.detection_efficiency_by_draw,
        vectorized["detection_efficiencies"],
        rtol=1e-10,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        loop_result.selected_log_likelihood_by_draw,
        vectorized["selected_log_likelihoods"],
        rtol=1e-10,
        atol=1e-10,
    )


def test_batch_mcse_reports_low_effective_batches_and_insufficient_draws():
    from dark_siren_robustness.finite_validation_h0_marginalization import (
        marginalize_selected_draws_with_mcse,
    )

    log_likelihood = np.zeros((100, 11))
    concentrated_weights = np.full(100, -20.0)
    concentrated_weights[0] = 0.0
    low_ess = marginalize_selected_draws_with_mcse(
        log_likelihood, concentrated_weights, n_batches=10, min_batch_ess=10
    )
    assert low_ess["status"] == "LOW_BATCH_ESS"
    assert low_ess["batch_ess"].shape == (10,)
    assert low_ess["log_likelihood_mcse"].shape == (11,)

    too_few = marginalize_selected_draws_with_mcse(
        np.zeros((2, 11)), np.zeros(2), n_batches=10
    )
    assert too_few["status"] == "INSUFFICIENT_DRAWS"
    assert too_few["log_likelihood_mcse"] is None


def test_end_to_end_parent_frame_validation_to_h0_posterior_smoke():
    from scipy.integrate import cumulative_trapezoid, trapezoid
    from scipy.optimize import brentq
    from scipy.special import expit

    from dark_siren_robustness.bayesian_inclusion import importance_sample_inclusion_posterior
    from dark_siren_robustness.finite_validation_density_recovery import parent_intensity
    from dark_siren_robustness.finite_validation_h0 import normalize_h0_posterior
    from dark_siren_robustness.selection_toy import C_KM_S
    from dark_siren_robustness.validation_design import split_disjoint_validation

    seed = 20261108
    n_parent = 1500
    n_validation = 200
    gamma = 1.2
    target_retention = 0.25
    h0_true = 70.0
    sigma = 140.0
    threshold = 250.0

    parent_grid = np.linspace(0.01, 0.40, 20001)
    intensity = parent_intensity(parent_grid)
    parent_cdf = np.r_[0.0, cumulative_trapezoid(intensity, parent_grid)]
    parent_cdf /= parent_cdf[-1]
    def stream(stream_id):
        return np.random.default_rng(np.random.SeedSequence([seed, 0, stream_id]))

    parent_rng = stream(1)
    parent_z = np.interp(parent_rng.random(n_parent), parent_cdf, parent_grid)
    x_parent = (parent_z - 0.01) / 0.39
    curvature_center = trapezoid(intensity * ((parent_grid - 0.01) / 0.39) ** 2,
                                 parent_grid) / trapezoid(intensity, parent_grid)

    def mean_retention(alpha):
        probability = expit(alpha - (parent_grid - 0.01) / 0.39
                            + gamma * (((parent_grid - 0.01) / 0.39) ** 2
                                       - curvature_center))
        return trapezoid(intensity * probability, parent_grid) / trapezoid(intensity, parent_grid)

    alpha = brentq(lambda value: mean_retention(value) - target_retention, -30.0, 30.0)
    true_probability = expit(alpha - x_parent + gamma * (x_parent**2 - curvature_center))
    retention_rng = stream(3)
    retention_uniform = retention_rng.random(n_parent)
    retained = (retention_uniform < true_probability).astype(int)
    parent = pd.DataFrame({
        "object_id": np.arange(n_parent),
        "redshift": parent_z,
        "retained": retained,
        "true_inclusion_probability": true_probability,
    })
    validation_seed = int(stream(2).integers(0, 2**32 - 1))
    split = split_disjoint_validation(parent, n_validation=n_validation, seed=validation_seed)
    assert np.intersect1d(
        split.validation_object_ids, split.analysis_object_ids
    ).size == 0
    assert split.validation.retained.nunique() == 2
    assert len(split.analysis_catalogue) >= 20

    inclusion_posterior = importance_sample_inclusion_posterior(
        split.validation.redshift.to_numpy(),
        split.validation.retained.to_numpy(),
        degree=2,
        prior_sd=1.0,
        n_draws=1000,
        seed=int(stream(5).integers(0, 2**32 - 1)),
        min_ess=50,
        min_ess_fraction=0.01,
    )
    assert inclusion_posterior.status == "SUCCESS", inclusion_posterior.message

    host_grid = np.linspace(0.01, 0.40, 20001)
    merger_parent = 1.0 + 1.5 * (host_grid - 0.01) / 0.39
    host_intensity = parent_intensity(host_grid) * merger_parent
    host_cdf = np.r_[0.0, cumulative_trapezoid(host_intensity, host_grid)]
    host_cdf /= host_cdf[-1]
    event_rng = stream(4)
    event_distances = []
    while len(event_distances) < 5:
        event_z = float(np.interp(event_rng.random(), host_cdf, host_grid))
        observed_distance = float(
            C_KM_S * event_z / h0_true + event_rng.normal(0.0, sigma)
        )
        if observed_distance >= threshold:
            event_distances.append(observed_distance)

    z_grid = np.linspace(0.01, 0.40, 101)
    h0_grid = np.linspace(55.0, 85.0, 61)
    merger = 1.0 + 1.5 * (z_grid - 0.01) / 0.39
    result = propagate_inclusion_posterior_to_h0(
        inclusion_posterior,
        split.analysis_catalogue.redshift.to_numpy(),
        z_grid,
        h0_grid,
        np.asarray(event_distances),
        merger_weight=merger,
        scott_multiplier=0.15,
        sigma=sigma,
        detection_threshold=threshold,
    )
    h0_posterior = normalize_h0_posterior(
        h0_grid,
        result.marginalized_selected_log_likelihood,
        prior_density=np.ones_like(h0_grid),
    )
    assert result.host_density_draws.shape == (1000, 101)
    assert result.selected_log_likelihood_by_draw.shape == (1000, 61)
    assert result.marginalization_mcse_status == "PASS"
    assert result.marginalization_log_likelihood_mcse.shape == (61,)
    assert result.marginalization_batch_ess.shape == (10,)
    assert np.all(result.marginalization_batch_ess >= 20.0)
    assert np.all(np.isfinite(result.detection_efficiency_by_draw))
    assert np.all((result.detection_efficiency_by_draw > 0)
                  & (result.detection_efficiency_by_draw <= 1.0 + 1e-12))
    assert abs(h0_posterior["posterior_integral"] - 1.0) < 1e-10
    assert h0_posterior["ci95_width"] > h0_posterior["ci68_width"]
    assert 0.0 <= h0_posterior["boundary_mass_fraction"] <= 1.0
