import numpy as np
import pytest
from scipy.special import logsumexp

from dark_siren_robustness.inclusion_posterior import (
    fit_inclusion_posterior_quadrature,
    log_marginalize_shared_parameter,
    posterior_predictive_inclusion,
)


def test_proper_posterior_quadrature_normalizes_and_predicts():
    rng = np.random.default_rng(4182)
    z = rng.uniform(0.01, 0.40, 500)
    x = (z - z.mean()) / z.std()
    probability = 1.0 / (1.0 + np.exp(-(-0.35 + 0.85 * x + 0.25 * x * x)))
    retained = rng.binomial(1, probability)

    posterior = fit_inclusion_posterior_quadrature(
        z, retained, degree=2, order=13
    )
    assert posterior.nodes.shape == (13**3, 3)
    assert posterior.weights.shape == (13**3,)
    assert np.all(np.isfinite(posterior.nodes))
    assert np.all(np.isfinite(posterior.weights))
    assert np.all(posterior.weights >= 0)
    assert np.isclose(posterior.weights.sum(), 1.0)
    assert 1.0 <= posterior.quadrature_ess <= posterior.weights.size
    assert 0 < posterior.maximum_normalized_weight <= 1
    assert np.all(np.linalg.eigvalsh(posterior.posterior_covariance) > 0)

    predictive = posterior_predictive_inclusion(posterior, np.linspace(0.01, 0.40, 51))
    assert predictive.shape == (51,)
    assert np.all(np.isfinite(predictive))
    assert np.all((predictive > 0) & (predictive < 1))


def test_quadrature_order_refinement_is_stable_for_regular_validation_sample():
    rng = np.random.default_rng(9102)
    z = rng.uniform(0.01, 0.40, 800)
    x = (z - z.mean()) / z.std()
    p = 1.0 / (1.0 + np.exp(-(-0.2 + 0.65*x - 0.18*x*x)))
    retained = rng.binomial(1, p)
    low = fit_inclusion_posterior_quadrature(z, retained, degree=2, order=9)
    high = fit_inclusion_posterior_quadrature(z, retained, degree=2, order=15)
    assert np.max(np.abs(low.posterior_mean - high.posterior_mean)) < 0.04
    grid = np.linspace(0.01, 0.40, 101)
    assert np.max(np.abs(
        posterior_predictive_inclusion(low, grid)
        - posterior_predictive_inclusion(high, grid)
    )) < 0.015


def test_proper_prior_keeps_one_class_validation_posterior_finite():
    z = np.linspace(0.01, 0.40, 80)
    retained = np.ones(z.size, dtype=int)
    posterior = fit_inclusion_posterior_quadrature(
        z, retained, degree=1, order=11
    )
    assert np.all(np.isfinite(posterior.nodes))
    assert np.isclose(posterior.weights.sum(), 1.0)
    assert np.all(np.isfinite(posterior.posterior_mean))


def test_rejects_rank_deficient_or_invalid_validation_data():
    with pytest.raises(ValueError, match="positive finite spread"):
        fit_inclusion_posterior_quadrature(np.ones(20), np.zeros(20), degree=1)
    z = np.linspace(0.01, 0.40, 20)
    bad = np.zeros(20)
    bad[0] = 2
    with pytest.raises(ValueError, match="binary"):
        fit_inclusion_posterior_quadrature(z, bad, degree=1)
    with pytest.raises(ValueError, match="prior_sd"):
        fit_inclusion_posterior_quadrature(
            z, np.arange(20) % 2, degree=1, prior_sd=[5.0, 0.0]
        )


def test_shared_parameter_marginalization_is_not_eventwise_factorization():
    # Two events prefer opposite values of the same shared inclusion parameter.
    # E_theta[L1(theta)L2(theta)] = 0.5*(1*9 + 9*1) = 9,
    # while E[L1]E[L2] = 5*5 = 25.
    event_log_likelihoods = np.log(np.asarray([[1.0, 9.0], [9.0, 1.0]]))
    joint_by_node = event_log_likelihoods.sum(axis=1)
    marginalized_joint = log_marginalize_shared_parameter(
        joint_by_node, np.asarray([0.5, 0.5])
    )
    eventwise = log_marginalize_shared_parameter(
        event_log_likelihoods, np.asarray([0.5, 0.5])
    )
    assert np.isclose(np.exp(marginalized_joint), 9.0)
    assert np.isclose(np.prod(np.exp(eventwise)), 25.0)
    assert not np.isclose(marginalized_joint, np.sum(eventwise))


def test_log_marginalization_is_stable_to_common_likelihood_offset():
    log_likelihood = np.asarray([[-10000.0, -10002.0], [-10001.0, -9999.0]])
    weights = np.asarray([0.3, 0.7])
    result = log_marginalize_shared_parameter(log_likelihood, weights)
    shifted = log_marginalize_shared_parameter(log_likelihood + 12000.0, weights)
    assert np.allclose(shifted - 12000.0, result)


def test_log_marginalization_validates_weights_and_log_likelihoods():
    with pytest.raises(ValueError, match="sum to one"):
        log_marginalize_shared_parameter(np.zeros(2), np.asarray([0.2, 0.2]))
    with pytest.raises(ValueError, match="NaN"):
        log_marginalize_shared_parameter(np.asarray([0.0, np.nan]), np.asarray([0.5, 0.5]))
    with pytest.raises(ValueError, match="match posterior weights"):
        log_marginalize_shared_parameter(np.zeros((3, 2)), np.asarray([0.5, 0.5]))


def test_quadrature_ladder_reports_convergence_and_nonconvergence_explicitly():
    from dark_siren_robustness.inclusion_posterior import (
        fit_inclusion_posterior_quadrature_ladder,
    )

    rng = np.random.default_rng(22061)
    z = rng.uniform(0.01, 0.40, 350)
    x = (z - z.mean()) / z.std()
    p = 1.0 / (1.0 + np.exp(-(-0.25 + 0.7 * x)))
    retained = rng.binomial(1, p)

    permissive = fit_inclusion_posterior_quadrature_ladder(
        z, retained, degree=1, orders=(5, 7, 9),
        mean_tolerance=1.0, covariance_tolerance=1.0,
        predictive_tolerance=1.0, consecutive_passes=1,
    )
    assert permissive.status == "CONVERGED"
    assert len(permissive.fits) >= 2
    assert len(permissive.mean_max_abs_change) == len(permissive.fits) - 1
    assert len(permissive.log_evidence_abs_change) == len(permissive.fits) - 1
    assert np.all(np.isfinite(permissive.log_evidence_abs_change))

    strict = fit_inclusion_posterior_quadrature_ladder(
        z, retained, degree=1, orders=(5, 7, 9),
        mean_tolerance=1e-14, covariance_tolerance=1e-14,
        predictive_tolerance=1e-14, consecutive_passes=2,
    )
    assert strict.status == "NON_CONVERGED"
    assert len(strict.fits) == 3



def test_quadrature_ladder_requires_log_evidence_stability():
    from dark_siren_robustness.inclusion_posterior import (
        fit_inclusion_posterior_quadrature_ladder,
    )

    rng = np.random.default_rng(22062)
    z = rng.uniform(0.01, 0.40, 300)
    x = (z - z.mean()) / z.std()
    retained = rng.binomial(1, 1.0 / (1.0 + np.exp(-(-0.2 + 0.55*x))))

    result = fit_inclusion_posterior_quadrature_ladder(
        z, retained, degree=1, orders=(5, 7, 9),
        mean_tolerance=10.0, covariance_tolerance=10.0,
        predictive_tolerance=10.0, log_evidence_tolerance=1e-15,
        consecutive_passes=1,
    )
    assert result.status == "NON_CONVERGED"
    assert len(result.log_evidence_abs_change) == 2
    assert any(change > 1e-15 for change in result.log_evidence_abs_change)


def test_quadrature_ladder_rejects_invalid_order_schedule():
    from dark_siren_robustness.inclusion_posterior import (
        fit_inclusion_posterior_quadrature_ladder,
    )

    with pytest.raises(ValueError, match="strictly increasing"):
        fit_inclusion_posterior_quadrature_ladder(
            np.linspace(0.01, 0.4, 30), np.arange(30) % 2,
            orders=(7, 7, 11),
        )


def test_log_target_gradient_and_hessian_match_finite_differences():
    from dark_siren_robustness.inclusion_posterior import (
        _log_target, _prior_vectors,
    )

    rng = np.random.default_rng(7712)
    z = rng.uniform(0.01, 0.40, 70)
    x = (z - z.mean()) / z.std()
    design = np.column_stack([np.ones(z.size), x, x*x])
    y = rng.binomial(1, 1.0 / (1.0 + np.exp(-(-0.3 + 0.5*x - 0.2*x*x))))
    mu, sd = _prior_vectors(3, None, None)
    precision = 1.0 / sd**2
    theta = np.asarray([0.2, -0.4, 0.15])
    eps = 2e-5

    def objective(b):
        return -_log_target(b, design, y, mu, sd)

    def gradient(b):
        from scipy.special import expit
        return design.T @ (expit(design @ b) - y) + (b-mu)*precision

    p = 3
    numerical_gradient = np.empty(p)
    numerical_hessian = np.empty((p, p))
    for j in range(p):
        step = np.zeros(p)
        step[j] = eps
        numerical_gradient[j] = (objective(theta+step)-objective(theta-step))/(2*eps)
        numerical_hessian[:, j] = (gradient(theta+step)-gradient(theta-step))/(2*eps)

    from scipy.special import expit
    probability = expit(design @ theta)
    analytic_hessian = design.T @ (
        (probability*(1-probability))[:, None] * design
    ) + np.diag(precision)
    assert np.allclose(numerical_gradient, gradient(theta), rtol=2e-6, atol=2e-6)
    assert np.allclose(numerical_hessian, analytic_hessian, rtol=2e-6, atol=2e-6)


def test_shared_node_marginalization_constant_likelihood_is_exact():
    log_likelihood = np.full((7, 4), -3.75)
    weights = np.asarray([0.02, 0.08, 0.1, 0.2, 0.25, 0.15, 0.2])
    result = log_marginalize_shared_parameter(log_likelihood, weights)
    assert np.allclose(result, -3.75, rtol=0, atol=1e-14)


def test_node_specific_host_weights_are_normalized_and_report_diagnostics():
    from dark_siren_robustness.inclusion_posterior import (
        build_node_specific_host_weights,
    )

    rng = np.random.default_rng(3107)
    z_validation = rng.uniform(0.01, 0.40, 180)
    x = (z_validation-z_validation.mean())/z_validation.std()
    p = 1.0/(1.0+np.exp(-(-0.3+0.65*x)))
    retained = rng.binomial(1, p)
    posterior = fit_inclusion_posterior_quadrature(
        z_validation, retained, degree=1, order=9,
    )
    catalogue = np.linspace(0.012, 0.398, 250)
    host = build_node_specific_host_weights(posterior, catalogue)

    assert host.inclusion_probability.shape == (9**2, 250)
    assert host.raw_inverse_weight.shape == (9**2, 250)
    assert host.normalized_weight.shape == (9**2, 250)
    assert np.all((host.inclusion_probability > 0) & (host.inclusion_probability < 1))
    assert np.all(host.raw_inverse_weight >= 1)
    assert np.allclose(host.normalized_weight.sum(axis=1), 1.0)
    assert np.all((host.kish_ess >= 1) & (host.kish_ess <= 250))
    assert np.all((host.entropy_ess >= 1) & (host.entropy_ess <= 250))
    assert np.all(host.maximum_normalized_weight <= 1)
    assert np.all(host.top_one_percent_mass <= host.top_five_percent_mass + 1e-14)
    assert np.all(host.minimum_inclusion_probability > 0)


def test_node_specific_host_weights_refuse_underflow_without_clipping():
    from dataclasses import replace
    from dark_siren_robustness.inclusion_posterior import (
        build_node_specific_host_weights,
    )

    z = np.linspace(0.01, 0.40, 40)
    posterior = fit_inclusion_posterior_quadrature(
        z, np.zeros(z.size, dtype=int), degree=1, order=5,
    )
    extreme_nodes = posterior.nodes.copy()
    extreme_nodes[:, 0] = -1000.0
    extreme = replace(posterior, nodes=extreme_nodes)
    with pytest.raises(ArithmeticError, match="underflows"):
        build_node_specific_host_weights(extreme, np.asarray([0.2]))



def test_laplace_gh_quadrature_matches_independent_prior_centered_gh_reference():
    """Cross-check moments and evidence using a different quadrature measure."""
    from scipy.special import expit, roots_hermitenorm

    rng = np.random.default_rng(20261108)
    z = rng.uniform(0.01, 0.40, 20)
    center, scale = z.mean(), z.std()
    x = (z - center) / scale
    design = np.column_stack([np.ones(z.size), x])
    probability = expit(-0.45 + 0.72 * x)
    retained = rng.binomial(1, probability)
    prior_mean = np.zeros(2)
    prior_sd = np.asarray([1.5, 0.8])

    adaptive = fit_inclusion_posterior_quadrature(
        z, retained, degree=1, prior_mean=prior_mean,
        prior_sd=prior_sd, order=15,
    )

    # Independent tensor GH directly under the Gaussian prior, rather than
    # the implementation's Laplace-centered Gaussian importance proposal.
    nodes_1d, weights_1d = roots_hermitenorm(100)
    standard_weights = weights_1d / np.sqrt(2.0 * np.pi)
    theta0, theta1 = np.meshgrid(
        prior_mean[0] + prior_sd[0] * nodes_1d,
        prior_mean[1] + prior_sd[1] * nodes_1d,
        indexing="ij",
    )
    theta = np.column_stack([theta0.ravel(), theta1.ravel()])
    log_prior_quadrature = np.log(
        np.outer(standard_weights, standard_weights).ravel()
    )
    eta = design @ theta.T
    log_likelihood = np.sum(
        retained[:, None] * eta - np.logaddexp(0.0, eta), axis=0
    )
    log_mass = log_prior_quadrature + log_likelihood
    log_evidence_reference = float(logsumexp(log_mass))
    reference_weights = np.exp(log_mass - log_evidence_reference)
    reference_mean = reference_weights @ theta
    centered = theta - reference_mean
    reference_covariance = (
        centered * reference_weights[:, None]
    ).T @ centered

    np.testing.assert_allclose(
        adaptive.posterior_mean, reference_mean, rtol=2e-5, atol=2e-5
    )
    np.testing.assert_allclose(
        adaptive.posterior_covariance, reference_covariance,
        rtol=5e-5, atol=5e-5,
    )
    assert abs(
        adaptive.log_evidence_unnormalized - log_evidence_reference
    ) < 5e-5


def test_fixed_domain_coordinates_are_invariant_to_validation_sample_moments():
    from dark_siren_robustness.inclusion_posterior import (
        fit_inclusion_posterior_quadrature,
    )

    rng = np.random.default_rng(20261108)
    z_a = rng.uniform(0.01, 0.40, 80)
    z_b = rng.uniform(0.07, 0.32, 40)
    y_a = rng.binomial(1, 1.0 / (1.0 + np.exp(-(-0.4 + 0.7 * (z_a - 0.205) / 0.195))))
    y_b = rng.binomial(1, 1.0 / (1.0 + np.exp(-(-0.4 + 0.7 * (z_b - 0.205) / 0.195))))
    assert not np.isclose(np.mean(z_a), np.mean(z_b))
    assert not np.isclose(np.std(z_a), np.std(z_b))

    fit_a = fit_inclusion_posterior_quadrature(
        z_a, y_a, degree=1, order=5,
        coordinate_center=0.205, coordinate_scale=0.195,
    )
    fit_b = fit_inclusion_posterior_quadrature(
        z_b, y_b, degree=1, order=5,
        coordinate_center=0.205, coordinate_scale=0.195,
    )
    assert fit_a.center == fit_b.center == 0.205
    assert fit_a.scale == fit_b.scale == 0.195

    prediction_z = np.linspace(0.01, 0.40, 51)
    fixed_draw = np.asarray([-0.4, 0.7])
    curve_a = 1.0 / (1.0 + np.exp(-(((prediction_z - fit_a.center) / fit_a.scale) * fixed_draw[1] + fixed_draw[0])))
    curve_b = 1.0 / (1.0 + np.exp(-(((prediction_z - fit_b.center) / fit_b.scale) * fixed_draw[1] + fixed_draw[0])))
    np.testing.assert_array_equal(curve_a, curve_b)


def test_fixed_coordinate_center_and_scale_must_be_supplied_together():
    from dark_siren_robustness.inclusion_posterior import (
        fit_inclusion_posterior_quadrature,
    )

    z = np.linspace(0.01, 0.40, 20)
    y = np.tile([0, 1], 10)
    with pytest.raises(ValueError, match="must be supplied together"):
        fit_inclusion_posterior_quadrature(
            z, y, degree=1, order=3, coordinate_center=0.205
        )

def test_quadrature_ladder_propagates_fixed_coefficient_coordinates():
    from dark_siren_robustness.inclusion_posterior import (
        fit_inclusion_posterior_quadrature_ladder,
    )

    z = np.linspace(0.01, 0.40, 90)
    retained = (np.arange(z.size) % 3 != 0).astype(int)
    ladder = fit_inclusion_posterior_quadrature_ladder(
        z, retained, degree=1, orders=(3, 5),
        coordinate_center=0.205, coordinate_scale=0.195,
        predictive_redshift=np.linspace(0.01, 0.40, 11),
        mean_tolerance=10.0, covariance_tolerance=10.0,
        predictive_tolerance=10.0, log_evidence_tolerance=10.0,
        consecutive_passes=1,
    )
    assert len(ladder.fits) == 2
    assert all(fit.center == 0.205 for fit in ladder.fits)
    assert all(fit.scale == 0.195 for fit in ladder.fits)


def test_default_degree_one_prior_matches_degree_two_centered_marginal():
    from dark_siren_robustness.inclusion_posterior import _prior_vectors

    mean_d1, sd_d1 = _prior_vectors(2, None, None)
    mean_d2, sd_d2 = _prior_vectors(3, None, None)

    np.testing.assert_array_equal(mean_d1, mean_d2[:2])
    np.testing.assert_array_equal(sd_d1, sd_d2[:2])
    np.testing.assert_array_equal(mean_d1, np.zeros(2))
    np.testing.assert_array_equal(sd_d1, np.asarray([5.0, 2.5]))


def test_centered_quadratic_basis_matches_generator_parameterization():
    from dark_siren_robustness.inclusion_posterior import _validation_inputs

    z = np.linspace(0.01, 0.40, 101)
    coordinate_center = 0.205
    coordinate_scale = 0.195
    c_g = 0.2845005579191241
    alpha, beta, gamma = -0.4, -3.9, 1.2
    _, _, design, _, _, returned_center = _validation_inputs(
        z, np.arange(z.size) % 2, 2,
        coordinate_center=coordinate_center,
        coordinate_scale=coordinate_scale,
        quadratic_center=c_g,
    )
    x = (z - coordinate_center) / coordinate_scale
    expected_eta = alpha + beta * x + gamma * (x * x - c_g)
    assert returned_center == pytest.approx(c_g)
    assert np.allclose(design @ np.asarray([alpha, beta, gamma]), expected_eta, rtol=0, atol=1e-14)


def test_centered_quadratic_center_is_validated():
    z = np.linspace(0.01, 0.40, 30)
    with pytest.raises(ValueError, match="quadratic_center"):
        fit_inclusion_posterior_quadrature(
            z, np.arange(z.size) % 2, degree=2, quadratic_center=np.nan, order=5
        )


def test_degree_two_map_refines_past_trust_exact_stagnation_case():
    """Regression for a real validation frame where trust-exact stopped at ~8e-9."""
    from scipy.special import expit

    z = np.asarray([0.204243956167295,0.078799076227147,0.323433162125299,0.384786235178502,0.031368086617152,0.257224364621394,0.311200021781614,0.2981594265316,0.096100601036989,0.102030666362966,0.30414669696271,0.130736266448311,0.278959307081609,0.142863695323556,0.327275169680881,0.36016293866641,0.12222135300735,0.279491250552181,0.110671868005239,0.270500107137441,0.09951719426675,0.096247235965451,0.25615143338092,0.225019574579263,0.065469904783305,0.318992928083608,0.13530873070727,0.074099634697719,0.120191142457842,0.055702443460497,0.091614494705879,0.276016655074223,0.277099566554171,0.138489198441801,0.317608697968678,0.086935532338777,0.114800770480209,0.083338847128375,0.287520565381811,0.384648636096427,0.237514522270354,0.14631362568244,0.204094680525688,0.13237993351005,0.121013177335119,0.324837651777765,0.181991863033787,0.254280624895468,0.293393087938743,0.237532466707367,0.044841262427127,0.251879978281662,0.258922163954626,0.276101155886185,0.330011287505283,0.117408180380185,0.250874056252566,0.206094009740936,0.272642199462038,0.113162015594812,0.256303914384715,0.162447230356023,0.321417983223736,0.272342713531255,0.089475329304586,0.09018760390125,0.387493090736846,0.238438253100436,0.160726514200788,0.27436346258253,0.340634794450516,0.028589659333903,0.249199342415785,0.025223557720507,0.122949337204189,0.320978239227486,0.303346670272685,0.168699522778013,0.398232043479284,0.017911481160736,0.066156392318619,0.366703139479695,0.249149088257972,0.091905569307109,0.223127122642021,0.109305205147906,0.103891093901193,0.255800793406196,0.131255925286943,0.263159652686803,0.238154270752714,0.106402873936288,0.147683789473106,0.050546241083711,0.0684551980519,0.19754807339843,0.386712683014908,0.266176954970337,0.278448582214318,0.085609630504441], dtype=float)
    retained = np.asarray([0,1,0,0,1,0,1,0,1,1,1,1,0,1,0,0,0,1,1,0,1,1,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,1,1,1,0,1,0,1,1,0,0,1,1,1,0,0,1,0,0,1,1,0,0,0,1,1,1,0,1,1,1,0,1,0,1,1,0,0,1,0,1,1,0,1,1,0,0,0,0,1,0,0,1,0,1,1,1,0,0,0,1], dtype=int)
    center, scale, quadratic_center = 0.205, 0.195, 0.2845005579191241
    fit = fit_inclusion_posterior_quadrature(
        z, retained, degree=2, order=5,
        coordinate_center=center, coordinate_scale=scale,
        quadratic_center=quadratic_center,
    )
    x = (z - center) / scale
    design = np.column_stack([np.ones_like(x), x, x*x - quadratic_center])
    probability = expit(design @ fit.mode)
    score = design.T @ (probability - retained) + fit.mode / np.asarray([5.0, 2.5, 2.5])**2
    assert np.linalg.norm(score, ord=np.inf) < 1e-9
