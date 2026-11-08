import numpy as np

from dark_siren_robustness.inclusion_posterior import sample_inclusion_posterior


def _synthetic(seed=123, n=240):
    rng = np.random.default_rng(seed)
    z = rng.uniform(0.01, 0.40, n)
    u = (z - 0.01) / 0.39
    probability = 1.0 / (1.0 + np.exp(-(-1.0 + 0.45*u - 0.35*u*u)))
    return z, rng.binomial(1, probability)


def test_reproducible_multichain_posterior_and_prediction():
    z, y = _synthetic()
    first = sample_inclusion_posterior(
        z, y, degree=2, warmup=600, draws_per_chain=600, seed=19
    )
    second = sample_inclusion_posterior(
        z, y, degree=2, warmup=600, draws_per_chain=600, seed=19
    )
    assert first.status == "SUCCESS", first.message
    assert first.draws.shape == (4, 600, 3)
    assert np.array_equal(first.draws, second.draws)
    assert np.all((first.acceptance_rates > 0.05) & (first.acceptance_rates < 0.8))
    assert np.max(first.rhat) < 1.05
    assert np.min(first.ess) >= 100
    predicted = first.predict(np.array([0.02, 0.10, 0.20, 0.39]))
    assert predicted.shape == (4, 600, 4)
    assert np.all(np.isfinite(predicted))
    assert np.all((predicted > 0) & (predicted < 1))


def test_proper_prior_regularizes_complete_separation():
    z = np.linspace(0.01, 0.40, 240)
    y = (z > np.median(z)).astype(int)
    fit = sample_inclusion_posterior(
        z, y, degree=1, warmup=800, draws_per_chain=800, seed=20
    )
    assert fit.map_coefficients is not None
    assert np.all(np.isfinite(fit.map_coefficients))
    assert np.max(np.abs(fit.map_coefficients)) < 25
    assert fit.draws is not None


def test_one_class_validation_is_non_evaluable():
    z = np.linspace(0.01, 0.40, 100)
    fit = sample_inclusion_posterior(
        z, np.ones(100), degree=1, warmup=500, draws_per_chain=500
    )
    assert fit.status == "NON_EVALUABLE"
    assert fit.draws is None


def test_rank_deficient_validation_is_non_evaluable():
    z = np.full(100, 0.2)
    y = np.tile([0, 1], 50)
    fit = sample_inclusion_posterior(
        z, y, degree=1, warmup=500, draws_per_chain=500
    )
    assert fit.status == "NON_EVALUABLE"
    assert fit.draws is None


def test_prior_basis_uses_fixed_domain_coordinates():
    z, y = _synthetic(n=300)
    fit = sample_inclusion_posterior(
        z, y, degree=2, warmup=600, draws_per_chain=600, seed=29
    )
    assert fit.z_min == 0.01
    assert fit.z_max == 0.40
    assert fit.status == "SUCCESS", fit.message
    physical_z = np.array([0.05, 0.20, 0.35])
    u = (physical_z - 0.01) / 0.39
    design = np.column_stack([u**k for k in range(3)])
    expected = 1.0 / (1.0 + np.exp(-np.einsum("cdp,np->cdn", fit.draws, design)))
    assert np.allclose(fit.predict(physical_z), expected)
