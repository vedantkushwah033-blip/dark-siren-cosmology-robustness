import numpy as np
import pytest
from scipy.special import expit

from dark_siren_robustness.independent_audit import generate_independent_inclusion_audit


def test_independent_audit_is_reproducible_and_has_separate_observable_truth_tables():
    kwargs = dict(
        seed=20261031, n_audit=300, z_min=0.01, z_max=0.40,
        intercept=0.25, slope=-1.2, quadratic=1.2,
    )
    observed_a, truth_a = generate_independent_inclusion_audit(**kwargs)
    observed_b, truth_b = generate_independent_inclusion_audit(**kwargs)

    assert observed_a.equals(observed_b)
    assert truth_a.equals(truth_b)
    assert list(observed_a.columns) == ["audit_object_id", "redshift", "included"]
    assert list(truth_a.columns) == ["audit_object_id", "true_inclusion_probability"]
    assert observed_a["audit_object_id"].is_unique
    assert observed_a["audit_object_id"].equals(truth_a["audit_object_id"])
    assert observed_a["redshift"].between(0.01, 0.40).all()
    assert set(observed_a["included"].unique()).issubset({0, 1})
    assert truth_a["true_inclusion_probability"].between(0, 1, inclusive="neither").all()


def test_independent_audit_ids_are_disjoint_from_integer_source_frame_ids():
    observed, _truth = generate_independent_inclusion_audit(
        seed=20261032, n_audit=100, z_min=0.01, z_max=0.40,
        intercept=0.25, slope=-1.2,
    )
    source_ids = set(range(4000))
    assert source_ids.isdisjoint(set(observed["audit_object_id"]))


def test_independent_audit_probabilities_follow_declared_quadratic_logit():
    observed, truth = generate_independent_inclusion_audit(
        seed=20261033, n_audit=100, z_min=0.01, z_max=0.40,
        intercept=0.25, slope=-1.2, quadratic=1.2,
    )
    z = observed["redshift"].to_numpy()
    x = (z - 0.205) / 0.195
    expected = expit(0.25 - 1.2 * x + 1.2 * (x**2 - 1.0 / 3.0))
    np.testing.assert_allclose(truth["true_inclusion_probability"], expected, rtol=0, atol=1e-14)


@pytest.mark.parametrize("overrides", [
    {"seed": -1},
    {"n_audit": 9},
    {"z_min": 0.4, "z_max": 0.01},
    {"intercept": float("nan")},
    {"grid_size": 100},
])
def test_independent_audit_rejects_invalid_configuration(overrides):
    kwargs = dict(
        seed=1, n_audit=100, z_min=0.01, z_max=0.40,
        intercept=0.25, slope=-1.2, quadratic=0.0, grid_size=20001,
    )
    kwargs.update(overrides)
    with pytest.raises(ValueError):
        generate_independent_inclusion_audit(**kwargs)
