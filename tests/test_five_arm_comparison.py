import numpy as np
import pytest

from dark_siren_robustness.five_arm_comparison import (
    build_five_arm_weights,
    weight_concentration_diagnostics,
)
from dark_siren_robustness.generative_selection import (
    catalogue_inclusion_probability,
    generate_quadratic_logit_inclusion_samples,
)


def _fixture():
    rng = np.random.default_rng(20261113)
    parent = np.sort(rng.uniform(0.01, 0.40, 600))
    probability = catalogue_inclusion_probability(parent, floor=0.10)
    catalogue = parent[rng.random(parent.size) < probability]
    validation = np.sort(rng.uniform(0.01, 0.40, 240))
    retained = rng.binomial(
        1, catalogue_inclusion_probability(validation, floor=0.10)
    )
    return parent, catalogue, validation, retained


def test_five_arm_weights_have_matched_shapes_and_normalized_quadrature():
    parent, catalogue, validation, retained = _fixture()
    result = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.10, quadrature_order=7,
    )
    assert result.arm_a_catalogue_only.shape == catalogue.shape
    assert result.arm_b_parent_oracle.shape == parent.shape
    assert result.arm_c_known_probability.shape == catalogue.shape
    assert result.arm_d_map_plugin.shape == catalogue.shape
    assert result.arm_e_node_weights.shape == (
        result.inclusion_posterior.nodes.shape[0], catalogue.size
    )
    assert result.arm_e_node_inclusion_probability.shape == result.arm_e_node_weights.shape
    assert np.isclose(result.arm_e_quadrature_weights.sum(), 1.0)
    assert np.all(result.arm_c_known_probability >= 1.0)
    assert np.all(result.arm_d_map_plugin >= 1.0)
    assert np.all(result.arm_e_node_weights >= 1.0)


def test_arm_a_and_b_are_equal_weight_oracles_on_their_respective_samples():
    parent, catalogue, validation, retained = _fixture()
    result = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.10, quadrature_order=5,
    )
    assert np.array_equal(result.arm_a_catalogue_only, np.ones(catalogue.size))
    assert np.array_equal(result.arm_b_parent_oracle, np.ones(parent.size))


def test_weight_diagnostics_are_finite_and_report_concentration():
    result = weight_concentration_diagnostics([1.0, 1.0, 2.0, 6.0])
    assert result["kish_ess"] == pytest.approx(100.0 / 42.0)
    assert result["kish_ess_fraction"] == pytest.approx((100.0 / 42.0) / 4.0)
    p = np.array([0.1, 0.1, 0.2, 0.6])
    assert result["entropy_ess"] == pytest.approx(np.exp(-np.sum(p * np.log(p))))
    assert result["maximum_normalized_weight"] == pytest.approx(0.6)
    assert result["top_1_percent_mass"] == pytest.approx(0.6)
    assert result["top_5_percent_mass"] == pytest.approx(0.6)


def test_weight_diagnostics_reject_nonpositive_or_nonfinite_weights():
    with pytest.raises(ValueError):
        weight_concentration_diagnostics([1.0, 0.0])
    with pytest.raises(ValueError):
        weight_concentration_diagnostics([1.0, np.inf])


def test_five_arm_builder_rejects_mismatched_validation_labels():
    parent, catalogue, validation, retained = _fixture()
    with pytest.raises(ValueError, match="match validation_redshift"):
        build_five_arm_weights(
            parent, catalogue, validation, retained[:-1],
            catalogue_floor=0.10,
        )

def test_weight_diagnostics_reject_positive_mass_underflow():
    with pytest.raises(ArithmeticError, match="underflowed positive mass"):
        weight_concentration_diagnostics([1e-320, 1e308])


def test_five_arm_weights_apply_common_base_population_mass_to_all_arms():
    parent, catalogue, validation, retained = _fixture()
    parent_base = 1.0 + 1.5 * (parent - 0.01) / 0.39
    catalogue_base = 1.0 + 1.5 * (catalogue - 0.01) / 0.39
    result = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.10,
        parent_base_weights=parent_base,
        catalogue_base_weights=catalogue_base,
        quadrature_order=5,
    )
    true_probability = catalogue_inclusion_probability(catalogue, floor=0.10)
    assert np.allclose(result.arm_a_catalogue_only, catalogue_base)
    assert np.allclose(result.arm_b_parent_oracle, parent_base)
    assert np.allclose(result.arm_c_known_probability, catalogue_base / true_probability)
    assert np.allclose(
        result.arm_d_map_plugin,
        catalogue_base / result.arm_d_map_inclusion_probability,
    )
    assert np.allclose(
        result.arm_e_node_weights,
        catalogue_base[None, :] / result.arm_e_node_inclusion_probability,
    )


@pytest.mark.parametrize(
    "parent_base,catalogue_base,match",
    [
        (np.ones(599), None, "parent_base_weights"),
        (None, np.ones(2), "catalogue_base_weights"),
        (np.r_[np.ones(599), np.nan], None, "parent_base_weights"),
        (None, np.r_[np.ones(2), 0.0], "catalogue_base_weights"),
    ],
)
def test_five_arm_builder_uses_explicit_truth_probabilities_for_nondefault_dgp():
    parent, catalogue, validation, retained = _fixture()
    true_probability = 0.35 + 0.25 * (catalogue - 0.01) / 0.39
    result = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.10,
        true_catalogue_inclusion_probability=true_probability,
        quadrature_order=5,
    )
    assert np.array_equal(
        result.arm_c_true_inclusion_probability, true_probability
    )
    assert np.allclose(
        result.arm_c_known_probability, 1.0 / true_probability
    )
    # The validation-only D/E arms must remain fitted from observable labels,
    # rather than being overwritten by the supplied simulation truth.
    assert np.allclose(
        result.arm_d_map_plugin,
        1.0 / result.arm_d_map_inclusion_probability,
    )
    assert np.allclose(
        result.arm_e_node_weights,
        1.0 / result.arm_e_node_inclusion_probability,
    )


@pytest.mark.parametrize(
    "kind",
    ["wrong_shape", "zero", "above_one", "nan"],
)
def test_quadratic_logit_dgp_known_arm_uses_truth_but_estimated_arms_use_validation():
    sample = generate_quadratic_logit_inclusion_samples(
        population_seed=20261122,
        validation_seed=20261123,
        n_parent=700,
        n_validation_frame=300,
        n_validation=240,
        gamma=1.0,
    )
    parent = sample["parent_truth"]
    catalogue = sample["catalogue_observed"]
    validation = sample["validation_observed"]
    true_probability = parent.loc[
        parent.retained.astype(bool), "inclusion_probability"
    ].to_numpy()
    result = build_five_arm_weights(
        parent.redshift.to_numpy(),
        catalogue.redshift.to_numpy(),
        validation.redshift.to_numpy(),
        validation.retained.to_numpy(),
        catalogue_floor=0.10,
        true_catalogue_inclusion_probability=true_probability,
        quadrature_order=5,
    )
    assert np.allclose(result.arm_c_true_inclusion_probability, true_probability)
    assert np.allclose(result.arm_c_known_probability, 1.0 / true_probability)
    assert not np.allclose(
        result.arm_d_map_inclusion_probability, true_probability,
        rtol=1e-3, atol=1e-3,
    )
    assert np.allclose(
        result.arm_d_map_plugin,
        1.0 / result.arm_d_map_inclusion_probability,
    )


def test_five_arm_builder_rejects_invalid_explicit_truth_probabilities(kind):
    parent, catalogue, validation, retained = _fixture()
    if kind == "wrong_shape":
        probability = np.ones(catalogue.size + 1)
    elif kind == "zero":
        probability = np.zeros(catalogue.size)
    elif kind == "above_one":
        probability = np.full(catalogue.size, 1.01)
    else:
        probability = np.full(catalogue.size, np.nan)
    with pytest.raises(ValueError, match="true_catalogue_inclusion_probability"):
        build_five_arm_weights(
            parent, catalogue, validation, retained,
            catalogue_floor=0.10,
            true_catalogue_inclusion_probability=probability,
            quadrature_order=5,
        )


def test_five_arm_builder_rejects_invalid_base_population_weights(
    parent_base, catalogue_base, match
):
    parent, catalogue, validation, retained = _fixture()
    with pytest.raises(ValueError, match=match):
        build_five_arm_weights(
            parent, catalogue, validation, retained,
            catalogue_floor=0.10,
            parent_base_weights=parent_base,
            catalogue_base_weights=catalogue_base,
            quadrature_order=5,
        )
