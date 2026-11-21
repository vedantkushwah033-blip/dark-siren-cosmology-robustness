import numpy as np
import pytest

import dark_siren_robustness.arm_cell_execution as execution
from dark_siren_robustness.five_arm_comparison import build_five_arm_weights
from dark_siren_robustness.five_arm_kde_likelihood import matched_five_arm_kde_selected_log_likelihood


def _fixture():
    parent = np.linspace(0.0101, 0.3999, 600)
    catalogue = parent[::2]
    validation = np.linspace(0.0101, 0.3999, 180)
    retained = (np.arange(validation.size) % 3 != 0).astype(int)
    parent_base = 1.0 + 1.5 * (parent - 0.01) / 0.39
    catalogue_base = 1.0 + 1.5 * (catalogue - 0.01) / 0.39
    arms = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.1, quadrature_order=5,
        parent_base_weights=parent_base,
        catalogue_base_weights=catalogue_base,
    )
    events = [
        {"gw_distance": 275.0, "gw_sigma": 140.0},
        {"gw_distance": 340.0, "gw_sigma": 140.0},
    ]
    return parent, catalogue, validation, retained, arms, events, parent_base, catalogue_base


def _evaluate(arm, fixture):
    parent, catalogue, validation, retained, weights, events, parent_base, catalogue_base = fixture
    return execution.evaluate_arm_cell(
        arm, parent_redshift=parent, catalogue_redshift=catalogue,
        validation_redshift=validation, validation_retained=retained,
        true_catalogue_inclusion_probability=weights.arm_c_true_inclusion_probability,
        events=events, H0_grid=np.linspace(55.0, 85.0, 41),
        redshift_grid=np.linspace(0.01, 0.40, 501),
        merger_weight=np.ones(501), h0_true=70.0,
        h0_prior_min=55.0, h0_prior_max=85.0,
        detection_threshold_mpc=250.0, scott_multiplier=0.7,
        quadrature_order=5, inclusion_degree=1,
        parent_base_weights=parent_base, catalogue_base_weights=catalogue_base,
    )


def test_arm_cell_results_match_existing_joint_five_arm_implementation():
    fixture = _fixture()
    parent, catalogue, validation, retained, arms, events, _, _ = fixture
    h0 = np.linspace(55.0, 85.0, 41)
    z = np.linspace(0.01, 0.40, 501)
    joint = matched_five_arm_kde_selected_log_likelihood(
        h0, arms, events, detection_threshold_mpc=250.0,
        redshift_grid=z, merger_weight=np.ones_like(z), scott_multiplier=0.7,
    )
    expected = {
        "A_catalogue_only": joint.arm_a_log_likelihood,
        "B_parent_oracle": joint.arm_b_log_likelihood,
        "C_known_probability_ipw": joint.arm_c_log_likelihood,
        "D_map_plugin_ipw": joint.arm_d_log_likelihood,
        "E_uncertainty_integrated_ipw": joint.arm_e_marginal_log_likelihood,
    }
    for arm, target in expected.items():
        result = _evaluate(arm, fixture)
        assert result.status == "SUCCESS"
        np.testing.assert_allclose(result.marginal_log_likelihood, target, rtol=1e-11, atol=1e-11)
        assert np.isfinite(result.posterior_summary["posterior_median"])
        assert result.component_density.shape[0] == (arms.arm_e_node_weights.shape[0] if arm.startswith("E_") else 1)


@pytest.mark.parametrize("arm", ["A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw"])
def test_reference_arms_do_not_fit_validation_inclusion_model(monkeypatch, arm):
    fixture = _fixture()
    def forbidden(*args, **kwargs):
        raise AssertionError("A/B/C must not fit the validation inclusion model")
    monkeypatch.setattr(execution, "fit_inclusion_posterior_quadrature", forbidden)
    result = _evaluate(arm, fixture)
    assert result.status == "SUCCESS"


def test_failure_in_one_arm_does_not_change_another_arm_result():
    fixture = _fixture()
    baseline = _evaluate("A_catalogue_only", fixture)
    with pytest.raises(ValueError, match="true inclusion probabilities"):
        execution.evaluate_arm_cell(
            "C_known_probability_ipw",
            parent_redshift=fixture[0], catalogue_redshift=fixture[1],
            validation_redshift=fixture[2], validation_retained=fixture[3],
            true_catalogue_inclusion_probability=np.full(fixture[1].size, np.nan),
            events=fixture[5], H0_grid=np.linspace(55.0, 85.0, 41),
            redshift_grid=np.linspace(0.01, 0.40, 501), merger_weight=np.ones(501),
            h0_true=70.0, h0_prior_min=55.0, h0_prior_max=85.0,
            detection_threshold_mpc=250.0, scott_multiplier=0.7,
            parent_base_weights=fixture[6], catalogue_base_weights=fixture[7],
        )
    replay = _evaluate("A_catalogue_only", fixture)
    np.testing.assert_array_equal(baseline.marginal_log_likelihood, replay.marginal_log_likelihood)


def test_arm_cell_rejects_unknown_arm_and_invalid_validation_labels():
    fixture = _fixture()
    with pytest.raises(ValueError, match="unknown arm"):
        _evaluate("ALL_ARMS", fixture)
    args = dict(
        parent_redshift=fixture[0], catalogue_redshift=fixture[1],
        validation_redshift=fixture[2], validation_retained=np.full(fixture[2].size, 2),
        true_catalogue_inclusion_probability=fixture[4].arm_c_true_inclusion_probability,
        events=fixture[5], H0_grid=np.linspace(55.0, 85.0, 41),
        redshift_grid=np.linspace(0.01, 0.40, 501), merger_weight=np.ones(501),
        h0_true=70.0, h0_prior_min=55.0, h0_prior_max=85.0,
        detection_threshold_mpc=250.0, scott_multiplier=0.7,
        parent_base_weights=fixture[6], catalogue_base_weights=fixture[7],
    )
    with pytest.raises(ValueError, match="validation_retained"):
        execution.evaluate_arm_cell("D_map_plugin_ipw", **args)

@pytest.mark.parametrize("arm", ["A_catalogue_only", "B_parent_oracle"])
def test_ab_do_not_require_validation_or_true_inclusion_truth(arm):
    fixture = _fixture()
    result = execution.evaluate_arm_cell(
        arm, parent_redshift=fixture[0], catalogue_redshift=fixture[1],
        events=fixture[5], H0_grid=np.linspace(55.0, 85.0, 41),
        redshift_grid=np.linspace(0.01, 0.40, 501),
        merger_weight=np.ones(501), h0_true=70.0,
        h0_prior_min=55.0, h0_prior_max=85.0,
        detection_threshold_mpc=250.0, scott_multiplier=0.7,
        validation_redshift=None, validation_retained=None,
        true_catalogue_inclusion_probability=None,
        parent_base_weights=fixture[6], catalogue_base_weights=fixture[7],
    )
    assert result.status == "SUCCESS"


def test_c_does_not_require_validation_but_d_and_e_do():
    fixture = _fixture()
    common = dict(
        parent_redshift=fixture[0], catalogue_redshift=fixture[1],
        events=fixture[5], H0_grid=np.linspace(55.0, 85.0, 41),
        redshift_grid=np.linspace(0.01, 0.40, 501),
        merger_weight=np.ones(501), h0_true=70.0,
        h0_prior_min=55.0, h0_prior_max=85.0,
        detection_threshold_mpc=250.0, scott_multiplier=0.7,
        validation_redshift=None, validation_retained=None,
        parent_base_weights=fixture[6], catalogue_base_weights=fixture[7],
    )
    c = execution.evaluate_arm_cell(
        "C_known_probability_ipw",
        true_catalogue_inclusion_probability=fixture[4].arm_c_true_inclusion_probability,
        **common,
    )
    assert c.status == "SUCCESS"
    for arm in ("D_map_plugin_ipw", "E_uncertainty_integrated_ipw"):
        with pytest.raises(ValueError, match="D/E require validation"):
            execution.evaluate_arm_cell(
                arm, true_catalogue_inclusion_probability=None, **common,
            )

def test_reference_arms_ignore_other_arms_base_weight_inputs():
    fixture = _fixture()
    common = dict(
        parent_redshift=fixture[0], catalogue_redshift=fixture[1],
        events=fixture[5], H0_grid=np.linspace(55.0, 85.0, 41),
        redshift_grid=np.linspace(0.01, 0.40, 501),
        merger_weight=np.ones(501), h0_true=70.0,
        h0_prior_min=55.0, h0_prior_max=85.0,
        detection_threshold_mpc=250.0, scott_multiplier=0.7,
    )
    a = execution.evaluate_arm_cell(
        "A_catalogue_only", parent_base_weights=np.full(fixture[0].size, np.nan),
        **common,
    )
    b = execution.evaluate_arm_cell(
        "B_parent_oracle", catalogue_base_weights=np.full(fixture[1].size, np.nan),
        **common,
    )
    assert a.status == "SUCCESS"
    assert b.status == "SUCCESS"

def test_a_and_b_require_only_their_own_host_population():
    fixture = _fixture()
    common = dict(
        events=fixture[5], H0_grid=np.linspace(55.0, 85.0, 41),
        redshift_grid=np.linspace(0.01, 0.40, 501),
        merger_weight=np.ones(501), h0_true=70.0,
        h0_prior_min=55.0, h0_prior_max=85.0,
        detection_threshold_mpc=250.0, scott_multiplier=0.7,
    )
    a = execution.evaluate_arm_cell(
        "A_catalogue_only", catalogue_redshift=fixture[1], **common,
    )
    b = execution.evaluate_arm_cell(
        "B_parent_oracle", parent_redshift=fixture[0], **common,
    )
    assert a.status == "SUCCESS"
    assert b.status == "SUCCESS"

