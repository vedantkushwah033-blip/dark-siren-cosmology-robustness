from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.special import expit

from dark_siren_robustness.generative_selection import (
    calibrate_quadratic_logit_intercept,
    generate_quadratic_logit_inclusion_samples,
    quadratic_logit_inclusion_probability,
)


@pytest.mark.parametrize("gamma", [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0])
def test_quadratic_logit_intercept_preserves_parent_average_rate(gamma: float) -> None:
    result = calibrate_quadratic_logit_intercept(gamma=gamma)
    assert result["x2_center"] == pytest.approx(0.2845005579191241, abs=2e-10)
    assert result["target_rate"] == pytest.approx(0.5101919962291073, abs=2e-10)
    assert result["achieved_rate"] == pytest.approx(result["target_rate"], abs=2e-12)
    assert np.isfinite(result["alpha"])


def test_quadratic_logit_probability_uses_declared_centered_polynomial() -> None:
    z = np.linspace(0.01, 0.40, 101)
    gamma = 1.0
    calibration = calibrate_quadratic_logit_intercept(gamma=gamma)
    x = (z - 0.205) / 0.195
    expected = expit(
        calibration["alpha"]
        - 3.9 * x
        + gamma * (x**2 - calibration["x2_center"])
    )
    actual = quadratic_logit_inclusion_probability(
        z, gamma=gamma, beta=-3.9, alpha=calibration["alpha"],
        target_rate=calibration["target_rate"],
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-14)
    assert np.all((actual > 0.0) & (actual < 1.0))


def test_quadratic_logit_generator_is_reproducible_nested_and_truth_separated() -> None:
    common = {
        "population_seed": 2026111903,
        "validation_seed": 2026111904,
        "n_parent": 1200,
        "n_validation_frame": 500,
        "gamma": 0.5,
    }
    small = generate_quadratic_logit_inclusion_samples(
        **common, n_validation=100
    )
    medium = generate_quadratic_logit_inclusion_samples(
        **common, n_validation=300
    )
    large = generate_quadratic_logit_inclusion_samples(
        **common, n_validation=500
    )

    # The same population seed must give the same parent catalogue realization.
    for result in (medium, large):
        pd.testing.assert_frame_equal(
            small["parent_truth"], result["parent_truth"]
        )
        pd.testing.assert_frame_equal(
            small["catalogue_observed"], result["catalogue_observed"]
        )

    # The validation frame and retention outcomes are identical; only the
    # permutation-prefix sample size changes.
    truth_columns = [
        "frame_id", "redshift", "inclusion_probability", "retained"
    ]
    for result in (medium, large):
        pd.testing.assert_frame_equal(
            small["validation_frame_truth"][truth_columns],
            result["validation_frame_truth"][truth_columns],
        )
    ids_small = small["validation_observed"]["validation_id"].to_numpy()
    ids_medium = medium["validation_observed"]["validation_id"].to_numpy()
    ids_large = large["validation_observed"]["validation_id"].to_numpy()
    np.testing.assert_array_equal(ids_small, ids_medium[:100])
    np.testing.assert_array_equal(ids_medium, ids_large[:300])

    assert list(small["catalogue_observed"].columns) == ["object_id", "redshift"]
    assert list(small["validation_observed"].columns) == [
        "validation_id", "redshift", "retained"
    ]
    assert "inclusion_probability" in small["parent_truth"]
    assert "inclusion_probability" in small["validation_frame_truth"]
    assert "inclusion_probability" not in small["validation_observed"]
    assert "inclusion_probability" not in small["catalogue_observed"]



def test_h0_auditor_rejects_stale_catalogue_count_metadata():
    from scripts.audit_quadratic_logit_inclusion_h0_smoke import (
        _validate_population_counts,
    )

    parent = pd.DataFrame({"object_id": [0, 1, 2]})
    catalogue = pd.DataFrame({"object_id": [0, 2]})
    validation_truth = pd.DataFrame({"frame_id": [0, 1, 2, 3]})
    validation = pd.DataFrame({"validation_id": [0, 2]})
    analysis_config = {"n_parent": 3, "n_catalogue_observed": 2}
    validation_config = {"n_validation_frame": 4, "n_validation": 2}
    _validate_population_counts(
        analysis_config, parent, catalogue,
        validation_config, validation_truth, validation,
    )

    analysis_config["n_catalogue_observed"] = 3
    with pytest.raises(AssertionError, match="n_catalogue_observed"):
        _validate_population_counts(
            analysis_config, parent, catalogue,
            validation_config, validation_truth, validation,
        )

    analysis_config["n_catalogue_observed"] = 2
    validation_config["n_validation"] = 3
    with pytest.raises(AssertionError, match="n_validation does not match"):
        _validate_population_counts(
            analysis_config, parent, catalogue,
            validation_config, validation_truth, validation,
        )
