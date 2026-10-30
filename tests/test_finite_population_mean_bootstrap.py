from __future__ import annotations

import numpy as np
import pytest

from dark_siren_robustness.finite_population_mean_bootstrap import (
    studentized_pseudopopulation_srswor_bootstrap,
)


def test_studentized_srswor_bootstrap_is_deterministic_and_accounted():
    sample = np.array([0.2, 0.5, 0.7, 1.1, 1.4, 2.0, 2.2, 3.1])
    first = studentized_pseudopopulation_srswor_bootstrap(
        sample, 80, n_bootstrap=300, seed=20261031
    )
    second = studentized_pseudopopulation_srswor_bootstrap(
        sample, 80, n_bootstrap=300, seed=20261031
    )
    for key in (
        "bootstrap_t_statistics", "bootstrap_pseudopopulation_means",
        "bootstrap_sample_means", "bootstrap_design_ses", "success_mask",
    ):
        np.testing.assert_array_equal(first[key], second[key])
    assert first["schema"] == "studentized_finite_population_mean_srswor_v1"
    assert first["n_attempted"] == 300
    assert first["n_successful"] + first["n_failed"] == 300
    assert first["status"] == "ok"
    assert first["finite_population_correction_variance"] == pytest.approx(0.9)


def test_studentized_statistic_and_interval_are_reconstructible():
    sample = np.array([1.0, 1.5, 2.0, 2.5, 4.0, 5.0])
    result = studentized_pseudopopulation_srswor_bootstrap(
        sample, 30, n_bootstrap=200, seed=71, confidence=0.90
    )
    success = result["success_mask"]
    reconstructed = (
        result["bootstrap_sample_means"][success]
        - result["bootstrap_pseudopopulation_means"][success]
    ) / result["bootstrap_design_ses"][success]
    np.testing.assert_allclose(
        result["bootstrap_t_statistics"][success], reconstructed,
        rtol=1e-13, atol=1e-13,
    )
    q_low, q_high = np.quantile(reconstructed, [0.05, 0.95])
    se = result["estimated_design_se"]
    estimate = np.mean(sample)
    np.testing.assert_allclose(
        result["interval"], [estimate - q_high * se, estimate - q_low * se]
    )


def test_studentized_srswor_uses_sample_variance_and_fpc():
    sample = np.array([1.0, 2.0, 4.0, 8.0])
    result = studentized_pseudopopulation_srswor_bootstrap(
        sample, 20, n_bootstrap=20, seed=9
    )
    expected_variance = np.var(sample, ddof=1)
    expected_fpc = 1.0 - len(sample) / 20
    assert result["sample_variance"] == pytest.approx(expected_variance)
    assert result["finite_population_correction_variance"] == pytest.approx(expected_fpc)
    assert result["estimated_design_se"] == pytest.approx(
        np.sqrt(expected_fpc * expected_variance / len(sample))
    )


@pytest.mark.parametrize(
    "sample,population_size,kwargs",
    [
        ([1.0], 10, {}),
        ([1.0, np.nan], 10, {}),
        ([1.0, 2.0], 2, {}),
        ([1.0, 2.0], 10, {"n_bootstrap": 1}),
        ([1.0, 2.0], 10, {"seed": -1}),
        ([1.0, 2.0], 10, {"confidence": 1.0}),
        ([1.0, 2.0], 10, {"confidence": 0.0}),
    ],
)
def test_studentized_srswor_rejects_invalid_inputs(sample, population_size, kwargs):
    with pytest.raises((TypeError, ValueError)):
        studentized_pseudopopulation_srswor_bootstrap(
            sample, population_size, **kwargs
        )


def test_studentized_srswor_rejects_boolean_and_fractional_design_integers():
    with pytest.raises(TypeError):
        studentized_pseudopopulation_srswor_bootstrap([1.0, 2.0], True)
    with pytest.raises(TypeError):
        studentized_pseudopopulation_srswor_bootstrap([1.0, 2.0], 10.5)
