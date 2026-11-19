import numpy as np

from scripts.analyze_five_arm_weight_error_association import spearman_bootstrap


def test_spearman_bootstrap_detects_monotone_association_deterministically():
    x = np.arange(100, dtype=float)
    y = x ** 2
    first = spearman_bootstrap(x, y, resamples=500, seed=42)
    second = spearman_bootstrap(x, y, resamples=500, seed=42)
    assert first == second
    assert first["status"] == "PASS"
    assert first["n"] == 100
    assert first["rho"] > 0.99
    assert first["bootstrap_low"] > 0.95


def test_spearman_bootstrap_fails_closed_for_constant_diagnostic():
    x = np.ones(30)
    y = np.arange(30, dtype=float)
    result = spearman_bootstrap(x, y, resamples=200, seed=7)
    assert result["status"] == "NON_EVALUABLE"
    assert np.isnan(result["rho"])
