import importlib.util
from pathlib import Path

import numpy as np
import pytest

PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_noise_scale_consistency_audit.py"
SPEC = importlib.util.spec_from_file_location("noise_scale_audit", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_posterior_is_normalized_and_finite_for_both_modes():
    grid = np.linspace(30.0, 120.0, 1801)
    for mode in ("fixed_reported_sigma", "sigma_depends_on_h0"):
        posterior = module.posterior_for_event(
            grid, 0.035, 150.0, 20.0, 0.15, mode
        )
        assert np.all(np.isfinite(posterior))
        assert np.all(posterior >= 0)
        assert np.trapezoid(posterior, grid) == pytest.approx(1.0, abs=1e-6)


def test_h0_dependent_noise_mode_includes_sigma_normalization():
    grid = np.linspace(30.0, 120.0, 1801)
    p = module.posterior_for_event(
        grid, 0.035, 150.0, 20.0, 0.15, "sigma_depends_on_h0"
    )
    # The H0-dependent likelihood is not just the fixed-sigma residual term.
    assert not np.allclose(p, module.posterior_for_event(
        grid, 0.035, 150.0, 20.0, 0.15, "fixed_reported_sigma"
    ))


def test_rejects_unknown_likelihood_mode():
    with pytest.raises(ValueError):
        module.posterior_for_event(
            np.linspace(30, 120, 100), 0.035, 150.0, 20.0, 0.15, "bad"
        )
