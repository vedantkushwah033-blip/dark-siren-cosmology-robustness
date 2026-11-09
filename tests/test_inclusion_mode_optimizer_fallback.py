from types import SimpleNamespace

import numpy as np
from scipy.optimize import minimize as scipy_minimize

import dark_siren_robustness.inclusion_posterior as inclusion_posterior


def test_inclusion_mode_retries_line_search_after_trust_region_stagnation(monkeypatch):
    rng = np.random.default_rng(99173)
    redshift = rng.uniform(0.01, 0.40, 240)
    x = (redshift - 0.205) / 0.195
    probability = 1.0 / (1.0 + np.exp(-(-0.35 + 0.8 * x)))
    retained = rng.binomial(1, probability)

    calls = []

    def controlled_minimize(fun, x0, **kwargs):
        calls.append(kwargs.get("method"))
        if kwargs.get("method") == "trust-exact":
            # Simulate a trust-region termination away from the mode. The
            # fallback must refine this point rather than immediately failing.
            return SimpleNamespace(
                x=np.asarray(x0, dtype=float),
                message="simulated trust-region stagnation",
            )
        return scipy_minimize(fun, x0, **kwargs)

    monkeypatch.setattr(inclusion_posterior, "minimize", controlled_minimize)
    fit = inclusion_posterior.fit_inclusion_posterior_quadrature(
        redshift,
        retained,
        degree=1,
        coordinate_center=0.205,
        coordinate_scale=0.195,
        order=5,
    )

    assert calls == ["trust-exact", "L-BFGS-B"]
    assert np.all(np.isfinite(fit.mode))
    assert np.all(np.isfinite(fit.posterior_mean))
    assert np.isclose(fit.weights.sum(), 1.0)
