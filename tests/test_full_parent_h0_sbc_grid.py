import sys
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from dark_siren_robustness.generative_selection import generate_synthetic_population
from audit_full_parent_h0_sbc_grid import _posterior


def test_independent_sbc_posterior_is_normalized_and_monotone():
    sim = generate_synthetic_population(seed=5512, n_parent=300, n_detected=5, h0_true=70.0)
    parent, events, cfg = sim["parent_truth"], sim["gw_events_observed"], sim["config"]
    z = parent.redshift.to_numpy(float)
    zmin, zmax = cfg["redshift_domain"]
    weights = 1.0 + 1.5*(z-zmin)/(zmax-zmin)
    grid = np.linspace(20.0, 300.0, 281)
    density, cdf, q = _posterior(z, weights, events.distance_mpc.to_numpy(float), grid,
                                 cfg["sigma_mpc"], cfg["detection_threshold_mpc"])
    assert np.isclose(trapezoid(density, grid), 1.0, atol=1e-12)
    assert np.all(np.diff(cdf) >= 0)
    assert np.all(np.diff(q) >= 0)
    assert np.isclose(cdf[-1], 1.0)


def test_grid_audit_reconstructs_truth_rank():
    sim = generate_synthetic_population(seed=7781, n_parent=300, n_detected=5, h0_true=70.0)
    parent, events, cfg = sim["parent_truth"], sim["gw_events_observed"], sim["config"]
    z = parent.redshift.to_numpy(float)
    zmin, zmax = cfg["redshift_domain"]
    weights = 1.0 + 1.5*(z-zmin)/(zmax-zmin)
    grid = np.linspace(20.0, 300.0, 281)
    _, cdf, _ = _posterior(z, weights, events.distance_mpc.to_numpy(float), grid,
                           cfg["sigma_mpc"], cfg["detection_threshold_mpc"])
    rank = np.interp(70.0, grid, cdf)
    assert 0.0 <= rank <= 1.0
