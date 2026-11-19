import sys
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from dark_siren_robustness.generative_selection import generate_synthetic_population
from audit_full_parent_event_count_grid_refinement import _posterior


def test_independent_posterior_is_normalized_and_quantiles_ordered():
    sim = generate_synthetic_population(seed=4567, n_parent=300, n_detected=30, h0_true=70.0)
    parent, events, cfg = sim["parent_truth"], sim["gw_events_observed"], sim["config"]
    z = parent.redshift.to_numpy(float)
    zmin, zmax = cfg["redshift_domain"]
    weights = 1.0 + 1.5*(z-zmin)/(zmax-zmin)
    grid = np.linspace(20.0, 300.0, 281)
    density, q = _posterior(
        z, weights, events.distance_mpc.to_numpy(float)[:4], grid,
        cfg["sigma_mpc"], cfg["detection_threshold_mpc"],
    )
    assert np.isclose(trapezoid(density, grid), 1.0, atol=1e-12)
    values = [q[k] for k in ("q025", "q050", "q160", "q500", "q840", "q950", "q975")]
    assert np.all(np.diff(values) >= 0)
    assert all(np.isfinite(values))


def test_finer_grid_quantiles_agree_on_smoke_realization():
    sim = generate_synthetic_population(seed=9876, n_parent=300, n_detected=30, h0_true=70.0)
    parent, events, cfg = sim["parent_truth"], sim["gw_events_observed"], sim["config"]
    z = parent.redshift.to_numpy(float)
    zmin, zmax = cfg["redshift_domain"]
    weights = 1.0 + 1.5*(z-zmin)/(zmax-zmin)
    coarse = np.linspace(20.0, 300.0, 141)
    fine = np.linspace(20.0, 300.0, 281)
    _, qc = _posterior(z, weights, events.distance_mpc.to_numpy(float)[:30], coarse,
                       cfg["sigma_mpc"], cfg["detection_threshold_mpc"])
    _, qf = _posterior(z, weights, events.distance_mpc.to_numpy(float)[:30], fine,
                       cfg["sigma_mpc"], cfg["detection_threshold_mpc"])
    assert max(abs(qc[k]-qf[k]) for k in qc) < 0.05
