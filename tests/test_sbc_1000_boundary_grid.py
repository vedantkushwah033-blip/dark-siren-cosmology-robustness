import sys
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from audit_sbc_1000_boundary_grid import density_l1, grid_pairs, posterior, seed
from dark_siren_robustness.generative_selection import generate_synthetic_population


def test_targeted_grid_posterior_is_normalized_and_monotone():
    sim = generate_synthetic_population(seed=77121, n_parent=300, n_detected=5, h0_true=70.0)
    z = sim["parent_truth"].redshift.to_numpy(float)
    cfg = sim["config"]
    zmin, zmax = cfg["redshift_domain"]
    weights = 1.0 + 1.5 * (z - zmin) / (zmax - zmin)
    distances = sim["gw_events_observed"].distance_mpc.to_numpy(float)
    grid = np.linspace(20.0, 300.0, 281)
    density, cdf, quantiles = posterior(
        z, weights, distances, grid, cfg["sigma_mpc"], cfg["detection_threshold_mpc"]
    )
    assert np.isclose(trapezoid(density, grid), 1.0, atol=1e-12)
    assert np.all(np.diff(cdf) >= 0)
    assert np.all(np.diff(quantiles) >= 0)
    assert np.isclose(cdf[-1], 1.0)


def test_density_l1_is_small_for_same_smooth_density():
    coarse = np.linspace(20.0, 300.0, 561)
    fine = np.linspace(20.0, 300.0, 1121)
    dc = np.exp(-0.5 * ((coarse - 70.0) / 8.0) ** 2)
    df = np.exp(-0.5 * ((fine - 70.0) / 8.0) ** 2)
    dc /= trapezoid(dc, coarse)
    df /= trapezoid(df, fine)
    assert density_l1(coarse, dc, fine, df) < 0.001


def test_seed_stream_is_deterministic_and_separated():
    assert seed(2026111919, 208, 0) == seed(2026111919, 208, 0)
    assert seed(2026111919, 208, 0) != seed(2026111919, 208, 1)


def test_audit_includes_all_requested_grid_pairs():
    assert grid_pairs([1121, 2241, 4481, 8961, 17921]) == [
        (1121, 2241), (2241, 4481), (4481, 8961), (8961, 17921)
    ]
