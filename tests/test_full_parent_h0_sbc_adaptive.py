import sys
from pathlib import Path

import numpy as np
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_full_parent_h0_sbc_adaptive import _adaptive_reference


def test_adaptive_cdf_matches_truncated_normal_reference():
    lower, upper, mean, sd, truth = 20.0, 300.0, 22.0, 0.7, 22.4
    log_density = lambda x: -0.5 * ((x - mean) / sd) ** 2
    grid = np.linspace(lower, upper, 2241)
    result = _adaptive_reference(log_density, lower, upper, truth, grid)
    norm_const = norm.cdf((upper - mean) / sd) - norm.cdf((lower - mean) / sd)
    expected_rank = (norm.cdf((truth - mean) / sd) - norm.cdf((lower - mean) / sd)) / norm_const
    assert abs(result["adaptive_rank"] - expected_rank) < 2e-7
    assert np.all(np.diff(result["adaptive_quantiles"]) > 0)
    assert result["adaptive_relative_quad_error_bound"] < 1e-7


def test_adaptive_cdf_handles_monotone_boundary_posterior():
    lower, upper, truth = 20.0, 300.0, 20.4
    log_density = lambda x: -(x - lower) / 0.3
    grid = np.linspace(lower, upper, 2241)
    result = _adaptive_reference(log_density, lower, upper, truth, grid)
    expected = 1.0 - np.exp(-(truth - lower) / 0.3)
    assert abs(result["adaptive_rank"] - expected) < 2e-6
    assert result["adaptive_quantiles"][5] < 21.0


def test_seed_helper_returns_deterministic_python_integer():
    from audit_full_parent_h0_sbc_adaptive import _seed
    first = _seed(2026111907, 100, 0)
    assert isinstance(first, int)
    assert first == _seed(2026111907, 100, 0)
    assert first != _seed(2026111907, 100, 1)


def test_lowest_q975_selection_is_deterministic():
    import pandas as pd
    from audit_full_parent_h0_sbc_adaptive import _select_replicates
    frame = pd.DataFrame({"replicate": [0, 1, 2, 3], "q975": [40.0, 22.0, 30.0, 21.0]})
    assert _select_replicates(frame, lowest_q975=2) == (3, 1)
    assert _select_replicates(frame, replicates=[2, 0]) == (2, 0)
