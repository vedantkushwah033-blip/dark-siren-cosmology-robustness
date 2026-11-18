"""Unit tests for density-level quadrature diagnostics."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_curvature_density_quadrature import density_total_variation

def test_density_total_variation_is_zero_for_identical_curves():
    grid = np.linspace(0.0, 1.0, 101)
    density = np.ones_like(grid)
    assert density_total_variation(density, density, grid) == pytest.approx(0.0)

def test_density_total_variation_is_symmetric_and_positive_for_different_curves():
    grid = np.linspace(0.0, 1.0, 101)
    left = 2.0 * grid
    right = 2.0 * (1.0 - grid)
    assert density_total_variation(left, right, grid) == pytest.approx(
        density_total_variation(right, left, grid)
    )
    assert density_total_variation(left, right, grid) > 0.0

def test_density_total_variation_rejects_mismatched_shapes():
    with pytest.raises(ValueError, match="matching one-dimensional"):
        density_total_variation(np.ones(3), np.ones(2), np.linspace(0, 1, 2))
