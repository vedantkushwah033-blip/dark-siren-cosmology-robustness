"""Regression tests for the cross-gamma pilot's complete-cell gate."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_quadratic_logit_curvature_pilot import _require_complete_replicates


def test_complete_gamma_cell_is_accepted():
    _require_complete_replicates(
        {"planned_replicates": 20, "successful_replicates": 20, "failed_replicates": 0},
        gamma=0.0,
    )


@pytest.mark.parametrize(
    "report",
    [
        {"planned_replicates": 20, "successful_replicates": 19, "failed_replicates": 1},
        {"planned_replicates": 20, "successful_replicates": 18, "failed_replicates": 2},
    ],
)
def test_incomplete_gamma_cell_is_rejected(report):
    with pytest.raises(AssertionError, match="incomplete gamma=1.0 cell"):
        _require_complete_replicates(report, gamma=1.0)
