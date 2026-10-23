import importlib.util
from pathlib import Path

import numpy as np
import pytest

PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_sbc_diagnostic.py"
SPEC = importlib.util.spec_from_file_location("sbc_diagnostic", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_central_interval_is_ordered_and_symmetric_for_uniform_density():
    grid = np.linspace(50.0, 90.0, 4001)
    cdf = (grid - 50.0) / 40.0
    low, high = module.central_interval(grid, cdf, 0.68)
    assert low == pytest.approx(56.4, abs=0.02)
    assert high == pytest.approx(83.6, abs=0.02)
    assert low < high


@pytest.mark.parametrize("level", [0.5, 0.68, 0.9, 0.95])
def test_central_interval_narrows_for_lower_credible_level(level):
    grid = np.linspace(0.0, 1.0, 1001)
    cdf = grid**2
    low, high = module.central_interval(grid, cdf, level)
    assert 0 <= low < high <= 1


def test_central_interval_rejects_invalid_level():
    grid = np.linspace(0.0, 1.0, 11)
    with pytest.raises(ValueError):
        module.central_interval(grid, grid, 1.0)
