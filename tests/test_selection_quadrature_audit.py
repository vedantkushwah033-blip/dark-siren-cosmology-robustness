import importlib.util
from pathlib import Path

PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_selection_quadrature_audit.py"
SPEC = importlib.util.spec_from_file_location("selection_quadrature_audit", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)

def test_adaptive_quadrature_agrees_with_grid_normalization():
    rows = module.run_audit(grid_size=4001)
    assert len(rows) == 3
    for row in rows:
        assert row["grid_normalization_constant"] > 0
        assert row["adaptive_quadrature_normalization_constant"] > 0
        assert row["quadrature_estimated_absolute_error"] >= 0
        assert row["relative_absolute_difference"] < 1e-4
        assert row["pass_relative_tolerance_1e-4"] is True

def test_quadrature_grid_requires_prior_endpoints_on_grid():
    try:
        module.run_audit(grid_size=1002)
    except ValueError as exc:
        assert "prior endpoints" in str(exc)
    else:
        raise AssertionError("Expected invalid grid-size rejection.")
