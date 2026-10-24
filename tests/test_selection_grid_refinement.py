import importlib.util
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts" / "run_selection_grid_refinement.py"
SPEC = importlib.util.spec_from_file_location("selection_grid_refinement", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)

def test_grid_refinement_output_is_finite_and_complete():
    rows = module.run_audit((501, 1001, 2001))
    assert len(rows) == 9
    assert {r["observed_distance_mpc"] for r in rows} == {70.0, 100.0, 140.0}
    assert all(np.isfinite(float(v)) for r in rows for v in r.values())
    assert all(r["abs_median_diff_vs_finest"] >= 0 for r in rows)
    assert all(r["max_abs_quantile_diff_vs_finest"] >= 0 for r in rows)

def test_grid_refinement_rejects_duplicate_sizes():
    import pytest
    with pytest.raises(ValueError, match="unique and increasing"):
        module.run_audit((501, 501, 1001))
