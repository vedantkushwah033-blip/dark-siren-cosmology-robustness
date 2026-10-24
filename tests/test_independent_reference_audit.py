import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_independent_reference_audit.py"
SPEC = importlib.util.spec_from_file_location("independent_reference_audit", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_adaptive_quadrature_reference_is_normalized_and_ordered():
    ref = module.reference_summary([0.02, 0.04], 90.0, 3.0, [0.7, 0.3])
    assert 50.0 < ref["lower_95"] < ref["lower_68"] < ref["median"]
    assert ref["median"] < ref["upper_68"] < ref["upper_95"] < 90.0
    assert ref["normalization_integral"] > 0
    assert np.isfinite(ref["quadrature_error_estimate"])


def test_grid_audit_compares_all_cases_and_grid_sizes():
    rows = module.run_independent_reference_audit()
    assert len(rows) == len(module.CASES) * 2 * len(module.GRID_SIZES)
    assert all(np.isfinite(row["abs_error_median"]) for row in rows)
    # The finest grid should be at least competitive with the coarsest in aggregate.
    med = lambda size: np.median([r["abs_error_median"] for r in rows if r["grid_size"] == size])
    assert med(16001) <= med(1001)
