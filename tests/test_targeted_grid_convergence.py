from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_targeted_grid_convergence.py"
SPEC = spec_from_file_location("targeted_grid_convergence", MODULE_PATH)
module = module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_targeted_convergence_is_finite_and_complete():
    rows = module.calculate_convergence()
    assert [row["grid_size"] for row in rows] == list(module.GRID_SIZES)
    fields = ("posterior_median", "lower_68", "upper_68", "lower_95", "upper_95")
    assert all(np.isfinite(float(row[key])) for row in rows for key in fields)
    assert all(np.isfinite(float(row["posterior_area_on_prior_support"])) for row in rows)
    assert all(abs(float(row["posterior_area_on_prior_support"]) - 1.0) <= 1e-8 for row in rows)
    assert all(np.isfinite(float(row["posterior_area_simpson_on_prior_support"])) for row in rows)
    assert all(float(row["posterior_area_simpson_on_prior_support"]) > 0.0 for row in rows)
    assert all(np.isfinite(float(row["abs_area_rule_difference"])) for row in rows)
    # Simpson's rule is a cross-method diagnostic, not required to match the
    # trapezoid-normalized area to machine precision on every sampled grid.
    assert all(row["lower_95"] <= row["lower_68"] <= row["posterior_median"]
               <= row["upper_68"] <= row["upper_95"] for row in rows)


def test_targeted_upper95_successive_changes_shrink():
    rows = module.calculate_convergence()
    changes = [float(row["abs_upper95_change_from_previous_grid"])
               for row in rows[1:]]
    assert all(np.isfinite(changes))
    # A fixed-case regression check: finer grids should reduce successive
    # endpoint changes for this known tail-sensitive synthetic observation.
    assert all(later < earlier for earlier, later in zip(changes, changes[1:]))
    assert changes[-1] < 0.001
