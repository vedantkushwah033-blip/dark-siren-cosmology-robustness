import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_randomized_numerical_stress.py"
SPEC = importlib.util.spec_from_file_location("randomized_numerical_stress", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_stress_audit_is_reproducible_and_complete():
    first = module.run_stress_audit(cases=30, seed=12)
    second = module.run_stress_audit(cases=30, seed=12)
    assert len(first) == len(second) == 60
    assert first == second
    assert all(np.isfinite(row["abs_error_median"]) for row in first)


def test_stress_audit_refines_grid_in_typical_median_error():
    rows = module.run_stress_audit(cases=30, seed=20261024)
    errors = {
        size: np.median([r["abs_error_median"] for r in rows if r["grid_size"] == size])
        for size in (1001, 4001)
    }
    assert errors[4001] < errors[1001]
