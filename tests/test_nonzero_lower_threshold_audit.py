import importlib.util
from pathlib import Path
import numpy as np
import pytest

PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_nonzero_lower_threshold_audit.py"
SPEC = importlib.util.spec_from_file_location("nonzero_lower_threshold_audit", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)

def test_nonzero_lower_threshold_audit_is_finite():
    rows = module.run_audit(1001)
    assert len(rows) == 6
    assert {r["observed_distance_mpc"] for r in rows} == {30.0, 70.0, 140.0}
    assert {r["method"] for r in rows} == {"uncorrected", "selection_normalized"}
    assert all(float(r["detection_min_mpc"]) == 20.0 for r in rows)
    assert all(float(r["detection_max_mpc"]) == 180.0 for r in rows)
    numeric_keys = ("observed_distance_mpc", "grid_size", "detection_min_mpc", "detection_max_mpc", "beta_min_on_grid", "beta_max_on_grid", "mean", "median", "posterior_sd", "lower_68", "upper_68", "lower_95", "upper_95")
    assert all(np.isfinite(float(r[key])) for r in rows for key in numeric_keys)

def test_observation_below_lower_detection_threshold_rejected():
    from dark_siren_robustness.selection import posterior_h0_selected
    with pytest.raises(ValueError, match="detection interval"):
        posterior_h0_selected(np.linspace(25, 125, 101), [0.02, 0.04], 10.0, 20.0, 180.0,
                              weights=[0.7, 0.3], detection_distance_min=20.0)
