import importlib.util
from pathlib import Path
import numpy as np

PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_lower_threshold_sensitivity.py"
SPEC = importlib.util.spec_from_file_location("lower_threshold_sensitivity", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)

def test_lower_threshold_sensitivity_returns_finite_paired_rows():
    rows = module.run_audit(1001)
    assert len(rows) == 24
    assert {r["detection_min_mpc"] for r in rows} == {0.0, 10.0, 20.0, 40.0}
    assert {r["observed_distance_mpc"] for r in rows} == {50.0, 90.0, 140.0}
    assert {r["method"] for r in rows} == {"uncorrected", "selection_normalized"}
    numeric_keys = ("observed_distance_mpc", "grid_size", "detection_min_mpc",
                    "detection_max_mpc", "beta_min_prior_support", "beta_max_prior_support",
                    "mean", "median", "posterior_sd", "lower_68", "upper_68", "lower_95", "upper_95")
    assert all(np.isfinite(float(row[key])) for row in rows for key in numeric_keys)
    assert all(0 < float(row["beta_min_prior_support"]) <= float(row["beta_max_prior_support"]) <= 1
               for row in rows)
