import importlib.util
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("calibration_uncertainty",ROOT/"scripts"/"analyze_posterior_calibration_uncertainty.py")
MODULE=importlib.util.module_from_spec(spec); spec.loader.exec_module(MODULE)

def test_stratified_bootstrap_is_deterministic_and_preserves_cell_structure():
    values=np.vstack([np.full(20,i/12) for i in range(12)])
    a=MODULE.stratified_bootstrap(values,seed=77,resamples=300,batch_size=30)
    b=MODULE.stratified_bootstrap(values,seed=77,resamples=300,batch_size=30)
    assert a==b
    assert a["mean_ci95_low"]==pytest.approx(values.mean())
    assert a["mean_ci95_high"]==pytest.approx(values.mean())

def test_stratified_bootstrap_rejects_wrong_design_shape():
    with pytest.raises(ValueError,match="12-cell array"):
        MODULE.stratified_bootstrap(np.ones((3,10)),seed=1,resamples=100)
