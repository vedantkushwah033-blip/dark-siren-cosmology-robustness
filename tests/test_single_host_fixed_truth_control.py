import importlib.util
from pathlib import Path
import numpy as np
import pytest
PATH=Path(__file__).resolve().parents[1]/"scripts"/"run_single_host_fixed_truth_control.py"
SPEC=importlib.util.spec_from_file_location("single_host_fixed_truth_control",PATH)
module=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(module)
def test_sampler_respects_detection_interval_and_is_reproducible():
    a=np.random.default_rng(21);b=np.random.default_rng(21)
    x=[module.draw(70.,a) for _ in range(100)]
    y=[module.draw(70.,b) for _ in range(100)]
    assert x==y and all(0<d<module.DMAX for d in x)
def test_small_run_pairing_metrics_and_manifest():
    rows,summaries,manifest=module.run(replicates=100,seed=92,grid_size=1001,bootstrap_draws=1000)
    assert len(rows)==len(module.TRUTHS)*100
    assert len(summaries)==len(module.TRUTHS)*3
    assert manifest["raw_rows"]==len(rows) and manifest["summary_rows"]==len(summaries)
    for truth in module.TRUTHS:
        cell=[r for r in rows if r["true_H0"]==truth]
        assert len(cell)==100
        assert all(0<r["observed_distance_mpc"]<module.DMAX for r in cell)
        for r in cell:
            assert r["corrected_lower_68"]<=r["corrected_upper_68"]
            assert r["uncorrected_lower_95"]<=r["uncorrected_upper_95"]
    for r in summaries:
        if r["method"]=="paired_corrected_minus_uncorrected_absolute_error":
            assert r["bootstrap_95_low"]<=r["bootstrap_95_high"]
        elif "coverage_68" in r:
            assert 0<=r["coverage_68"]<=1 and r["coverage_68_wilson_low"]<=r["coverage_68"]<=r["coverage_68_wilson_high"]
def test_rejects_invalid_design_sizes():
    with pytest.raises(ValueError,match="at least 100"):
        module.run(replicates=50)
