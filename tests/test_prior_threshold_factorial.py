import importlib.util
from pathlib import Path
import numpy as np

PATH=Path(__file__).resolve().parents[1]/"scripts"/"run_prior_threshold_factorial.py"
SPEC=importlib.util.spec_from_file_location("prior_threshold_factorial",PATH)
module=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(module)

def test_factorial_rows_pair_same_observations_and_have_complete_cells():
    rows=module.run(replicates=100,seed=44,grid_size=1001)
    expected=len(module.TRUTHS)*100*len(module.PRIORS)*len(module.ASSUMED_MAXIMA)*len(module.METHODS)
    assert len(rows)==expected
    for truth in module.TRUTHS:
        for rep in range(100):
            cell=[r for r in rows if r["true_H0"]==truth and r["replicate"]==rep]
            assert len({r["observed_distance_mpc"] for r in cell})==1
            assert len({r["host_index"] for r in cell})==1
            assert len(cell)==len(module.PRIORS)*len(module.ASSUMED_MAXIMA)*len(module.METHODS)
    assert all(r["true_detection_max_mpc"]==120.0 for r in rows)

def test_uncorrected_posterior_is_identical_across_threshold_assumptions():
    rows=module.run(replicates=100,seed=45,grid_size=1001)
    for truth in module.TRUTHS:
        for rep in range(100):
            cell=[r for r in rows if r["true_H0"]==truth and r["replicate"]==rep and r["prior_min"]==50 and r["prior_max"]==90 and r["method"]=="uncorrected"]
            assert len(cell)==len(module.ASSUMED_MAXIMA)
            assert len({r["posterior_median"] for r in cell})==1

def test_summaries_and_contrasts_are_complete_and_finite():
    rows=module.run(replicates=100,seed=46,grid_size=1001)
    summaries=module.summarize(rows)
    assert len(summaries)==len(module.TRUTHS)*len(module.PRIORS)*len(module.ASSUMED_MAXIMA)*len(module.METHODS)
    assert all(np.isfinite(r["rmse"]) and 0<=r["coverage_68"]<=1 for r in summaries)
    contrasts=module.paired_contrasts(rows,bootstrap_draws=1000,seed=47)
    assert len(contrasts)==len(module.TRUTHS)*((len(module.PRIORS)*len(module.ASSUMED_MAXIMA))-1)*len(module.METHODS)*6
    assert all(r["bootstrap_95_low"]<=r["bootstrap_95_high"] for r in contrasts)
