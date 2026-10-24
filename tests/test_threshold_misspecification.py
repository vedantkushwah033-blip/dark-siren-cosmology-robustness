import importlib.util
from pathlib import Path
import numpy as np

PATH=Path(__file__).resolve().parents[1]/"scripts"/"run_threshold_misspecification.py"
SPEC=importlib.util.spec_from_file_location("threshold_misspecification",PATH)
module=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(module)

def test_threshold_misspecification_design_is_paired_and_complete():
    rows,summaries=module.run(replicates=100,seed=123,grid_size=1001)
    assert len(rows)==len(module.TRUTHS)*100*len(module.INFERRED_MAXIMA)
    assert len(summaries)==len(module.TRUTHS)*len(module.INFERRED_MAXIMA)
    for truth in module.TRUTHS:
        for rep in range(100):
            cell=[r for r in rows if r["true_H0"]==truth and r["replicate"]==rep]
            assert len(cell)==len(module.INFERRED_MAXIMA)
            assert len({r["observed_distance_mpc"] for r in cell})==1
            assert len({r["host_index"] for r in cell})==1
            assert all(r["true_detection_max_mpc"]==module.TRUE_MAX for r in cell)
    assert {r["inferred_detection_max_mpc"] for r in summaries}==set(module.INFERRED_MAXIMA)
    for row in summaries:
        assert np.isfinite(row["bias"]) and np.isfinite(row["rmse"])
        assert row["rmse"]>=0 and 0<=row["coverage_68"]<=1 and 0<=row["coverage_95"]<=1

def test_correct_threshold_row_matches_direct_posterior_computation():
    from dark_siren_robustness.hardened_inference import posterior_summary
    from dark_siren_robustness.selection import posterior_h0_selected
    rows,_=module.run(replicates=100,seed=123,grid_size=1001)
    row=next(r for r in rows if r["true_H0"]==52 and r["replicate"]==0 and r["inferred_detection_max_mpc"]==module.TRUE_MAX)
    grid=np.linspace(25,125,1001)
    posterior=posterior_h0_selected(grid,module.HOST_REDSHIFTS,row["observed_distance_mpc"],module.SIGMA,
        module.TRUE_MAX,weights=module.HOST_WEIGHTS,prior_min=50,prior_max=90,detection_distance_min=module.TRUE_MIN)
    summary=posterior_summary(grid,posterior)
    assert row["corrected_median"]==summary["median"]
    assert row["corrected_lower_68"]==summary["lower_68"]
    assert row["corrected_upper_68"]==summary["upper_68"]


def test_paired_threshold_contrasts_have_expected_cells_and_finite_intervals():
    rows,_=module.run(replicates=100,seed=124,grid_size=1001)
    contrasts=module.paired_threshold_contrasts(rows,seed=125,bootstrap_draws=1000)
    assert len(contrasts)==len(module.TRUTHS)*(len(module.INFERRED_MAXIMA)-1)*7
    assert {r["metric"] for r in contrasts}=={
        "posterior_median","absolute_error","squared_error","width_68","width_95",
        "coverage_68_indicator","coverage_95_indicator"
    }
    for row in contrasts:
        assert row["n_pairs"]==100
        assert row["bootstrap_95_low"]<=row["bootstrap_95_high"]
        assert np.isfinite(row["mean_wrong_minus_correct"])
        assert row["bootstrap_draws"]==1000
