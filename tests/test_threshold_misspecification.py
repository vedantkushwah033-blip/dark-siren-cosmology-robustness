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
    assert {r["inferred_detection_max_mpc"] for r in summaries}==set(module.INFERRED_MAXIMA)
    for row in summaries:
        assert np.isfinite(row["bias"]) and np.isfinite(row["rmse"])
        assert row["rmse"]>=0 and 0<=row["coverage_68"]<=1 and 0<=row["coverage_95"]<=1

def test_correct_threshold_matches_recomputed_posterior_median():
    from dark_siren_robustness.hardened_inference import posterior_summary
    from dark_siren_robustness.selection import posterior_h0_selected
    from dark_siren_robustness.hardened_inference import posterior_h0
    grid=np.linspace(25,125,1001); truth=70.; host,distance=module.draw_detected(truth,np.random.default_rng(1))
    posterior=posterior_h0_selected(grid,module.HOST_REDSHIFTS,distance,module.SIGMA,
        module.TRUE_MAX,weights=module.HOST_WEIGHTS,prior_min=50,prior_max=90,detection_distance_min=0)
    summary=posterior_summary(grid,posterior)
    rows,_=module.run(replicates=100,seed=123,grid_size=1001)
    row=next(r for r in rows if r["true_H0"]==52 and r["replicate"]==0 and r["inferred_detection_max_mpc"]==module.TRUE_MAX)
    assert row["corrected_median"]>=50 and row["corrected_median"]<=90
    assert 50<=summary["median"]<=90
