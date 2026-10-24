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


def test_uncorrected_interval_endpoints_are_invariant_to_threshold():
    rows=module.run(replicates=100,seed=48,grid_size=1001)
    fields=("posterior_median","lower_68","upper_68","lower_95","upper_95")
    for truth in module.TRUTHS:
        for rep in range(100):
            for pmin,pmax in module.PRIORS:
                cell=[r for r in rows if r["true_H0"]==truth and r["replicate"]==rep
                      and r["prior_min"]==pmin and r["prior_max"]==pmax
                      and r["method"]=="uncorrected"]
                assert len(cell)==len(module.ASSUMED_MAXIMA)
                for field in fields:
                    assert np.ptp([r[field] for r in cell])==0.0


def test_paired_bootstrap_is_deterministic_and_preserves_pairing():
    # Tiny complete fixture: two matched replicates per truth and every
    # prior/threshold/method cell. This isolates the contrast/bootstrap logic
    # without rerunning the computationally heavier posterior simulation.
    rows=[]
    for truth in module.TRUTHS:
        for rep in range(2):
            for pmin,pmax in module.PRIORS:
                for threshold in module.ASSUMED_MAXIMA:
                    for method in module.METHODS:
                        shift=(pmin-50.0)*0.01+(threshold-120.0)*0.001
                        if method=="corrected":
                            shift += 0.15
                        median=truth+(-1.0 if rep==0 else 1.0)+shift
                        width=4.0+0.01*(pmax-pmin)+0.001*(threshold-120.0)
                        rows.append({
                            "true_H0":truth,"replicate":rep,"prior_min":pmin,"prior_max":pmax,
                            "assumed_detection_max_mpc":threshold,"method":method,
                            "posterior_median":median,"lower_68":median-width/2,
                            "upper_68":median+width/2,"lower_95":median-width,
                            "upper_95":median+width,
                        })
    first=module.paired_contrasts(rows,bootstrap_draws=257,seed=90210)
    second=module.paired_contrasts(rows,bootstrap_draws=257,seed=90210)
    assert first==second
    assert len(first)==len(module.TRUTHS)*((len(module.PRIORS)*len(module.ASSUMED_MAXIMA))-1)*len(module.METHODS)*6
    assert all(np.isfinite(r["mean_difference_vs_baseline"]) for r in first)
    assert all(r["bootstrap_95_low"]<=r["bootstrap_95_high"] for r in first)
