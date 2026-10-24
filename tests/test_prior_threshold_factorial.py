import importlib.util
from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import trapezoid

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


def test_detection_sampler_observations_obey_declared_selection_interval():
    rng=np.random.default_rng(20261024)
    observations=[module.draw_detected(70.0,rng) for _ in range(500)]
    hosts=np.array([host for host,_ in observations])
    distances=np.array([distance for _,distance in observations])
    assert np.all((distances > module.TRUE_MIN) & (distances < module.TRUE_MAX))
    assert set(np.unique(hosts)).issubset({0,1})
    assert np.all(np.isfinite(distances))


def test_selection_efficiency_increases_with_upper_detection_threshold():
    from dark_siren_robustness.selection import detection_efficiency_h0
    grid=np.linspace(30.0,120.0,901)
    beta_120=detection_efficiency_h0(
        grid,module.HOST_REDSHIFTS,module.SIGMA,120.0,
        weights=module.HOST_WEIGHTS,detection_distance_min=module.TRUE_MIN)
    beta_150=detection_efficiency_h0(
        grid,module.HOST_REDSHIFTS,module.SIGMA,150.0,
        weights=module.HOST_WEIGHTS,detection_distance_min=module.TRUE_MIN)
    assert np.all(np.isfinite(beta_120)) and np.all(np.isfinite(beta_150))
    assert np.all((beta_120 > 0) & (beta_120 <= 1))
    assert np.all((beta_150 > 0) & (beta_150 <= 1))
    assert np.all(beta_150 >= beta_120)
    assert np.any(beta_150 > beta_120)


def test_selected_posterior_is_normalized_inside_prior_support_for_multiple_thresholds():
    grid=np.linspace(25.0,125.0,4001)
    support=(grid >= 50.0) & (grid <= 90.0)
    for threshold in (120.0,150.0,180.0,210.0):
        posterior=module.posterior_h0_selected(
            grid,module.HOST_REDSHIFTS,100.0,module.SIGMA,threshold,
            weights=module.HOST_WEIGHTS,prior_min=50.0,prior_max=90.0,
            detection_distance_min=module.TRUE_MIN)
        assert np.all(np.isfinite(posterior))
        assert np.all(posterior >= 0)
        assert trapezoid(posterior[support],grid[support]) == pytest.approx(1.0,abs=1e-10)



def test_experiment_manifest_records_design_and_artifact_counts():
    manifest=module.build_manifest(
        seed=12345,replicates=100,grid_size=2001,
        rows=[{"row":1},{"row":2}],summaries=[{"summary":1}],
        contrasts=[{"contrast":1},{"contrast":2},{"contrast":3}],
    )
    assert manifest["experiment"]=="prior_support_x_assumed_detection_threshold_factorial"
    assert manifest["scope"].startswith("synthetic simulation")
    assert manifest["seed"]==12345
    assert manifest["replicates_per_true_H0"]==100
    assert manifest["posterior_grid_points"]==2001
    assert manifest["raw_rows"]==2
    assert manifest["summary_rows"]==1
    assert manifest["paired_contrast_rows"]==3
    assert manifest["true_detection_interval_mpc"]==[module.TRUE_MIN,module.TRUE_MAX]
    assert manifest["assumed_detection_upper_thresholds_mpc"]==list(module.ASSUMED_MAXIMA)
    assert manifest["numpy_version"]==np.__version__
