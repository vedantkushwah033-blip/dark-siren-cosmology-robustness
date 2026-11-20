"""Focused tests for quadratic five-arm replicate aggregation."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_quadratic_five_arm_replicates import (
    ARMS, _aggregate, _h0_truth, _paired_contrasts, _seed, run as run_replicates,
)
from audit_quadratic_five_arm_replicates import audit as audit_replicates
from audit_quadratic_five_arm_smoke import _independent_inclusion_quadrature
from dark_siren_robustness.generative_selection import generate_quadratic_logit_inclusion_samples


def test_child_seed_schedule_is_deterministic_and_distinct():
    first=[(_seed(2026111929,r,0),_seed(2026111929,r,1)) for r in range(20)]
    second=[(_seed(2026111929,r,0),_seed(2026111929,r,1)) for r in range(20)]
    assert first==second
    flattened=[x for pair in first for x in pair]
    assert len(set(flattened))==40


def test_quadratic_five_arm_aggregate_and_paired_contrasts():
    summaries=[]
    ledger=[]
    errors={
        "A_catalogue_only":-4.0,
        "B_parent_oracle":1.0,
        "C_known_probability_ipw":-1.0,
        "D_map_plugin_ipw":2.0,
        "E_uncertainty_integrated_ipw":1.5,
    }
    for arm in ARMS:
        e=errors[arm]
        summaries.append({
            "replicate":0,"gamma":0.0,"validation_size":100,"arm":arm,
            "status":"SUCCESS","h0_true":70.0,"median_error":e,
            "posterior_median":70.0+e,"q16":65.0,"q84":75.0,
        })
        ledger.append({"replicate":0,"gamma":0.0,"validation_size":100,
            "arm":arm,"status":"SUCCESS"})
        ledger.append({"replicate":1,"gamma":0.0,"validation_size":100,
            "arm":arm,"status":"FAILED"})
    summary=pd.DataFrame(summaries)
    lifecycle=pd.DataFrame(ledger)
    contrasts=_paired_contrasts(summary)
    assert len(contrasts)==8
    assert set(contrasts.reference_arm)=={"A_catalogue_only","C_known_probability_ipw"}
    aggregate=_aggregate(summary,lifecycle,replicates=2,max_abs_error=20.0)
    row=aggregate.loc[aggregate.arm=="A_catalogue_only"].iloc[0]
    assert row.successful_replicates==1
    assert row.failed_or_non_evaluable==1
    assert row.success_only_bias==-4.0
    assert row.success_only_rmse==4.0
    assert row.failure_inclusive_rmse_lower==np.sqrt(8.0)
    assert row.failure_inclusive_rmse_upper==np.sqrt(208.0)
    assert row.failure_inclusive_coverage_lower==0.5
    assert row.failure_inclusive_coverage_upper==1.0
    assert row.failure_inclusive_width_lower==5.0
    assert row.failure_inclusive_width_upper==25.0


def test_all_failed_cells_remain_in_aggregate_denominator():
    ledger=pd.DataFrame([
        {"replicate":r,"gamma":g,"validation_size":n,"arm":a,"status":"FAILED"}
        for r in (0,1) for g in (-1.0,0.0) for n in (100,300) for a in ARMS
    ])
    summary=pd.DataFrame(columns=[
        "replicate","gamma","validation_size","arm","status","h0_true",
        "median_error","q16","q84",
    ])
    aggregate=_aggregate(summary,ledger,replicates=2,max_abs_error=20.0)
    assert len(aggregate)==2*2*len(ARMS)
    assert (aggregate.successful_replicates==0).all()
    assert (aggregate.failed_or_non_evaluable==2).all()
    assert aggregate.success_only_bias.isna().all()



def test_independent_quadrature_solver_replays_previously_non_evaluable_seed():
    generated=generate_quadratic_logit_inclusion_samples(
        population_seed=1067343507,
        validation_seed=1671174138,
        n_parent=1200,
        n_validation_frame=500,
        n_validation=500,
        gamma=-1.0,
        beta=-3.9,
        z_min=0.01,
        z_max=0.40,
    )
    validation=generated["validation_observed"]
    for n in (100,300,500):
        sample=validation.iloc[:n]
        mode,covariance,nodes,weights=_independent_inclusion_quadrature(
            sample.redshift.to_numpy(float),
            sample.retained.to_numpy(int),
            order=5,
            center=0.205,
            scale=0.195,
        )
        assert np.all(np.isfinite(mode))
        assert np.all(np.isfinite(covariance))
        assert nodes.shape==(25,2)
        assert weights.shape==(25,)
        assert np.isclose(weights.sum(),1.0,atol=1e-12)


def test_expanded_prior_failure_bounds_use_full_support_width():
    summaries = pd.DataFrame(columns=[
        "replicate", "gamma", "validation_size", "arm", "status",
        "h0_true", "median_error", "q16", "q84",
    ])
    ledger = pd.DataFrame([
        {"replicate": r, "gamma": 0.0, "validation_size": 100,
         "arm": arm, "status": "FAILED"}
        for r in (0, 1) for arm in ARMS
    ])
    aggregate = _aggregate(
        summaries, ledger, replicates=2, max_abs_error=40.0,
        max_interval_width=80.0,
    )
    assert (aggregate.failure_inclusive_rmse_upper == 40.0).all()
    assert (aggregate.failure_inclusive_width_upper == 80.0).all()


def test_prior_support_sensitivity_quantiles_use_same_likelihood():
    from audit_five_arm_prior_support_sensitivity import posterior_on_support

    grid = np.linspace(30.0, 110.0, 1121)
    log_likelihood = -0.5 * ((grid - 44.0) / 3.0) ** 2
    full = posterior_on_support(grid, log_likelihood, 30.0, 110.0)
    restricted = posterior_on_support(grid, log_likelihood, 50.0, 90.0)
    assert full["median"] < 50.0
    assert restricted["median"] > 50.0
    assert restricted["q025"] >= 50.0
    assert np.isclose(full["support_mass"], 1.0, atol=1e-12)
    assert np.isclose(restricted["support_mass"], 1.0, atol=1e-12)


def test_degree_two_replicate_orchestrator_roundtrip(tmp_path: Path):
    root = tmp_path / "degree_two_replicates"
    args = type("Args", (), {
        "master_seed": 20261120,
        "replicates": 2,
        "gammas": [0.0],
        "n_parent": 300,
        "n_detected": 4,
        "validation_frame": 100,
        "validation_sizes": [60, 90],
        "curvature_beta": -3.9,
        "inclusion_degree": 2,
        "h0_true": 70.0,
        "h0_prior_min": 30.0,
        "h0_prior_max": 110.0,
        "h0_points": 41,
        "scott_multiplier": 0.7,
        "quadrature_order": 3,
        "sigma_mpc": 140.0,
        "detection_threshold_mpc": 250.0,
        "output_dir": root,
    })()
    result = run_replicates(args)
    assert result["status"] == "PASS"
    assert result["planned_cells"] == 20
    config = __import__("json").loads((root / "config.json").read_text())
    assert config["inclusion_model_degree"] == 2
    report = audit_replicates(root)
    assert report["status"] == "PASS"
    assert report["successful_child_runs"] == 2
    assert report["planned_cells"] == 20



def test_prior_predictive_h0_replicate_orchestrator_roundtrip(tmp_path: Path):
    root = tmp_path / "prior_predictive_replicates"
    args = type("Args", (), {
        "master_seed": 20261121,
        "replicates": 2,
        "gammas": [-1.0, 0.0, 1.0],
        "n_parent": 300,
        "n_detected": 4,
        "validation_frame": 100,
        "validation_sizes": [60, 90],
        "curvature_beta": -3.9,
        "inclusion_degree": 1,
        "h0_true": 70.0,
        "h0_prior_min": 30.0,
        "h0_prior_max": 110.0,
        "prior_predictive_h0": True,
        "h0_points": 41,
        "scott_multiplier": 0.7,
        "quadrature_order": 3,
        "sigma_mpc": 140.0,
        "detection_threshold_mpc": 250.0,
        "output_dir": root,
    })()
    result = run_replicates(args)
    assert result["status"] == "PASS"
    assert result["planned_cells"] == 2 * 3 * 2 * len(ARMS)
    config = __import__("json").loads((root / "config.json").read_text())
    schedule = config["h0_truth_schedule"]
    assert config["h0_truth_mode"] == "prior_predictive_uniform"
    assert len(schedule) == 2
    assert len({row["h0_true"] for row in schedule}) == 2
    assert all(30.0 <= row["h0_true"] < 110.0 for row in schedule)
    seeds = config["seed_schedule"]
    assert len(seeds) == 2
    assert len({row["population_seed"] for row in seeds}) == 2
    assert len({row["event_seed"] for row in seeds}) == 2
    assert len({row["validation_frame_seed"] for row in seeds}) == 2
    assert len({row["validation_sample_seed"] for row in seeds}) == 2
    assert len({row["h0_truth_seed"] for row in seeds}) == 2
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    for rep in (0, 1):
        values = ledger.loc[ledger.replicate == rep, "h0_true"].unique()
        assert len(values) == 1
        assert values[0] == schedule[rep]["h0_true"]
    report = audit_replicates(root)
    assert report["status"] == "PASS"
    assert report["successful_child_runs"] == 6
    assert report["planned_cells"] == 2 * 3 * 2 * len(ARMS)


def test_h0_truth_seed_helper_is_deterministic_and_validates_support():
    first = [_h0_truth(20261121, rep, 30.0, 110.0) for rep in range(20)]
    second = [_h0_truth(20261121, rep, 30.0, 110.0) for rep in range(20)]
    assert first == second
    assert len({seed for seed, _ in first}) == 20
    assert len({truth for _, truth in first}) == 20
    assert all(30.0 <= truth < 110.0 for _, truth in first)
    import pytest
    with pytest.raises(ValueError):
        _h0_truth(1, 0, 30.0, 30.0)
