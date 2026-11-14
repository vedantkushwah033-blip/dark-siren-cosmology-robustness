from __future__ import annotations
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import numpy as np
import pytest


ROOT=Path(__file__).resolve().parents[1]

def _load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/"scripts"/path)
    module=importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module

def test_absolute_bandwidth_density_smoke_independent_reconstruction_and_corruption(tmp_path):
    runner=_load("absolute_density_runner","run_absolute_bandwidth_density_smoke.py")
    auditor=_load("absolute_density_auditor","audit_absolute_bandwidth_density_smoke.py")
    output=tmp_path/"absolute_density"
    args=SimpleNamespace(
        output_dir=output, seed=20261114, validation_seed=20261115,
        n_parent=500, validation_frame=120, validation_sizes=(40,80),
        catalogue_floor=0.10, scott_multiplier=0.7, quadrature_order=3, grid_points=101,
        absolute_widths=(0.012,0.024),
    )
    result=runner.run(args)
    assert result["status"]=="PASS"
    report=auditor.audit(output)
    assert report["status"]=="PASS"
    assert report["metric_rows_reconstructed"]==len(pd.read_csv(output/"results"/"density_metrics.csv"))
    assert report["scott_width_rows_reconstructed"]==len(pd.read_csv(output/"results"/"scott_widths.csv"))

    metrics_path=output/"results"/"density_metrics.csv"
    metrics=pd.read_csv(metrics_path)
    metrics.loc[0,"iae"] += 0.01
    metrics.to_csv(metrics_path,index=False)
    with pytest.raises(AssertionError,match="manifest verification failed"):
        auditor.audit(output)


def test_quadratic_logit_density_smoke_independently_audits_true_probability_arm(tmp_path):
    import json

    runner=_load("quadratic_density_runner","run_absolute_bandwidth_density_smoke.py")
    auditor=_load("quadratic_density_auditor","audit_absolute_bandwidth_density_smoke.py")
    output=tmp_path/"quadratic_density"
    args=SimpleNamespace(
        output_dir=output, seed=20261124, validation_seed=20261125,
        n_parent=700, validation_frame=300, validation_sizes=(80,160),
        gamma=1.0, beta=-3.9, scott_multiplier=0.7,
        quadrature_order=3, grid_points=101, absolute_widths=(0.012,0.024),
    )
    result=runner.run(args)
    assert result["status"]=="PASS"
    cfg=json.loads((output/"config.json").read_text())
    assert cfg["catalogue_selection_model"]=="quadratic_logit"
    assert cfg["catalogue_floor"] is None
    assert cfg["quadratic_logit_calibration"]["gamma"]==pytest.approx(1.0)

    report=auditor.audit(output)
    assert report["status"]=="PASS"
    assert report["metric_rows_reconstructed"]==len(
        pd.read_csv(output/"results"/"density_metrics.csv")
    )

    # A refreshed checksum must not allow a false parent-measure centering constant.
    import hashlib
    config_path=output/"config.json"
    cfg=json.loads(config_path.read_text())
    cfg["quadratic_logit_calibration"]["x2_center"] += 0.01
    config_path.write_text(json.dumps(cfg,indent=2,sort_keys=True)+"\n")
    manifest_path=output/"manifest.json"
    manifest=json.loads(manifest_path.read_text())
    payload=config_path.read_bytes()
    manifest["files"]["config.json"]={"bytes":len(payload),"sha256":hashlib.sha256(payload).hexdigest()}
    manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True))
    with pytest.raises(AssertionError,match="quadratic-logit parent-measure x2 centering mismatch"):
        auditor.audit(output)


def test_absolute_bandwidth_auditor_rejects_corrupt_weight_diagnostics_after_manifest_refresh(tmp_path):
    import hashlib
    import json

    runner=_load("absolute_density_runner_weight_diag","run_absolute_bandwidth_density_smoke.py")
    auditor=_load("absolute_density_auditor_weight_diag","audit_absolute_bandwidth_density_smoke.py")
    output=tmp_path/"absolute_density_weight_diag"
    args=SimpleNamespace(
        output_dir=output, seed=20261120, validation_seed=20261121,
        n_parent=500, validation_frame=120, validation_sizes=(40,80),
        catalogue_floor=0.10, scott_multiplier=0.7, quadrature_order=3,
        grid_points=101, absolute_widths=(0.012,0.024),
    )
    runner.run(args)
    path=output/"results"/"weight_diagnostics.csv"
    frame=pd.read_csv(path)
    assert "entropy_ess" in frame
    frame.loc[0,"entropy_ess"] += 0.25
    frame.to_csv(path,index=False)

    manifest_path=output/"manifest.json"
    manifest=json.loads(manifest_path.read_text())
    rel=str(path.relative_to(output))
    payload=path.read_bytes()
    manifest["files"][rel]={"bytes":len(payload),"sha256":hashlib.sha256(payload).hexdigest()}
    manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True))

    with pytest.raises(AssertionError,match="independently reconstructed weight diagnostic entropy_ess"):
        auditor.audit(output)


def test_absolute_bandwidth_auditor_rejects_non_nested_validation_ids(tmp_path):
    import hashlib
    import json

    runner=_load("absolute_density_runner_ids","run_absolute_bandwidth_density_smoke.py")
    auditor=_load("absolute_density_auditor_ids","audit_absolute_bandwidth_density_smoke.py")
    output=tmp_path/"absolute_density_ids"
    args=SimpleNamespace(
        output_dir=output, seed=20261116, validation_seed=20261117,
        n_parent=500, validation_frame=120, validation_sizes=(40,80),
        catalogue_floor=0.10, scott_multiplier=0.7, quadrature_order=3, grid_points=101,
        absolute_widths=(0.012,0.024),
    )
    runner.run(args)
    path=output/"observable"/"validation_n80.csv"
    frame=pd.read_csv(path)
    frame.loc[0,"validation_id"] = frame.loc[1,"validation_id"]
    frame.to_csv(path,index=False)

    manifest_path=output/"manifest.json"
    manifest=json.loads(manifest_path.read_text())
    rel=str(path.relative_to(output))
    payload=path.read_bytes()
    manifest["files"][rel]={"bytes":len(payload),"sha256":hashlib.sha256(payload).hexdigest()}
    manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True))

    with pytest.raises(AssertionError,match="validation IDs are not the declared nested prefix"):
        auditor.audit(output)

def test_absolute_bandwidth_auditor_rejects_modified_scott_width_even_with_updated_manifest(tmp_path):
    import hashlib
    import json

    runner=_load("absolute_density_runner_scott","run_absolute_bandwidth_density_smoke.py")
    auditor=_load("absolute_density_auditor_scott","audit_absolute_bandwidth_density_smoke.py")
    output=tmp_path/"absolute_density_scott"
    args=SimpleNamespace(
        output_dir=output, seed=20261118, validation_seed=20261119,
        n_parent=500, validation_frame=120, validation_sizes=(40,80),
        catalogue_floor=0.10, scott_multiplier=0.7, quadrature_order=3,
        grid_points=101, absolute_widths=(0.012,0.024),
    )
    runner.run(args)
    path=output/"results"/"scott_widths.csv"
    frame=pd.read_csv(path)
    frame.loc[0,"realized_bandwidth_redshift"] += 0.01
    frame.to_csv(path,index=False)
    manifest_path=output/"manifest.json"
    manifest=json.loads(manifest_path.read_text())
    rel=str(path.relative_to(output))
    payload=path.read_bytes()
    manifest["files"][rel]={"bytes":len(payload),"sha256":hashlib.sha256(payload).hexdigest()}
    manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True))
    with pytest.raises(AssertionError,match="independently reconstructed realized Scott widths"):
        auditor.audit(output)

def test_posterior_mean_density_metric_is_distinct_and_independently_audited(tmp_path):
    import hashlib
    import json

    runner=_load("absolute_density_runner_mean","run_absolute_bandwidth_density_smoke.py")
    auditor=_load("absolute_density_auditor_mean","audit_absolute_bandwidth_density_smoke.py")
    output=tmp_path/"absolute_density_mean"
    args=SimpleNamespace(
        output_dir=output, seed=20261231, validation_seed=20261232,
        n_parent=500, validation_frame=120, validation_sizes=(40,80),
        catalogue_floor=0.10, scott_multiplier=0.7, quadrature_order=3,
        grid_points=101, absolute_widths=(0.012,0.024),
    )
    runner.run(args)
    metrics_path=output/"results"/"density_metrics.csv"
    metrics=pd.read_csv(metrics_path)
    mean_rows=metrics[metrics.arm=="E_posterior_mean_density"]
    node_rows=metrics[metrics.arm=="E_posterior_weighted_node_metric"]
    e_node_rows=metrics[metrics.arm=="E_uncertainty_integrated_ipw"]
    original_metrics=metrics.copy(deep=True)
    assert len(mean_rows)==4
    assert len(node_rows)==4
    # Convexity: loss of the posterior-mean density cannot exceed
    # posterior expected node-specific loss for these convex metrics.
    joined=mean_rows.merge(node_rows,on=["replicate","validation_size","absolute_bandwidth_sd"],suffixes=("_mean","_nodes"))
    for metric in ("iae","ise","boundary_iae"):
        assert (joined[f"{metric}_mean"] <= joined[f"{metric}_nodes"] + 1e-10).all()
    assert mean_rows.kish_ess.isna().all()
    assert node_rows.kish_ess.isna().all()
    assert e_node_rows.kish_ess.notna().all()
    assert (e_node_rows.kish_ess > 0).all()
    assert auditor.audit(output)["status"]=="PASS"

    # Refresh the manifest after corrupting the primary posterior-mean metric.
    metrics.loc[metrics.arm=="E_posterior_mean_density","iae"] += 0.01
    metrics.to_csv(metrics_path,index=False)
    manifest_path=output/"manifest.json"
    manifest=json.loads(manifest_path.read_text())
    rel=str(metrics_path.relative_to(output))
    payload=metrics_path.read_bytes()
    manifest["files"][rel]={"bytes":len(payload),"sha256":hashlib.sha256(payload).hexdigest()}
    manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True))
    with pytest.raises(AssertionError,match="metric iae"):
        auditor.audit(output)

    # Restore the original metrics, then corrupt a node-specific ESS while
    # refreshing the manifest; the auditor must reconstruct ESS from raw weights.
    metrics=original_metrics.copy(deep=True)
    metrics.to_csv(metrics_path,index=False)
    payload=metrics_path.read_bytes()
    manifest["files"][rel]={"bytes":len(payload),"sha256":hashlib.sha256(payload).hexdigest()}
    manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True))
    metrics.loc[metrics.arm=="E_uncertainty_integrated_ipw","kish_ess"] += 5.0
    metrics.to_csv(metrics_path,index=False)
    payload=metrics_path.read_bytes()
    manifest["files"][rel]={"bytes":len(payload),"sha256":hashlib.sha256(payload).hexdigest()}
    manifest_path.write_text(json.dumps(manifest,indent=2,sort_keys=True))
    with pytest.raises(AssertionError,match="metric Kish ESS"):
        auditor.audit(output)

def test_floor_025_failed_validation_posterior_uses_strict_newton_refinement():
    from scipy.special import expit
    from dark_siren_robustness.generative_selection import generate_independent_validation_sample
    from dark_siren_robustness.inclusion_posterior import fit_inclusion_posterior_quadrature

    sample=generate_independent_validation_sample(
        seed=20261266, n_validation_frame=1500, n_validation=1000,
        z_min=0.01, z_max=0.40, catalogue_floor=0.25,
        catalogue_midpoint=0.22, catalogue_width=0.05,
    )
    validation=sample["validation_observed"]
    z=validation.redshift.to_numpy(float)
    y=validation.retained.to_numpy(int)
    fit=fit_inclusion_posterior_quadrature(
        z,y,degree=1,prior_mean=(0.0,0.0),prior_sd=(5.0,2.5),
        coordinate_center=0.205,coordinate_scale=0.195,order=5,
    )
    design=np.column_stack((np.ones(len(z)),(z-0.205)/0.195))
    score=design.T@(expit(design@fit.mode)-y)+fit.mode/np.square((5.0,2.5))
    assert np.linalg.norm(score,ord=np.inf)<=1e-6
    assert np.all(np.linalg.eigvalsh(fit.laplace_covariance)>0)

def test_quadratic_logit_curvature_dgp_preserves_parent_mean_and_separates_truth(tmp_path):
    from dark_siren_robustness.generative_selection import (
        calibrate_quadratic_logit_intercept,
        generate_quadratic_logit_inclusion_samples,
        quadratic_logit_inclusion_probability,
        quadratic_logit_parent_measure,
    )

    measure=quadratic_logit_parent_measure()
    assert measure["x2_center"] == pytest.approx(0.284500557,abs=2e-8)
    gammas=(-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0)
    calibrations=[calibrate_quadratic_logit_intercept(gamma=g) for g in gammas]
    rates=np.array([x["achieved_rate"] for x in calibrations])
    assert np.ptp(rates) < 2e-12
    assert all(x["achieved_rate"] == pytest.approx(x["target_rate"],abs=2e-12) for x in calibrations)
    assert calibrations[3]["alpha"] == pytest.approx(0.0,abs=1e-12)

    z=np.linspace(0.01,0.40,501)
    p0=quadratic_logit_inclusion_probability(z,gamma=0.0)
    expected=1.0/(1.0+np.exp(3.9*(z-0.205)/0.195))
    assert np.allclose(p0,expected,rtol=0,atol=2e-12)
    p_curved=quadratic_logit_inclusion_probability(z,gamma=2.0)
    assert not np.allclose(p_curved,p0)

    sample=generate_quadratic_logit_inclusion_samples(
        population_seed=2026111401,validation_seed=2026111402,
        n_parent=1000,n_validation_frame=500,n_validation=300,gamma=1.0,
    )
    assert set(sample["validation_observed"].columns)=={"validation_id","redshift","retained"}
    assert "inclusion_probability" not in sample["validation_observed"].columns
    assert len(sample["parent_truth"])==1000
    assert len(sample["catalogue_observed"])==int(sample["parent_truth"].retained.sum())
    assert len(sample["validation_observed"])==300
    assert sample["config"]["catalogue_selection_model"]=="quadratic_logit"
    repeat=generate_quadratic_logit_inclusion_samples(
        population_seed=2026111401,validation_seed=2026111402,
        n_parent=1000,n_validation_frame=500,n_validation=300,gamma=1.0,
    )
    pd.testing.assert_frame_equal(sample["parent_truth"],repeat["parent_truth"])
    pd.testing.assert_frame_equal(sample["validation_observed"],repeat["validation_observed"])



def test_quadratic_logit_parent_measure_rejects_nonfinite_or_invalid_grid_inputs():
    from dark_siren_robustness.generative_selection import quadratic_logit_parent_measure

    invalid_kwargs = (
        {"z_min": float("nan")},
        {"z_max": float("inf")},
        {"coordinate_center": float("nan")},
        {"coordinate_scale": float("inf")},
        {"coordinate_scale": 0.0},
        {"grid_size": 1000},
        {"grid_size": 1001.5},
        {"grid_size": True},
    )
    for kwargs in invalid_kwargs:
        with pytest.raises(ValueError):
            quadratic_logit_parent_measure(**kwargs)
