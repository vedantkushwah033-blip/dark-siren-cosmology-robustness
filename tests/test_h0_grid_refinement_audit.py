from __future__ import annotations
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

MODULE_PATH=Path(__file__).resolve().parents[1]/"scripts"/"audit_h0_grid_refinement.py"
SPEC=importlib.util.spec_from_file_location("h0_grid_refinement",MODULE_PATH)
MODULE=importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def make_grid(root: Path, points: int, shift: float=0.0):
    root.mkdir(parents=True)
    h0=np.linspace(20,300,points)
    density=np.exp(-0.5*((h0-(100+shift))/20)**2)
    density/=np.trapezoid(density,h0)
    pd.DataFrame({"H0":h0,"density_order_5":density}).to_csv(
        root/"h0_posterior_density_by_order.csv",index=False)
    (root/"config.json").write_text(json.dumps({
        "inference":{"h0_grid_points":points,"uniform_h0_prior":[20.0,300.0]},
        "generator":{"seed":123,"gamma":1.0},
    }))
    (root/"quadrature_convergence.json").write_text(json.dumps({
        "status":"CONVERGED","orders_evaluated":[3,5],
        "posterior_by_order":[{"order":5,"h0_median":100+shift,
          "h0_q16":80+shift,"h0_q84":120+shift,"h0_q025":60+shift,
          "h0_q975":140+shift}]}))
    (root/"observable").mkdir()
    pd.DataFrame({"object_id":[1,2],"redshift":[.1,.2]}).to_csv(
        root/"observable/analysis_catalogue.csv",index=False)
    pd.DataFrame({"event_id":[1],"distance_mpc":[100.0],"detected":[1]}).to_csv(
        root/"observable/gw_events_observed.csv",index=False)
    pd.DataFrame({"validation_id":[1,2,3],"redshift":[.1,.2,.3],
                  "retained":[1,0,1]}).to_csv(
        root/"observable/validation_observed.csv",index=False)


def test_refinement_rejects_material_observable_difference(tmp_path):
    for n in (141,281,561):
        make_grid(tmp_path/f"grid_{n}",n)
    path=tmp_path/"grid_141"/"observable/validation_observed.csv"
    frame=pd.read_csv(path)
    frame.loc[0,"redshift"] += 1e-6
    frame.to_csv(path,index=False)
    records=[MODULE.summarize(tmp_path/f"grid_{n}",n) for n in (141,281,561)]
    with pytest.raises(AssertionError,match="beyond atol"):
        MODULE.compare(records,tmp_path/"summary.json")


def test_refinement_tolerates_machine_roundoff_and_row_order(tmp_path):
    for n in (141,281,561):
        make_grid(tmp_path/f"grid_{n}",n)
    path=tmp_path/"grid_281"/"observable/validation_observed.csv"
    frame=pd.read_csv(path)
    frame.loc[0,"redshift"] += 1e-16
    frame=frame.iloc[::-1]
    frame.to_csv(path,index=False)
    records=[MODULE.summarize(tmp_path/f"grid_{n}",n) for n in (141,281,561)]
    report=MODULE.compare(records,tmp_path/"summary.json")
    assert report["same_observable_inputs"] is True
    assert report["observable_comparison_atol"]==1e-12
    assert report["observable_comparison_diagnostics"]["281"]["observable/validation_observed.csv"]["max_abs_float_difference"]<1e-12


def test_refinement_compares_coarse_and_fine_grids(tmp_path):
    for n in (141,281,561):
        make_grid(tmp_path/f"grid_{n}",n)
    records=[MODULE.summarize(tmp_path/f"grid_{n}",n) for n in (141,281,561)]
    report=MODULE.compare(records,tmp_path/"summary.json")
    assert report["status"]=="PASS"
    assert report["same_observable_inputs"] is True
    assert report["reference_grid_points"]==561
    assert len(report["grid_summaries"])==3
    assert report["grid_summaries"][-1]["total_variation_vs_finest"]<1e-12


def test_refinement_rejects_model_configuration_drift(tmp_path):
    for n in (141,281,561):
        make_grid(tmp_path/f"grid_{n}",n)
    config_path=tmp_path/"grid_281"/"config.json"
    config=json.loads(config_path.read_text())
    config["generator"]["gamma"]=2.0
    config_path.write_text(json.dumps(config))
    records=[MODULE.summarize(tmp_path/f"grid_{n}",n) for n in (141,281,561)]
    with pytest.raises(AssertionError,match="configuration differs"):
        MODULE.compare(records,tmp_path/"summary.json")
