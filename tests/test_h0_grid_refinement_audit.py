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
    density/=np.trapz(density,h0)
    pd.DataFrame({"H0":h0,"density_order_5":density}).to_csv(
        root/"h0_posterior_density_by_order.csv",index=False)
    (root/"quadrature_convergence.json").write_text(json.dumps({
        "status":"CONVERGED","orders_evaluated":[3,5],
        "posterior_by_order":[{"order":5,"h0_median":100+shift,
          "h0_q16":80+shift,"h0_q84":120+shift,"h0_q025":60+shift,
          "h0_q975":140+shift}]}))
    (root/"observable").mkdir()
    for name in ("analysis_catalogue.csv","gw_events_observed.csv","validation_observed.csv"):
        (root/"observable"/name).write_text("same inputs\n")


def test_refinement_requires_identical_observables(tmp_path):
    for n in (141,281,561):
        make_grid(tmp_path/f"grid_{n}",n,shift=(1 if n==141 else 0))
    records=[MODULE.summarize(tmp_path/f"grid_{n}",n) for n in (141,281,561)]
    with pytest.raises(AssertionError,match="observable inputs differ"):
        MODULE.compare(records,tmp_path/"summary.json")


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
