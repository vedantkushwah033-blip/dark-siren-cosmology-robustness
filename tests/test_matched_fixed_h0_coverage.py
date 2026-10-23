import importlib.util
from pathlib import Path
import numpy as np
import pytest

PATH=Path(__file__).resolve().parents[1]/"scripts"/"run_matched_fixed_h0_coverage.py"
SPEC=importlib.util.spec_from_file_location("matched_cov",PATH)
m=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(m)

def test_cdf_normalizes_and_interval_contains_median():
    grid=np.linspace(50,90,1001)
    density=np.exp(-0.5*((grid-70)/5)**2)
    cdf=m.posterior_cdf(grid,density)
    lo,hi=m.central_interval(grid,cdf,.68)
    med=float(np.interp(.5,cdf,grid))
    assert cdf[0]==pytest.approx(0)
    assert cdf[-1]==pytest.approx(1)
    assert lo < med < hi

def test_summary_reports_coverage_and_mcse():
    rows=[]
    for i in range(100):
        rows.append({"median_error":1.0,"covers_50":int(i<50),"covers_68":int(i<68),
          "covers_90":int(i<90),"covers_95":int(i<95),
          "width_50":10.,"width_68":15.,"width_90":20.,"width_95":25.})
    result=m.summarize(rows,"known_host",70.)
    assert result["coverage_68"]==pytest.approx(.68)
    assert result["coverage_68_mc_se"]>0
    assert result["rmse_median"]==pytest.approx(1.)
