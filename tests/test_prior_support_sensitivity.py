import importlib.util
from pathlib import Path
PATH=Path(__file__).resolve().parents[1]/"scripts"/"run_prior_support_sensitivity.py"
SPEC=importlib.util.spec_from_file_location("prior_support_sensitivity",PATH)
module=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(module)

def test_prior_support_sensitivity_has_complete_paired_design():
    rows,summaries=module.run_audit(replicates=100,seed=17,grid_size=1001)
    assert len(rows)==len(module.TRUTHS)*len(module.PRIORS)*100
    assert len(summaries)==len(module.TRUTHS)*len(module.PRIORS)*len(module.METHODS)
    contrasts=module.paired_prior_differences(rows)
    assert len(contrasts)==len(module.TRUTHS)*(len(module.PRIORS)-1)*100*len(module.METHODS)
    assert all(r["absolute_median_difference"]>=0 for r in contrasts)
    assert all(0<=r["coverage_68"]<=1 and 0<=r["coverage_95"]<=1 for r in summaries)
    assert all(r["rmse"]>=0 and r["coverage_68_wilson_low"]<=r["coverage_68"]<=r["coverage_68_wilson_high"] for r in summaries)

def test_wilson_boundary_cases():
    assert module.wilson(0,100)[0]==0.0
    assert module.wilson(100,100)[1]==1.0
