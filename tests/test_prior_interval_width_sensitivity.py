import importlib.util
from pathlib import Path
PATH=Path(__file__).resolve().parents[1]/"scripts"/"run_prior_interval_width_sensitivity.py"
SPEC=importlib.util.spec_from_file_location("interval_width_audit",PATH)
module=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(module)

def test_width_summary_and_paired_differences_are_complete():
    rows,_=module.source.run_audit(replicates=100,seed=29,grid_size=1001)
    summary,paired=module.summarize_widths(rows)
    assert len(summary)==len(module.source.TRUTHS)*len(module.source.PRIORS)*len(module.source.METHODS)
    assert len(paired)==len(module.source.TRUTHS)*(len(module.source.PRIORS)-1)*100*len(module.source.METHODS)
    assert all(r["n"]==100 for r in summary)
    assert all(r["mean_width_68"]>0 and r["mean_width_95"]>=r["mean_width_68"] for r in summary)
    assert all(r["absolute_width_difference_68"]>=0 and r["absolute_width_difference_95"]>=0 for r in paired)
