import importlib.util
from pathlib import Path
PATH=Path(__file__).resolve().parents[1]/"scripts"/"run_joint_prior_scorecard.py"
SPEC=importlib.util.spec_from_file_location("joint_scorecard",PATH)
module=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(module)

def test_joint_scorecard_has_complete_keys_and_metrics():
    raw, summaries=module.prior.run_audit(replicates=100,seed=37,grid_size=1001)
    width_summary,_=module.width.summarize_widths(raw)
    score=module.build_scorecard(summaries,width_summary)
    assert len(score)==18
    assert len({(r["true_H0"],r["prior_min"],r["prior_max"],r["method"]) for r in score})==18
    assert all(0<=r["abs_coverage_68_error_from_nominal"]<=1 for r in score)
    assert all(r["rmse"]>=0 and r["mean_width_68"]>0 and r["mean_width_95"]>=r["mean_width_68"] for r in score)
