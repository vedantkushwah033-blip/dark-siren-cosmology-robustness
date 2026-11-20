from pathlib import Path
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_prior_predictive_sbc import analyze, _mean_ci

ARMS = [
    "A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
    "D_map_plugin_ipw", "E_uncertainty_integrated_ipw",
]

def test_mean_ci_handles_one_and_multiple_observations():
    assert _mean_ci(np.array([2.0])) == (2.0, np.nan, np.nan, np.nan)
    mean, se, lo, hi = _mean_ci(np.array([1.0, 3.0, 5.0]))
    assert mean == 3.0
    assert se > 0 and lo < mean < hi

def test_detailed_readout_requires_complete_audited_paired_cells(tmp_path):
    root = tmp_path / "input"
    root.mkdir()
    ledger, summaries, diagnostics = [], [], []
    offsets = {
        "A_catalogue_only": 0.0,
        "B_parent_oracle": 0.2,
        "C_known_probability_ipw": -0.1,
        "D_map_plugin_ipw": 0.3,
        "E_uncertainty_integrated_ipw": 0.1,
    }
    for rep, truth in [(0, 65.0), (1, 75.0), (2, 85.0)]:
        for arm in ARMS:
            error = offsets[arm] + (rep - 1) * 0.5
            ledger.append({
                "replicate": rep, "gamma": 0.0, "validation_size": 100,
                "arm": arm, "status": "SUCCESS", "audit_status": "PASS",
            })
            summaries.append({
                "replicate": rep, "gamma": 0.0, "validation_size": 100,
                "arm": arm, "status": "SUCCESS", "h0_true": truth,
                "median_error": error, "posterior_rank": (rep + 0.5) / 3,
                "q16": truth - 2.0, "q84": truth + 2.0,
            })
            diagnostics.append({
                "replicate": rep, "gamma": 0.0, "validation_size": 100,
                "arm": arm, "kish_ess": 25.0 + rep,
                "minimum_node_kish_ess": 15.0 + rep,
            })
    pd.DataFrame(ledger).to_csv(root / "replicate_ledger.csv", index=False)
    pd.DataFrame(summaries).to_csv(root / "replicate_posterior_summaries.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(root / "replicate_weight_diagnostics.csv", index=False)
    out = tmp_path / "out"
    result = analyze(root, out)
    assert result["status"] == "PASS"
    assert result["planned_cells"] == 15
    strata = pd.read_csv(out / "prior_predictive_detailed_strata.csv")
    assert len(strata) == 5
    paired = pd.read_csv(out / "prior_predictive_paired_contrasts.csv")
    assert len(paired) == 15
    e_vs_c = paired[
        (paired.arm == "E_uncertainty_integrated_ipw")
        & (paired.reference_arm == "C_known_probability_ipw")
        & (paired.metric == "squared_median_error")
    ].iloc[0]
    assert e_vs_c.replicates == 3
