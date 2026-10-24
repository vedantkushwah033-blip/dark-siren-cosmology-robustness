import importlib.util
from pathlib import Path
import numpy as np

PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_fixed_truth_coverage_map.py"
SPEC = importlib.util.spec_from_file_location("fixed_truth_coverage_map", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)

def test_fixed_truth_coverage_map_is_finite_and_reports_mc_uncertainty():
    rows, summaries = module.run_audit(replicates_per_truth=50, seed=12, grid_size=1001)
    assert len(rows) == 7 * 50
    assert len(summaries) == 7 * 3
    assert {r["true_H0"] for r in summaries} == set(module.TRUTHS)
    coverage_rows = [r for r in summaries if r["method"] in {"uncorrected", "corrected"}]
    for row in coverage_rows:
        for key in ("bias", "rmse", "coverage_68", "coverage_68_wilson_low",
                    "coverage_68_wilson_high", "coverage_95", "coverage_95_wilson_low",
                    "coverage_95_wilson_high", "mean_width_68", "mean_width_95"):
            assert np.isfinite(float(row[key]))
        assert 0 <= row["coverage_68"] <= 1
        assert row["coverage_68_wilson_low"] <= row["coverage_68"] <= row["coverage_68_wilson_high"]
        assert row["coverage_95_wilson_low"] <= row["coverage_95"] <= row["coverage_95_wilson_high"]
