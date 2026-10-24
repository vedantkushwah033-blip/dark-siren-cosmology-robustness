import importlib.util
import math
from pathlib import Path

PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_joint_prior_scorecard.py"
SPEC = importlib.util.spec_from_file_location("joint_scorecard", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_joint_scorecard_has_complete_keys_and_metrics():
    raw, summaries = module.prior.run_audit(replicates=100, seed=37, grid_size=1001)
    width_summary, _ = module.width.summarize_widths(raw)
    score = module.build_scorecard(summaries, width_summary)

    expected_keys = {
        (truth, pmin, pmax, method)
        for truth in module.prior.TRUTHS
        for pmin, pmax in module.prior.PRIORS
        for method in module.prior.METHODS
    }
    observed_keys = {
        (r["true_H0"], r["prior_min"], r["prior_max"], r["method"])
        for r in score
    }
    assert len(score) == 18
    assert observed_keys == expected_keys
    assert len(observed_keys) == len(score)

    for row in score:
        assert 0 <= row["abs_coverage_68_error_from_nominal"] <= 1
        assert 0 <= row["abs_coverage_95_error_from_nominal"] <= 1
        assert row["rmse"] >= 0
        assert row["mean_width_68"] > 0
        assert row["mean_width_95"] >= row["mean_width_68"]
        for metric in ("bias", "rmse", "coverage_68", "coverage_95",
                       "mean_width_68", "median_width_68", "p90_width_68",
                       "mean_width_95", "median_width_95", "p90_width_95"):
            assert math.isfinite(float(row[metric]))
        for coverage, low_key, high_key in (
            ("coverage_68", "coverage_68_wilson_low", "coverage_68_wilson_high"),
            ("coverage_95", "coverage_95_wilson_low", "coverage_95_wilson_high"),
        ):
            assert 0 <= row[low_key] <= row[coverage] <= row[high_key] <= 1
