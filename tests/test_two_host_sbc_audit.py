import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_two_host_sbc_audit.py"
SPEC = importlib.util.spec_from_file_location("two_host_sbc_audit", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_two_host_generator_uses_valid_host_support_and_weights():
    rows = module.generate_replicates(200, 7)
    assert len(rows) == 400
    for row in rows:
        assert row["host_index"] in (0, 1)
        assert row["host_redshift"] in (0.02, 0.04)
        assert module.PRIOR_MIN <= row["true_H0"] <= module.PRIOR_MAX
        assert row["observed_distance_mpc"] > 0
        assert row["sigma_mpc"] > 0


def test_two_host_sbc_and_grid_refinement_are_finite():
    rows, summaries, diffs = module.run_audit(
        replicates=200, seed=91, grid_sizes=(501, 1001)
    )
    assert len(rows) == 2 * 200 * 2
    assert len(summaries) == 4
    assert len(diffs) == 2
    assert all(np.isfinite(row["posterior_cdf_rank"]) for row in rows)
    assert all(0 <= row["posterior_cdf_rank"] <= 1 for row in rows)
    assert all(0 <= row["ks_pvalue"] <= 1 for row in summaries)
    assert all(row["max_abs_rank_difference"] >= 0 for row in diffs)
