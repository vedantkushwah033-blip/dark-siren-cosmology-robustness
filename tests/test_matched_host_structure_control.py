import importlib.util
from pathlib import Path
import pytest

PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_matched_host_structure_control.py"
SPEC = importlib.util.spec_from_file_location("matched_host_structure_control", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_generator_is_reproducible_and_selected():
    a = __import__("numpy").random.default_rng(31)
    b = __import__("numpy").random.default_rng(31)
    x = [module.draw_detected(70.0, a) for _ in range(100)]
    y = [module.draw_detected(70.0, b) for _ in range(100)]
    assert x == y
    assert all(host in (0, 1) and 0 < distance < module.DETECTION_MAX_MPC for host, distance in x)


def test_matched_run_uses_same_observation_for_all_methods():
    rows, summaries, manifest = module.run(replicates=100, seed=44, grid_size=1001, bootstrap_draws=1000)
    assert len(rows) == len(module.TRUTHS) * 100
    assert len(summaries) == len(module.TRUTHS) * 6
    assert manifest["raw_rows"] == len(rows)
    assert manifest["summary_rows"] == len(summaries)
    for row in rows:
        assert 0 < row["observed_distance_mpc"] < module.DETECTION_MAX_MPC
        for method in module.METHODS:
            assert row[f"{method}_lower_68"] <= row[f"{method}_upper_68"]
            assert row[f"{method}_lower_95"] <= row[f"{method}_upper_95"]
    for row in summaries:
        if row["method"].startswith("paired_"):
            assert row["bootstrap_95_low"] <= row["bootstrap_95_high"]
        else:
            assert 0 <= row["coverage_68"] <= 1
            assert row["coverage_68_wilson_low"] <= row["coverage_68"] <= row["coverage_68_wilson_high"]


def test_rejects_small_or_invalid_design():
    with pytest.raises(ValueError, match="at least 100"):
        module.run(replicates=50)
    with pytest.raises(ValueError, match="at least 1000"):
        module.run(replicates=100, bootstrap_draws=100)
