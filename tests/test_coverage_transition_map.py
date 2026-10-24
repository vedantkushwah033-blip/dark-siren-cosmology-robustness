import importlib.util
from pathlib import Path

PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_coverage_transition_map.py"
SPEC = importlib.util.spec_from_file_location("coverage_transition", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_transition_map_is_complete_and_uses_paired_design():
    rows, summaries = module.run(replicates=100, seed=42, grid_size=1001)
    assert len(rows) == len(module.TRUTH_GRID) * len(module.PRIOR_GRID) * 100
    assert len(summaries) == len(module.TRUTH_GRID) * len(module.PRIOR_GRID) * 2
    assert {r["true_H0"] for r in rows} == set(module.TRUTH_GRID)
    assert {(r["prior_min"], r["prior_max"]) for r in rows} == set(module.PRIOR_GRID)
    assert {r["method"] for r in summaries} == {"corrected", "uncorrected"}
    # Each truth/replicate shares the same generated distance across prior supports.
    for truth in module.TRUTH_GRID:
        for rep in range(100):
            cell = [r for r in rows if r["true_H0"] == truth and r["replicate"] == rep]
            assert len(cell) == len(module.PRIOR_GRID)
            assert len({r["observed_distance_mpc"] for r in cell}) == 1
            assert len({r["host_index"] for r in cell}) == 1
    for row in summaries:
        assert 0 <= row["coverage_68"] <= 1
        assert 0 <= row["coverage_95"] <= 1
        assert 0 <= row["coverage_68_wilson_low"] <= row["coverage_68"] <= row["coverage_68_wilson_high"] <= 1
        assert 0 <= row["coverage_95_wilson_low"] <= row["coverage_95"] <= row["coverage_95_wilson_high"] <= 1


def test_paired_method_contrasts_are_complete_and_finite():
    rows, _ = module.run(replicates=100, seed=43, grid_size=1001)
    contrasts = module.paired_method_contrasts(rows, seed=44, bootstrap_draws=1000)
    expected = len(module.TRUTH_GRID) * len(module.PRIOR_GRID) * 7
    assert len(contrasts) == expected
    assert {r["metric"] for r in contrasts} == {
        "posterior_median", "absolute_error", "squared_error",
        "width_68", "width_95", "coverage_68_indicator", "coverage_95_indicator"
    }
    for row in contrasts:
        assert row["n_pairs"] == 100
        assert row["bootstrap_95_low"] <= row["mean_corrected_minus_uncorrected"] <= row["bootstrap_95_high"] or (
            row["bootstrap_95_low"] <= row["bootstrap_95_high"]
        )
        assert row["bootstrap_95_low"] <= row["bootstrap_95_high"]
        assert row["bootstrap_draws"] == 1000
