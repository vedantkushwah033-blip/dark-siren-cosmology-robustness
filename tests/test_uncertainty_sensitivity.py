import numpy as np
import pytest

from dark_siren_robustness.uncertainty_sensitivity import (
    ARMS, paired_ab_contrasts, posterior_summary, run_experiment,
)


def test_grid_posterior_summary_is_ordered_and_within_support():
    grid = np.linspace(50.0, 90.0, 1601)
    result = posterior_summary(grid, -0.5 * ((grid - 70.0) / 2.0) ** 2)
    assert 50.0 <= result["lower_95"] <= result["lower_68"] <= result["median"]
    assert result["median"] <= result["upper_68"] <= result["upper_95"] <= 90.0
    assert result["posterior_sd"] > 0


def test_experiment_is_deterministic_and_a_b_observations_are_paired():
    raw1, summaries1 = run_experiment(3, 123, 201, truths=(55.0, 70.0), fractional_errors=(0.1,))
    raw2, summaries2 = run_experiment(3, 123, 201, truths=(55.0, 70.0), fractional_errors=(0.1,))
    assert raw1 == raw2
    assert summaries1 == summaries2
    by_key = {(r["true_h0"], r["fractional_error"], r["replicate"], r["arm"]): r for r in raw1}
    for truth in (55.0, 70.0):
        for rep in range(3):
            a = by_key[(truth, 0.1, rep, ARMS[0])]
            b = by_key[(truth, 0.1, rep, ARMS[1])]
            assert a["observed_distance_mpc"] == b["observed_distance_mpc"]
            assert a["standard_normal"] == b["standard_normal"]
            assert a["generated_sigma_mpc"] == b["generated_sigma_mpc"]
    assert len(raw1) == 2 * 1 * 3 * 3
    assert len(summaries1) == 2 * 1 * 3
    assert all("pull_median" in row for row in raw1)
    assert all(0 <= row["coverage_68_wilson_low"] <= row["coverage_68"] <= row["coverage_68_wilson_high"] <= 1
               for row in summaries1)
    assert all(0 <= row["coverage_95_wilson_low"] <= row["coverage_95"] <= row["coverage_95_wilson_high"] <= 1
               for row in summaries1)


def test_paired_contrasts_are_reproducible_and_have_bootstrap_intervals():
    raw, _ = run_experiment(3, 123, 201, truths=(55.0, 70.0), fractional_errors=(0.1,))
    contrasts = paired_ab_contrasts(raw, seed=123, bootstrap_resamples=100)
    assert len(contrasts) == 2 * 1 * 5
    for item in contrasts:
        assert item["direction"] == "A_minus_B"
        assert item["n_paired_replicates"] == 3
        assert item["bootstrap_resamples"] == 100
        assert item["bootstrap_95_low"] <= item["bootstrap_95_high"]
    assert contrasts == paired_ab_contrasts(raw, seed=123, bootstrap_resamples=100)


@pytest.mark.parametrize("kwargs", [
    {"replicates": 1}, {"seed": -1}, {"grid_points": 100},
    {"truths": (49.0,)}, {"fractional_errors": (0.0,)},
])
def test_experiment_rejects_invalid_design(kwargs):
    with pytest.raises(ValueError):
        run_experiment(**kwargs)
