import math

import pytest

from scripts.run_paired_seed_sensitivity import mean_se


def test_mean_se_known_values():
    mean, se = mean_se([1.0, 2.0, 3.0])
    assert mean == pytest.approx(2.0)
    assert se == pytest.approx(math.sqrt(1.0 / 3.0))


def test_mean_se_rejects_empty_input():
    with pytest.raises(ValueError, match="empty"):
        mean_se([])


def test_mean_se_rejects_non_finite_values():
    with pytest.raises(ValueError, match="finite"):
        mean_se([1.0, float("nan")])


def test_mean_se_single_value_has_nan_standard_error():
    mean, se = mean_se([4.0])
    assert mean == pytest.approx(4.0)
    assert math.isnan(se)



def test_paired_differences_uses_b_minus_a():
    from scripts.run_paired_seed_sensitivity import paired_differences

    a = [{"run": 0, "median": 70.0}, {"run": 1, "median": 71.0}]
    b = [{"run": 0, "median": 72.0}, {"run": 1, "median": 70.5}]
    assert paired_differences(a, b, "median") == pytest.approx([2.0, -0.5])


def test_paired_differences_rejects_different_row_counts():
    from scripts.run_paired_seed_sensitivity import paired_differences

    with pytest.raises(ValueError, match="different row counts"):
        paired_differences([{"run": 0, "median": 1.0}], [], "median")


def test_paired_differences_rejects_misaligned_run_ids():
    from scripts.run_paired_seed_sensitivity import paired_differences

    a = [{"run": 0, "median": 70.0}, {"run": 1, "median": 71.0}]
    b = [{"run": 1, "median": 72.0}, {"run": 0, "median": 70.5}]
    with pytest.raises(ValueError, match="Pair alignment mismatch"):
        paired_differences(a, b, "median")


def test_paired_differences_rejects_non_finite_metrics():
    from scripts.run_paired_seed_sensitivity import paired_differences

    a = [{"run": 0, "median": 70.0}]
    b = [{"run": 0, "median": float("nan")}]
    with pytest.raises(ValueError, match="Non-finite"):
        paired_differences(a, b, "median")


def test_paired_differences_supports_coverage_indicators():
    from scripts.run_paired_seed_sensitivity import paired_differences

    complete = [
        {"run": 0, "covered_68": 1},
        {"run": 1, "covered_68": 1},
        {"run": 2, "covered_68": 0},
    ]
    incomplete = [
        {"run": 0, "covered_68": 0},
        {"run": 1, "covered_68": 1},
        {"run": 2, "covered_68": 1},
    ]
    assert paired_differences(complete, incomplete, "covered_68") == pytest.approx(
        [-1.0, 0.0, 1.0]
    )
