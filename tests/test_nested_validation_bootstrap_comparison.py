from __future__ import annotations

import pytest

from scripts.run_nested_validation_bootstrap_comparison import _seed, run_comparison


def test_nested_comparison_seed_derivation_is_reproducible_and_keyed():
    first = _seed(20261030, 2, 4, 100, 1)
    assert first == _seed(20261030, 2, 4, 100, 1)
    assert first != _seed(20261030, 2, 4, 100, 2)


def test_nested_comparison_rejects_duplicate_validation_sizes_before_simulation():
    with pytest.raises(ValueError, match="must not contain duplicates"):
        run_comparison(
            scenarios=1, draws=2, n_parent=100, n_detected=2,
            validation_sizes=(20, 20), bootstrap_replicates=2, seed=7,
        )


def test_nested_comparison_rejects_invalid_bootstrap_count_before_simulation():
    with pytest.raises(ValueError, match="bootstrap_replicates"):
        run_comparison(
            scenarios=1, draws=2, n_parent=100, n_detected=2,
            validation_sizes=(20,), bootstrap_replicates=1, seed=7,
        )


@pytest.mark.parametrize("sizes", [(20.5,), (True,), (20, 30.5)])
def test_nested_comparison_rejects_noninteger_validation_sizes(sizes):
    with pytest.raises(TypeError, match="validation sizes must be integers"):
        run_comparison(
            scenarios=1, draws=2, n_parent=100, n_detected=2,
            validation_sizes=sizes, bootstrap_replicates=2, seed=7,
        )


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"scenarios": True}, "scenarios must be an integer"),
        ({"draws": 2.5}, "draws must be an integer"),
        ({"h0_true": float("nan")}, "h0_true must be finite and positive"),
        ({"sigma_mpc": 0.0}, "sigma_mpc must be finite and positive"),
        ({"threshold_mpc": float("inf")}, "threshold_mpc must be finite and positive"),
    ],
)
def test_nested_comparison_rejects_malformed_parameters_before_simulation(kwargs, message):
    parameters = dict(
        scenarios=1, draws=2, n_parent=100, n_detected=2,
        validation_sizes=(20,), bootstrap_replicates=2, seed=7,
    )
    parameters.update(kwargs)
    with pytest.raises((TypeError, ValueError), match=message):
        run_comparison(**parameters)
