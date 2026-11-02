from __future__ import annotations

import numpy as np
import pytest

from dark_siren_robustness.multiframe_seed_schedule import (
    frame_seed,
    inner_attempt_seed,
    method_seed,
    outer_sample_seed,
    seed_schedule_metadata,
)


def test_seed_schedule_is_deterministic_and_domain_separated():
    args = (20261031, 1, 7)
    assert frame_seed(*args) == frame_seed(*args)
    assert frame_seed(*args) != outer_sample_seed(20261031, 1, 7, 100, 0)
    assert outer_sample_seed(20261031, 1, 7, 100, 0) != method_seed(
        20261031, 1, 7, 100, 0, 1
    )


def test_frame_and_sample_seeds_change_with_each_coordinate():
    base = frame_seed(9, 0, 0)
    assert len({base, frame_seed(9, 0, 1), frame_seed(9, 1, 0), frame_seed(10, 0, 0)}) == 4
    sample = outer_sample_seed(9, 0, 0, 100, 0)
    assert len({
        sample,
        outer_sample_seed(9, 0, 0, 100, 1),
        outer_sample_seed(9, 0, 0, 200, 0),
        outer_sample_seed(9, 0, 1, 100, 0),
    }) == 4


def test_method_stream_is_independent_of_method_execution_order():
    method_ids = (11, 22, 33)
    forward = {m: method_seed(42, 1, 2, 100, 3, m) for m in method_ids}
    reverse = {m: method_seed(42, 1, 2, 100, 3, m) for m in reversed(method_ids)}
    assert forward == reverse
    assert len(set(forward.values())) == len(method_ids)


def test_inner_attempt_seeds_are_distinct_and_replayable():
    values = [
        inner_attempt_seed(42, 1, 2, 100, 3, 11, attempt)
        for attempt in range(20)
    ]
    assert len(set(values)) == 20
    assert values[5] == inner_attempt_seed(42, 1, 2, 100, 3, 11, 5)


def test_seed_metadata_identifies_schedule():
    metadata = seed_schedule_metadata()
    assert metadata["seed_schedule"] == "multiframe_seed_schedule_v1"
    assert metadata["bit_generator"] == "numpy.random.PCG64 via default_rng"


@pytest.mark.parametrize("value", [True, -1, 2**32])
def test_seed_coordinates_reject_invalid_values(value):
    with pytest.raises((TypeError, ValueError)):
        frame_seed(1, value, 0)


def test_seed_rejects_boolean_root_seed():
    with pytest.raises(TypeError, match="root_seed"):
        frame_seed(True, 0, 0)


def test_seeds_generate_usable_numpy_generators():
    seed = outer_sample_seed(20261031, 0, 4, 100, 9)
    first = np.random.default_rng(seed).choice(1000, size=100, replace=False)
    second = np.random.default_rng(seed).choice(1000, size=100, replace=False)
    assert np.array_equal(first, second)
    assert np.unique(first).size == 100
