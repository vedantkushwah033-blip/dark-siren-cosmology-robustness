"""Stable, order-independent seed derivation for multi-frame experiments.

The functions in this module are deliberately pure.  They do not mutate a
Generator, depend on loop order, or share state between methods.
"""
from __future__ import annotations

from numbers import Integral

import numpy as np

_SEED_SCHEMA = "multiframe_seed_schedule_v1"
_DOMAIN_FRAME = 101
_DOMAIN_OUTER_SAMPLE = 211
_DOMAIN_METHOD = 307
_DOMAIN_INNER = 401


def _nonnegative_integer(name: str, value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be an integer, not a boolean")
    result = int(value)
    if result < 0:
        raise ValueError(f"{name} must be non-negative")
    if result > 2**32 - 1:
        raise ValueError(f"{name} must be no larger than 2**32 - 1")
    return result


def _seed(root_seed: int, domain: int, *coordinates: int) -> int:
    root = _nonnegative_integer("root_seed", root_seed)
    values = [_nonnegative_integer("seed coordinate", x) for x in coordinates]
    sequence = np.random.SeedSequence([root, domain, *values])
    return int(sequence.generate_state(1, dtype=np.uint64)[0])


def seed_schedule_metadata() -> dict[str, str]:
    """Return immutable metadata identifying the seed derivation contract."""
    return {
        "seed_schedule": _SEED_SCHEMA,
        "frame_domain": str(_DOMAIN_FRAME),
        "outer_sample_domain": str(_DOMAIN_OUTER_SAMPLE),
        "method_domain": str(_DOMAIN_METHOD),
        "inner_domain": str(_DOMAIN_INNER),
        "bit_generator": "numpy.random.PCG64 via default_rng",
    }


def frame_seed(root_seed: int, scenario_id: int, frame_id: int) -> int:
    """Seed one independently generated frame."""
    return _seed(root_seed, _DOMAIN_FRAME, scenario_id, frame_id)


def outer_sample_seed(
    root_seed: int, scenario_id: int, frame_id: int,
    sample_size: int, replicate_id: int,
) -> int:
    """Seed one SRSWOR sample, shared across all methods for that key."""
    return _seed(
        root_seed, _DOMAIN_OUTER_SAMPLE,
        scenario_id, frame_id, sample_size, replicate_id,
    )


def method_seed(
    root_seed: int, scenario_id: int, frame_id: int,
    sample_size: int, replicate_id: int, method_id: int,
) -> int:
    """Seed one method-specific stochastic stream, independent of call order."""
    return _seed(
        root_seed, _DOMAIN_METHOD,
        scenario_id, frame_id, sample_size, replicate_id, method_id,
    )


def inner_attempt_seed(
    root_seed: int, scenario_id: int, frame_id: int,
    sample_size: int, replicate_id: int, method_id: int, attempt_id: int,
) -> int:
    """Seed an individual inner attempt when attempt-level replay is needed."""
    return _seed(
        root_seed, _DOMAIN_INNER,
        scenario_id, frame_id, sample_size, replicate_id, method_id, attempt_id,
    )
