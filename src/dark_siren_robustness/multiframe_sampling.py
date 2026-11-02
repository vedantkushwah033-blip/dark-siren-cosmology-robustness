"""Paired SRSWOR outer-sample generation for multi-frame v3 experiments."""
from __future__ import annotations

import hashlib
from collections.abc import Mapping
from numbers import Integral

import numpy as np

from dark_siren_robustness.multiframe_design import generate_finite_frame, validate_multiframe_manifest
from dark_siren_robustness.multiframe_seed_schedule import outer_sample_seed

OUTER_SAMPLE_SCHEMA = "linear_finite_population_outer_sample_v3"


def _coordinate(name: str, value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be an integer, not a boolean")
    value = int(value)
    if value < 0:
        raise ValueError(f"{name} must be non-negative")
    return value


def _indices_sha256(indices: np.ndarray) -> str:
    canonical = np.ascontiguousarray(np.asarray(indices, dtype="<i8"))
    digest = hashlib.sha256()
    digest.update(str(canonical.shape).encode("ascii"))
    digest.update(canonical.tobytes())
    return digest.hexdigest()


def generate_outer_sample(
    manifest: Mapping[str, object],
    scenario_id: int,
    frame_id: int,
    sample_size: int,
    replicate_id: int,
) -> dict[str, object]:
    """Generate one deterministic SRSWOR sample shared by every method.

    Method identifiers intentionally do not enter the seed. The returned
    immutable indices and values are the common data input for paired method
    comparisons.
    """
    design = validate_multiframe_manifest(manifest)
    sid = _coordinate("scenario_id", scenario_id)
    fid = _coordinate("frame_id", frame_id)
    rid = _coordinate("replicate_id", replicate_id)
    n = _coordinate("sample_size", sample_size)
    if fid >= design["frames_per_scenario"]:
        raise ValueError("frame_id is outside the manifest frame grid")
    if n not in design["sample_sizes"]:
        raise ValueError("sample_size is not declared in the manifest")
    if rid >= design["outer_repetitions"]:
        raise ValueError("replicate_id is outside the manifest repetition grid")

    frame = generate_finite_frame(design, sid, fid)
    seed = outer_sample_seed(design["root_seed"], sid, fid, n, rid)
    rng = np.random.default_rng(seed)
    indices = np.asarray(
        rng.choice(design["population_size"], size=n, replace=False),
        dtype=np.int64,
    )
    values = np.asarray(frame["values"][indices], dtype=np.float64)
    indices.setflags(write=False)
    values.setflags(write=False)
    return {
        "schema": OUTER_SAMPLE_SCHEMA,
        "scenario_id": sid,
        "frame_id": fid,
        "sample_size": n,
        "replicate_id": rid,
        "seed": seed,
        "indices_sha256": _indices_sha256(indices),
        "indices": indices,
        "values": values,
        "sample_mean": float(np.mean(values)),
        "sample_variance": float(np.var(values, ddof=1)),
        "finite_population_mean": frame["population_mean"],
        "frame_seed": frame["seed"],
        "frame_sha256": hashlib.sha256(
            np.ascontiguousarray(frame["values"], dtype="<f8").tobytes()
        ).hexdigest(),
    }


def generate_outer_sample_grid(
    manifest: Mapping[str, object],
) -> list[dict[str, object]]:
    """Generate all declared outer samples in canonical manifest order."""
    design = validate_multiframe_manifest(manifest)
    samples = []
    for scenario in design["scenarios"]:
        for frame_id in range(design["frames_per_scenario"]):
            for sample_size in design["sample_sizes"]:
                for replicate_id in range(design["outer_repetitions"]):
                    samples.append(generate_outer_sample(
                        design, scenario["scenario_id"], frame_id,
                        sample_size, replicate_id,
                    ))
    return samples
