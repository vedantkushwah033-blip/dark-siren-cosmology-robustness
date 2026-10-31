from __future__ import annotations

import numpy as np
import pytest

from dark_siren_robustness.multiframe_design import (
    FRAME_SCHEMA,
    MANIFEST_SCHEMA,
    generate_finite_frame,
    generate_frame_grid,
    validate_multiframe_manifest,
)


@pytest.fixture
def manifest():
    return {
        "schema": MANIFEST_SCHEMA,
        "root_seed": 20261031,
        "population_size": 500,
        "frames_per_scenario": 4,
        "outer_repetitions": 5,
        "bootstrap_replicates": 19,
        "nominal_level": 0.95,
        "sample_sizes": [50, 200],
        "methods": [
            {"method_id": "pseudo_studentized"},
            {"method_id": "wald_srswor"},
            {"method_id": "ordinary_percentile"},
        ],
        "scenarios": [
            {"scenario_id": 0, "name": "Normal", "family": "normal",
             "params": {"loc": 0.0, "scale": 1.0}},
            {"scenario_id": 1, "name": "Heavy right tail", "family": "gamma",
             "params": {"shape": 0.7, "scale": 1.0}},
            {"scenario_id": 2, "name": "Lognormal", "family": "lognormal",
             "params": {"meanlog": 0.0, "sigmalog": 0.8}},
        ],
    }


def test_manifest_normalizes_and_preserves_explicit_scenarios(manifest):
    result = validate_multiframe_manifest(manifest)
    assert result["schema"] == MANIFEST_SCHEMA
    assert [s["scenario_id"] for s in result["scenarios"]] == [0, 1, 2]
    assert result["sample_sizes"] == [50, 200]


def test_manifest_rejects_duplicate_scenario_and_method_ids(manifest):
    duplicate = {**manifest, "scenarios": manifest["scenarios"] + [manifest["scenarios"][0]]}
    with pytest.raises(ValueError, match="scenario_id"):
        validate_multiframe_manifest(duplicate)
    duplicate_methods = {**manifest, "methods": [manifest["methods"][0], manifest["methods"][0]]}
    with pytest.raises(ValueError, match="method_id"):
        validate_multiframe_manifest(duplicate_methods)


def test_manifest_rejects_sample_size_not_below_population(manifest):
    bad = {**manifest, "sample_sizes": [500]}
    with pytest.raises(ValueError, match="smaller"):
        validate_multiframe_manifest(bad)


@pytest.mark.parametrize("bad", [True, -1, 1.5])
def test_manifest_rejects_invalid_integer_fields(manifest, bad):
    with pytest.raises((TypeError, ValueError)):
        validate_multiframe_manifest({**manifest, "frames_per_scenario": bad})


def test_frame_generation_is_reproducible_and_immutable(manifest):
    a = generate_finite_frame(manifest, 1, 2)
    b = generate_finite_frame(manifest, 1, 2)
    assert a["schema"] == FRAME_SCHEMA
    assert a["seed"] == b["seed"]
    assert np.array_equal(a["values"], b["values"])
    assert a["population_mean"] == b["population_mean"]
    with pytest.raises(ValueError):
        a["values"][0] = 100.0


def test_frame_coordinates_generate_distinct_frames(manifest):
    frames = [generate_finite_frame(manifest, 0, i) for i in range(4)]
    assert len({frame["seed"] for frame in frames}) == 4
    assert len({frame["population_mean"] for frame in frames}) == 4


def test_frame_grid_has_balanced_explicit_hierarchy(manifest):
    frames = generate_frame_grid(manifest)
    assert len(frames) == 3 * 4
    assert [(x["scenario_id"], x["frame_id"]) for x in frames] == [
        (sid, fid) for sid in range(3) for fid in range(4)
    ]


def test_unknown_scenario_and_negative_frame_are_rejected(manifest):
    with pytest.raises(ValueError, match="unknown"):
        generate_finite_frame(manifest, 99, 0)
    with pytest.raises(ValueError):
        generate_finite_frame(manifest, 0, -1)


def test_generation_does_not_mutate_manifest_or_global_rng(manifest):
    before = np.random.default_rng(123).random(5)
    generate_frame_grid(manifest)
    after = np.random.default_rng(123).random(5)
    assert np.array_equal(before, after)
    assert manifest["frames_per_scenario"] == 4


@pytest.mark.parametrize("bad", [0, 1, True, 1.5])
def test_manifest_requires_at_least_two_bootstrap_replicates(manifest, bad):
    with pytest.raises((TypeError, ValueError), match="bootstrap_replicates"):
        validate_multiframe_manifest({**manifest, "bootstrap_replicates": bad})
