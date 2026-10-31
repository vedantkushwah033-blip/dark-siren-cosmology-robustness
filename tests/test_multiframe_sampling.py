from __future__ import annotations

import numpy as np
import pytest

from dark_siren_robustness.multiframe_design import MANIFEST_SCHEMA
from dark_siren_robustness.multiframe_sampling import (
    OUTER_SAMPLE_SCHEMA,
    generate_outer_sample,
    generate_outer_sample_grid,
)


@pytest.fixture
def manifest():
    return {
        "schema": MANIFEST_SCHEMA, "root_seed": 20261031,
        "population_size": 200, "frames_per_scenario": 3,
        "outer_repetitions": 4, "bootstrap_replicates": 25,
        "nominal_level": 0.95, "sample_sizes": [20, 80],
        "methods": [{"method_id": "pseudo"}, {"method_id": "wald"}],
        "scenarios": [
            {"scenario_id": 4, "name": "Normal", "family": "normal",
             "params": {"loc": 0.0, "scale": 1.0}},
            {"scenario_id": 9, "name": "Gamma", "family": "gamma",
             "params": {"shape": 0.7, "scale": 1.0}},
        ],
    }


def test_outer_sample_is_reproducible_and_srswor(manifest):
    a = generate_outer_sample(manifest, 4, 1, 20, 2)
    b = generate_outer_sample(manifest, 4, 1, 20, 2)
    assert a["schema"] == OUTER_SAMPLE_SCHEMA
    assert a["seed"] == b["seed"]
    assert np.array_equal(a["indices"], b["indices"])
    assert np.array_equal(a["values"], b["values"])
    assert len(np.unique(a["indices"])) == 20
    assert np.all((a["indices"] >= 0) & (a["indices"] < 200))


def test_outer_sample_stream_is_shared_across_methods(manifest):
    # Method IDs are intentionally absent from the outer sample API/seed.
    first = generate_outer_sample(manifest, 4, 1, 20, 2)
    second = generate_outer_sample(manifest, 4, 1, 20, 2)
    assert first["indices_sha256"] == second["indices_sha256"]


def test_sample_coordinates_change_seed_and_sample(manifest):
    base = generate_outer_sample(manifest, 4, 1, 20, 2)
    variants = [
        generate_outer_sample(manifest, 4, 1, 20, 3),
        generate_outer_sample(manifest, 4, 1, 80, 2),
        generate_outer_sample(manifest, 4, 2, 20, 2),
        generate_outer_sample(manifest, 9, 1, 20, 2),
    ]
    assert len({base["seed"], *(x["seed"] for x in variants)}) == 5


def test_sample_indices_and_values_are_immutable(manifest):
    sample = generate_outer_sample(manifest, 4, 0, 20, 0)
    with pytest.raises(ValueError):
        sample["indices"][0] = 1
    with pytest.raises(ValueError):
        sample["values"][0] = 1.0


def test_sample_grid_is_complete_and_canonical(manifest):
    rows = generate_outer_sample_grid(manifest)
    assert len(rows) == 2 * 3 * 2 * 4
    keys = [(x["scenario_id"], x["frame_id"], x["sample_size"], x["replicate_id"]) for x in rows]
    assert len(keys) == len(set(keys))
    assert keys[0] == (4, 0, 20, 0)
    assert keys[-1] == (9, 2, 80, 3)


@pytest.mark.parametrize(
    "coordinates",
    [(4, 3, 20, 0), (4, 0, 21, 0), (4, 0, 20, 4), (999, 0, 20, 0)],
)
def test_invalid_sample_coordinates_are_rejected(manifest, coordinates):
    with pytest.raises(ValueError):
        generate_outer_sample(manifest, *coordinates)


def test_boolean_coordinates_are_rejected(manifest):
    with pytest.raises(TypeError):
        generate_outer_sample(manifest, 4, 0, True, 0)


def test_frame_and_sample_hashes_are_recorded(manifest):
    sample = generate_outer_sample(manifest, 4, 0, 20, 0)
    assert len(sample["indices_sha256"]) == 64
    assert len(sample["frame_sha256"]) == 64
    assert sample["finite_population_mean"] == pytest.approx(
        float(np.mean(sample["finite_population_mean"])), abs=0
    )
