from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np\nimport pytest

from dark_siren_robustness.shared_input_bundle import load_shared_input_manifest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "generate_shared_inputs.py"
SPEC = importlib.util.spec_from_file_location("generate_shared_inputs", SCRIPT)
GENERATOR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(GENERATOR)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def test_shared_input_only_generator_freezes_complete_bundle_without_inference(tmp_path):
    output = tmp_path / "shared"
    result = GENERATOR.run(SimpleNamespace(
        output_dir=output,
        seed=12001,
        event_seed=12002,
        validation_seed=12003,
        validation_sample_seed=12004,
        n_parent=600,
        n_detected=2,
        validation_frame=60,
        validation_sizes=[20, 40],
        inclusion_degree=1,
        curvature_gamma=0.0,
        curvature_beta=-3.9,
        h0_true=70.0,
        h0_prior_min=50.0,
        h0_prior_max=90.0,
        h0_points=81,
        scott_multiplier=0.7,
        quadrature_order=3,
        sigma_mpc=140.0,
        detection_threshold_mpc=250.0,
    ))
    assert result["status"] == "PASS"
    assert result["inference_performed"] is False
    assert not (output / "results").exists()
    config = json.loads((output / "config.json").read_text())
    assert config["status"] == "shared_inputs_only_no_inference"
    assert config["validation_sizes"] == [20, 40]

    shared = load_shared_input_manifest(output, expected_validation_sizes=[20, 40])
    assert shared["status"] == "FROZEN"
    assert len(shared["files"]) == 8

    outer = json.loads((output / "manifest.json").read_text())
    assert outer["status"] == "PASS"
    for relative, metadata in outer["files"].items():
        path = output / relative
        assert path.stat().st_size == metadata["bytes"]
        assert _sha256(path) == metadata["sha256"]


def test_shared_input_loader_rejects_modified_validation_slice(tmp_path):
    output = tmp_path / "shared"
    GENERATOR.run(SimpleNamespace(
        output_dir=output,
        seed=13001,
        event_seed=13002,
        validation_seed=13003,
        validation_sample_seed=13004,
        n_parent=600,
        n_detected=2,
        validation_frame=30,
        validation_sizes=[20],
        inclusion_degree=1,
        curvature_gamma=0.0,
        curvature_beta=-3.9,
        h0_true=70.0,
        h0_prior_min=50.0,
        h0_prior_max=90.0,
        h0_points=41,
        scott_multiplier=0.7,
        quadrature_order=3,
        sigma_mpc=140.0,
        detection_threshold_mpc=250.0,
    ))
    with (output / "observable" / "validation_n20.csv").open("a") as stream:
        stream.write("\n")
    with pytest.raises(ValueError, match="SHA-256"):
        load_shared_input_manifest(output, expected_validation_sizes=[20])
