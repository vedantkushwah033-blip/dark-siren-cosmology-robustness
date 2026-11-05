from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts import run_estimated_inclusion_curvature_screen_from_config as configured


def valid_config(output_dir: Path) -> dict:
    return {
        "schema_version": 1,
        "replicates": 2,
        "n_parent": 1000,
        "n_detected": 5,
        "validation_sizes": [100, 300],
        "seed": 20261105,
        "h0_points": 31,
        "z_points": 101,
        "curvature_levels": [-2, 0, 2],
        "output_dir": str(output_dir),
    }


def test_load_config_validates_and_hashes_canonical_manifest(tmp_path):
    path = tmp_path / "config.json"
    config = valid_config(tmp_path / "artifacts")
    path.write_text(json.dumps(config), encoding="utf-8")

    loaded, digest = configured.load_config(path)

    assert loaded == config
    assert len(digest) == 64
    assert configured.load_config(path)[1] == digest


@pytest.mark.parametrize(
    "mutate",
    [
        lambda c: c.update({"unexpected": True}),
        lambda c: c.update({"validation_sizes": [100, 100]}),
        lambda c: c.update({"replicates": True}),
        lambda c: c.update({"curvature_levels": []}),
        lambda c: c.update({"schema_version": 99}),
    ],
)
def test_load_config_rejects_invalid_manifests(tmp_path, mutate):
    config = valid_config(tmp_path / "artifacts")
    mutate(config)
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(config), encoding="utf-8")

    with pytest.raises(ValueError):
        configured.load_config(path)


def test_configured_runner_records_manifest_and_hash(tmp_path, monkeypatch):
    output = tmp_path / "artifacts"
    config = valid_config(output)
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    monkeypatch.setattr("sys.argv", ["runner", "--config", str(config_path)])

    def fake_run(command, check):
        assert check is True
        assert "--replicates" in command and "2" in command
        output.mkdir(parents=True)
        (output / "screen_report.json").write_text(
            json.dumps({"status": "SCREEN_NOT_CONFIRMATORY"}), encoding="utf-8"
        )
        (output / "raw.csv").write_text("replicate,status\\n1,success\\n", encoding="utf-8")

    monkeypatch.setattr(configured.subprocess, "run", fake_run)
    configured.main()

    report = json.loads((output / "screen_report.json").read_text(encoding="utf-8"))
    canonical = json.loads((output / "experiment_config.canonical.json").read_text(encoding="utf-8"))
    assert report["configuration"] == config
    assert canonical == config
    assert report["configuration_sha256"] == configured.load_config(config_path)[1]
    assert report["configuration_schema_version"] == 1
    provenance = json.loads((output / "provenance.json").read_text(encoding="utf-8"))
    assert provenance["artifact_sha256"]["experiment_config.canonical.json"]
    assert provenance["artifact_sha256"]["raw.csv"]
    assert "source_revision" in provenance and "package_versions" in provenance
