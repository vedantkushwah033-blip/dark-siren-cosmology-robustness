from pathlib import Path

import pytest

from dark_siren_robustness.shared_input_bundle import write_shared_input_manifest


def _populate(root: Path) -> None:
    for relative in (
        "observable/catalogue_observed.csv",
        "observable/gw_events_observed.csv",
        "observable/validation_n10.csv",
        "observable/validation_n20.csv",
        "truth/parent_truth.csv",
        "truth/gw_event_truth.csv",
        "truth/gw_generation_trace_truth.csv",
        "truth/validation_frame_truth.csv",
    ):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("id,value\n1,2\n", encoding="utf-8")


def test_shared_input_manifest_freezes_hashes_and_provenance(tmp_path: Path):
    _populate(tmp_path)
    result = write_shared_input_manifest(
        tmp_path,
        validation_sizes=[10, 20],
        provenance={"seed": 1, "curvature_gamma": 0.5},
    )
    assert result["status"] == "FROZEN"
    assert len(result["files"]) == 8
    assert (tmp_path / "shared_input_manifest.json").is_file()


def test_shared_input_manifest_refuses_overwrite(tmp_path: Path):
    _populate(tmp_path)
    kwargs = {"validation_sizes": [10, 20], "provenance": {"seed": 1}}
    write_shared_input_manifest(tmp_path, **kwargs)
    with pytest.raises(FileExistsError, match="immutable"):
        write_shared_input_manifest(tmp_path, **kwargs)


def test_shared_input_manifest_rejects_missing_validation_slice(tmp_path: Path):
    _populate(tmp_path)
    (tmp_path / "observable/validation_n20.csv").unlink()
    with pytest.raises(FileNotFoundError, match="validation_n20.csv"):
        write_shared_input_manifest(
            tmp_path, validation_sizes=[10, 20], provenance={"seed": 1}
        )


def test_stale_legacy_temp_file_does_not_block_manifest_creation(tmp_path: Path):
    _populate(tmp_path)
    stale = tmp_path / "shared_input_manifest.json.tmp"
    stale.write_text("incomplete prior attempt", encoding="utf-8")
    result = write_shared_input_manifest(
        tmp_path,
        validation_sizes=[10, 20],
        provenance={"seed": 7},
    )
    assert result["status"] == "FROZEN"
    assert (tmp_path / "shared_input_manifest.json").is_file()
    assert stale.read_text(encoding="utf-8") == "incomplete prior attempt"
