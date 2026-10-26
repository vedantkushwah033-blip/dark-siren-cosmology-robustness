import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "chimera_fixture_downloader", ROOT / "scripts/prepare_chimera_v2_benchmark_data.py"
)
DOWNLOADER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DOWNLOADER)


class FakeResponse:
    def __init__(self, payload):
        self.payload = payload
        self.offset = 0

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def read(self, size):
        chunk = self.payload[self.offset:self.offset + size]
        self.offset += len(chunk)
        return chunk


def test_download_verified_writes_only_after_checksum_pass(tmp_path, monkeypatch):
    payload = b"small public-data fixture for checksum test"
    expected = hashlib.md5(payload).hexdigest()
    monkeypatch.setattr(DOWNLOADER.urllib.request, "urlopen", lambda url, timeout: FakeResponse(payload))
    destination = tmp_path / "fixture.bin"
    result = DOWNLOADER.download_verified("https://example.invalid/fixture.bin", destination, expected)
    assert destination.read_bytes() == payload
    assert result["size_bytes"] == len(payload)
    assert result["md5"] == expected
    assert result["sha256"] == hashlib.sha256(payload).hexdigest()
    assert not list(tmp_path.glob("*.part"))


def test_download_verified_rejects_corruption_and_cleans_partial(tmp_path, monkeypatch):
    payload = b"corrupted payload"
    monkeypatch.setattr(DOWNLOADER.urllib.request, "urlopen", lambda url, timeout: FakeResponse(payload))
    destination = tmp_path / "fixture.bin"
    with pytest.raises(ValueError, match="MD5 mismatch"):
        DOWNLOADER.download_verified("https://example.invalid/fixture.bin", destination, "0" * 32)
    assert not destination.exists()
    assert not list(tmp_path.glob("*.part"))


def test_prepare_writes_manifest_for_already_verified_file(tmp_path, monkeypatch):
    payload = b"manifest fixture"
    md5 = hashlib.md5(payload).hexdigest()
    monkeypatch.setitem(DOWNLOADER.FILES, "tiny", {
        "filename": "tiny.bin", "url": "https://example.invalid/tiny.bin",
        "md5": md5, "size_bytes": len(payload), "role": "test fixture",
    })
    (tmp_path / "tiny.bin").write_bytes(payload)
    manifest = DOWNLOADER.prepare(tmp_path, selected=("tiny",))
    saved = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert saved["files"]["tiny"]["status"] == "already_verified"
    assert saved["files"]["tiny"]["sha256"] == hashlib.sha256(payload).hexdigest()
    assert saved["license_and_redistribution_reviewed"] is False
    assert manifest["files"]["tiny"]["size_bytes"] == len(payload)
