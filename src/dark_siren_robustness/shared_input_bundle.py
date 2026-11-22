"""Frozen, hash-addressed common inputs for paired five-arm experiments.

The manifest is written before model fitting so arm-level workers can later
consume the same immutable replicate/gamma input bundle.
"""
from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Sequence

_REQUIRED = (
    "observable/catalogue_observed.csv",
    "observable/gw_events_observed.csv",
    "truth/parent_truth.csv",
    "truth/gw_event_truth.csv",
    "truth/gw_generation_trace_truth.csv",
    "truth/validation_frame_truth.csv",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_shared_input_manifest(
    root: Path,
    *,
    validation_sizes: Sequence[int],
    provenance: dict[str, Any],
) -> dict[str, Any]:
    """Freeze hashes for all observable/truth inputs before inference starts."""
    root = Path(root).resolve()
    sizes = [int(value) for value in validation_sizes]
    if not sizes or sizes != sorted(set(sizes)) or any(value < 1 for value in sizes):
        raise ValueError("validation sizes must be unique, positive, and increasing")
    paths = [root / relative for relative in _REQUIRED]
    paths.extend(root / "observable" / f"validation_n{value}.csv" for value in sizes)
    entries: dict[str, dict[str, Any]] = {}
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"required shared input is missing: {path.relative_to(root)}")
        relative = str(path.relative_to(root))
        entries[relative] = {"bytes": path.stat().st_size, "sha256": _sha256(path)}
    manifest = {
        "schema_version": 1,
        "status": "FROZEN",
        "validation_sizes": sizes,
        "provenance": dict(provenance),
        "files": entries,
    }
    destination = root / "shared_input_manifest.json"
    if destination.exists():
        raise FileExistsError("shared-input manifest already exists; frozen bundles are immutable")
    fd, temporary_name = tempfile.mkstemp(
        prefix=".shared_input_manifest.", suffix=".tmp", dir=root
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as stream:
            json.dump(manifest, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        # Atomic no-replacement publication: a concurrent writer cannot silently
        # replace a bundle manifest after another worker has frozen it.
        os.link(temporary, destination)
        directory_fd = os.open(root, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
    return manifest


def load_shared_input_manifest(
    root: Path,
    *,
    expected_validation_sizes: Sequence[int] | None = None,
    expected_provenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Load and verify a frozen bundle before a cell worker consumes it.

    This verifies bytes and provenance without deserializing every CSV. Cell
    workers can then read only the inputs they need from the verified root.
    """
    root = Path(root).resolve()
    manifest_path = root / "shared_input_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError("shared-input manifest is missing")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("shared-input manifest is not valid UTF-8 JSON") from exc
    if not isinstance(manifest, dict) or manifest.get("schema_version") != 1:
        raise ValueError("unsupported shared-input manifest schema")
    if manifest.get("status") != "FROZEN":
        raise ValueError("shared-input bundle is not frozen")

    sizes = manifest.get("validation_sizes")
    if (not isinstance(sizes, list) or not sizes
            or any(type(value) is not int or value < 1 for value in sizes)
            or sizes != sorted(set(sizes))):
        raise ValueError("shared-input validation-size plan is invalid")
    if expected_validation_sizes is not None:
        expected = [int(value) for value in expected_validation_sizes]
        if sizes != expected:
            raise ValueError("shared-input validation-size plan differs from requested plan")

    provenance = manifest.get("provenance")
    if not isinstance(provenance, dict):
        raise ValueError("shared-input provenance must be an object")
    if expected_provenance is not None:
        for key, value in expected_provenance.items():
            if provenance.get(key) != value:
                raise ValueError(f"shared-input provenance mismatch: {key}")

    files = manifest.get("files")
    expected_paths = set(_REQUIRED) | {
        f"observable/validation_n{value}.csv" for value in sizes
    }
    if not isinstance(files, dict) or set(files) != expected_paths:
        raise ValueError("shared-input manifest file set differs from its declared plan")
    for relative, metadata in files.items():
        if (not isinstance(relative, str) or not isinstance(metadata, dict)
                or type(metadata.get("bytes")) is not int or metadata["bytes"] < 0
                or not isinstance(metadata.get("sha256"), str)
                or len(metadata["sha256"]) != 64):
            raise ValueError(f"invalid shared-input file metadata: {relative}")
        candidate = (root / relative).resolve()
        if not candidate.is_relative_to(root):
            raise ValueError(f"shared-input path escapes bundle root: {relative}")
        if not candidate.is_file():
            raise FileNotFoundError(f"shared input is missing: {relative}")
        if candidate.stat().st_size != metadata["bytes"]:
            raise ValueError(f"shared-input byte count mismatch: {relative}")
        if _sha256(candidate) != metadata["sha256"]:
            raise ValueError(f"shared-input SHA-256 mismatch: {relative}")
    return manifest
