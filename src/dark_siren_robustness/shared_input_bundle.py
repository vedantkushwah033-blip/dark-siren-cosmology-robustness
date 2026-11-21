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
