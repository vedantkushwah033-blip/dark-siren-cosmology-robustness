"""Atomic, restart-aware persistence for the cell lifecycle ledger.

The plan and ledger share one checkpoint file so a process interruption cannot
leave a newly-created plan without its matching initial NOT_RUN ledger.
"""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Mapping

from .lifecycle_ledger import (
    atomic_write_json,
    initial_ledger,
    set_cell_outcome,
    validate_ledger,
)

SCHEMA_VERSION = 1


def _canonical_hash(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _cell_key(row: Mapping) -> tuple:
    return (row["replicate"], row["gamma"], row["validation_size"], row["arm"])


@contextmanager
def _checkpoint_lock(path: str | Path):
    """Serialize checkpoint read/modify/write transactions across processes.

    The advisory lock is held on a stable sidecar inode (not the atomically
    replaced checkpoint inode), and is automatically released if a process exits.
    """
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    lock_path = target.with_name(target.name + ".lock")
    with lock_path.open("a", encoding="utf-8") as lock_stream:
        fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)


def _atomic_create_json(path: str | Path, value) -> None:
    """Atomically create a JSON file without replacing an existing run checkpoint."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(value, indent=2, sort_keys=True) + "\n"
    fd, temporary = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".tmp", dir=target.parent
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        # A hard link is an atomic no-replace operation on the same filesystem.
        os.link(temporary, target)
        os.unlink(temporary)
        directory_fd = os.open(target.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def initialize_checkpoint(path: str | Path, planned_rows: Iterable[Mapping]) -> dict:
    """Atomically persist immutable plan + complete NOT_RUN ledger; never overwrite."""
    target = Path(path)
    plan = [dict(row) for row in planned_rows]
    ledger = initial_ledger(plan)
    checkpoint = {
        "schema_version": SCHEMA_VERSION,
        "plan_sha256": _canonical_hash(plan),
        "planned_cells": plan,
        "cell_ledger": ledger,
    }
    with _checkpoint_lock(target):
        _atomic_create_json(target, checkpoint)
    return checkpoint


def read_checkpoint(path: str | Path) -> dict:
    """Load and validate a checkpoint, including immutable plan provenance."""
    checkpoint = json.loads(Path(path).read_text(encoding="utf-8"))
    if checkpoint.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("unsupported cell-ledger checkpoint schema")
    plan = checkpoint.get("planned_cells")
    ledger = checkpoint.get("cell_ledger")
    if not isinstance(plan, list) or not isinstance(ledger, list):
        raise ValueError("checkpoint must contain planned_cells and cell_ledger lists")
    if checkpoint.get("plan_sha256") != _canonical_hash(plan):
        raise ValueError("immutable planned-cell hash mismatch")
    validate_ledger(plan, ledger)
    plan_by_key = {_cell_key(row): row for row in plan}
    for row in ledger:
        expected = plan_by_key[_cell_key(row)]
        for field, value in expected.items():
            if row.get(field) != value:
                raise ValueError(f"planned provenance mismatch for {field!r}")
    return checkpoint


def _write_updated(path: str | Path, checkpoint: dict) -> dict:
    validate_ledger(checkpoint["planned_cells"], checkpoint["cell_ledger"])
    atomic_write_json(path, checkpoint)
    return checkpoint


def mark_cell_started(path: str | Path, key: tuple, started_at: str) -> dict:
    """Record an attempted cell while retaining NOT_RUN until audited completion."""
    if not isinstance(started_at, str) or not started_at.strip():
        raise ValueError("mark_cell_started requires a non-empty ISO-8601 timestamp")
    with _checkpoint_lock(path):
        checkpoint = read_checkpoint(path)
        matches = [i for i, row in enumerate(checkpoint["cell_ledger"])
                   if _cell_key(row) == tuple(key)]
        if len(matches) != 1:
            raise KeyError(f"planned cell key must match exactly once: {key}")
        index = matches[0]
        row = checkpoint["cell_ledger"][index]
        if row["status"] != "NOT_RUN":
            raise ValueError("cell has already reached a terminal outcome")
        if row.get("started_at"):
            if row.get("reason_code") != "RUN_INTERRUPTED":
                raise ValueError("cell is already active or is not marked recoverable")
            # Retry an interrupted cell without rewriting its immutable first-start time.
            checkpoint["cell_ledger"][index] = set_cell_outcome(
                row, status="NOT_RUN", audit_status="NOT_RUN",
                reason_code="RUN_NOT_STARTED",
                reason_detail="Cell resumed after an interrupted attempt; first-start timestamp retained.",
                started_at=row["started_at"],
            )
        else:
            checkpoint["cell_ledger"][index] = set_cell_outcome(
                row, status="NOT_RUN", audit_status="NOT_RUN",
                reason_code="RUN_NOT_STARTED",
                reason_detail="Execution is active; no audited terminal outcome has been persisted.",
                started_at=started_at,
            )
        return _write_updated(path, checkpoint)


def record_cell_outcome(path: str | Path, key: tuple, **outcome) -> dict:
    """Persist one audited terminal cell outcome without changing plan provenance."""
    with _checkpoint_lock(path):
        checkpoint = read_checkpoint(path)
        matches = [i for i, row in enumerate(checkpoint["cell_ledger"])
                   if _cell_key(row) == tuple(key)]
        if len(matches) != 1:
            raise KeyError(f"planned cell key must match exactly once: {key}")
        index = matches[0]
        current = checkpoint["cell_ledger"][index]
        if current["status"] != "NOT_RUN":
            raise ValueError("terminal cell outcome cannot be overwritten")
        if outcome.get("status") in {"SUCCESS", "NON_EVALUABLE", "FAILED"} and not current.get("started_at"):
            raise ValueError("terminal outcomes require a persisted start transition")
        persisted_start = current.get("started_at")
        supplied_start = outcome.get("started_at")
        if persisted_start and supplied_start not in (None, persisted_start):
            raise ValueError("persisted cell start timestamp cannot be changed")
        outcome["started_at"] = persisted_start or supplied_start
        checkpoint["cell_ledger"][index] = set_cell_outcome(current, **outcome)
        return _write_updated(path, checkpoint)


def finalize_interrupted_cells(path: str | Path, reason_detail: str = "") -> dict:
    """Mark every still-unfinished planned cell as interrupted, preserving started_at."""
    with _checkpoint_lock(path):
        checkpoint = read_checkpoint(path)
        for index, row in enumerate(checkpoint["cell_ledger"]):
            if row["status"] == "NOT_RUN" and row.get("reason_code") == "RUN_NOT_STARTED":
                was_started = bool(row.get("started_at"))
                checkpoint["cell_ledger"][index] = set_cell_outcome(
                    row, status="NOT_RUN", audit_status="NOT_RUN",
                    reason_code="RUN_INTERRUPTED" if was_started else "NOT_SCHEDULED",
                    reason_detail=reason_detail or (
                        "Run ended while this cell was active."
                        if was_started else "Run ended before this planned cell was scheduled."
                    ),
                    started_at=row.get("started_at"),
                )
        return _write_updated(path, checkpoint)
