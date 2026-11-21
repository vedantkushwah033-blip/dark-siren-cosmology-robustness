"""Atomic, restart-aware persistence for the cell lifecycle ledger.

The plan and ledger share one checkpoint file so a process interruption cannot
leave a newly-created plan without its matching initial NOT_RUN ledger.
Read-modify-write transitions are protected by an OS advisory lock: atomic
replacement alone prevents torn JSON but does not prevent lost concurrent
updates from independent workers.
"""
from __future__ import annotations

import errno
import hashlib
import json
import os
import tempfile
import time
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
    """Serialize checkpoint mutations across threads and independent processes.

    The lock file is intentionally persistent. The operating system releases
    the advisory lock automatically if a process exits or crashes, so a stale
    lock *file* cannot strand future workers.
    """
    target = Path(path)
    lock_path = target.with_name(target.name + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+b") as stream:
        if os.name == "nt":
            import msvcrt

            stream.seek(0, os.SEEK_END)
            if stream.tell() == 0:
                stream.write(b"\\0")
                stream.flush()
            while True:
                try:
                    stream.seek(0)
                    msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
                    break
                except OSError as exc:
                    if exc.errno not in (errno.EACCES, errno.EDEADLK, errno.EAGAIN):
                        raise
                    time.sleep(0.05)
            try:
                yield
            finally:
                stream.seek(0)
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl

            fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


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
    if target.exists():
        raise FileExistsError(f"refusing to overwrite existing run checkpoint: {target}")
    plan = [dict(row) for row in planned_rows]
    ledger = initial_ledger(plan)
    checkpoint = {
        "schema_version": SCHEMA_VERSION,
        "plan_sha256": _canonical_hash(plan),
        "planned_cells": plan,
        "cell_ledger": ledger,
    }
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


def _mutate_checkpoint(path: str | Path, mutation) -> dict:
    """Apply one read-modify-write transition under the checkpoint lock."""
    with _checkpoint_lock(path):
        checkpoint = read_checkpoint(path)
        mutation(checkpoint)
        return _write_updated(path, checkpoint)


def mark_cell_started(path: str | Path, key: tuple, started_at: str) -> dict:
    """Record an attempted cell while retaining NOT_RUN until audited completion."""
    def mutation(checkpoint: dict) -> None:
        matches = [i for i, row in enumerate(checkpoint["cell_ledger"])
                   if _cell_key(row) == tuple(key)]
        if len(matches) != 1:
            raise KeyError(f"planned cell key must match exactly once: {key}")
        index = matches[0]
        row = checkpoint["cell_ledger"][index]
        if row["status"] != "NOT_RUN" or row.get("started_at"):
            raise ValueError("cell has already been started or completed")
        checkpoint["cell_ledger"][index] = set_cell_outcome(
            row, status="NOT_RUN", audit_status="NOT_RUN",
            reason_code="RUN_INTERRUPTED",
            reason_detail="Execution started; no audited terminal outcome has been persisted.",
            started_at=started_at,
        )

    return _mutate_checkpoint(path, mutation)


def record_cell_outcome(path: str | Path, key: tuple, **outcome) -> dict:
    """Persist one audited terminal cell outcome without changing plan provenance."""
    def mutation(checkpoint: dict) -> None:
        matches = [i for i, row in enumerate(checkpoint["cell_ledger"])
                   if _cell_key(row) == tuple(key)]
        if len(matches) != 1:
            raise KeyError(f"planned cell key must match exactly once: {key}")
        index = matches[0]
        current = checkpoint["cell_ledger"][index]
        if current["status"] != "NOT_RUN":
            raise ValueError("terminal cell outcome cannot be overwritten")
        if not current.get("started_at") and outcome.get("status") in {"SUCCESS", "NON_EVALUABLE"}:
            raise ValueError("successful outcomes require a persisted start transition")
        persisted_start = current.get("started_at")
        supplied_start = outcome.get("started_at")
        if persisted_start and supplied_start not in (None, persisted_start):
            raise ValueError("persisted cell start timestamp cannot be changed")
        cell_outcome = dict(outcome)
        cell_outcome["started_at"] = persisted_start or supplied_start
        checkpoint["cell_ledger"][index] = set_cell_outcome(current, **cell_outcome)

    return _mutate_checkpoint(path, mutation)


def finalize_interrupted_cells(path: str | Path, reason_detail: str = "") -> dict:
    """Mark every still-unfinished planned cell as interrupted, preserving started_at."""
    def mutation(checkpoint: dict) -> None:
        for index, row in enumerate(checkpoint["cell_ledger"]):
            if row["status"] == "NOT_RUN" and row.get("reason_code") == "RUN_NOT_STARTED":
                checkpoint["cell_ledger"][index] = set_cell_outcome(
                    row, status="NOT_RUN", audit_status="NOT_RUN",
                    reason_code="RUN_INTERRUPTED",
                    reason_detail=reason_detail or "Run ended before this planned cell completed.",
                )

    return _mutate_checkpoint(path, mutation)
