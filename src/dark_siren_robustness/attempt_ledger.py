"""Reusable planned-attempt ledger primitives for nested simulation campaigns.

The ledger deliberately separates planned attempts from completed outcomes.
It is a foundation for benchmark runners; callers must still integrate it with
scenario generation, reference fitting, and artifact publication.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
from typing import Any

import pandas as pd

ATTEMPT_STATUSES = ("SUCCESS", "FAILED", "NON_EVALUABLE", "NOT_RUN")
TERMINAL_STATUSES = frozenset(("SUCCESS", "FAILED", "NON_EVALUABLE"))
LEDGER_SCHEMA = "planned_attempt_ledger_v1"
KEY_COLUMNS = ("scenario", "draw", "validation_size")


def initialize_attempt_ledger(
    *,
    n_scenarios: int,
    n_draws: int,
    validation_sizes: tuple[int, ...] | list[int],
    metadata: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Preallocate every planned scenario/draw/size key as NOT_RUN."""
    for name, value, minimum in (
        ("n_scenarios", n_scenarios, 1),
        ("n_draws", n_draws, 1),
    ):
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(f"{name} must be an integer")
        if value < minimum:
            raise ValueError(f"{name} must be >= {minimum}")
    sizes = tuple(validation_sizes)
    if not sizes:
        raise ValueError("validation_sizes must not be empty")
    if any(isinstance(v, bool) or not isinstance(v, int) for v in sizes):
        raise TypeError("validation_sizes must contain integers")
    if any(v < 1 for v in sizes):
        raise ValueError("validation sizes must be positive")
    if len(set(sizes)) != len(sizes):
        raise ValueError("validation_sizes must be unique")

    rows = []
    for scenario in range(n_scenarios):
        for draw in range(n_draws):
            for size in sorted(sizes):
                rows.append({
                    "scenario": scenario,
                    "draw": draw,
                    "validation_size": size,
                    "ledger_schema": LEDGER_SCHEMA,
                    "status": "NOT_RUN",
                    "stage": "not_started",
                    "failure_code": "",
                    "failure_reason": "",
                    "transition_count": 0,
                    **(metadata or {}),
                })
    return pd.DataFrame(rows)


def transition_attempt(
    ledger: pd.DataFrame,
    *,
    scenario: int,
    draw: int,
    validation_size: int,
    status: str,
    stage: str,
    failure_code: str = "",
    failure_reason: str = "",
    payload: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Transition one unique planned key from NOT_RUN to a terminal status."""
    if status not in TERMINAL_STATUSES:
        raise ValueError(f"status must be one of {sorted(TERMINAL_STATUSES)}")
    if not isinstance(stage, str) or not stage.strip():
        raise ValueError("stage must be a non-empty string")
    if status in {"FAILED", "NON_EVALUABLE"}:
        if not failure_code.strip() or not failure_reason.strip():
            raise ValueError(f"{status} requires failure_code and failure_reason")
    elif failure_code or failure_reason:
        raise ValueError("SUCCESS must not carry failure metadata")

    matches = (
        ledger["scenario"].eq(scenario)
        & ledger["draw"].eq(draw)
        & ledger["validation_size"].eq(validation_size)
    )
    positions = ledger.index[matches].tolist()
    if len(positions) != 1:
        raise ValueError(f"planned attempt key must match exactly one row; found {len(positions)}")
    index = positions[0]
    if ledger.at[index, "status"] != "NOT_RUN":
        raise ValueError("attempt has already transitioned from NOT_RUN")
    ledger.at[index, "status"] = status
    ledger.at[index, "stage"] = stage
    ledger.at[index, "failure_code"] = failure_code
    ledger.at[index, "failure_reason"] = failure_reason
    ledger.at[index, "transition_count"] = int(ledger.at[index, "transition_count"]) + 1
    for key, value in (payload or {}).items():
        if key in KEY_COLUMNS or key in {"status", "ledger_schema", "transition_count"}:
            raise ValueError(f"payload cannot overwrite reserved ledger field {key!r}")
        ledger.at[index, key] = value
    return ledger


def reconcile_attempt_ledger(
    ledger: pd.DataFrame,
    *,
    n_scenarios: int,
    n_draws: int,
    validation_sizes: tuple[int, ...] | list[int],
    require_terminal: bool = True,
) -> dict[str, int]:
    """Validate key completeness, uniqueness, states, and planned denominator."""
    required = set(KEY_COLUMNS) | {
        "ledger_schema", "status", "stage", "failure_code",
        "failure_reason", "transition_count",
    }
    missing = required.difference(ledger.columns)
    if missing:
        raise ValueError(f"ledger missing columns: {sorted(missing)}")
    if not ledger["ledger_schema"].eq(LEDGER_SCHEMA).all():
        raise ValueError("ledger schema mismatch")
    if ledger.duplicated(list(KEY_COLUMNS)).any():
        raise ValueError("duplicate planned attempt keys")
    expected = {
        (scenario, draw, size)
        for scenario in range(n_scenarios)
        for draw in range(n_draws)
        for size in validation_sizes
    }
    observed = set(map(tuple, ledger[list(KEY_COLUMNS)].itertuples(index=False, name=None)))
    if observed != expected:
        raise ValueError(
            f"planned attempt keys do not reconcile: missing={len(expected-observed)}, "
            f"unexpected={len(observed-expected)}"
        )
    if not ledger["status"].isin(ATTEMPT_STATUSES).all():
        raise ValueError("ledger contains an unknown attempt status")
    if require_terminal and ledger["status"].eq("NOT_RUN").any():
        raise ValueError("ledger contains planned attempts that remain NOT_RUN")
    if not ledger["transition_count"].isin((0, 1)).all():
        raise ValueError("transition_count must be zero or one")
    terminal = ledger["status"].isin(TERMINAL_STATUSES)
    if not ledger.loc[terminal, "transition_count"].eq(1).all():
        raise ValueError("terminal attempt must have exactly one transition")
    if not ledger.loc[~terminal, "transition_count"].eq(0).all():
        raise ValueError("NOT_RUN attempt must have zero transitions")
    for _, row in ledger.loc[terminal].iterrows():
        if row["status"] == "SUCCESS":
            if str(row["failure_code"]) or str(row["failure_reason"]):
                raise ValueError("SUCCESS row carries failure metadata")
        elif not str(row["failure_code"]).strip() or not str(row["failure_reason"]).strip():
            raise ValueError("failed/non-evaluable row lacks failure metadata")
    return {
        "n_planned": int(len(ledger)),
        "n_success": int(ledger["status"].eq("SUCCESS").sum()),
        "n_failed": int(ledger["status"].eq("FAILED").sum()),
        "n_non_evaluable": int(ledger["status"].eq("NON_EVALUABLE").sum()),
        "n_not_run": int(ledger["status"].eq("NOT_RUN").sum()),
    }


def _atomic_replace(path: Path, writer) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer(handle)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def atomic_write_ledger(ledger: pd.DataFrame, path: str | Path) -> None:
    """Atomically checkpoint the current ledger as CSV in its destination directory."""
    destination = Path(path)
    _atomic_replace(destination, lambda handle: ledger.to_csv(handle, index=False))


def atomic_write_run_state(state: dict[str, Any], path: str | Path) -> None:
    """Atomically checkpoint scenario/run lifecycle metadata as JSON."""
    destination = Path(path)
    _atomic_replace(
        destination,
        lambda handle: handle.write(json.dumps(state, indent=2, sort_keys=True) + "\n"),
    )
