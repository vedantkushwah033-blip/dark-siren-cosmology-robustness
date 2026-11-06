"""Failure-aware orchestration for scenario-nested simulation campaigns.

This module is a tested orchestration layer. It does not itself define a
scientific simulator: callers provide one scenario executor returning one row
per planned draw/validation-size pair. Planned rows are checkpointed before
execution, so an interrupted process leaves explicit NOT_RUN records.
"""
from __future__ import annotations

from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any
import json
import os
import tempfile

import pandas as pd

from dark_siren_robustness.attempt_ledger import (
    atomic_write_ledger,
    atomic_write_run_state,
    initialize_attempt_ledger,
    reconcile_attempt_ledger,
    transition_attempt,
)

SCENARIO_STATES = (
    "NOT_STARTED",
    "RUNNING",
    "SUCCESS",
    "PARTIAL_FAILURE",
    "FAILED",
)
SCENARIO_FAILURE_STAGES = {
    "scenario_generation",
    "oracle_reference",
    "catalogue_reference",
    "finite_parent_reference",
    "scenario_execution",
    "result_reconciliation",
}


class ScenarioLifecycleError(RuntimeError):
    """Expected scenario-level failure with an explicit scientific stage."""

    def __init__(self, stage: str, failure_code: str, message: str) -> None:
        if stage not in SCENARIO_FAILURE_STAGES:
            raise ValueError(f"unsupported scenario failure stage {stage!r}")
        if not failure_code.strip() or not message.strip():
            raise ValueError("scenario failure code and message must be non-empty")
        super().__init__(message)
        self.stage = stage
        self.failure_code = failure_code


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    """Write JSON atomically without destroying the previous checkpoint."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def _write_checkpoint(root: Path, ledger: pd.DataFrame, state: dict[str, Any]) -> None:
    """Publish sidecars and then one authoritative atomic two-artifact snapshot."""
    atomic_write_ledger(ledger, root / "attempt_ledger.csv")
    atomic_write_run_state(state, root / "run_state.json")
    snapshot = {
        "schema": "scenario_lifecycle_checkpoint_v1",
        "ledger_records": json.loads(ledger.to_json(orient="records")),
        "run_state": state,
    }
    _atomic_json(root / "checkpoint.json", snapshot)


def _scenario_rows(ledger: pd.DataFrame, scenario: int) -> pd.DataFrame:
    return ledger.loc[ledger["scenario"].eq(scenario)]


def _mark_scenario_non_evaluable(
    ledger: pd.DataFrame,
    *,
    scenario: int,
    stage: str,
    code: str,
    reason: str,
) -> None:
    pending = _scenario_rows(ledger, scenario)
    for row in pending.loc[pending["status"].eq("NOT_RUN")].itertuples(index=False):
        transition_attempt(
            ledger,
            scenario=int(row.scenario),
            draw=int(row.draw),
            validation_size=int(row.validation_size),
            status="NON_EVALUABLE",
            stage=stage,
            failure_code=code,
            failure_reason=reason,
        )


def _apply_scenario_records(
    ledger: pd.DataFrame,
    *,
    scenario: int,
    records: pd.DataFrame,
) -> None:
    """Validate and apply a scenario's draw-level outcomes to its planned rows."""
    required = {"draw", "validation_size", "status"}
    if not required.issubset(records.columns):
        raise ValueError(f"scenario result missing columns: {sorted(required-set(records.columns))}")
    if records.duplicated(["draw", "validation_size"]).any():
        raise ValueError("scenario result contains duplicate draw/validation-size keys")
    expected = set(
        map(
            tuple,
            _scenario_rows(ledger, scenario)[["draw", "validation_size"]].itertuples(
                index=False, name=None
            ),
        )
    )
    observed = set(map(tuple, records[["draw", "validation_size"]].itertuples(index=False, name=None)))
    unexpected = observed - expected
    if unexpected:
        raise ValueError(f"scenario result contains unexpected attempt keys: {sorted(unexpected)}")

    by_key = {
        (int(row.draw), int(row.validation_size)): row
        for row in records.itertuples(index=False)
    }
    for planned in _scenario_rows(ledger, scenario).itertuples(index=False):
        key = (int(planned.draw), int(planned.validation_size))
        row = by_key.get(key)
        if row is None:
            transition_attempt(
                ledger,
                scenario=scenario,
                draw=key[0],
                validation_size=key[1],
                status="NON_EVALUABLE",
                stage="result_reconciliation",
                failure_code="missing_attempt_record",
                failure_reason="scenario executor omitted a planned attempt record",
            )
            continue
        raw_status = str(row.status).lower()
        mapping = {
            "success": "SUCCESS",
            "failed": "FAILED",
            "non_evaluable": "NON_EVALUABLE",
        }
        if raw_status not in mapping:
            raise ValueError(f"unsupported scenario attempt status {raw_status!r}")
        terminal = mapping[raw_status]
        stage = str(getattr(row, "failure_stage", "") or "validation_attempt")
        code = str(getattr(row, "failure_code", "") or "")
        reason = str(getattr(row, "failure_reason", "") or "")
        if terminal in {"FAILED", "NON_EVALUABLE"}:
            if not code:
                code = "draw_level_failure"
            if not reason:
                reason = "scenario executor reported a non-success attempt"
        transition_attempt(
            ledger,
            scenario=scenario,
            draw=key[0],
            validation_size=key[1],
            status=terminal,
            stage=stage,
            failure_code=code,
            failure_reason=reason,
        )


def run_scenario_campaign(
    *,
    n_scenarios: int,
    n_draws: int,
    validation_sizes: Sequence[int],
    executor: Callable[[int], pd.DataFrame],
    checkpoint_dir: str | Path,
    run_metadata: dict[str, Any] | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Execute scenarios with complete planned-attempt accounting.

    The executor is called once per scenario and must return rows keyed by
    (draw, validation_size), with status values success/failed/non_evaluable.
    Any exception escaping the executor is recorded as a scenario-level
    failure; all not-yet-completed attempts for that scenario become
    NON_EVALUABLE. Attempts for scenarios not reached remain NOT_RUN.
    """
    root = Path(checkpoint_dir)
    root.mkdir(parents=True, exist_ok=True)
    raw_sizes = tuple(validation_sizes)
    ledger = initialize_attempt_ledger(
        n_scenarios=n_scenarios,
        n_draws=n_draws,
        validation_sizes=raw_sizes,
        metadata=run_metadata,
    )
    sizes = tuple(int(value) for value in raw_sizes)
    # Fail before creating any artifacts if run metadata cannot be checkpointed.
    json.dumps(run_metadata or {})
    state: dict[str, Any] = {
        "schema": "scenario_lifecycle_state_v1",
        "run_status": "RUNNING",
        "n_scenarios": int(n_scenarios),
        "n_draws": int(n_draws),
        "validation_sizes": list(sizes),
        "scenario_states": {str(i): "NOT_STARTED" for i in range(n_scenarios)},
        "scenario_failures": {},
        "metadata": run_metadata or {},
        "checkpoint_sequence": 0,
    }
    # Persist the complete denominator before starting any scientific work.
    _write_checkpoint(root, ledger, state)

    for scenario in range(n_scenarios):
        state["scenario_states"][str(scenario)] = "RUNNING"
        state["checkpoint_sequence"] += 1
        _write_checkpoint(root, ledger, state)
        try:
            records = executor(scenario)
            if not isinstance(records, pd.DataFrame):
                raise ScenarioLifecycleError(
                    "result_reconciliation", "invalid_executor_result",
                    "scenario executor must return a pandas DataFrame",
                )
            try:
                _apply_scenario_records(ledger, scenario=scenario, records=records)
            except Exception as exc:
                raise ScenarioLifecycleError(
                    "result_reconciliation", "scenario_result_reconciliation_failure",
                    f"{type(exc).__name__}: {exc}",
                ) from exc
            statuses = set(_scenario_rows(ledger, scenario)["status"])
            if statuses == {"SUCCESS"}:
                state["scenario_states"][str(scenario)] = "SUCCESS"
            elif statuses.issubset({"SUCCESS", "FAILED"}) and "FAILED" in statuses:
                state["scenario_states"][str(scenario)] = "PARTIAL_FAILURE"
            else:
                state["scenario_states"][str(scenario)] = "PARTIAL_FAILURE"
        except Exception as exc:
            stage = exc.stage if isinstance(exc, ScenarioLifecycleError) else "scenario_execution"
            code = exc.failure_code if isinstance(exc, ScenarioLifecycleError) else f"scenario_{type(exc).__name__.lower()}"
            reason = f"{type(exc).__name__}: {exc}"
            _mark_scenario_non_evaluable(
                ledger, scenario=scenario, stage=stage, code=code, reason=reason,
            )
            state["scenario_states"][str(scenario)] = "FAILED"
            state["scenario_failures"][str(scenario)] = {
                "stage": stage, "failure_code": code, "failure_reason": reason,
            }
        # Each scenario boundary is a durable recovery point. A crash during
        # the next scenario leaves its planned rows explicitly NOT_RUN.
        state["checkpoint_sequence"] += 1
        _write_checkpoint(root, ledger, state)

    counts = reconcile_attempt_ledger(
        ledger,
        n_scenarios=n_scenarios,
        n_draws=n_draws,
        validation_sizes=sizes,
        require_terminal=True,
    )
    state["run_status"] = (
        "SUCCESS" if counts["n_failed"] == 0 and counts["n_non_evaluable"] == 0
        else "COMPLETED_WITH_FAILURES"
    )
    state["attempt_counts"] = counts
    state["checkpoint_sequence"] += 1
    _write_checkpoint(root, ledger, state)
    return ledger, state
