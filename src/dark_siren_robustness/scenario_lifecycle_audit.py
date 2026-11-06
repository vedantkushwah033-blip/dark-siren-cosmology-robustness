"""Independent structural audit for scenario-lifecycle checkpoint artifacts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

STATUSES = {"SUCCESS", "FAILED", "NON_EVALUABLE", "NOT_RUN"}
TERMINAL = {"SUCCESS", "FAILED", "NON_EVALUABLE"}
SCENARIO_STATES = {"NOT_STARTED", "RUNNING", "SUCCESS", "PARTIAL_FAILURE", "FAILED"}


def audit_scenario_lifecycle(artifact_dir: Path) -> dict:
    """Independently reconcile ledger keys, lifecycle state, and denominators."""
    snapshot_path = artifact_dir / "checkpoint.json"
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
    if snapshot.get("schema") != "scenario_lifecycle_checkpoint_v1":
        raise AssertionError("unsupported atomic checkpoint snapshot schema")
    if not isinstance(snapshot.get("ledger_records"), list):
        raise AssertionError("checkpoint snapshot has no ledger record list")
    ledger = pd.DataFrame(snapshot["ledger_records"])
    state = snapshot.get("run_state")
    if not isinstance(state, dict):
        raise AssertionError("checkpoint snapshot has no run-state object")

    required = {
        "scenario", "draw", "validation_size", "ledger_schema", "status",
        "stage", "failure_code", "failure_reason", "transition_count",
    }
    missing = required - set(ledger.columns)
    if missing:
        raise AssertionError(f"ledger missing required fields: {sorted(missing)}")
    if state.get("schema") != "scenario_lifecycle_state_v1":
        raise AssertionError("unsupported scenario lifecycle state schema")
    if not isinstance(state.get("checkpoint_sequence"), int) or state["checkpoint_sequence"] < 0:
        raise AssertionError("invalid checkpoint sequence")
    n_scenarios = state.get("n_scenarios")
    n_draws = state.get("n_draws")
    sizes = state.get("validation_sizes")
    if isinstance(n_scenarios, bool) or not isinstance(n_scenarios, int) or n_scenarios < 1:
        raise AssertionError("invalid scenario count in run state")
    if isinstance(n_draws, bool) or not isinstance(n_draws, int) or n_draws < 1:
        raise AssertionError("invalid draw count in run state")
    if not isinstance(sizes, list) or not sizes or any(
        isinstance(x, bool) or not isinstance(x, int) or x < 1 for x in sizes
    ) or len(set(sizes)) != len(sizes):
        raise AssertionError("invalid validation-size plan in run state")
    if not ledger["ledger_schema"].eq("planned_attempt_ledger_v1").all():
        raise AssertionError("ledger schema mismatch")
    if ledger.duplicated(["scenario", "draw", "validation_size"]).any():
        raise AssertionError("duplicate planned attempt keys")
    expected = {
        (scenario, draw, size)
        for scenario in range(n_scenarios)
        for draw in range(n_draws)
        for size in sizes
    }
    observed = set(map(tuple, ledger[
        ["scenario", "draw", "validation_size"]
    ].itertuples(index=False, name=None)))
    if expected != observed:
        raise AssertionError(
            f"planned denominator mismatch: missing={len(expected-observed)}, "
            f"unexpected={len(observed-expected)}"
        )
    if not ledger["status"].isin(STATUSES).all():
        raise AssertionError("unknown attempt status")
    if not ledger["transition_count"].isin((0, 1)).all():
        raise AssertionError("invalid transition count")
    for row in ledger.itertuples(index=False):
        if row.status in TERMINAL:
            if row.transition_count != 1:
                raise AssertionError("terminal attempt must have exactly one transition")
            if not str(row.stage).strip():
                raise AssertionError("terminal attempt has no stage")
            if row.status == "SUCCESS":
                if row.failure_code or row.failure_reason:
                    raise AssertionError("successful attempt carries failure metadata")
            elif not str(row.failure_code).strip() or not str(row.failure_reason).strip():
                raise AssertionError("failed attempt lacks failure metadata")
        elif row.transition_count != 0:
            raise AssertionError("NOT_RUN attempt must have zero transitions")

    counts = {
        "n_planned": int(len(ledger)),
        "n_success": int(ledger["status"].eq("SUCCESS").sum()),
        "n_failed": int(ledger["status"].eq("FAILED").sum()),
        "n_non_evaluable": int(ledger["status"].eq("NON_EVALUABLE").sum()),
        "n_not_run": int(ledger["status"].eq("NOT_RUN").sum()),
    }
    if sum(counts[k] for k in (
        "n_success", "n_failed", "n_non_evaluable", "n_not_run"
    )) != counts["n_planned"]:
        raise AssertionError("attempt-state denominator does not reconcile")
    if state.get("attempt_counts") is not None and state["attempt_counts"] != counts:
        raise AssertionError("run-state attempt counts disagree with independently derived counts")

    scenario_states = state.get("scenario_states")
    if not isinstance(scenario_states, dict) or set(scenario_states) != {
        str(i) for i in range(n_scenarios)
    }:
        raise AssertionError("scenario state map does not cover planned scenarios")
    for scenario in range(n_scenarios):
        group = ledger.loc[ledger["scenario"].eq(scenario)]
        statuses = set(group["status"])
        state_value = scenario_states[str(scenario)]
        if state_value not in SCENARIO_STATES:
            raise AssertionError("unknown scenario lifecycle state")
        scenario_failure_recorded = str(scenario) in state.get("scenario_failures", {})
        if "NOT_RUN" in statuses:
            if state_value not in {"RUNNING", "NOT_STARTED"}:
                raise AssertionError("unfinished scenario has an incompatible lifecycle state")
        elif scenario_failure_recorded:
            if state_value != "FAILED":
                raise AssertionError("scenario-level failure record requires FAILED lifecycle state")
        elif statuses == {"SUCCESS"}:
            if state_value != "SUCCESS":
                raise AssertionError("all-success scenario has an inconsistent lifecycle state")
        elif state_value != "PARTIAL_FAILURE":
            raise AssertionError("draw-level or reconciliation outcomes require PARTIAL_FAILURE state")

    run_status = state.get("run_status")
    if counts["n_not_run"]:
        if run_status != "RUNNING":
            raise AssertionError("incomplete ledger must retain RUNNING run status")
    else:
        expected_run_status = (
            "SUCCESS" if counts["n_failed"] == 0 and counts["n_non_evaluable"] == 0
            else "COMPLETED_WITH_FAILURES"
        )
        if run_status != expected_run_status:
            raise AssertionError("run status disagrees with terminal attempt states")
    return {
        "status": "SCENARIO_LIFECYCLE_AUDIT_PASS",
        **counts,
        "n_scenarios": n_scenarios,
        "validation_sizes": sizes,
        "run_status": run_status,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit_scenario_lifecycle(args.artifact_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
