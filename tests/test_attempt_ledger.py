from __future__ import annotations
import json
import pandas as pd
import pytest
from dark_siren_robustness.attempt_ledger import (
    initialize_attempt_ledger, transition_attempt, reconcile_attempt_ledger,
    atomic_write_ledger, atomic_write_run_state,
)

def make_ledger():
    return initialize_attempt_ledger(n_scenarios=2, n_draws=2,
        validation_sizes=(100, 200), metadata={"run_id": "test"})

def test_preallocation_and_denominator():
    ledger = make_ledger()
    assert len(ledger) == 8
    assert set(ledger.status) == {"NOT_RUN"}
    counts = reconcile_attempt_ledger(ledger, n_scenarios=2, n_draws=2,
        validation_sizes=(100, 200), require_terminal=False)
    assert counts == {"n_planned":8,"n_success":0,"n_failed":0,
        "n_non_evaluable":0,"n_not_run":8}
    with pytest.raises(ValueError, match="NOT_RUN"):
        reconcile_attempt_ledger(ledger, n_scenarios=2, n_draws=2,
            validation_sizes=(100, 200))

def test_terminal_transitions_and_failure_states():
    ledger = make_ledger()
    transition_attempt(ledger, scenario=0, draw=0, validation_size=100,
        status="SUCCESS", stage="posterior")
    transition_attempt(ledger, scenario=0, draw=0, validation_size=200,
        status="FAILED", stage="reference_fit", failure_code="fit",
        failure_reason="optimizer failed")
    transition_attempt(ledger, scenario=0, draw=1, validation_size=100,
        status="NON_EVALUABLE", stage="scenario_generation",
        failure_code="generation", failure_reason="frame unavailable")
    for s in (0,1):
        for d in (0,1):
            for size in (100,200):
                mask=(ledger.scenario==s)&(ledger.draw==d)&(ledger.validation_size==size)
                if ledger.loc[mask,"status"].iloc[0] == "NOT_RUN":
                    transition_attempt(ledger, scenario=s, draw=d, validation_size=size,
                        status="SUCCESS", stage="posterior")
    counts = reconcile_attempt_ledger(ledger, n_scenarios=2, n_draws=2,
        validation_sizes=(100,200))
    assert counts == {"n_planned":8,"n_success":6,"n_failed":1,
        "n_non_evaluable":1,"n_not_run":0}
    with pytest.raises(ValueError, match="already transitioned"):
        transition_attempt(ledger, scenario=0, draw=0, validation_size=100,
            status="SUCCESS", stage="again")

@pytest.mark.parametrize("status,code,reason", [
    ("FAILED","","reason"), ("FAILED","code",""), ("NON_EVALUABLE","","reason"),
    ("SUCCESS","","bad"), ("NOT_RUN","",""),
])
def test_invalid_transitions_rejected(status,code,reason):
    ledger=make_ledger()
    with pytest.raises(ValueError):
        transition_attempt(ledger,scenario=0,draw=0,validation_size=100,
            status=status,stage="stage",failure_code=code,failure_reason=reason)

def test_plan_validation():
    with pytest.raises(ValueError):
        initialize_attempt_ledger(n_scenarios=0,n_draws=2,validation_sizes=(10,))
    with pytest.raises(TypeError):
        initialize_attempt_ledger(n_scenarios=1.5,n_draws=2,validation_sizes=(10,))
    with pytest.raises(ValueError):
        initialize_attempt_ledger(n_scenarios=1,n_draws=2,validation_sizes=(10,10))

def test_reconciler_rejects_duplicate_and_missing_keys():
    ledger=make_ledger()
    with pytest.raises(ValueError,match="duplicate"):
        reconcile_attempt_ledger(pd.concat([ledger,ledger.iloc[[0]]]),
            n_scenarios=2,n_draws=2,validation_sizes=(100,200),require_terminal=False)
    with pytest.raises(ValueError,match="missing"):
        reconcile_attempt_ledger(ledger.iloc[:-1],
            n_scenarios=2,n_draws=2,validation_sizes=(100,200),require_terminal=False)

def test_atomic_checkpoint_files(tmp_path):
    ledger=make_ledger()
    csv=tmp_path/"nested"/"attempts.csv"
    js=tmp_path/"nested"/"run.json"
    atomic_write_ledger(ledger,csv)
    atomic_write_run_state({"stage":"running"},js)
    assert len(pd.read_csv(csv))==8
    assert json.loads(js.read_text())["stage"]=="running"
    assert not list(csv.parent.glob(".*.tmp"))

def test_metadata_cannot_overwrite_attempt_identity_or_state():
    with pytest.raises(ValueError, match="reserved fields"):
        initialize_attempt_ledger(n_scenarios=1,n_draws=1,
            validation_sizes=(10,),metadata={"status":"SUCCESS"})

def test_reserved_payload_rejection_does_not_partially_transition():
    ledger = make_ledger()
    with pytest.raises(ValueError, match="reserved ledger fields"):
        transition_attempt(ledger, scenario=0, draw=0, validation_size=100,
            status="SUCCESS", stage="posterior", payload={"status":"FAILED"})
    row = ledger.iloc[0]
    assert row["status"] == "NOT_RUN"
    assert row["transition_count"] == 0


def test_failed_atomic_serialization_preserves_previous_checkpoint(tmp_path):
    path = tmp_path / "run.json"
    atomic_write_run_state({"stage":"stable"}, path)
    with pytest.raises(TypeError):
        atomic_write_run_state({"bad": object()}, path)
    assert json.loads(path.read_text(encoding="utf-8")) == {"stage":"stable"}
    assert not list(tmp_path.glob(".*.tmp"))
