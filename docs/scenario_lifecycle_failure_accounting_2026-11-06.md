# Scenario-level failure accounting and checkpoint protocol

**Status:** implementation protocol / non-confirmatory infrastructure.  
**Date:** 2026-11-06.

## Purpose

The finite-parent SRSWOR experiment has a nested design: parent scenarios are
the top-level units, and validation draws and validation sizes are nested
within each scenario. A scenario-level reference failure can occur before any
validation draw is evaluated. If the runner aborts at that point, an ordinary
draw-level records file can silently omit the planned attempts.

The lifecycle layer makes the planned denominator explicit before scientific
execution and separates upstream scenario failures from individual
validation-attempt failures.

## Attempt-level states

Every planned `(scenario, draw, validation_size)` key is initialized as
`NOT_RUN`. It may transition once to one terminal state:

- `SUCCESS`: the draw-level estimator and downstream posterior completed.
- `FAILED`: the draw was evaluated but its estimator or downstream numerical
  calculation failed.
- `NON_EVALUABLE`: an upstream scenario/reference failure or missing
  attempt-level output prevented a valid draw-level evaluation.
- `NOT_RUN`: the attempt was planned but had not reached a terminal outcome
  at the last durable checkpoint.

The denominator is always the number of planned Cartesian keys, not the
number of returned scientific records. The independent audit checks key
uniqueness, completeness, legal states, one-way transition counts, failure
metadata, and exact reconciliation of all four states.

## Scenario-level states and failure stages

Scenario lifecycle state is tracked separately:

- `NOT_STARTED`
- `RUNNING`
- `SUCCESS`
- `PARTIAL_FAILURE`
- `FAILED`

Upstream failures can be classified at:

- `scenario_generation`
- `oracle_reference`
- `catalogue_reference`
- `finite_parent_reference`
- `validation_sample_generation`
- `result_reconciliation`
- `scenario_execution` (fallback for otherwise uncategorized executor errors)

When a scenario fails before its draw records are available, its still-pending
planned attempts are terminalized as `NON_EVALUABLE`. Draw-level failures
remain `FAILED` and are not conflated with scenario-level failures. Missing
planned records are made explicit as `NON_EVALUABLE`, while duplicate or
unexpected keys cause a scenario reconciliation failure.

## Checkpoint protocol

The checkpointed runner writes the complete planned ledger before invoking
the first scenario. It writes a new checkpoint before each scenario and after
each scenario boundary.

Files:

- `attempt_ledger.csv`: human-readable attempt-ledger sidecar.
- `run_state.json`: human-readable scenario-state sidecar.
- `checkpoint.json`: authoritative atomic snapshot containing both the
  ledger records and run-state object.

The single-file `checkpoint.json` snapshot is the consistency boundary for
the independent auditor. It is atomically replaced only after the sidecars
have been written. If the process stops while a scenario is running, the last
snapshot preserves all planned keys and records the active scenario as
`RUNNING`; the current scenario's attempts remain `NOT_RUN` unless the
scenario-boundary checkpoint was completed.

This is durable accounting, **not resume support**. The runner does not yet
reload a prior checkpoint and continue from it. Use a fresh output directory
for a new execution; do not treat the sidecar CSV or JSON as a resumable state
store.

## Finite-parent benchmark integration

The opt-in entry point is:

```python
from pathlib import Path
from dark_siren_robustness.finite_parent_sampling import (
    run_finite_parent_benchmark_checkpointed,
)

records, ledger, state = run_finite_parent_benchmark_checkpointed(
    checkpoint_dir=Path("artifacts/finite_parent_lifecycle_smoke"),
    n_scenarios=2,
    n_validation_draws=2,
    n_parent=1200,
    n_detected=6,
    validation_sizes=(100,),
    seed=20261110,
    h0_points=31,
    z_points=101,
)
```

Each scenario is run independently, but its random stream is keyed by the
original global scenario index. This is intended to preserve the same
simulation, thinning, and validation seeds as the corresponding scenario in
the multi-scenario runner. A dedicated regression test compares those
identities and generated artifacts.

`records` contains scientific draw-level outputs only for scenarios whose
executor returned. `ledger` contains the complete planned denominator,
including upstream failures. These outputs answer different questions and
must not be substituted for one another.

Run the independent structural audit:

```bash
python scripts/audit_scenario_lifecycle.py artifacts/finite_parent_lifecycle_smoke
```

The lifecycle auditor verifies accounting and checkpoint consistency. It
does **not** independently reconstruct KDEs, logistic fits, likelihoods, or
posteriors. The pre-existing finite-parent numerical auditor remains a
separate scientific validation gate; it has not yet been upgraded to consume
the new lifecycle snapshot for partially evaluable campaigns.

## Scientific boundaries

This infrastructure does not change the finite-parent estimand, the
synthetic data-generating process, the fixed detected-event-count likelihood,
or the scientific interpretation of the experiment. It does not validate a
real gravitational-wave catalogue, provide a measured value of (H_0), or
establish superiority of an estimator. A confirmatory campaign remains
blocked until the lifecycle outputs are fully integrated with the independent
scientific artifact audit and the remaining numerical/model-misspecification
gates are closed.
