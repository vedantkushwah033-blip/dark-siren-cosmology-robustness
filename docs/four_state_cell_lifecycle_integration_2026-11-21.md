# Four-state cell lifecycle integration (2026-11-21)

## Purpose and scope

The quadratic-logit five-arm H0 replicate runner now persists an immutable plan and a four-state lifecycle record for every replicate × curvature × validation-size × arm cell. The integration is intended to make execution accounting auditable and restartable. It does not change the statistical estimators or establish calibration.

## Lifecycle states

- `SUCCESS`: the cell has a completed posterior summary and the independent child audit passed.
- `NON_EVALUABLE`: reserved for a predeclared numerical/estimand failure that is still independently auditable. The current child runner does not emit this state yet.
- `FAILED`: execution or independent audit failed. The reason code distinguishes child-process errors from audit failures.
- `NOT_RUN`: the cell has not reached a terminal outcome. Its reason distinguishes not-started, interrupted, upstream-failure and not-scheduled cases.

Successful and non-evaluable rows require an independent audit PASS, an audit-report path, and offset-aware start/completion timestamps. Failed rows cannot claim audit PASS. Terminal rows cannot be overwritten.

## Checkpoint and transaction model

The runner writes `cell_lifecycle_checkpoint.json` before starting any child. It contains:
- schema version;
- a SHA-256 hash of the canonical immutable planned-cell list;
- the full planned-cell records;
- one lifecycle record per planned cell.

The checkpoint store serializes read/modify/write transactions with a stable sidecar lock (POSIX `flock` or Windows byte-range locking, plus an in-process thread lock) and uses atomic file replacement. A complete replicate × gamma child produces all validation-size × arm outcomes. Those outcomes are validated and committed in one batch transaction, so an interruption cannot leave only some cells in that child block terminal.

## Resume behavior

Use `--resume` with the same output directory and exactly the same simulation design. The runner compares the complete planned-cell records before continuing; changed seeds, curvature, prior support, sample sizes, grid settings, quadrature order or other recorded design parameters are rejected.

At resume startup, active cells are marked `RUN_INTERRUPTED`; cells never reached are marked `NOT_SCHEDULED`. Unfinished child blocks are rerun. Their first-start timestamps are retained. Already terminal child blocks are skipped and their audited saved summaries/diagnostics are reloaded.

**Execution granularity limitation:** the current child process computes all validation sizes and five arms together. If it is interrupted, resume replays that entire replicate × gamma child block. It does not continue halfway through an individual arm's numerical calculation.

## Independent audit

The aggregate auditor does not import the checkpoint store. It independently:
- verifies the outer artifact manifest;
- recomputes the canonical plan hash;
- reconciles checkpoint plans and lifecycle rows against their CSV exports;
- checks allowed status/reason/audit combinations and timestamp offsets/order;
- rejects unaccounted `NOT_RUN` cells in a finalized artifact;
- binds successful cells to the expected child audit report and corresponding child posterior-summary row;
- reconstructs the existing child audit, seed schedule, pairing, aggregate metrics and paired contrasts.

## Current limitations

1. The present child smoke currently emits only `SUCCESS`. The orchestrator fails closed if it sees another child-summary status. Arm-specific `NON_EVALUABLE` decision rules and their independently auditable evidence still need a frozen protocol and implementation.
2. Child blocks are currently all-success or all-failed in ordinary runs; the ledger schema can represent four states, but the analysis path is not yet designed to make arm-specific evaluability decisions.
3. The integration is validated using small synthetic smoke configurations. It is not confirmatory H0 calibration, observational inference, or evidence of estimator superiority.
4. This integration branch is stacked on the prior-predictive runner branch. Reconcile the branch stack before any merge; neither branch should be merged independently on the assumption that it is based on `main`.
