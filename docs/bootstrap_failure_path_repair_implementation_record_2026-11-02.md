# Bootstrap failure-path repair branch: implementation record

**Base:** PR #17 head `7c612efb1b075ea2df0e665b93091fa78c77fe84`  
**Branch:** `research/bootstrap-failure-path-repairs-2026-11-02`  
**Date:** 2026-11-02  
**Scope:** narrow failure-path integrity repairs; exploratory synthetic software only.

## Changes made

### 1. Preserve wrapper-level exception metadata

In `scripts/run_nested_validation_bootstrap_comparison.py`, the exception handler around the bootstrap wrapper now stores:
- `wrapper_failure_class = type(exc).__name__`;
- `wrapper_failure_message = str(exc)`.

These values flow into `attempt_audit.jsonl` through the existing serialization fields. Successful wrapper calls continue to set both fields to null. This makes the wrapper-failure replay branch's expected metadata available for actual wrapper exceptions.

A regression test in `tests/test_nested_validation_bootstrap_comparison.py` forces a wrapper exception and checks the saved attempt rows contain the exact class and message.

### 2. Make index generation a prerequisite for an inner attempt

In `src/dark_siren_robustness/posterior_validation_bootstrap.py`, paired-bootstrap and design-aware index generation now occur before the per-attempt posterior-fit exception handler.

Rationale: an attempt is only defined after its resampling indices have been generated. If index generation itself raises, the function now propagates the exception to the runner, where it is recorded as a wrapper-level failure. It no longer returns a nominally complete attempt trace with missing index arrays. Failures after index generation remain attempt-level failures and retain their generated indices.

Tests in `tests/test_posterior_validation_bootstrap.py` inject an index-generation exception for each method and assert it propagates as a wrapper-level failure.

## Additional repairs

### 3. Preserve scenario-level target failures as explicit partial-run artifacts

In `scripts/run_nested_validation_bootstrap_comparison.py`:
- Finite-parent target precondition failures and expected posterior-computation exceptions are captured as scenario failures instead of aborting the entire run before artifacts are written.
- Planned outer-draw/validation-size/method rows are still emitted. Their outer status is failed and bootstrap status is not_run when the scenario target is unavailable.
- Failure reasons are persisted in both `calibration_records.csv` and `attempt_audit.jsonl`.
- A separate `scenario_failure_records.csv` ledger records every scenario, its complete/failed status, and its failure message.
- The v2 manifest records planned/completed/failed scenario counts, planned/completed row counts, and scenario IDs included in summary tables.
- Scenario summaries are calculated only over complete scenarios; failed scenarios remain visible in raw records and the failure ledger.

A regression test injects a finite-parent target failure and checks that the run continues with explicit failed/not-run records.

### 4. Reject incomplete artifacts and recompute summaries from raw rows

In `scripts/audit_nested_validation_bootstrap_replay.py`:
- The auditor requires the v2 manifest and rejects partial/incomplete runs for a complete replay-audit pass.
- It checks unique and exact scenario/draw/validation-size/method keys against the planned design.
- It reconciles planned/completed method-row counts and scenario denominators.
- It validates the per-scenario failure ledger and its consistency with manifest counts.
- It recomputes both summary CSVs from raw calibration records using the same summary implementation, then compares columns, rows, and values. This is a tamper/reconciliation check, not an independent second implementation of the summary statistics.

A regression test confirms that a partial manifest is rejected before expensive reconstruction.

## Important limits and remaining work

1. **No tests were executed in this review environment.** New and existing tests are authored but remain unverified until CI or a local test run.
2. The replay auditor still reruns the production wrapper to reproduce wrapper-level failures. This remains deterministic failure reproduction, not independent reconstruction.
3. The scenario-level failure schema currently handles failures during finite-parent target definition. Unexpected exceptions during parent/GW generation or before attempt records are emitted can still abort the run before a complete artifact set is written.
4. These software repairs do not establish bootstrap calibration, inferential validity, estimator superiority, real-survey validity, or detector realism.

## Verification state

- Changes were committed to a separate branch based directly on PR #17's inspected head.
- The branch is maintained separately from PR #17; its live ahead/behind count should be taken from GitHub rather than this document.
- No PR was opened, no existing PR was modified, and nothing was merged.
- No tests, simulations, or artifact replays were executed.
- A small CI run is required before treating these repairs as verified.
