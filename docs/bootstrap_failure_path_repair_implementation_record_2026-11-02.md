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

## Important limits and remaining work

1. **No tests were executed in this review environment.** The new tests are authored but unverified until CI runs on this branch or a PR.
2. The scenario-level finite-parent target checks in `run_comparison` still raise before the runner writes scenario/draw records. Planned-versus-completed scenario accounting and a formal partial-run artifact schema remain unresolved.
3. The replay auditor still reruns the production wrapper to reproduce wrapper-level failures. The newly persisted exception metadata fixes the serialization mismatch, but this remains deterministic failure reproduction, not independent reconstruction.
4. The artifact auditor still needs stronger manifest/summary reconciliation and explicit rejection of incomplete or duplicated scenario/draw/method keys.
5. These software repairs do not establish bootstrap calibration, inferential validity, estimator superiority, real-survey validity, or detector realism.

## Verification state

- Changes were committed to a separate branch based directly on PR #17's inspected head.
- The branch is two commits ahead of the PR #17 base at the point of this record; later commits may update that count.
- No PR was opened, no existing PR was modified, and nothing was merged.
- No tests, simulations, or artifact replays were executed.
- A small CI run is required before treating these repairs as verified.
