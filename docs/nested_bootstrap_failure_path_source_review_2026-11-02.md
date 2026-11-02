# Follow-up source review: nested-bootstrap failure-path integrity

**Reviewed source reference:** PR #17 head `7c612efb1b075ea2df0e665b93091fa78c77fe84`  
**Review branch:** `research/hostile-referee-audit-response-2026-11-02`  
**Review date:** 2026-11-02  
**Scope:** targeted static inspection of the nested comparison runner, production bootstrap wrapper, and replay auditor. No simulation or test suite was executed.

## Executive finding

The prior review identified scenario-level aborts and exceptional-path index alignment as risks. A more specific inspection exposes an additional, deterministic incompatibility between the runner's wrapper-failure serialization and the replay auditor: the runner does not persist the wrapper exception class/message, but the auditor requires those exact fields to independently reproduce and compare a wrapper-level failure. Therefore, any row produced by the runner's current outer `except` branch is not auditable through the current wrapper-failure branch as written.

This is a source-level finding, not a claim that a particular historical artifact contains such a failed row.

## 1. Wrapper-failure metadata is dropped

In `scripts/run_nested_validation_bootstrap_comparison.py::run_comparison`, the exception handler around `bootstrap_validation_h0_posterior(...)` creates a failed method result containing:
- `bootstrap_status = "failed"`;
- NaN interval/SE values;
- zero successful attempts and the requested attempted count;
- a generic `failure_records` entry containing only a formatted reason;
- the FPC factor.

It does **not** set `wrapper_failure_class` or `wrapper_failure_message`.

Later, when the runner builds `attempt_log`, it serializes those values using `result.get("wrapper_failure_class")` and `result.get("wrapper_failure_message")`. Both are therefore null for the wrapper-exception path.

In `scripts/audit_nested_validation_bootstrap_replay.py::audit`, a row whose `bootstrap_status` is `failed` enters the wrapper replay branch. The auditor converts the stored class/message to strings and compares them with the exception class/message from rerunning the production wrapper. With the current runner output, the stored class is effectively `"nan"`/null-like rather than the actual exception class, so the equality check cannot pass for a genuine wrapper exception.

**Required repair:** persist the original `type(exc).__name__` and `str(exc)` in dedicated wrapper-level fields, distinct from inner-attempt `failure_records`. Add a regression test that forces a wrapper-level failure, writes the artifact, and verifies that the replay auditor can classify the record. Do not infer wrapper failure class from a generic free-text reason.

## 2. Design-aware index-generation exceptions can break attempt alignment

In `src/dark_siren_robustness/posterior_validation_bootstrap.py`, the paired-bootstrap branch appends an empty `sample_indices` placeholder if an exception occurs before sample indices have been appended. The design-aware branch appends `pseudo_population_indices`, `sample_indices`, and `selected_pseudo_population_positions` only after `design_aware_srswor_bootstrap_indices(...)` returns; it has no corresponding placeholder logic if that call itself raises.

The independent auditor requires the design-aware index arrays to have one entry per attempted replicate. An exception at the index-generation boundary can therefore produce a result whose attempt-level status arrays are complete but whose saved index lists are shorter than `n_attempted`, making the artifact fail audit and losing exact attempt alignment.

**Required repair:** preserve one explicit index-record slot per attempt for both methods. Use a schema that distinguishes “index generation failed” from a valid empty index vector; do not use an empty vector as an ambiguous sentinel without an explicit status. Add injected-failure tests for both index generators and verify all attempt-aligned collections retain the planned length.

## 3. Artifact completeness checks should be made explicit

The replay auditor reconstructs expected scenario/draw/size coordinates and checks that each expected draw has three method rows. It also performs a one-to-one merge of the attempt trace with calibration status. These are useful checks, but the audit entry point does not define a formal partial-run schema with planned/completed counts or scenario-level failure rows. It also does not independently verify the scenario-summary and scenario-balanced-summary CSVs against the raw records.

Before any long run, the artifact contract should explicitly state:
- planned scenario count;
- completed and failed scenario counts;
- planned outer-draw × validation-size coordinates;
- planned method rows and completed method rows;
- planned and completed inner attempts per method;
- a canonical list of scenario/draw/size/method keys;
- whether partial artifacts are permitted and how they are labelled.

The auditor should reject missing, duplicated, and unexpected keys; reconcile every denominator; and independently recompute summary tables from raw rows. A run interrupted before all scenarios finish must be represented as partial/incomplete, never as a complete benchmark.

## 4. Recommended implementation order

1. Repair wrapper exception metadata and add a forced wrapper-failure serialization/replay test.
2. Repair design-aware attempt-index alignment and add injected index-generation failure tests.
3. Introduce explicit scenario-level failure records and planned/completed counts in a versioned artifact schema.
4. Harden the auditor against missing/duplicate/unexpected keys and independently reconstruct summary CSVs.
5. Only then run a small adversarial pilot and audit it end-to-end. Do not launch a confirmatory simulation before the pilot passes.

## 5. Current execution evidence

At the time of this review:
- PR #17 is open, draft, unmerged, and its inspected head is `7c612efb1b075ea2df0e665b93091fa78c77fe84`.
- A PR-triggered GitHub Actions workflow run associated with that exact head completed successfully (run `41439822334`). This is evidence for the CI workflow at that commit, not evidence that the newly identified failure paths are covered.
- PR #15 remains open, draft, unmerged, and not mergeable at inspected head `4f700fcf120ce547de598e6199d3829bc88c8498`; no workflow run was returned for that exact head by the queried workflow-runs endpoint.
- No local tests, simulations, or artifact replay were executed during this review.
- No source code, PR, branch history, or simulation artifact was changed by this review.

## Scientific interpretation

These are software/artifact-integrity blockers. Resolving them would make failure accounting and auditability more reliable; it would not establish that either bootstrap procedure is calibrated for the nonlinear posterior-mean functional. The paired bootstrap, pseudo-population comparator, and scalar FPC adjustment must remain exploratory until repeated-sampling performance is measured against the finite-parent target with scenario-level uncertainty.
