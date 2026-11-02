# Independent audit boundary review — 2026-11-02

## Review target

Branch: `research/bootstrap-failure-path-repairs-2026-11-02`

This note records a code-inspection review of the nested validation-bootstrap runner and replay auditor. It is not a test result, a simulation result, or a claim that the branch has passed CI.

## Findings

### 1. Attempt-level failure accounting is materially clearer

In `posterior_validation_bootstrap.py`, random-index generation occurs before the per-attempt `try` block. This is a defensible boundary: if the index-generation step fails, no complete resampling attempt has been defined, so the exception propagates to the comparison runner as a wrapper-level failure. Once indices have been generated and recorded, model/posterior failures are caught and stored with the attempt index, exception class, stable failure code, and message.

This avoids manufacturing an attempt record with missing resampling provenance. It does, however, mean the artifact schema must preserve two distinct failure levels:

- attempt-level failures: an index exists and the posterior attempt can be independently replayed;
- wrapper-level failures: index generation or another pre-attempt operation fails before a complete attempt trace exists.

The runner currently records wrapper failure class/message separately, which is the right distinction.

### 2. Scenario failure ledger is reconciled against row-level artifacts

The replay auditor checks that the scenario ledger has exactly the planned scenario identifiers, that complete/failed counts reconcile with the manifest, and that both calibration rows and attempt traces carry the same scenario failure annotation. It also requires rows from failed scenarios to have failed outer status and `not_run` bootstrap status.

These checks address the earlier class of defect where a scenario-level failure artifact could be absent or internally inconsistent. A successful audit still depends on the complete artifact set being present and on the manifest declaring the expected schema/status.

### 3. The posterior reconstruction is substantially independent, but not end-to-end independent

The posterior-level auditor implements its own logistic objective/fit, weighted reflected KDE, selection-aware likelihood, posterior normalization, quantiles, diagnostics, and RNG-index reconstruction. It does not simply call the producer's posterior estimator for successful bootstrap attempts.

There are two important shared-code boundaries in the enclosing replay auditor:

- The stored scenario and aggregate summary CSVs were initially rebuilt through the production `summarize_bootstrap_calibration` function. This gap has now been addressed in code: `independent_calibration_audit.py` implements an auditor-owned reconstruction, and the replay auditor calls it instead of the production summarizer. A unit test compares its output with the production reference on a controlled fixture. Exact-head CI and end-to-end adversarial artifact tests remain outstanding.
- For a recorded wrapper-level bootstrap failure, the auditor reruns `bootstrap_validation_h0_posterior` itself with the same inputs and seed, then checks that the exception class and message match. This is deterministic replay of the producer's wrapper path, not an independent reconstruction of the wrapper failure.

Consequently, the status string `INDEPENDENT_SCENARIO_AND_POSTERIOR_REPLAY_PASS` should be interpreted narrowly: independent reconstruction of the generated scenario, finite-parent target, outer fits, and successful inner posterior resamples, alongside producer-code replay for wrapper-level failures and production-code summary reconstruction. It should not be described as a wholly independent implementation of every artifact transformation.

### 4. The current audit does not establish scientific validity of the resampling method

A successful artifact replay establishes that recorded computations are reproducible and internally consistent under the declared software environment. It does not show that ordinary paired bootstrap, pseudo-population/SRSWOR resampling, or the diagnostic FPC adjustment gives valid uncertainty for the nonlinear logistic-fit → inverse-weighted KDE → selection-aware posterior-mean functional.

The design-aware branch remains an exploratory multinomial-pseudo-population/SRSWOR construction. The finite-population-mean theory motivating it does not automatically establish coverage for this nonlinear posterior functional. The FPC-adjusted standard error remains diagnostic-only.

## Required next improvements before a stronger independence claim

1. Completed in code: implement an auditor-owned summary reconstruction. Remaining: add end-to-end tests that alter stored scenario and aggregate CSVs and prove the full replay audit rejects them.
2. Define a wrapper-failure replay contract that independently reconstructs the pre-attempt failure condition from inputs and RNG state. If exact independent reconstruction is not practical for a given exception, label that class as producer-replayed rather than independently audited.
3. Add adversarial tests that mutate, one at a time:
   - scenario failure status/message;
   - a missing scenario-ledger row;
   - wrapper failure class/message;
   - attempt-level failure index/code/reason;
   - a successful posterior summary;
   - a summary-table aggregate;
   - seed or index provenance.
   Each mutation must cause the auditor to reject the artifact for the intended reason.
4. Add a positive control for each failure class, demonstrating that an unmodified artifact still passes.
5. Keep audit outcomes separated in machine-readable output, for example:
   - `independent_reconstruction`;
   - `producer_replay_only`;
   - `artifact_consistency_only`;
   - `not_checked`.
6. Run the full Python-version CI matrix on the exact branch head after these changes. A green run on a parent commit or a different PR head is not exact-head verification.

## Release interpretation

The auditor can now be described as a **mixed independent-reconstruction and deterministic-replay auditor with independent summary reconstruction**. Do not call it fully independent end-to-end validation. Do not promote pilot artifacts to confirmatory evidence based on a replay PASS alone.
