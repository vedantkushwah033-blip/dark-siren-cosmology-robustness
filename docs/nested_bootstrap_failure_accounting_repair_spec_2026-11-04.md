# Nested-bootstrap failure-accounting repair specification

**Date:** 2026-11-04  
**Target branch:** `research/bootstrap-failure-path-repairs-2026-11-02`  
**Implementation commits:** `7cd508b`, `04133e1`, `127c132`, `81ff9b9`, `c150039`, `a270777`, `cab8b36`, `d024026`, `beed9f0`, `491c1fd`, `a2e37bb`, `675239b`, `8805412`, `0e42f13`, `e065db8`, `58d8aa6`  
**Status:** Runner, summary, replay-auditor, and focused-test changes committed. Tests have not been executed in this environment; release gate remains closed.

## Scope

This document records a source-level defect found in the nested outer-SRSWOR / inner-bootstrap runner and defines the accounting contract required before another confirmatory experiment. It does not claim that the repair has been implemented, tested, or scientifically validated.

The benchmark remains synthetic and exploratory. The ordinary paired bootstrap, the FPC-scaled diagnostic, and the multinomial-pseudo-population/SRSWOR construction must not be described as theoretically validated uncertainty procedures for the nonlinear logistic-fit → inverse-weighting → host-density → posterior-mean pipeline.

## Source-level findings

The runner currently calls `bootstrap_validation_h0_posterior`, reads its `n_successful` value, and raises a `RuntimeError` when fewer than two inner resamples succeeded. The surrounding exception handler then writes a generic failed record with zero successes and a replacement failure record. This discards information already returned by the bootstrap function, including the actual success mask, per-attempt posterior values, failure records, and actual attempt counts.

The same exception handler currently assigns `n_bootstrap_attempted = bootstrap_replicates` even when the wrapper raises before returning a result. In that case the runner has no returned attempt ledger from which to verify that the requested number of inner resamples was actually reached. The requested count must not be presented as an observed count.

The calibration summarizer currently derives `n_bootstrap_attempted` from the number of successful outer fits. That is a count of bootstrap procedure opportunities, not the number of individual inner resampling attempts. These two denominators are distinct and must be named and reported separately.

The `paired_bootstrap_fpc` row is a deterministic transformation of the paired-bootstrap result. It is a derived diagnostic arm and does not represent a second, independently executed set of inner resamples.

## Required status model

Use these statuses for each method-specific outer-draw procedure:

- `SUCCESS`: the wrapper returned an attempt ledger and at least two inner resamples succeeded, allowing the requested resampling summary to be calculated.
- `NON_EVALUABLE`: the wrapper returned an attempt ledger, but fewer than two inner resamples succeeded. Preserve every returned attempt and failure record; do not calculate an SE or percentile interval.
- `FAILED`: the wrapper raised before returning a usable result. Preserve exception class/message and classify the number of inner attempts as known only when independently established. Do not substitute the requested replicate count for an observed count.
- `NOT_RUN`: the outer validation fit failed, so the inner bootstrap procedure was not invoked.

Status values may be serialized in lowercase for consistency with the existing tables, but the meanings above are normative.

## Required count fields

Keep procedure-level and inner-attempt-level quantities separate.

1. `n_bootstrap_procedures_attempted`: number of outer-draw opportunities on which a given independently executed bootstrap method was invoked.
2. `n_bootstrap_procedures_successful`: number of those procedures with status `success`.
3. `n_bootstrap_procedures_non_evaluable`: number with status `non_evaluable`.
4. `n_bootstrap_procedures_failed`: number with status `failed`.
5. `n_bootstrap_procedures_not_run`: number with status `not_run`.
6. `n_inner_resamples_attempted`: sum of actual per-procedure inner attempts, not the requested count.
7. `n_inner_resamples_successful`: sum of true entries in the returned success masks.
8. `n_inner_resamples_failed`: sum of false entries in the returned success masks, reconciled against the failure ledger.

For a wrapper failure with no returned attempt ledger, an inner-attempt count must be represented as unknown/null unless the code path proves that zero attempts occurred. A null/unknown count is not zero.

The derived `paired_bootstrap_fpc` arm must carry a reference to the paired-bootstrap source procedure and must not increment independently executed procedure or inner-resample totals. If it remains in method-level tables, summaries must explicitly label it as derived.

## Required rates and denominators

Report at least these distinct quantities:

- **Procedure non-success fraction:** (non-evaluable procedures + failed procedures) / procedures attempted. Exclude `not_run` procedures from this denominator and report them separately.
- **Inner-resample failure fraction:** failed inner resamples / observed inner resamples attempted. Calculate only from returned attempt ledgers with known denominators.
- **Wrapper-failure fraction:** wrapper-failed procedures / procedures attempted.
- **Non-evaluable fraction:** non-evaluable procedures / procedures attempted.

Do not call the procedure non-success fraction an inner-bootstrap replicate failure rate.

## Required persistence and invariants

For each independently executed method/outer-draw key, persist:

- procedure status;
- requested inner-resample count;
- actual inner-resample count, or explicit unknown;
- actual success and failure counts;
- complete success mask when returned;
- complete per-attempt failure records;
- per-attempt random indices and diagnostics when returned;
- wrapper exception class and message when applicable;
- input fingerprints and seed when available;
- explicit provenance linking the FPC-derived row to its paired-bootstrap source.

Validate the following invariants:

- `n_inner_resamples_attempted == len(success_mask)` whenever a success mask exists.
- `n_inner_resamples_successful == sum(success_mask)`.
- `n_inner_resamples_failed == len(success_mask) - sum(success_mask)`.
- The failure ledger indexes exactly the false entries in the success mask, with no missing or duplicate attempt identifiers.
- `success` requires at least two successful inner resamples and finite summary metrics.
- `non_evaluable` requires a returned ledger and fewer than two successes; summary metrics must be null/NaN.
- `failed` requires wrapper exception metadata; no fabricated inner attempts or failure records.
- `not_run` requires a failed outer fit and zero bootstrap procedure attempts.
- The three status counts sum to the procedure-attempt denominator, with `not_run` tracked separately.
- Derived FPC rows must not inflate independent procedure or inner-resample denominators.

## Required tests before execution

1. All inner resamples succeed.
2. Mixed success/failure with at least two successes.
3. Exactly one successful resample: `non_evaluable`, preserving the one success and all failures.
4. Zero successful resamples with a returned ledger: `non_evaluable`, preserving all attempts.
5. Wrapper raises before returning: `failed`, preserving exception metadata and leaving unobserved inner-attempt count unknown rather than assigning the requested count.
6. Outer fit fails: `not_run`, with no bootstrap invocation.
7. Tamper independently with each count, success mask, failure ledger, and procedure status; the replay auditor must reject each inconsistency.
8. Confirm that paired-bootstrap FPC-derived rows do not double-count inner resamples.
9. Confirm summary denominators by hand on a small fixture containing all four statuses.
10. Round-trip the artifact files through serialization and independently reconstruct the same counts and rates.

## Release gate

Do not launch the large confirmatory experiment until the runner, calibration summarizer, artifact schema, and independent replay auditor implement this contract; the focused tests pass on the supported Python/dependency matrix; exact-head CI is available; and a small persisted smoke artifact passes independent audit.

Even after this gate, the study remains a synthetic robustness experiment. A passing software audit does not establish bootstrap validity, real-survey transportability, realistic detector selection, or a new observational H0 result.


## Implementation progress (2026-11-04)

The implementation now records returned bootstrap outcomes even when fewer than two resamples succeed. Such procedures receive `non_evaluable`, retain the returned success mask, attempt-level posterior values, failure ledger, and actual counts, and do not receive an SE or percentile interval.

Wrapper exceptions now retain exception class/message but leave inner attempt, success, and failure counts unknown (`null`) and do not fabricate attempt-level failure records. The FPC diagnostic row is explicitly marked as derived from `paired_bootstrap`; calibration denominators exclude it.

The calibration summary and auditor-owned independent reconstruction now use matching procedure/inner-resample denominators. A parity regression test compares both reconstructed tables under a mixed success, non-evaluable, wrapper-failure, not-run, and derived-FPC fixture. Count-ledger validation rejects incomplete counts, status/count disagreement, fabricated wrapper-failure counts, and nonzero counts for not-run procedures. The independent replay auditor also checks the zero-failure invariant for not-run rows.

The calibration summary now also reports the number of procedures whose inner-attempt count is unknown after a wrapper exception. It distinguishes:
- bootstrap procedure attempts, successes, non-evaluable procedures, failures, and not-run procedures;
- observed inner resample attempts, successes, and failures;
- procedure-level non-success fraction;
- inner-resample failure fraction;
- wrapper-failure fraction.

The replay auditor now checks count/mask/failure-ledger reconciliation, rejects fabricated counts on wrapper failure, checks the derived FPC provenance label, and deterministically replays non-evaluable ledgers. Its report separates successful and non-evaluable replay totals.

Focused tests were added for partial outcomes, unknown wrapper-failure counts, denominator separation, FPC exclusion, independent reconstruction parity, status/count disagreement, and inconsistent count rejection. The CI workflow was updated to include these tests and the nested-runner regression suite.

**Verification limitation:** These commits were written through the GitHub contents API. No local test execution or exact-head CI run has yet been obtained. Therefore the changes are implemented but not yet verified as passing. Do not describe the repair as complete or launch the confirmatory experiment.

## Additional source-level correction (2026-11-04)

A follow-up review found that the runner's exception handler covered both the bootstrap wrapper call and validation of the result it returned. Consequently, an internally inconsistent returned ledger could be caught and mislabeled as a wrapper failure, with its known attempt counts and returned evidence discarded as unknown. This conflated two distinct cases:

- Wrapper exception before a result is returned: inner-attempt counts may be unknown.
- Malformed result returned by the wrapper: the returned ledger exists but fails its integrity checks; the run must stop rather than relabel the result as an unobserved wrapper failure.

The runner now catches exceptions only around the wrapper invocation. Count/mask reconciliation occurs after that exception boundary, so an invalid returned ledger raises an error and cannot silently enter the calibration table as a generic wrapper failure. A regression test supplies a deliberately inconsistent returned ledger and requires the runner to reject it.

Implementation commits for this correction:
- 79ea5e3 — initial targeted runner edit; superseded by the complete block correction below.
- f127af1 — isolate wrapper exception handling from returned-ledger validation.
- c113b22 — add the malformed-returned-ledger regression test.

Verification status: The updated source and regression test were fetched back from the repair branch and inspected. No local test execution was performed. The GitHub connector returned no workflow run or status-check records for commit c113b22; exact-head CI therefore remains unverified. The release gate remains closed.

## Independent replay boundary correction (2026-11-04)

A further audit found that `scripts/audit_nested_validation_bootstrap_replay.py` called the production `bootstrap_validation_h0_posterior` wrapper to reproduce wrapper failures and non-evaluable outcomes. That made those branches same-implementation replays, despite successful-attempt posterior calculations being reconstructed by the separate auditor implementation.

The replay auditor now:

- Contains no import or call to the production bootstrap wrapper.
- Reconstructs successful and non-evaluable returned ledgers through `audit_posterior_validation_bootstrap`, which independently regenerates the resampling indices and recomputes the posterior attempts.
- For wrapper-level failures, where the producer returned no attempt ledger and actual inner-attempt counts are unknown, checks only that exception class/message metadata are present. It does **not** claim to independently reproduce those failures.
- Reports wrapper-level failures as `wrapper_level_bootstrap_failures_metadata_checked_not_replayed`, not as independently replayed failures.

A static regression test prevents the replay auditor from importing or calling the production wrapper. The existing fabricated-wrapper-failure test continues to require rejection when failure metadata is absent.

**Verification status:** Source and test edits were fetched back from the repair branch and inspected. Tests have not been executed locally, and exact-head GitHub Actions status is not yet available. The release gate remains closed. The wrapper-failure branch remains a documented audit limitation until an independently implemented preflight/index-generation replay is added; no returned attempt ledger exists to reconstruct for that path.

A final replay-integrity refinement also ensures `_reconstruct_result` passes the persisted `n_bootstrap_failed` count to the independent posterior-bootstrap auditor. It no longer derives that value as `attempted - successful`, which could otherwise conceal a corrupted stored failure count. A regression guard checks that the persisted field is used directly. This source/test change is committed as `253a126` and `4adf11d`; exact-head execution remains pending.
