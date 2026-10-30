# Internal bootstrap-attempt audit gap and remediation protocol

**Status:** required follow-up before the linear benchmark is treated as a fully attempt-audited implementation.  
**Observed on:** source commit `74e6dacd6690875fdeb6c6d6fc0a1fcaa40fde65`.  
**Scope:** studentized pseudo-population SRSWOR method in the linear finite-population mean benchmark.

## 1. Why this matters

The bootstrap primitive records an outcome for every inner replicate: its success mask and stable failure code, along with replay arrays. The outer benchmark currently stores the resulting interval or an interval-level failure, but does not serialize the inner attempt, success, and failure counts or the inner failure-code histogram into its replicate-level CSV. The independent benchmark auditor independently reconstructs the interval endpoints, but its helper discards the inner failure accounting.

Consequently, the current pilot's independent audit establishes reconstruction of the outer sample/method attempt grid and interval results. It does not establish that the artifact itself preserves an independently checkable accounting of every internal bootstrap attempt. This is an auditability limitation, even if the observed pilot happened to have no inner failures.

## 2. Required replicate-level schema extension

For `studentized_pseudopop_srswor`, record:

- `inner_bootstrap_attempted` (must equal the manifest's bootstrap replicate count);
- `inner_bootstrap_successful`;
- `inner_bootstrap_failed`;
- `inner_bootstrap_failure_counts_json` (canonical JSON mapping stable failure code to count);
- optionally, `inner_bootstrap_success_mask_sha256` and `inner_bootstrap_failure_codes_sha256` for compact integrity checks, while retaining replayable deterministic seeds.

For Wald and ordinary percentile rows, use explicit not-applicable values, not zero: those methods do not run this inner studentized bootstrap.

For every studentized row, enforce:
`attempted = successful + failed`, and
`sum(failure-code counts) = failed`.
The number of successes must be sufficient for interval construction whenever the row is marked successful.

## 3. Summary-level accounting

For each scenario × sample-size × method cell, include:
- total outer attempts;
- successful and failed outer intervals;
- total inner bootstrap attempts;
- total inner bootstrap successes and failures;
- inner failure fraction with its denominator;
- per-code failure totals.

Do not combine inner-bootstrap failures with outer interval-construction failures. Report both layers separately. A row that fails interval construction after some inner replicates succeeded still contributes all of its inner attempts to the inner accounting.

## 4. Independent reconstruction requirements

The auditor must independently replay the studentized bootstrap using the recorded deterministic seed and reconstruct:

1. all multinomial pseudo-population count vectors;
2. all selected SRSWOR positions;
3. each bootstrap replicate's variance and studentized statistic;
4. the success mask and stable failure code for every inner attempt;
5. the exact attempt/success/failure totals and per-code histogram;
6. the final quantiles and interval endpoints.

Compare the reconstructed values with the serialized row-level accounting and reject any mismatch. Recompute summary totals from replicate rows and reject summary discrepancies.

## 5. Adversarial tests

Add tests that:
- modify inner failure counts and refresh the CSV checksum;
- alter a failure-code histogram while preserving its total;
- alter the inner attempt count;
- change a success count but preserve attempted = success + failure;
- alter summary inner-failure totals;
- create a controlled degenerate bootstrap draw to exercise a genuine zero/nonfinite-variance failure path.

Each tampering test must be rejected by the independent auditor even after the relevant file checksum in the manifest is refreshed. The controlled failure test must verify that the failure is retained and counted rather than silently dropped.

## 6. Versioning and compatibility

Bump the benchmark artifact schema to `linear_finite_population_bootstrap_benchmark_v2`. Update runner, auditor, tests, documentation, and workflow together. Do not silently reinterpret v1 artifacts as v2. Preserve the existing v1 pilot unchanged as a historical pilot artifact.

## 7. Release gate

The v2 benchmark is ready for a new pilot only when:
- all unit and adversarial tests pass;
- the exact-head Python matrix passes;
- the v2 pilot completes;
- the independent auditor reconstructs all outer and inner attempts;
- raw artifacts, manifest, audit report, and checksums are uploaded;
- the manuscript describes v1 as an audited outer-interval pilot with this inner-attempt serialization limitation, and v2 as the first attempt-level audited release.

Passing this linear release gate still does not validate the nonlinear dark-siren posterior functional. That remains a separate scientific validation problem.
