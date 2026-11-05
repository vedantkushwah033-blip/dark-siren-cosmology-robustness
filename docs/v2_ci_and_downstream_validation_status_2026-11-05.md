# V2 cross-version CI and downstream validation status

**Audit date:** 2026-11-05  
**Research branch:** `research/logistic-inclusion-v2-newton-2026-11-05`  
**Current branch head observed:** `8380001095a879dfde63e1392868f0434e586d9c`  
**Full workflow run:** `42027804117`  
**Code commit tested by that run:** `46e96b9d0d90b8b4c6bfd1e6411a27ae07e68caf`

## 1. Verified CI result

The full `numerical-diagnostics` matrix completed successfully on all three configured Python versions:

| Python | Job ID | Conclusion |
|---|---:|---|
| 3.10 | 127246328200 | success |
| 3.11 | 127246327686 | success |
| 3.12 | 127246328335 | success |

The following relevant steps completed successfully in all three jobs:
- Versioned damped-Newton logistic inclusion estimator V2 tests.
- Estimated-inclusion pilot runner tests.
- Persisted validation-bootstrap artifact generation and independent refit audit tests.
- Integrated nested validation-bootstrap comparison smoke experiment.
- Independent audit of nested validation-bootstrap smoke artifacts.

The Python 3.12 job also completed the persisted validation-bootstrap sensitivity smoke and audit successfully.

Artifacts were uploaded for all three Python jobs. This establishes successful execution of the configured workflow on the tested code commit. It does not establish that every optional Python-3.12-only workflow step ran on Python 3.10 or 3.11.

## 2. Commit-provenance boundary

The workflow tested code commit `46e96b9d0d90b8b4c6bfd1e6411a27ae07e68caf`. The observed research-branch head is `8380001095a879dfde63e1392868f0434e586d9c`, a later documentation update. Therefore the CI result supports the code at the tested commit; the later head must not be described as having its own separate workflow run unless a run for that exact SHA is found.

No merge, main-branch modification, PR closure, or release declaration was made.

## 3. V2 integration status

Source inspection confirms an opt-in `solver_version="v1"|"v2"` path in the validation-to-host-density wrapper, pilot function, and paired validation-bootstrap path. V1 remains the default. Solver identity is persisted in pilot records and inference audit payloads, and the independent pilot auditor accepts both versions and independently refits the same bounded logistic objective using L-BFGS-B.

The CI result verifies the relevant unit, integration, smoke, and artifact-audit steps listed above. It does **not** demonstrate a repeated-sampling H0 performance advantage for V2.

## 4. V1/V2 parity evidence already present

The branch contains same-data fitter-level tests for:
- Linear model under ordinary conditions.
- Quadratic model under ordinary conditions.
- Linear and quadratic complete separation.
- Quadratic near-separation with one discordant label.
- Rank-deficient quadratic design.
- Analytic gradient and Hessian finite-difference checks.
- Hessian rank and conditioning diagnostics.

These tests compare fitted probability curves and objective values, not coefficients alone where parameters may be non-identifiable.

A remaining distinction is important: fitter-level parity tests are not the same as a paired end-to-end H0 comparison across independently generated replicates. The latter must explicitly verify that V1 and V2 receive identical parent populations, catalogue retention, validation records, and GW distances.

## 5. D-01 status

The historical D-01 root cause remains **OPEN / UNKNOWN**. The inspected branch documentation does not contain the original raw validation records or optimizer trace needed to classify the historical event as separation, rank deficiency, boundary limitation, numerical line-search failure, or environment failure.

The successful V2 tests and current CI do not retrospectively identify the cause of the historical V1 failures. Do not close D-01 as a resolved statistical defect without recovering and replaying the original records. If the records cannot be recovered, close only the replayability investigation as unavailable and preserve the root cause as unknown.

## 6. Remaining research gates

1. Recover the original D-01 inputs, optimizer status, objective/gradient trace, package versions, and environment if available.
2. Run a paired V1/V2 end-to-end comparison with identical synthetic inputs and persisted artifacts; report failures in both arms.
3. Compare host-density error, inclusion-probability prediction error, weight concentration, H0 posterior means, interval widths, and coverage over repeated sampling.
4. Calibrate any weight-stability warning thresholds using simulation; do not treat ESS or a fitted-probability floor as a universal validity certificate.
5. Freeze estimands, truth strata, prior/domain policy, failure denominators, and Monte Carlo precision before a confirmatory campaign.

## 7. Evidence classification

- **VERIFIED:** Workflow run 42027804117 completed successfully for Python 3.10, 3.11, and 3.12 on commit 46e96b9d0d90b8b4c6bfd1e6411a27ae07e68caf; listed V2 integration and artifact-audit steps passed.
- **VERIFIED:** The observed research-branch head was 8380001095a879dfde63e1392868f0434e586d9c.
- **VERIFIED:** V2 is opt-in in the inspected source; V1 remains the default.
- **UNRESOLVED:** Historical D-01 root cause.
- **NOT YET ESTABLISHED:** V2 improves host reconstruction or H0 posterior calibration relative to V1.
- **NOT AUTHORIZED:** Confirmatory simulation campaign or release declaration.

This is a CI and source-state update, not a new scientific result or a publication-readiness declaration.
