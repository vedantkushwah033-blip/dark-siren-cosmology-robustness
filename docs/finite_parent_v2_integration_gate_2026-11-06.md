# V2 integration gate for the finite-parent SRSWOR benchmark

**Status:** opt-in V2 dispatch and a separate SLSQP audit path have been implemented on the research branch; no V2 finite-parent result is yet claimed as validated because tests have not been observed running.

## Repository inspection

The finite-parent runner now accepts an opt-in `solver_version="v1"|"v2"` argument and passes it to both the full-parent finite-frame reference and each repeated validation-sample fit. The default remains V1. The CLI exposes `--solver-version`; the output manifest uses schema `finite_parent_srswor_smoke_v3` and records the selected version, while every record carries the same version.

The auditor now requires the record and manifest versions, solver revision and model family to agree. For V2 it independently solves the bounded linear-logistic likelihood using SciPy SLSQP, then reconstructs the fitted probabilities, inverse weights, reflected KDE and downstream posterior. This is algorithmically separate from the production damped-Newton iterations; it is not a high-precision oracle. The benchmark currently fixes the inclusion model to linear, so the model family is not a selectable dimension in this experiment. V2 records carry revision `2.1-armijo-noise-slack`; V1 records carry the explicit legacy label `v1-lbfgsb-kkt`.

The separate estimated-inclusion pilot already supports `--solver-version v2` and its CI workflow runs a small V2 smoke experiment. That does not replace the finite-parent repeated-SRSWOR benchmark: the pilot and finite-parent benchmark have different sampling estimands and artifact schemas.

## Required implementation contract

1. Implemented: explicit `solver_version` argument in the finite-parent runner; V1 remains default and only V1/V2 are accepted.
2. Implemented: selected version is passed to both the per-draw and full-parent estimated-density calls.
3. Partially implemented: manifest and records carry `solver_version`; solver revision and model-family fields are not yet persisted per record.
4. Existing scenario seeds and nested validation indices are unchanged by solver selection. A paired V1/V2 comparison artifact and cross-version pairing auditor are still outstanding.
5. Partially implemented: V2 audit uses an independently coded likelihood/gradient with SLSQP and an independent box-KKT check. A high-precision oracle and a full dedicated V2 stress suite remain outstanding.
6. Existing failure retention is preserved; exact V2 failure replay still needs runtime confirmation.
7. Existing finite-frame, generating-truth and oracle references remain distinct.
8. Added regression coverage for V2 dispatch, solver-version validation, an independently audited V2 smoke fixture, rejection of altered solver revision/model-family metadata, and paired V1/V2 comparisons. CI is configured to run V1 and V2 finite-parent smoke experiments, independently audit each arm, and execute the paired comparison regression. Runtime results remain pending observation.

## Numerical gate before benchmark interpretation

The V2 solver is not ready to be treated as independently established solely from V1 agreement or the current float64 SLSQP tests. The historical D-01 evidence package should be inspected and reproduced where available. At minimum, independently verify ordinary interior optima, lower- and upper-bound optima, near separation, severe but full-rank conditioning, and line-search termination. Record objective gaps and projected-KKT residuals using an implementation that does not call the production solver's private objective or derivative functions.

The eight-ULP Armijo allowance is a line-search acceptance tolerance, not a convergence certificate. The final projected-KKT check remains necessary, and boundary-limited fits must remain distinguishable from interior convergence.

## Statistical interpretation gate

The experiment estimates repeated-validation behavior conditional on each fixed simulated parent frame, catalogue, and GW-distance sample. Validation draws within a scenario are nested and paired; they are not independent top-level scenarios. Scenario-balanced summaries and uncertainty calculations must respect that hierarchy.

A solver-reliability result is not a calibration result. Even if V2 reproduces the constrained MLE reliably, this alone does not establish unbiased host-density reconstruction, calibrated H0 intervals, robustness to inclusion-model misspecification, or superiority over V1.

## Current conclusion

V2 is now an opt-in finite-parent benchmark path with an independent SLSQP reconstruction path in the auditor, explicit solver provenance, and a paired V1/V2 smoke-comparison runner. This closes the implementation and pairing-code gaps, but it does **not** establish validated V2 results: exact-head CI has not been observed running, no paired comparison artifact has yet been produced and audited in a live run, and the high-precision D-01 oracle evidence remains unverified. Until those gates close, label any V2 output **SMOKE ONLY — AUDIT/EXECUTION PENDING**.
