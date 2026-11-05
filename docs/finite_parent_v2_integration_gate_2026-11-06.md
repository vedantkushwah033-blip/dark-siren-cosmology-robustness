# V2 integration gate for the finite-parent SRSWOR benchmark

**Status:** design/audit note; no V2 finite-parent result is claimed.

## Repository inspection

The finite-parent runner in `src/dark_siren_robustness/finite_parent_sampling.py` currently calls `estimate_host_density_from_validation` without a `solver_version` argument. The wrapper in `src/dark_siren_robustness/estimated_inclusion.py` already accepts `solver_version="v1"|"v2"` and dispatches to the corresponding linear/quadratic fitters and probability predictor. Therefore V2 integration is a small production-code change, but not yet an independently audited experiment.

The existing finite-parent auditor, `scripts/audit_finite_parent_validation_benchmark.py`, independently reconstructs the V1 fitting path. It must not be treated as a V2 auditor merely because the production wrapper can select V2.

The separate estimated-inclusion pilot already supports `--solver-version v2` and its CI workflow runs a small V2 smoke experiment. That does not replace the finite-parent repeated-SRSWOR benchmark: the pilot and finite-parent benchmark have different sampling estimands and artifact schemas.

## Required implementation contract

1. Add an explicit `solver_version` argument to the finite-parent runner and CLI. Preserve `v1` as the default for backwards compatibility; accept only `v1` and `v2`.
2. Pass the selected version to both the per-draw estimated-density call and the full-parent finite-frame target call. The two fits must use the same declared model type and solver version.
3. Record `solver_version`, `solver_revision`, and the fitted model family in every success and failure record and in the summary manifest. Do not infer version from a branch name.
4. Preserve the same simulated parent frame, catalogue thinning, GW distances, random seeds, nested validation indices, and overlap identities when comparing V1 and V2. The comparison must be paired at the scenario/draw/validation-size level.
5. Extend the independent auditor with a separately implemented V2 objective, gradient, box-KKT calculation, and probability reconstruction. Do not import production V2 fitting functions into the auditor. Recompute the full-parent target and every successful row independently.
6. Retain all failures. A V2 optimizer failure must not cause the corresponding V1 row or scenario to be dropped.
7. Keep the finite-frame target distinct from the generating H0 truth and from the oracle-propensity diagnostic. Report errors against each target separately.
8. Add deterministic tests for V1 backward compatibility, V2 dispatch, paired input identity, failure retention, and auditor rejection of altered solver metadata or fit outputs.

## Numerical gate before benchmark interpretation

The V2 solver is not ready to be treated as independently established solely from V1 agreement or the current float64 SLSQP tests. The historical D-01 evidence package should be inspected and reproduced where available. At minimum, independently verify ordinary interior optima, lower- and upper-bound optima, near separation, severe but full-rank conditioning, and line-search termination. Record objective gaps and projected-KKT residuals using an implementation that does not call the production solver's private objective or derivative functions.

The eight-ULP Armijo allowance is a line-search acceptance tolerance, not a convergence certificate. The final projected-KKT check remains necessary, and boundary-limited fits must remain distinguishable from interior convergence.

## Statistical interpretation gate

The experiment estimates repeated-validation behavior conditional on each fixed simulated parent frame, catalogue, and GW-distance sample. Validation draws within a scenario are nested and paired; they are not independent top-level scenarios. Scenario-balanced summaries and uncertainty calculations must respect that hierarchy.

A solver-reliability result is not a calibration result. Even if V2 reproduces the constrained MLE reliably, this alone does not establish unbiased host-density reconstruction, calibrated H0 intervals, robustness to inclusion-model misspecification, or superiority over V1.

## Current conclusion

V2 integration into the production finite-parent runner is technically supported by an existing wrapper argument, but the finite-parent benchmark and its independent auditor have not yet been updated as a validated V2 experiment. Until the auditor is independently extended and exact-head tests are observed, report V2 finite-parent results as **not yet audited**.
