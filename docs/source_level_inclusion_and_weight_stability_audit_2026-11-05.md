# Source-level inclusion-estimator and host-reconstruction audit

**Date:** 5 November 2026  
**Branch:** `research/logistic-inclusion-v2-newton-2026-11-05`  
**Status:** source inspection and a new regression test committed; test execution pending.

## Executive findings

1. **D-01 remains historically unresolved.** A repository commit search for `D-01` returned no matches. The historical input rows, traceback, optimizer result object, and environment record were not recovered. Prospective tests cannot retrospectively diagnose the old failure.
2. **V1/V2 parity was incomplete.** A quadratic same-data parity test existed; a linear equivalent did not. The missing linear test was added in commit `ad8dbff38130bfe105589171812380a5352a3181`. It has not yet been executed.
3. **The coefficient boxes permit extremely small fitted probabilities.** On (x\in[-1,1]), the linear model can reach logit -75 and inverse probability about (3.73\times10^{32}); the quadratic model can reach -125 and inverse probability about (1.94\times10^{54}). These are worst-case feasible-box calculations, not observed pilot values.
4. **The KDE rescales weights by their maximum**, protecting arithmetic from a common large scale. This does not prevent a few catalogue objects dominating the normalized weighted empirical distribution.
5. **No hidden-host mixture occurs in the inspected estimated-inclusion pilot path.** It builds one IPW catalogue density, applies the merger-rate multiplier, normalizes, and passes that density to the fixed-count selection-aware likelihood.
6. **V2 remains opt-in and is not wired into the pilot or its diagnostics.** The pilot imports V1 fit/predict/bootstrap functions; pre-posterior diagnostics accept only V1 model labels.

## 1. D-01 evidence boundary

The commit search for `D-01` returned no results. This does not rule out local, unpushed, archived, or external records. The historical root cause is therefore **OPEN / UNKNOWN**.

Keep these failure mechanisms separate:

| Mechanism | Diagnostic signature | Prospective discriminator |
|---|---|---|
| Complete separation | Separating direction; no finite unconstrained logistic MLE | Separable fixture; LP/separation diagnostic; coefficient path; constrained boundary/KKT check |
| Quasi-complete separation | Likelihood has a non-finite maximizing sequence with boundary observations | Boundary/tie fixture; compare LP diagnostic and coefficient path |
| Rank deficiency | Design matrix lacks full column rank; Hessian rank deficient | Constant/duplicated predictor columns; inspect design rank and Hessian eigenvalues |
| Near rank deficiency | Small positive eigenvalue, large condition number | Cluster predictor values and vary spread; inspect condition and sensitivity |
| Boundary-limited constrained optimum | Coefficient at box face with small projected-gradient residual | Strong-signal/separation fixture; verify active-constraint KKT signs |
| Optimizer defect | Non-finite/invalid step or failure on a well-conditioned interior problem | Derivative finite differences, known-answer convex fixture, independent solver comparison |
| Classification defect | Status, KKT residual, boundary flag and failure label disagree | Assert status invariants across all adversarial fixtures |

A boundary-limited constrained optimum is not an unconstrained MLE. V2's Hessian rank and boundary fields improve prospective diagnosis but do not identify the historical D-01 cause.

## 2. V1/V2 parity status

The inspected test file already had:
- V2 linear and quadratic curve-recovery tests;
- a separated-data boundary test;
- a rank-deficient quadratic-design test;
- Hessian rank/conditioning diagnostics;
- a same-data quadratic V1/V2 parity test.

It lacked the same-data linear parity test. The added test uses identical deterministic redshifts/outcomes, the same standardized linear design and coefficient box, and compares coefficients, grid predictions and negative log likelihood.

**Commit:** `ad8dbff38130bfe105589171812380a5352a3181`.  
**Execution:** pending; no test pass is claimed.

The parity tests still need adversarial comparisons under separation, near-separation, near-rank-deficiency and boundary limitation. For separated data, compare the bounded objective and KKT conditions rather than a nonexistent finite unconstrained MLE.

## 3. Feasible-box probability and weight range

The domain standardizes redshift to (x\in[-1,1]).

Linear:
[
\eta(x)=a+bx,\quad a\in[-25,25],\ b\in[-50,50].
]
The feasible minimum is -75, so:
[
\operatorname{expit}(-75)\approx2.68\times10^{-33},\quad
1/\pi_{min}\approx3.73\times10^{32}.
]

Quadratic:
[
\eta(x)=a+bx+qx^2,\quad q\in[-50,50].
]
The feasible minimum is -125, so:
[
\operatorname{expit}(-125)\approx5.17\times10^{-55},\quad
1/\pi_{min}\approx1.94\times10^{54}.
]

These remain representable as positive double-precision probabilities, but finiteness is not statistical usability, adequate overlap, or stable reconstruction. These are theoretical extrema of the permitted boxes, not estimates of fitted pilot behavior.

The richer `preposterior_diagnostics.py` reports probability minima/quantiles, Kish ESS, entropy ESS, ESS fraction, maximum normalized weight share, top 1%/5% weight mass, inverse-weight quantiles, and validation-support gaps/boundary proximity. The KDE implementation scales weights by their maximum before moments and normalization. This prevents avoidable common-scale overflow but cannot fix concentration or lack of support. No universal ESS cutoff should be invented; warning thresholds require predeclared, design-specific calibration. Do not silently clip probabilities or cap weights.

## 4. Host-density construction and mixture check

The inspected estimated-inclusion route:
1. Fits `P(included=1 | z)` from external validation redshifts and retention indicators.
2. Evaluates fitted inclusion probabilities at retained catalogue redshifts.
3. Forms (w_i=1/\hat\pi(z_i)).
4. Builds a bounded, reflected Gaussian KDE using those weights.
5. Multiplies the KDE by the positive grid-level merger-rate function.
6. Renormalizes the resulting host density to unit integral.
7. Passes that single density to the selection-aware fixed-count likelihood.

At population level, the retained density is proportional to (f(z)\pi(z)), and inverse weighting targets the parent (f(z)) only under appropriate validation-design identification, positivity, model adequacy and controlled estimation error. Merger-rate multiplication then changes the target to the merger-weighted host density.

No separate ((1-f)L_{catalogue}+fL_{missing}) term appears in this inspected pilot path. Thus double correction is **not demonstrated here**. This finding does not apply automatically to other likelihood modules or branches; any future mixture requires a complete joint factorization audit.

## 5. Selection-likelihood interpretation

The pilot supplies the same normalized host density to the event numerator and the detection-efficiency normalizer. Under the documented Gaussian-distance/hard-threshold toy model, this is internally coherent for a fixed detected count. It is not an event-rate likelihood and is not a realistic LVK selection function.

## 6. V2 integration gap

The pilot runner imports V1 fitting, prediction and bootstrap routines. The pre-posterior diagnostics API accepts only `logistic_inclusion_v1` and `quadratic_logistic_inclusion_v1` labels, whereas V2 emits separate labels. V2 therefore has unit-test coverage but is not the estimator used by the current pilot. Do not describe existing pilot results as V2 results.

Any future integration needs an explicit solver selector, preserved solver identity/diagnostics in output artifacts, and tests for identical-data parity and failure accounting before campaign use.

## 7. Next gates

1. Recover original D-01 records if available; otherwise keep the historical cause unknown and execute the prospective failure matrix.
2. Run exact-head tests including the new linear parity test and existing quadratic parity test.
3. Add feasible-box probability/weight stress tests and verify diagnostic behavior near the extremes.
4. Establish warning thresholds from repeated-sampling behavior, not arbitrary ESS cutoffs.
5. Keep V2 opt-in until pilot and artifact schemas explicitly support it.
6. Freeze the population target and complete likelihood factorization before confirmatory simulations.

## Evidence boundary

This is source inspection and a test-file change. It does not report a new simulation, test pass, CI pass, recovered D-01 traceback, calibrated estimator, or confirmatory result. No production estimator, likelihood, DGP, main-branch file, merge, or release declaration was changed.
