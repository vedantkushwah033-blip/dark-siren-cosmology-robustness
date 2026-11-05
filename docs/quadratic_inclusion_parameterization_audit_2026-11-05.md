# Quadratic inclusion-model parameterization and coefficient-bound audit

**Date:** 5 November 2026  
**Status:** independent mathematical review and proposed sensitivity-test specification; no estimator change or empirical result is claimed  
**Applies to:** estimated_inclusion.py and the nonlinear curvature-screen study

## 1. Finding

The synthetic data-generating process writes the quadratic inclusion model using a parent-law-centered basis,

    eta_c(x) = alpha + beta*x + gamma*(x^2 - c),

where c = E_parent[x^2] is computed from the declared parent-redshift intensity. The fitted quadratic logistic estimator currently uses the uncentered basis,

    eta_u(x) = a + b*x + q*x^2.

These are exactly the same unconstrained function family. The coefficient transformation is

    a = alpha - gamma*c
    b = beta
    q = gamma

and its inverse is

    alpha = a + q*c
    beta = b
    gamma = q.

However, the current fitter constrains the uncentered intercept to [-25, 25], slope to [-50, 50], and quadratic coefficient to [-50, 50]. Applying the same numerical intercept bounds to a centered-basis fit would generally define a different feasible set of probability curves.

For a fixed q, the current uncentered constraint -25 <= a <= 25 is equivalent to the q-dependent centered-intercept constraint

    -25 + q*c <= alpha <= 25 + q*c.

It is not equivalent to the fixed centered constraint -25 <= alpha <= 25. Comparing two parameterizations with identical-looking coefficient boxes would therefore confound basis choice with a change in the constrained model.

This is a potential sensitivity issue, not evidence that the existing fit is biased or defective. The synthetic truth uses modest curvature values, so the transformed intercept shift gamma*c may be small in the prespecified screen. Actual impact must be measured rather than presumed.

## 2. What must be compared

The audit should separate three questions.

### A. Algebraic prediction invariance

For arbitrary finite coefficient vectors and valid x, verify numerically that the transformed centered and uncentered linear predictors agree to floating-point tolerance. This is a deterministic algebra check, not a simulation result.

### B. Constrained-optimization equivalence

Fit the same validation records under:

1. The existing uncentered parameterization and existing coefficient box.
2. The centered parameterization with the transformed, q-dependent intercept bounds above and unchanged slope/quadratic bounds.

These two optimizations represent the same feasible set in function space. Their fitted linear predictors, probabilities, objective values, and boundary classifications should agree within tolerances justified by the optimizer.

A centered implementation that simply reuses the fixed intercept interval [-25, 25] is not an equivalent parameterization.

### C. Sensitivity to the declared coefficient box

Separately repeat the uncentered fit with a prespecified wider coefficient box. This asks whether the existing box is influential; it is not a basis-invariance test. Record any changes in optimizer status, KKT/projected-gradient residual, coefficient-bound hits, predicted probabilities, inverse-weight diagnostics, KDE, and downstream H0 summaries.

Do not widen production bounds silently. The sensitivity arm must retain the original fit as its named baseline.

## 3. Required deterministic regression tests

Add focused tests before using this sensitivity analysis in the scientific screen:

- **Coefficient transformation:** generate several finite (alpha, beta, gamma, c) tuples and verify equality of the centered and transformed uncentered predictors over a dense x grid.
- **Probability equality:** verify equality after applying the logistic inverse link.
- **Equivalent feasible bounds:** for several q values, verify that the centered q-dependent intercept interval maps exactly to the uncentered intercept interval.
- **Fit equivalence:** on a fixed, non-separated validation fixture with both outcome classes, fit both equivalent parameterizations and compare objective, predicted probabilities, and KKT conditions. Preserve any optimizer discrepancy as a test failure or explicitly justified tolerance; do not compare raw intercepts without transforming them.
- **Non-equivalent fixed-box demonstration:** include a small algebraic test showing that fixed centered and uncentered intercept boxes are not generally the same when q*c is nonzero. This guards against accidentally claiming equivalence where none exists.
- **Boundary case:** construct a deterministic fixture near an intercept or quadratic bound and check that boundary status is reported consistently under the equivalent transformed constraints.
- **Prediction API:** if the fitted-model representation is changed, retain backward compatibility or add explicit model-version dispatch and tests.

Tests must be deterministic and fast. Do not add stochastic Monte Carlo tests where an algebraic identity suffices.

## 4. Diagnostics to add to the sensitivity report

For each validation size and curvature condition, preserve:

- Fit status and optimizer message.
- Coefficient-bound hit flags and which coefficient hit.
- Projected-gradient/KKT residual.
- Maximum absolute difference between paired fitted probabilities.
- Maximum absolute difference between paired linear predictors.
- Difference in negative log likelihood.
- Difference in inverse-probability weight quantiles and maximum normalized weight share.
- Difference in ESS and entropy ESS.
- Difference in KDE density on the fixed redshift grid.
- Difference in event log likelihood and normalized H0 posterior summaries, where the likelihood contract permits that comparison.
- Attempted, successful, and failed fit counts by parameterization.

Report both conditional-on-joint-success comparisons and unconditional operational failure rates. Do not discard a fit failure in one parameterization merely because its paired counterpart succeeded.

## 5. Interpretation rules

- Exact algebraic equivalence does not establish numerical optimizer equivalence.
- Numerical optimizer equivalence under transformed bounds does not establish robustness to a wider coefficient box.
- A small change in fitted inclusion probabilities does not guarantee a small change in downstream H0 summaries; the inverse-weight/KDE/likelihood map is nonlinear.
- A difference between same-number fixed coefficient boxes is not evidence of a basis effect, because those boxes define different feasible function classes.
- A boundary-hit reduction under wider bounds is not automatically an improvement; assess probability calibration, overlap, weight concentration, and downstream performance together.
- Do not select a parameterization or coefficient box based on whichever produces the narrowest H0 posterior.

## 6. Release gate

Before the quadratic curvature screen is used to motivate confirmatory design:

1. Land deterministic algebra and equivalent-bound regression tests.
2. Verify the existing uncentered baseline remains unchanged.
3. Run the deterministic fit-equivalence and wider-box sensitivity tests.
4. Record exact source revision and test evidence.
5. If the downstream sensitivity is material, revise the exploratory interpretation and freeze the chosen parameterization/bounds before a confirmatory campaign.
6. Keep the full five-curvature exploratory screen labelled non-confirmatory until its actual execution and independent artifact audit are complete.

## 7. Evidence boundary

This note is a mathematical review and test specification only. No code, estimator, coefficient bound, configuration, simulation result, or manuscript result was changed by writing it. The full five-curvature screen remains unexecuted in the verified evidence record. All project claims remain restricted to synthetic toy-model studies.

## 8. Implementation update — 5 November 2026

A deterministic regression-test file has now been added at
`tests/test_quadratic_parameterization_audit.py`. It covers:

- Equality of centered and uncentered linear predictors over a dense fixed grid for several coefficient tuples.
- Equality after applying the logistic inverse link.
- Exact mapping of the q-dependent centered intercept interval to the existing uncentered intercept interval.
- A counterexample showing that reusing the same fixed numeric intercept box in both bases is not generally equivalent.
- The zero-centering special case.

The production estimator, its default coefficient bounds, the data-generating process, and all experiment configurations remain unchanged.

**Not yet implemented:** the paired constrained-optimizer equivalence fixture, boundary-case optimizer comparison, and wider-box downstream sensitivity arm. The algebra tests do not establish optimizer equivalence or downstream H0 robustness.

**CI status at the time of this update:** the new test commit has triggered workflow run
`41999055457`; it was pending when last queried. The preceding documentation-only commit's run
`41998117284` was still in progress. No pass is claimed for the new tests until the exact-head run completes.
