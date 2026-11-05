# Versioned logistic inclusion estimator v2

**Date:** 5 November 2026  
**Branch:** `research/logistic-inclusion-v2-newton-2026-11-05`  
**Status:** active-set Newton refinement and adversarial V1/V2 tests added; the dedicated curvature regression passed on the active-set code commit; full three-version test workflow for the latest test additions is pending.  
**Relationship to v1:** opt-in alternative. No v1 production call path, pilot result, or manuscript result is changed.

## Motivation

The v1 linear and quadratic inclusion fitters use bounded L-BFGS-B and then check a projected-gradient/KKT residual. The historical D-01 pilot failures show why optimizer status and KKT checks need to be distinguished. V2 adds a separately named damped-Newton solver so comparisons can preserve estimator provenance.

## Algorithm

For design matrix X, binary outcomes y, and coefficient vector theta, v2 minimizes the convex logistic negative log-likelihood

`ell(theta) = -sum_i [y_i log sigmoid(x_i' theta) + (1-y_i) log sigmoid(-x_i' theta)]`.

The quadratic model uses the same uncentered `x**2` basis as the production v1 fitter. This is essential because the intercept coefficient bounds are part of the constrained estimator; centering shifts the intercept by `quadratic/3` and does not preserve the same feasible coefficient box. It computes the analytic gradient `X' (p-y)` and Hessian `X' diag(p(1-p)) X`. At each iteration, coordinates on a box face whose gradient sign satisfies the KKT condition are treated as active; a damped Newton system is solved on the remaining free-coordinate subspace. This prevents an outward Newton component in an already-active coordinate from destroying descent after projection. If the reduced Hessian solve is singular or does not yield a descent direction, the method uses a projected-gradient direction. A projected Armijo backtracking search clips candidate coefficients to the declared v1-compatible box, with a final projected-gradient descent fallback.

Convergence is declared only when the infinity norm of the unit-step projected-gradient mapping is at most `1e-7`. A non-converged run raises an explicit RuntimeError; it is not relabelled as success. Coefficients touching a box face are returned as `boundary_limited`, including separated data where the finite box constrains the estimate.

## Public functions

- `fit_logistic_inclusion_model_v2`
- `fit_quadratic_logistic_inclusion_model_v2`
- `predict_inclusion_probability_v2`

Results identify their model version and retain iteration, accepted/rejected step, objective, projected-gradient, boundary, and validation-count diagnostics.

## Tests added

- Linear curve recovery on a well-conditioned synthetic sample.
- Quadratic curve recovery.
- Complete-separation fixture reported as boundary-limited.
- Identical-data V1/V2 comparison on a deterministic quadratic near-separation fixture with one discordant label; compare the fitted probability curve, objective, boundary status, and projected KKT residual.
- Identical-data V1/V2 parity under complete separation for both linear and quadratic models.
- Rank-deficient quadratic design and deterministic replay.
- Analytic gradient and Hessian checked against central finite differences.
- Final observed-Hessian rank, minimum/maximum eigenvalue, and condition-number diagnostics, including explicit infinite condition number for rank deficiency.
- Degenerate outcome and prediction-domain validation.

## Interpretation and limitations

1. The quadratic basis is deliberately production-v1-compatible: `intercept + slope*x + quadratic*x**2`. Although centering is an algebraic reparameterization of the same unconstrained curves, it shifts the intercept by `quadratic/3`; with fixed coefficient bounds, the feasible parameter sets differ. Therefore v1/v2 comparisons must retain the uncentered basis.
2. This is a convex optimization routine for the stated finite-dimensional logistic likelihood with box constraints. It does not establish that the logistic inclusion model is correctly specified.
3. The coefficient box is part of the estimator. Boundary-limited fits are scientifically important diagnostics, not successful unconstrained MLEs.
4. A small projected-gradient mapping establishes approximate first-order constrained stationarity. For a convex objective this is strong optimization evidence, but finite precision, ill-conditioning, and the chosen box still matter.
5. The Hessian can be singular under rank deficiency or near separation. Damping supports computation but does not create information; rank and conditioning diagnostics remain necessary for a confirmatory release.
6. The Hessian condition number is reported as infinity when the observed Hessian is numerically rank-deficient under a scale-aware floating-point threshold; this is a diagnostic, not a regularization or a proof of statistical identifiability.
7. V2 is not yet wired into the host-density estimator, pilot runner, bootstrap, or artifact schemas. Do not compare it with v1 downstream until an explicit versioned integration and independent artifact audit are added.
8. Unit tests and CI do not establish D-01 closure in the historical pilot. The historical failure cases must be replayed against both v1 and v2 using identical validation records and recorded outcomes.

## Release gate

D-01 remains open until:
- the derivative/conditioning and adversarial near-separation tests pass on supported Python versions;
- historical abnormal-line-search fixtures and the recorded 2/840 pilot failures are replayed;
- v1/v2 coefficients, objective values, predicted probabilities, and boundary diagnostics are independently compared;
- failure accounting is preserved;
- a downstream density/H0 comparison is performed only after the fitter-level gate passes.

No confirmatory simulation is authorized by this implementation alone.
