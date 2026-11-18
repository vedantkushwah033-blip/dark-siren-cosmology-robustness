# Full-support inclusion quadrature gate — implementation protocol

Date: 2026-11-18  
Implementation branch: `research/quadrature-gate-implementation-20261118`

## Purpose

The curvature density pilot estimates catalogue inclusion probabilities from a
finite validation sample and evaluates those probabilities over the observed
catalogue. Its configured order-5 tensor Gauss–Hermite approximation has 25
nodes for the two-parameter working logistic model. The new gate evaluates
whether the posterior summaries and downstream normalized host weights are
stable as quadrature order increases.

This gate is numerical. It does not establish that the degree-one working
model is correctly specified, that inverse-probability weighting is unbiased
under curvature misspecification, or that the resulting H0 posterior is
calibrated.

## Inputs and computation

For each persisted replicate and validation size, the audit:

1. Reads the saved validation redshifts and retention indicators, plus the
   saved observed catalogue redshifts.
2. Refits the proper degree-one logistic validation posterior at orders
   `(5, 7, 9, 11, 13)` using the same coefficient coordinate system
   (`center=0.205`, `scale=0.195`) as the density pilot.
3. Explicitly evaluates posterior-mean inclusion predictions and all
   quadrature-node inclusion probabilities over 401 points spanning
   `[0.01, 0.40]`.
4. Records posterior coefficient moments, log evidence, quadrature ESS,
   maximum normalized quadrature-node weight, and full-support inclusion
   ranges at every evaluated order.
5. Computes posterior-averaged normalized host weights for every catalogue
   object. These use the same linear merger-rate factor as the density pilot:
   `b(z)=1+1.5(z-0.01)/(0.40-0.01)`, multiplied by the node-specific inverse
   inclusion probability before normalization.
6. Records the fraction of catalogue objects outside the realized validation
   redshift interval and the posterior-averaged rate-weighted IPW mass in that
   extrapolation region.
7. Compares posterior-mixture host weights between adjacent orders and between
   the first and final evaluated orders.

The audit writes JSON and CSV diagnostics outside the raw replicate artifact
directories. This avoids changing the original replicate manifests.

## Gate logic

The workflow currently fixes the order schedule to `(5, 7, 9, 11, 13)`,
requires two consecutive passes in the existing quadrature-ladder function,
and uses these initial numerical tolerances:

- coefficient posterior mean maximum absolute change: `1e-3`;
- coefficient posterior covariance maximum absolute change: `1e-3`;
- full-support posterior-mean inclusion maximum absolute change: `1e-3`;
- log-evidence absolute change: `1e-3`;
- posterior-mixture rate-weighted host-weight maximum absolute change:
  `1e-3`;
- posterior-mixture host-weight total variation distance: `5e-3`.

These thresholds are **provisional engineering tolerances**, not theoretically
derived error bounds and not yet empirically justified for H0 inference.
They must be reviewed against numerical sensitivity and downstream likelihood
requirements before confirmatory use. A gate pass means only that the listed
numerical summaries are stable under this finite order ladder and these
thresholds.

The audit fails closed if a planned replicate/validation-size case is missing,
the quadrature ladder does not converge, a fit raises an exception, or
posterior-mixture host weights exceed the declared stability tolerances. Raw
replicate artifacts remain preserved even when the gate fails.

## Interpretation boundaries

- A pass does not correct functional-form misspecification in the degree-one
  logistic inclusion model.
- A pass does not show that the validation sample covers the catalogue support;
  extrapolation fraction and extrapolation-region weight mass are reported as
  diagnostics, not acceptance criteria in this first implementation.
- The gate checks the posterior-mixture normalized host weights used by the
  fixed-bandwidth density construction. It does not validate all density
  metrics, bandwidth selection, or statistical uncertainty across replicates.
- The original density runner still computes its five-arm density outputs at
  the configured order 5. The gate checks whether that order is numerically
  stable relative to higher orders; it does not silently replace the original
  outputs with a different quadrature order.
- There is no H0 likelihood or posterior calculation in this workflow. A
  separate Q3 gate must compare the complete shared-coefficient-node selected
  likelihood and H0 posterior across quadrature orders.
- Numerical convergence is separate from simulation-based calibration and
  repeated-sampling coverage.

## Required next steps

1. Run CI on the implementation branch and repair any test or lint failures.
2. Execute a small, separately labelled local smoke fixture or GitHub Actions
   run to inspect the gate's artifact schema and computational cost.
3. Review the provisional tolerances and decide whether extrapolation-region
   weight mass needs a predeclared maximum or remains descriptive.
4. Run the seven-gamma density pilot only after the gate implementation is
   tested and its runtime is acceptable.
5. Implement Q3 on the integrated H0 likelihood path before making any H0
   calibration or astrophysical performance claim.

## Current execution status

The code and workflow integration are committed on a draft pull request.
The seven-gamma workflow is manual-dispatch-only and has not been dispatched
by this implementation. The project test suite has not been executed locally
because this environment cannot resolve GitHub to obtain a repository checkout.
Any CI result must be checked against the exact current head commit before the
implementation is considered tested.
