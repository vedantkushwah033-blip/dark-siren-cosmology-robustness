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

The posterior-averaged host-weight vector is a pre-KDE diagnostic. It is not
algebraically identical to the posterior-mean density, because each node KDE
is normalized separately before posterior averaging. A separate curved-DGP regression now also checks the posterior-mean Arm-E
density curve after node-wise KDE normalization; this remains a single-seed
numerical regression, not repeated-sampling density validation.

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
- The gate checks posterior-averaged normalized rate-weighted host weights as a
  pre-KDE diagnostic. It does not establish stability of the posterior-mean
  density curves (whose node-wise normalization matters), all density metrics,
  bandwidth selection, or statistical uncertainty across replicates.
- The original density runner still computes its five-arm density outputs at
  the configured order 5. The gate checks whether that order is numerically
  stable relative to higher orders; it does not silently replace the original
  outputs with a different quadrature order.
- There is no H0 likelihood or posterior calculation in this seven-gamma
  density workflow. The repository does already contain a narrow H0
  quadrature-regression test: one floor-logistic synthetic case with validation
  size 300, three toy events, and orders 7/9/11, comparing the order-9 and
  order-11 Arm-E H0 posteriors. A separate H0-grid test checks 81/161/321-point
  grids for one fixed synthetic case. These are useful existing tests, but they
  do not cover the quadratic-curvature DGP, all three validation sizes, the
  seven-gamma repeated-sampling design, or H0 posterior calibration. A Q3 gate
  for the curvature project must compare the complete shared-coefficient-node
  selected likelihood and H0 posterior across quadrature orders under the
  actual curvature design.
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
5. Add a density-level Q2 comparison of posterior-mean density curves across
   quadrature orders and the declared absolute-bandwidth grid.
6. Extend the existing H0 quadrature and H0-grid tests to the curvature DGP,
   then implement Q3 under the actual repeated-sampling design before making
   any H0 calibration or astrophysical performance claim.

## Current execution status

The code and workflow integration are committed on a draft pull request.
The seven-gamma workflow is manual-dispatch-only and has not been dispatched
by this implementation. The project test suite has not been executed locally
because this environment cannot resolve GitHub to obtain a repository checkout.
Any CI result must be checked against the exact current head commit before the
implementation is considered tested.

## Curved-DGP H0 quadrature regression added after initial protocol

A separate regression, `tests/test_curved_dgp_h0_quadrature_convergence.py`,
now constructs a quadratic-logit catalogue-inclusion DGP with
`gamma=1.2`, draws a finite validation sample, fits the deliberately
degree-one working logistic model, and compares the complete Arm-E
selected-likelihood H0 posterior at coefficient quadrature orders 9 and 11
(with order 7 retained as an additional ladder point).

This is a numerical regression only. It checks whether the shared-coefficient
node integration is stable for one deterministic curved-DGP realization. It
does **not** test repeated-sampling H0 bias, coverage, PIT/SBC calibration,
validation-size sensitivity, the seven-gamma design, or correction of the
working-model misspecification. It is now included in the project CI workflow;
its result must be checked against the exact commit that contains it.



## Posterior-mean density convergence regression

The curved-DGP regression now retains the Arm-E node-specific KDE densities
at coefficient quadrature orders 7, 9, and 11. At each order it forms the
posterior-mean density by averaging the separately normalized node densities
with their posterior quadrature weights, verifies unit integral on the
redshift grid, and compares the order-9 and order-11 curves using total
variation distance. The provisional numerical acceptance threshold is
TV <= 0.01 for this deterministic test fixture.

This closes the specific algebraic gap between convergence of posterior-
averaged raw host weights and convergence of a posterior mean of normalized
node-wise KDE curves. It does not establish density-reconstruction accuracy,
statistical uncertainty, or performance across repeated samples. The
threshold is an engineering regression tolerance and should be checked
against finer redshift grids and repeated seeds before confirmatory use.
