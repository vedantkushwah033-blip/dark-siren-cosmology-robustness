# Quadrature-convergence gate audit — 2026-11-18

## Status and scope

This is a static source audit of the current research branch, not a run of the
curvature experiment, not a test-suite result, and not evidence about estimator
performance. The inspected branch is
`research/nov18-project-audit-2026`, based on the November 18 audit branch.
The curvature workflow remains manual-dispatch-only; no workflow was dispatched
and no simulation artifacts were produced in this audit.

## Verified implementation path

The current curvature workflow
(`.github/workflows/quadratic-logit-curvature-density-pilot.yml`) invokes
`scripts/run_absolute_bandwidth_density_replicates.py` with
`--quadrature-order 5`. The multi-replicate runner forwards that setting to
the single-replicate smoke runner. The smoke runner passes it to
`build_five_arm_weights`; that function calls
`fit_inclusion_posterior_quadrature(..., order=quadrature_order)` directly.

Therefore, each validation-size fit in the configured workflow uses a single
tensor Gauss–Hermite order of 5 for a two-coefficient logistic model: 25
quadrature nodes. The configured seven-gamma design has 20 replicates per
gamma, three validation sizes, and 7 gamma cells, but no quadrature-order
refinement is performed in this workflow.

## Existing ladder and its limitations

The repository already implements
`fit_inclusion_posterior_quadrature_ladder` in
`src/dark_siren_robustness/inclusion_posterior.py`. It compares adjacent
orders on posterior coefficient means, coefficient covariance, posterior-mean
inclusion predictions, and log evidence, and requires a configurable number of
consecutive passes.

However:

1. If `predictive_redshift` is omitted, its prediction grid is constructed
   from the minimum and maximum redshifts in the validation sample. This does
   not test the full catalogue support `[0.01, 0.40]` when the validation
   sample covers a narrower interval.
2. The curvature workflow does not call the ladder at all; it uses one
   order-5 fit.
3. The ladder docstring explicitly limits its convergence claim: it does not
   establish convergence of a downstream H0 likelihood.
4. Existing ladder tests establish function behavior on generated examples;
   they do not establish convergence for the actual curvature design, its
   validation sizes, or the downstream dark-siren likelihood.

The concern is consequently confirmed at source level. It is a numerical
validation gap, not proof that order 5 is inaccurate in every replicate.

## Why this matters

The estimated inclusion model is fitted on validation redshifts, then its
coefficient nodes are evaluated at catalogue redshifts to construct inverse
inclusion weights. When catalogue support extends beyond the realized
validation range, the fitted model is extrapolated. Small changes in
coefficient nodes can have amplified effects on inverse probabilities and
weight concentration in those regions.

Coefficient-moment convergence and posterior-mean inclusion convergence are
necessary diagnostics, but they do not guarantee stability of:
- node-specific inverse weights;
- normalized host weights and their ESS/tail concentration;
- reconstructed host-redshift density;
- the selected event likelihood as a function of H0; or
- the marginalized H0 posterior and its calibration.

The current density-only pilot is not an H0 validation experiment. No H0
calibration conclusion should be drawn from its output.

## Required numerical gate before confirmatory use

### Gate Q1 — inclusion-posterior quadrature

For every validation sample used in a confirmatory analysis:

- Fit an order ladder such as `(5, 7, 9, 11, 13)`, extending it if the
  predeclared tolerances are not met.
- Supply an explicit fixed prediction grid spanning the entire declared
  catalogue support, e.g. `np.linspace(0.01, 0.40, 401)`; do not rely on the
  ladder's validation-range default.
- Compare posterior means/covariances, log evidence, and posterior-mean
  inclusion curves over the full support.
- Record node-weight ESS and maximum normalized quadrature weight at each
  order. A small number of dominant quadrature nodes is a separate warning
  even if summary differences appear small.
- Record validation minimum/maximum, the fraction of catalogue objects outside
  the validation interval, and the inverse-weight mass contributed by those
  objects for each posterior node.

The order schedule and numerical tolerances must be fixed before examining
performance contrasts. A non-converged ladder is a non-evaluable numerical
case, not permission to select a favorable order after looking at H0 results.

### Gate Q2 — downstream host-weight and density stability

Across the converged quadrature orders, compare:
- node-specific inclusion probabilities over all catalogue objects;
- normalized node-specific host weights;
- Kish and entropy ESS, maximum normalized weight, and top 1%/5% mass;
- density curves and density-risk metrics on the same fixed redshift grid.

This gate is distinct from Q1 because nonlinear inverse weighting can amplify
small inclusion-probability differences.

### Gate Q3 — downstream H0 stability

Before any confirmatory H0 statement, evaluate the complete shared-node
selected likelihood at each accepted quadrature order. Compare:
- the full H0 log-likelihood curve after a common additive normalization;
- marginalized posterior density on the identical H0 grid;
- posterior mean/median and interval endpoints;
- numerical normalization and grid-refinement diagnostics.

The same coefficient node must remain shared across the joint event likelihood.
Do not marginalize each event independently and then multiply those
marginalized likelihoods.

Quadrature stability is not posterior calibration. After Q3, repeated-sampling
simulation-based calibration and coverage remain separate required gates.

## Recommended implementation sequence

1. Add a test that the ladder is always called with an explicit full-support
   prediction grid in the curvature analysis path.
2. Add an order-refinement artifact per replicate and validation size, including
   order, node count, quadrature ESS, maximum node weight, coefficient moments,
   evidence, and full-support inclusion predictions.
3. Add a strict audit that fails closed on non-convergence or missing
   full-support coverage; preserve the complete failure record.
4. Add a separate H0 quadrature-ladder audit only after the density experiment
   and selected-likelihood implementation are integrated and independently
   validated.
5. Keep the existing order-5 curvature run labelled exploratory unless it is
   rerun with the numerical gate and all planned replicates pass.

## Decision

**Current order-5 curvature density workflow: numerical-convergence gate NOT
PASSED / NOT EVALUATED.** The source confirms that no quadrature ladder is
invoked. The existing ladder API is a useful starting point, but it defaults
to validation-range predictions and has no downstream H0 convergence gate.

No code or test was executed as part of this audit. No claim is made that the
configured workflow has run, that order 5 fails, or that any estimator wins.
