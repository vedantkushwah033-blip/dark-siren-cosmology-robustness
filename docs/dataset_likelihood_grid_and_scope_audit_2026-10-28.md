# Dataset-level likelihood review: prior-grid and model-scope checks

**Date:** 2026-10-28  
**Reviewed revision:** `80568631442d160b9b273df1ccbb10a568a2ba4f`  
**Status:** source review; latest workflow was pending when checked.

## Summary

The new `dataset_likelihood.py` is a useful minimal composition layer for conditionally independent events sharing a single H0. Its docstring correctly excludes selection conditioning, event-rate inference, and shared nuisance parameters. The accompanying tests cover several important algebraic limiting cases.

This is not yet a population-level dark-siren likelihood. It is a discrete-grid posterior prototype for a declared toy event likelihood.

## Finding 1 — grid support does not require prior boundaries to be represented

The function checks that at least two grid points lie within `[prior_min, prior_max]`, then normalizes the posterior using only those supported grid points and trapezoidal integration. It does not require the grid to include either prior boundary, nor does it insert the boundaries by interpolation.

This behavior is a legitimate discrete-grid convention if explicitly intended, but it can silently approximate a continuous uniform prior on a truncated grid interval rather than the declared full interval. Near a prior edge, the difference can affect posterior normalization and interval summaries.

**Recommendation:** choose and test one explicit policy:
1. Require grid points to include both prior boundaries within a documented tolerance; or
2. Interpolate log likelihood to the exact prior boundaries and integrate on the augmented grid; or
3. Clearly define the prior as the discrete grid restricted to supported grid nodes and document its quadrature convention.

For research-facing continuous-grid inference, option 2 is generally the more flexible choice, but it needs independent numerical tests. Do not silently alter the existing API without reviewing callers.

## Finding 2 — non-negative redshift validation is broader than the intended extragalactic domain

The new module accepts `z = 0`. The low-redshift approximation then predicts zero distance for that host at every positive H0. The observation itself must have positive distance, and the current positive-truncated likelihood has a well-defined formal expression at zero mean, but the host is not a meaningful extragalactic catalogue entry for this scientific application.

**Recommendation:** define whether this module is a generic mathematical primitive or an extragalactic host-catalogue API. For the latter, require strictly positive redshifts and add a regression test. If zero is retained for mathematical tests, label it explicitly as a boundary case, not a realistic host.

## Finding 3 — independent-event product requires conditional independence and common parameter semantics

The sum of event log likelihoods is mathematically appropriate only conditional on the assumptions in the docstring: events are conditionally independent given H0 and all event-specific inputs are treated as fixed. Shared calibration uncertainties, shared population parameters, event selection, and uncertain event counts cannot be added by simply multiplying these event terms.

**Recommendation:** keep the current function explicitly scoped as a conditional-independence prototype. A future hierarchical extension should place shared nuisance parameters inside a joint model and marginalize them once, rather than independently per event.

## Finding 4 — missing-population parameterization is per event

Each event mapping has its own `missing_weight` and missing-host support. This is computationally useful, but it does not enforce a common catalogue completeness or shared missing-host population across events. Different values could be supplied accidentally without an explicit warning.

**Recommendation:** the confirmatory experiment should define whether completeness is event-specific or shared, validate that design in its runner, and save the actual per-event values in raw records. If the intended model uses a common fraction, use a shared top-level parameter or assert equality in the experiment configuration.

## Required tests before confirmatory use

- Prior bounds exactly represented, and behavior when they are not.
- Convergence against a dense-grid or independent quadrature reference.
- Posterior summary stability for truths near each prior edge.
- Strictly positive host-redshift validation, or a documented zero-redshift boundary test.
- Multiple events with different uncertainties and deliberately shared versus event-specific nuisance settings.
- Explicit rejection or documented support for selection-conditioned inputs (currently unsupported).
- End-to-end raw-data provenance: every simulated dataset ID, event ID, truth, model arm, seed, and parameter configuration must be recoverable independently.

## Conclusion

No new runtime failure was demonstrated by this source inspection alone. The principal concern is a potentially silent mismatch between a declared continuous prior interval and a grid whose endpoints do not coincide with the prior boundaries. This should be resolved before using the prototype for confirmatory calibration or diagnostic-identifiability claims.
