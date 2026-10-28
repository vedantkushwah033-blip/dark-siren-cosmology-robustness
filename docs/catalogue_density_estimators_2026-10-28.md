# Bounded host-density estimation for the synthetic selection benchmark

**Status:** estimator prototype on an experimental branch; CI validation required before use in the main benchmark.  
**Scope:** synthetic-method development only; not a real survey model or H0 result.

## Intended benchmark arms

1. **Oracle:** parent-population redshifts from truth-only artifacts, with equal inclusion weights. This is a reference/scoring arm, not an observable-data analysis.
2. **Catalogue-only:** observed catalogue redshifts with equal weights. This deliberately ignores catalogue incompleteness.
3. **Known-completeness correction:** observed catalogue redshifts weighted by supplied per-object inclusion probabilities. This assumes probabilities are known; it does not estimate them.

All arms must apply the same declared merger-weight function and bandwidth sensitivity protocol. Non-oracle arms must never read `truth/` files.

## Estimator and assumptions

The estimator uses a Gaussian KDE with reflection at both redshift-domain boundaries, optional inverse-inclusion weights, and an optional positive merger-rate weight. It normalizes the final density on the requested grid so it can be supplied to the existing selection-aware toy likelihood.

- Inclusion probabilities are assumed known for the correction arm; their uncertainty or misspecification is not propagated.
- The inclusion probability must be positive throughout the domain. Zero-probability regions cannot be recovered by inverse weighting.
- Reflection reduces boundary leakage but does not eliminate bandwidth bias.
- Results are bandwidth-sensitive and should be included in robustness checks.
- The density is normalized over the specified finite interval.

This prototype is not yet wired into an H0 posterior runner, repeated-seed calibration, or coverage experiment. No inference result should be reported until those stages and independent numerical checks pass.
