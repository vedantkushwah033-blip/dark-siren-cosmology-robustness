# Bounded host-density estimation for the synthetic selection benchmark

**Status:** estimator implementation and direct kernel-sum regression test added on an experimental branch; exact-head CI and regenerated benchmark artifacts are required before interpreting scientific results.  
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
- Reflection reduces boundary leakage but does not eliminate bandwidth bias. The bandwidth is computed **before reflection**, from the original weighted observed sample; it must not be scaled using the spread of the reflected/tripled sample.
- If `bandwidth_factor=None`, the implementation uses Scott's one-dimensional rule, with effective sample size `n_eff = (sum w)^2 / sum(w^2)` and `h = SD_w * n_eff^(-1/5)`. If an explicit factor is supplied, the convention is `h = factor * SD_w`, where `SD_w` is the weighted population standard deviation of the original sample. This is intentionally not SciPy `gaussian_kde`'s `bw_method`-float convention.
- The reflected density is evaluated as an explicit weighted Gaussian kernel sum in chunks, then truncated and normalized over the specified finite interval. The direct kernel-sum test independently reconstructs this calculation for both Scott-default and explicit-factor cases.
- Results remain bandwidth-sensitive and must be tested over predeclared factors and independent simulation conditions.

The estimator is now used by the synthetic selection-function misspecification runner, but the scientific campaign must be regenerated after this bandwidth correction. Do not reuse earlier artifacts generated under the prior bandwidth semantics. No numerical performance conclusion should be reported until exact-head CI, regenerated artifacts, and independent audits pass.
