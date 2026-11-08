# Prior-coordinate audit for Bayesian inclusion inference

Date: 2026-11-08  
Scope: synthetic independent-validation, shared-inclusion H0 pipeline  
Audited head: `bdd84413bb1da8e380ce6cbad70d25ec0b88c362`

## Finding: default coefficient priors depend on the realized validation sample

The current `_validation_inputs` implementation in
`src/dark_siren_robustness/inclusion_posterior.py` defines

```python
center = float(np.mean(z))
scale = float(np.std(z, ddof=0))
x = (z - center) / scale
design = [1, x, x**2, ...]
```

The default prior is then applied to these sample-standardized coefficients:
`sd = [5.0, 2.5, 2.5]` (truncated to the selected model dimension).

For a degree-one model,
```
eta(z) = beta0 + beta1 * (z - mean(V)) / sd(V)
```
so the implied prior standard deviation of the physical-redshift slope is
`2.5 / sd(V)`. It therefore changes with the realized validation sample.
For degree two, the physical quadratic coefficient has prior scale
`2.5 / sd(V)^2`, and the sample-dependent center also changes the induced
joint prior on physical intercept, slope, and curvature.

This is not an algebraic defect in the posterior calculation. It is a
comparability and experimental-design issue: across repeated validation
replicates or validation sample sizes, the prior on the inclusion-probability
function itself changes along with the data. Consequently, an apparent
validation-size effect can partly reflect a changing prior rather than only
additional validation information.

The current single-realization smoke test does not reveal this issue because
it fits all quadrature orders to one fixed validation sample. A repeated-
sampling experiment would be affected.

## Required correction before Bayesian coverage/SBC claims

Use a fixed, declared coordinate system tied to the simulation domain, for
example
`x = (z - z_mid) / z_half_range`, where both constants are fixed by the
predeclared redshift domain, not estimated from each validation sample.
Apply the same coefficient prior in this fixed basis for every replicate,
validation size, and model fit.

An alternative is to retain sample standardization but transform the prior
mean and covariance from a single fixed physical-coordinate prior into each
sample-specific standardized basis. Merely reusing the same diagonal
`prior_sd` after changing center/scale does not preserve the same function
prior.

## Tests required

1. Verify that the fixed-coordinate design uses identical center and scale
   across validation samples and sample sizes.
2. Draw coefficient vectors from the declared prior and compare the induced
   inclusion-probability curves across validation samples; they must be
   identical for a fixed coefficient draw.
3. Verify the algebraic transformation between fixed-basis and
   sample-standardized coefficients if the transformed-prior route is used.
4. Add prior-scale sensitivity using at least a predeclared narrow, baseline,
   and broad prior, while keeping the coordinate system fixed.
5. In repeated-sampling work, report prior sensitivity separately from
   validation-size effects. Do not pool these effects.

## Interpretation boundary

This finding does not invalidate the current smoke artifact as a check of
the implemented calculation under its declared sample-dependent prior.
It does mean the current default prior is not yet suitable for clean
cross-replicate or cross-validation-size Bayesian comparisons without an
explicit prior-coordinate decision.

No calibration, simulation-based calibration, frequentist coverage, or
real-data conclusion follows from this audit.
