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


## Additional finding: validation samples are not nested across requested sizes

The current `generate_independent_validation_sample` implementation selects rows with

```python
selected = rng.choice(n_validation_frame, size=n_validation, replace=False)
```

Reinitializing the generator with the same seed but changing `size` does not
guarantee that the smaller selected set is a subset of the larger selected set.
A direct NumPy check with frame size 4000 and seed 123 produced different first
10 selected IDs for requested sizes 100, 300, and 1000. Therefore, the PR-body
claim that same-seed larger SRS samples retain all earlier validation
observations is not supported by this implementation.

Each requested size is individually an SRS, but comparisons across sizes are
not paired by nested validation information. This matters if the experiment
interprets size contrasts as adding observations to an otherwise fixed
validation sample.

**Correction:** generate one full random permutation (or one maximum-size SRS
ordering) from the validation-frame RNG and take prefixes for each requested
sample size. Persist the ordering/selection IDs and independently audit the
subset relation. The analysis frame, validation frame, and all other random
streams must remain fixed across the paired size conditions.

This is separate from the sample-dependent prior-coordinate issue above.
Both need resolution before attributing repeated-sampling differences cleanly
to validation sample size.


## Follow-up implementation status (2026-11-08)

The pipeline-level correction is now implemented on the inclusion-posterior branch:

- The H0 smoke runner sets coefficient coordinates from the fixed simulation domain, using center 0.205 and scale 0.195 for the [0.01, 0.40] redshift support.
- The same fixed coordinates are passed to the posterior fitter at every quadrature order and to every prior-scale sensitivity refit.
- The prior-sensitivity utility accepts explicit coordinate center/scale and forwards them to each fit.
- A posterior test checks that the same fixed coefficient draw induces the same inclusion curve in that fixed basis.
- A prior-sensitivity test now checks that the three prior-scale refits retain identical coordinate center/scale while only the slope prior scale changes.

**Important scope:** the low-level fitter still permits its historical sample-standardized default when coordinate arguments are omitted. Therefore, callers conducting repeated-sampling comparisons must explicitly supply fixed coordinates. This is not silently changed for all existing callers.

The validation sampling implementation now constructs one full random permutation of the independent validation frame and takes the requested prefix. A new test regenerates the frame with the same seed at sizes 100, 300, and 1000, and verifies identical truth frames and strict nested ID prefixes.

These tests establish the intended code-level invariants; they do not constitute a completed multi-replicate experiment. The CI run for the latest branch head must pass before treating these tests as verified.

The earlier paragraphs above document the original defects and remain as an audit trail. The nested-sample description refers to the pre-correction implementation, not the current permutation-prefix implementation.


## Follow-up status at current branch head (2026-11-08)

The two findings above were addressed in subsequent branch commits after the
initial audit:

- `fit_inclusion_posterior_quadrature` now accepts an explicit fixed
  `coordinate_center` and `coordinate_scale`. The independent-validation
  H0 smoke runner uses the midpoint and half-range of the declared redshift
  domain, records the coordinate definition and coefficient prior, and the
  independent artifact auditor checks that metadata and regenerates the
  posterior using the same fixed coordinate definition.
- `generate_independent_validation_sample` now selects a prefix of one
  seeded random permutation. Repeated calls with the same frame seed and
  frame size therefore produce nested validation samples. A regression test
  checks subset inclusion, identical validation-frame truth, and literal
  prefix ordering for n=100, 300, and 1000.

These are implementation-level resolutions in the current branch snapshot;
they are not yet evidence of passing CI because current-head workflow/check
results have not been returned. Also, the lower-level fitter still supports
sample-standardized coordinates when the optional fixed coordinates are
omitted. Repeated-sampling and validation-size experiments must explicitly
supply the fixed domain coordinates; callers that omit them retain the
sample-dependent prior-coordinate behavior.


## Implementation status update

The sample-dependent-prior finding above has been addressed in the independent
validation to H0 propagation path. The fitter now accepts an explicit fixed
coordinate center and scale, and the H0 runner sets them from the predeclared
redshift domain [0.01, 0.40]: center 0.205 and scale 0.195. The same coordinates
are propagated through slope-prior sensitivity, recorded in the run
configuration, and checked by the independent artifact auditor. The
quadrature-ladder API also accepts these explicit coordinates.

The fitter's historical sample-standardization behavior remains the default
for backward compatibility in unrelated analyses. It must not be used with a
fixed diagonal coefficient prior in cross-replicate Bayesian comparisons
unless that sample-dependent prior is itself the declared estimand.

A regression test now fits validation samples with different sizes and
different sample moments under the same fixed coordinate system, and checks
that the same physical coefficient draw induces the same inclusion curve.
This test and the end-to-end artifact-audit test have been committed, but
current-head CI has not yet been observed. Therefore the implementation is
committed, not yet CI-verified.

## H0 convergence-ladder coordinate propagation (2026-11-08)

A final API-path audit found that `run_inclusion_h0_quadrature_convergence`
previously called the posterior fitter without exposing its optional fixed
coordinate arguments. Thus a caller could use fixed coordinates in the main
smoke runner but accidentally fall back to sample-standardized coordinates
when using the standalone downstream convergence utility.

The convergence utility now accepts `coordinate_center` and
`coordinate_scale` and forwards both to every quadrature-order fit. A
regression test supplies the declared domain values (0.205, 0.195) and checks
that every fit in the ladder retains them.

This closes the identified coordinate-propagation gap in that API. The new
regression test is committed but has not yet been observed passing in CI.
