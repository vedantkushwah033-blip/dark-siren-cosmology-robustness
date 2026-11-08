# Audit: conditional H0 PIT is not joint SBC

**Status:** methodological review of the newly added multi-replicate calibration runner.  
**Scope:** synthetic low-redshift Euclidean dark-siren benchmark.  
**No calibration experiment is claimed to have run.**

## Finding

The runner's experiment labelled `SBC` currently randomizes the generating H0 from the inference H0 prior, but it holds the inclusion-function truth fixed within each catalogue-floor cell. That is a useful **H0 prior-predictive PIT conditional on a fixed inclusion DGP**, but it is not ordinary joint simulation-based calibration (SBC) for the declared inclusion-coefficient posterior.

The distinction matters because the inference procedure assigns a proper Gaussian prior to the inclusion coefficients, while the data generator does not draw the true coefficients from that prior.

## Exact floor-zero mismatch

The current floor-zero inclusion truth is

```
pi(z) = expit((z_mid - z) / width)
z_mid = 0.22
width = 0.05
```

The inference basis is fixed as

```
x = (z - 0.205) / 0.195
logit pi(z) = theta_0 + theta_1 x
```

Therefore the fixed generating coefficient vector is

```
theta_true = (0.30, -3.90)
```

The inference prior is independent Gaussian with standard deviations (5.0, 2.5) for intercept and slope. The true coefficient vector is not drawn from this prior; in particular, the true slope is fixed at -3.90. The positive-floor cases are additionally misspecified relative to the logistic-polynomial fitting family.

Drawing only H0 from its uniform prior does not turn this into joint SBC over (H0, theta). Nor does the floor-zero correctly specified functional form solve the coefficient-prior mismatch.

## What the current PIT does measure

For each fixed catalogue-floor cell, the runner draws

```
H0_true ~ Uniform(H0_min, H0_max)
```

and computes the fitted posterior CDF at that value. This can be retained as a **conditional H0 prior-predictive PIT diagnostic**, conditional on the fixed inclusion-truth parameters and the rest of the declared synthetic DGP.

Interpret it as a diagnostic of the H0 posterior pipeline under that fixed nuisance truth, not as a standard joint SBC rank test. Uniformity is not guaranteed merely because H0 was drawn from its prior: the nuisance truth is fixed while the inference procedure averages over a different coefficient prior, and the modular/cut posterior is not automatically the full joint posterior.

The fixed-truth repeated-sampling coverage experiment at H0=70 remains a separate frequentist operating-characteristic experiment. It must not be conflated with either PIT or SBC.

## Requirements for a genuine joint SBC study

Before using the label `SBC` for the inclusion-posterior model:

1. Draw H0 truth from the declared H0 prior.
2. Draw the true inclusion coefficients from the declared coefficient prior in the same fixed domain-standardized basis used by inference.
3. Generate both analysis-catalogue retention and independent validation retention from that same coefficient draw, with separate Bernoulli random streams.
4. Use a correctly specified logistic inclusion DGP for this reference arm (no positive floor).
5. Record the true coefficient vector in truth-only artifacts; never expose it to the fitter or likelihood.
6. Define precisely whether the target is the full joint posterior or a modular/cut posterior. Since the current method deliberately prevents GW data from updating inclusion coefficients, ordinary full-joint SBC uniformity must not be assumed for the cut target without a derivation.
7. If validating the cut posterior, predeclare a calibration construction appropriate to that modular target, and independently audit the resulting ranks/PITs.

The existing positive-floor scenarios should remain labelled as controlled misspecification stress tests. They can use conditional H0 PIT and fixed-truth coverage, but they are not a correctly specified joint-SBC reference.

## Recommended reporting terminology now

Until the coefficient-generating process and target posterior are aligned:

- Rename the current `SBC` arm to `H0_PIT_FIXED_INCLUSION_TRUTH` (or an equally explicit label).
- State that H0 is prior-predictive while inclusion truth is held fixed.
- Keep `FIXED_TRUTH` coverage as a separate experiment.
- Do not report the current PIT histogram as joint SBC or claim joint SBC calibration.

## Interpretation boundary

This is a methodological naming and estimand issue, not evidence that the implemented H0 likelihood is numerically wrong. The existing per-replicate runner and independent artifact audit can still support a useful conditional-PIT/coverage pilot once current-head tests pass. A proper joint or modular SBC experiment requires the additional generative and inferential alignment described above.
