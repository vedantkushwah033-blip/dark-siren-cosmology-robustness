# Smooth-selection misspecification experiment protocol

**Status:** implementation and small CI integration test added; confirmatory science run not yet performed.  
**Scope:** synthetic methodology only. This is not a realistic detector-response model and does not estimate the real Hubble constant.

## Question

Under a fixed-truth toy model, how do posterior summaries and interval coverage respond when the inference model uses a logistic event-selection curve whose midpoint differs from the curve used to generate the detected sample, and when the assumed pre-selection missing-host fraction is also misspecified?

## Generative model

The two toy host components have normalized redshift densities on the finite redshift grid:

- catalogue component: Gaussian-shaped density with mean 0.12 and standard deviation 0.05;
- missing-host component: Gaussian-shaped density with mean 0.25 and standard deviation 0.07.

The observed-distance model is Gaussian with mean (d=cz/H_0) and fixed standard deviation 140 Mpc. The selection probability is

[
s(d;d_{50},w)=operatorname{expit}left[-rac{d-d_{50}}{w}ight],
]

with fixed width (w=120) Mpc. The true pre-selection missing-host fraction is 0.30 and the fixed true value is (H_0=70) km s(^{-1}) Mpc(^{-1}).

The pilot crosses true and assumed selection midpoints (d_{50}in{700,1000,1300}) Mpc and assumed missing-host fractions (f_{mathrm{model}}in{0.1,0.3,0.5}). The matched selection midpoint and (f_{mathrm{model}}=0.3) form the correctly specified control.

## Likelihood contract

For each component (kin{C,M}),

[
q_k(dmid H_0)=int p(dmid z,H_0)g_k(z),dz,
]

and

[
alpha_k(H_0;d_{50},w)
=int g_k(z)left[int p(dmid z,H_0)s(d;d_{50},w),ddight]dz.
]

For the pre-selection mixture fraction (f), the fixed-count detected-event likelihood is

[
p(dmid H_0,f,mathrm{selected})
=
rac{s(d;d_{50},w)left[(1-f)q_C(dmid H_0)+fq_M(dmid H_0)ight]}
{(1-f)alpha_C(H_0;d_{50},w)+falpha_M(H_0;d_{50},w)}.
]

The implementation evaluates the inner expectation in the selection efficiency using Gauss-Hermite quadrature and the redshift integral using trapezoidal quadrature. A separate audit recomputes the selection efficiencies, likelihood grids, posterior summaries, and aggregate metrics without importing the production likelihood or runner.

## Sampling requirements

The generator must sample component identity using the **selected** mixture probability, sample redshift from the component's selected marginal, and then sample distance conditional on acceptance at that fixed redshift. Redrawing redshift after a rejected distance draw would apply the selection weighting twice and is prohibited.

The implementation should preserve latent component and true-redshift columns only as audit metadata. They must not enter the inference likelihood.

## Endpoints and diagnostics

For every true-midpoint × assumed-midpoint × assumed-fraction cell, report:

- bias and RMSE of posterior median relative to the fixed truth;
- empirical 68% and 95% interval coverage;
- mean 68% interval width;
- replicate count and Monte Carlo uncertainty;
- matched-versus-mismatched contrasts using the same generated dataset where applicable.

These are fixed-truth repeated-sampling diagnostics. They are not SBC results and do not establish model adequacy for real gravitational-wave data.

## Execution

A small plumbing run is exercised in CI. For a pilot outside CI, from a clean checkout:

```bash
python scripts/run_smooth_selection_misspecification.py \
  --replicates 50 --n-events 20 --seed 20261029 \
  --h0-points 121 --z-points 301 --output-dir artifacts/smooth_selection_pilot

python scripts/audit_smooth_selection_misspecification.py artifacts/smooth_selection_pilot
```

Do not call this a confirmatory run merely because it has 50 replicates. Inspect runtime, Monte Carlo uncertainty, and numerical stability first. Increase the replicate count only after defining a precision target and freezing the confirmatory configuration.

## Required gates before interpreting results

1. All supported-Python CI jobs pass.
2. The independent artifact audit returns PASS and reconstructs every posterior fit.
3. Selection-efficiency quadrature is checked against a higher-order calculation.
4. H0 and redshift grid refinement is checked for representative matched and mismatched cells.
5. The matched model is evaluated as a negative control.
6. Monte Carlo uncertainty is reported for coverage and paired contrasts.
7. All cells are retained and reported; no favorable cell is selected post hoc.
8. The results are described as synthetic sensitivity analysis, not a real detector or observational H0 result.
