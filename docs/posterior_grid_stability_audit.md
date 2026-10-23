# Posterior-grid stability audit

## Scope

This is a code-level audit of the pinned synthetic inference implementation at commit
`b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f` in
`vedantkushwah033-blip/bayesian-dark-siren-cosmology`. It is not a report of completed
grid-convergence experiments and does not establish the scientific performance of the model.

## Current numerical setup

The repeated-simulation routine constructs
`np.linspace(50.0, 90.0, 1601)`, giving a uniform grid spacing of 0.025 km/s/Mpc.
The posterior is normalized by trapezoidal integration over that grid. Posterior means,
variances, and the CDF used for quantiles are also calculated numerically on the grid.

The grid and integration choices therefore affect the numerical summaries, especially when
posterior probability approaches either boundary or becomes narrower than the grid can resolve.

## Required stability checks

For a fixed set of simulated observations, recompute the posterior with:

1. Baseline: 50–90 km/s/Mpc, 1601 grid points.
2. Wider range: 30–110 km/s/Mpc, with the same grid spacing.
3. Finer baseline grid: 50–90 km/s/Mpc, 3201 points.
4. Wider and finer grid: 30–110 km/s/Mpc, 6401 points.

Keep the observed distance, uncertainty, redshifts, catalogue weights, missing-host inputs,
and prior assumptions unchanged. A grid comparison that also changes the simulated data
cannot isolate numerical-grid effects.

Record for each case:

- trapezoidal posterior integral before normalization;
- posterior mass in a small band near each grid boundary;
- posterior mean, median, standard deviation, and 68%/95% interval endpoints;
- absolute differences in median and interval endpoints relative to the wider/finer reference.

## Interpretation rules

- A normalized posterior is not necessarily grid-stable: normalization can succeed even when
  material posterior mass is truncated by the selected range.
- A stable median alone is insufficient if interval endpoints or boundary mass change materially.
- Tolerances should be declared before reviewing the results, based on the precision needed for
  the scientific contrasts.
- If the wider range changes results materially, inspect the prior and likelihood before
  selecting a replacement grid. Do not simply widen the range until results look favorable.
- Report any excluded or failed replicates; do not silently discard them.

## Separate code-validation issue

The pinned `posterior_h0` routine assumes that the H0 grid and likelihood inputs are suitable
for numerical integration. Before treating it as a general-purpose function, add explicit checks
for a one-dimensional, finite, strictly increasing H0 grid; compatible shapes; finite,
non-negative redshifts and weights; and a finite, positive GW distance uncertainty. These checks
are separate from grid convergence and should be tested independently.

## Status

This document specifies an audit protocol only. The grid variants above have not yet been run,
and no claim is made that the current grid is inadequate. The main factorial experiment remains
a synthetic study and its scientific outputs must be inspected before interpretation.
