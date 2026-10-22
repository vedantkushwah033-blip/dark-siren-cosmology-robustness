# Calibration experiment protocol

## Purpose

This experiment is a small diagnostic of the simplified single-host inference model. It checks repeated-sampling coverage for equal-tail posterior intervals when the true H0 and host redshift are fixed.

It is not simulation-based calibration (SBC), a realistic dark-siren population analysis, or a reproduction of the original project's multi-host coverage result.

## Fixed design

- True H0: 70 km s^-1 Mpc^-1
- Host redshift: 0.02
- Distance model: d_L = c z / H0
- Prior: uniform over H0 = 50–90 km s^-1 Mpc^-1
- Interval: central/equal-tail 68% posterior interval
- Default replicates: 500
- Default random seed: 20261022
- Default Gaussian distance uncertainty: 10 Mpc

The script draws a Gaussian observed distance and resamples any non-positive draw. This makes the observation generator technically a Gaussian truncated at zero. With the default parameters, negative draws should be extremely rare, but the detail is documented because the generator and likelihood should be consistent in a precision calibration study.

## Primary outputs

- empirical coverage of the nominal 68% interval;
- Monte Carlo standard error for the estimated coverage;
- bias and RMSE of the posterior-grid mode;
- mean equal-tail interval width.

For a rough diagnostic, the Monte Carlo standard error is sqrt(p_hat * (1 - p_hat) / n), where p_hat is observed coverage and n is the number of independent replicates. Do not interpret a small difference from 0.68 as meaningful without considering Monte Carlo uncertainty.

## Interpretation rules

1. Do not expect this one-host experiment to reproduce the original multi-host project's reported 94.35% coverage.
2. Coverage far above 68% indicates conservative intervals under this design; coverage far below 68% indicates undercoverage. Neither outcome alone identifies the cause.
3. Inspect interval construction, prior boundaries, grid resolution, likelihood assumptions, and generator/likelihood consistency before modifying the model.
4. Report this experiment separately from any experiment involving catalogue incompleteness, host-weight misspecification, or GW detection selection.
5. A successful unit-test workflow means the tested code paths passed. It does not prove scientific validity or realistic modelling.

## Reproducible execution

Run the GitHub Actions workflow named `single-host-calibration` manually with explicit values for replicates, seed, and distance uncertainty. Download the `single-host-calibration-results` artifact and retain it with the run URL and commit SHA. Do not copy numerical results into the paper until the artifact has been inspected and the experiment has been independently checked.
