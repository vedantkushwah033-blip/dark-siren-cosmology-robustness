# Hidden-host fraction/redshift assumption sensitivity protocol

## Question

The previous hidden-host factorial used the true hidden-population fraction and redshift in its modelled-host inference arm. That arm is an oracle-style benchmark. This experiment asks how inference changes when those hidden-host assumptions are wrong, including when the true hidden fraction is zero but the analyst adds a spurious hidden component.

## Frozen design

- Injected H0: 55, 70, 85 km s^-1 Mpc^-1.
- Fractional distance-noise scale: 0.05 and 0.20 × c × 0.03/H0.
- True hidden-host fraction: 0, 0.15, 0.30, 0.45.
- True hidden-host redshift: 0.025, 0.035.
- Assumed hidden-host fractions: 0.15, 0.30, 0.45.
- Assumed hidden-host redshifts: 0.025, 0.030, 0.035.
- Inference arms: hidden component omitted; oracle hidden component using true fraction/redshift; all nine fixed assumed-fraction × assumed-redshift combinations.
- 48 true-data cells × 200 replicates = 9,600 shared observations; 11 inference arms per observation = 105,600 posterior fits.
- Uniform H0 prior [50,90], 4,001-point grid from 25 to 125.
- Positive-truncated Gaussian distance likelihood, matching the generator's redraw of nonpositive observations.
- Data seed 20261032; bootstrap seed 20261033; 1,000 paired bootstrap resamples per true-data cell.

## Estimands

For each truth/noise/fraction/redshift cell, calculate posterior-median bias, RMSE, 68%/95% coverage, and mean interval widths. Compare every fixed assumed hidden-component model with the oracle model using paired replicate-level bootstrap intervals. Also compare the omitted-host arm with the oracle. A positive RMSE contrast (assumed minus oracle) means the assumed model performs worse on that metric; negative coverage contrasts mean lower empirical coverage.

## Negative controls and audit

- When the true hidden fraction is zero, the oracle hidden-host model must exactly coincide with the omitted-host model for every replicate.
- Every inference arm is evaluated on the same observation in each replicate.
- The independent auditor replays seeded data generation, computes posteriors from a separately implemented likelihood equation, and reconstructs all summaries and paired-bootstrap intervals.
- CI verifies exact row counts, finite values, interval ordering, and the independent audit. An adversarial test corrupts a contrast and requires audit failure.

## Guardrails

This is a low-redshift, discrete-host Gaussian-distance toy model. Fixed assumed fraction/redshift combinations are deliberately misspecified stress cases, not a recommendation for real analyses. The oracle arm is not deployable without knowing hidden population parameters. No real GW/catalogue data, detector-selection effects, peculiar velocities, or full hierarchical population inference are included. Do not interpret the results as an observed H0 estimate or as proof of astrophysical validity.
