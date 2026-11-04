# Paired host-misspecification sensitivity surface protocol

## Research question

How do fixed-truth error and interval coverage change when the fraction and redshift of an unlisted host population vary, and how much of the result is attributable to reversing the relative weights of the two listed hosts? The earlier paired experiment used one truth, noise scale, and hidden-host configuration. This prespecified factorial expands those dimensions while keeping each model comparison paired on identical observations.

## Factorial design

- Injected H0: 60, 70, 80 km s^-1 Mpc^-1.
- Fractional distance-noise scale at reference redshift z=0.03: 0.05, 0.10, 0.20; sigma = noise_fraction × c × 0.03 / H0_true.
- Hidden-host population fraction: 0, 0.15, 0.30, 0.45.
- Hidden-host redshift: 0.025, 0.030, 0.035.
- Full design: 3 × 3 × 4 × 3 = 108 cells.
- Default: 200 observations per cell (21,600 unique observations), analyzed under four arms (86,400 posterior fits).
- Uniform prior [50, 90], posterior grid of 4,001 points on [25, 125].
- Base population redshifts 0.02 and 0.04 with relative weights 0.7 and 0.3. The hidden component has the cell-specific redshift and fraction f; true population weights are [(1-f)×0.7, (1-f)×0.3, f].
- Non-positive observed distances are redrawn. Each component likelihood is conditioned on positive distance.

## Paired inference arms

1. **base_weights**: two listed hosts weighted 0.7/0.3; hidden component omitted.
2. **reversed_weights**: two listed hosts weighted 0.3/0.7; hidden component omitted.
3. **equal_weights**: two listed hosts weighted 0.5/0.5; hidden component omitted.
4. **hidden_host_modelled**: listed hosts weighted 0.7/0.3 and hidden component assigned the generating fraction and redshift.

All four arms analyze the exact same observation in each replicate. The primary paired contrasts are reversed minus base weights and hidden-omitted base arm minus hidden-modelled arm. The equal-weight arm is a descriptive comparator. At f=0, the base and hidden-modelled likelihoods should coincide up to numerical precision; this is an explicit negative control.

## Outcomes and uncertainty

For each cell and arm, report posterior-median mean error, RMSE, empirical 68%/95% coverage, and mean interval widths. For each primary contrast, report point estimates and percentile-bootstrap 95% intervals, resampling replicate IDs jointly across arms within each cell. The bootstrap seed is fixed and recorded. These intervals describe Monte Carlo uncertainty conditional on each cell; they do not account for model uncertainty or imply the cells are representative of astrophysical conditions.

## Independent audit requirements

The independent auditor must regenerate each cell's host choices and distances from the declared RNG stream, independently reconstruct all posterior summaries from the likelihood equation, reproduce every per-cell arm summary, and regenerate all paired bootstrap contrasts. It must check cell completeness, replicate ordering, shared observations across arms, row counts, finite values, and credible-interval ordering. The f=0 negative control must be checked for posterior equivalence between **base_weights** and **hidden_host_modelled**.

## Scope

This is a toy-model sensitivity study, not real-data validation. It does not include realistic gravitational-wave likelihoods, detector selection, galaxy-catalogue completeness, peculiar velocities, host-population uncertainty beyond the chosen discrete mixture, or a full hierarchical cosmological model. Results must be interpreted as conditional diagnostics, not general claims about dark-siren inference.
