# Paired host-misspecification experiment protocol

## Question

When host-weight misspecification and omission of a possible host are evaluated on different random datasets, observed metric differences combine model effects with Monte Carlo sampling variation. This experiment removes that avoidable confounding by applying every inference arm to the exact same simulated distance observations.

## Fixed design

- Synthetic population redshifts: 0.02, 0.04, and 0.03.
- Population probabilities in that order: 0.35, 0.35, and 0.30.
- Injected H0: 70 km s^-1 Mpc^-1.
- Gaussian distance noise: sigma = 0.10 c (0.03) / 70 Mpc; non-positive observations are redrawn.
- Prior: uniform H0 on [50, 90].
- Posterior grid: 4,001 points on [25, 125].
- Default seed: 20261027; default 1,000 independently generated observations.
- Four inference arms per observation: correct two-host weights [0.7, 0.3] while omitting z=0.03; reversed two-host weights [0.3, 0.7] while omitting z=0.03; equal two-host weights [0.5, 0.5] while omitting z=0.03; and equal observed-host weights plus a 0.30 hidden-host component at z=0.03.
- Every Gaussian component likelihood is conditioned on positive observed distance, matching the rejection-sampling data generator.

The “correct weights” label means correct relative weights for the two listed hosts only; because the generating population also includes z=0.03, that arm still omits a real population component. It is not a fully specified oracle model. The hidden-host-modelled arm uses the generating host mixture.

## Estimands

The raw file records the same observed distance and true host for all four arms per replicate. For each arm, report median error, RMSE, empirical 68%/95% coverage, and mean credible-interval widths. The paired contrast file reports reversed-minus-correct two-host-weight contrasts and omitted-minus-modelled-hidden-host contrasts for bias, RMSE, coverage, and interval width, with percentile bootstrap 95% intervals resampling replicate IDs jointly across the compared arms.

The bootstrap interval quantifies Monte Carlo uncertainty conditional on this one fixed toy design. It does not include model uncertainty, variation across the injected truth/noise design, or astrophysical uncertainty. Contrasts in coverage are differences in proportions (multiply by 100 for percentage points).

## Independent audit

The auditor must independently regenerate the seeded observations, reconstruct every posterior from the likelihood equation without importing production inference or runner functions, reconstruct each summary and paired bootstrap interval, and reject changes in data, model outputs, row counts, pair alignment, or summaries.

## Interpretation boundaries

This is a controlled synthetic diagnostic, not a real gravitational-wave/cosmology analysis. One H0, one noise scale, one prior, one seed, and one host-redshift configuration cannot establish general robustness. The hidden-host-modelled arm is not a full hierarchical dark-siren model and does not include detector selection or catalogue incompleteness. The experiment's contrasts are paired and more statistically efficient than unpaired comparisons, but they do not establish universal scientific importance.
