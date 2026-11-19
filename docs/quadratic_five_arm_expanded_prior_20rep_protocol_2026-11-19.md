# Expanded common H0 prior: 20-replicate five-arm pilot protocol

**Frozen before inspecting this pilot's output.** This is a prior-support sensitivity follow-up to the completed 200-replicate [50,90] five-arm pilot. It is not a confirmatory calibration experiment.

## Motivation

The 200-replicate pilot showed that catalogue-only arm A's lower posterior quantiles repeatedly approach the [50,90] prior boundary. A separate expanded-support reconstruction of the earlier 10-replicate smoke also found non-negligible posterior mass below H0=50. The [50,90] support may therefore truncate A's lower tail and distort cross-arm comparisons.

The new experiment uses one common wider support for all five arms. It does not change the galaxy-inclusion model, the GW selection model, the generated events, or the estimators. Its primary purpose is to assess prior-boundary sensitivity under a shared support.

## Frozen configuration

- Master seed: 2026111940; 20 parent/event replicates.
- Parent population: 1,200 galaxies; 30 detected events.
- Generating H0: 70 km s^-1 Mpc^-1.
- Curvature values: gamma = -1, 0, +1, paired within each replicate.
- Validation frame: 500; nested validation sizes 100, 300, 500.
- Five matched arms A-E as in the prior 100-replicate experiment.
- Common H0 prior support: [30,110].
- H0 grid: 1,121 equally spaced points (spacing approximately 0.07143).
- Measurement sigma: 140 Mpc; detection threshold: 250 Mpc.
- Scott bandwidth multiplier: 0.7.
- Inclusion coefficient quadrature order: 5 per coefficient dimension (25 tensor-product nodes).

Planned denominator: 20 × 3 curvature values × 3 validation sizes × 5 arms = **900 lifecycle cells**, grouped into 60 independently audited replicate-by-curvature child runs.

## Primary readout

Report all gamma × validation-size × arm cells:
- bias and RMSE of posterior medians;
- empirical 68% equal-tail coverage with Wilson 95% intervals;
- mean 68% interval width;
- failure-inclusive bounds and explicit lifecycle statuses;
- paired replicate-level contrasts E-D, E-C, D-C, and B-C.

The independent unit is the parent/event replicate. Curvature conditions share parent/event seeds; validation sizes are nested. Do not treat rows across arms, gamma, or validation size as independent.

## Prior-support sensitivity estimand

The saved likelihood curves are evaluated under the same likelihood with two posterior supports:
1. the full [30,110] support;
2. the restricted [50,90] support, obtained by truncating and renormalizing the same likelihood curve.

This is a within-replicate prior-support comparison, not a different data realization. It should report the shift in posterior median and lower quantiles, as well as posterior mass below H0=50 under the wider support. The same support must be applied to every arm in each comparison.

## Interpretation safeguards

- The 20-replicate sample is exploratory and too small for precise coverage certification.
- The wider support is a sensitivity choice, not an assertion that [30,110] is the uniquely correct cosmological prior.
- A posterior that still presses against 30 would require a further support expansion; do not interpret the boundary as resolved merely because the lower edge changed.
- The 1,121-point grid is finer than the earlier 281-point grid but still requires a grid-refinement check on difficult cases.
- The fixed Scott multiplier is not holdout-selected; no KDE optimality claim is allowed.
- Oracle arm B remains diagnostic only.
- No LVK/GWTC/public event catalogue, real galaxy survey, or observational H0 estimate is used.

## Stop/go conditions

A result may be described as a successful computational pilot only if all 900 planned cells have explicit statuses, all 60 child runs pass independent reconstruction, and the top-level orchestration auditor passes. Any failed/non-evaluable cell remains in the planned denominator. A clean computational audit does not establish posterior calibration or estimator superiority.
