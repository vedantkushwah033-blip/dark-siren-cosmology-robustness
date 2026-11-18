# Paired quadratic-logit H0 calibration pilot — design note (2026-11-19)

## Question

How does finite validation size interact with quadratic misspecification in the catalogue-inclusion model when the resulting inverse-probability correction is propagated through the selected-event H0 likelihood?

This pilot is a software-validated, fixed-truth repeated-sampling experiment. It does not estimate the real Universe's H0.

## Design

- Generating H0: 70 km s^-1 Mpc^-1, fixed across all cells.
- Parent population: 1,200 low-redshift synthetic galaxies on [0.01, 0.40].
- GW process: four detected events per replicate; generated independently of catalogue retention.
- Catalogue inclusion: parent-measure-calibrated quadratic-logit law, beta=-3.9, gamma in {-1, 0, +1}.
- Inference: degree-one logistic inclusion model in all cells. Gamma=0 is the correctly specified reference; gamma=-1 and +1 are deliberate functional-form misspecification cases.
- Validation frame: 500 objects; nested SRS sizes 100, 300, 500.
- Replication: 10 independent replicate streams, shared across gamma and validation sizes for paired contrasts.
- H0 prior/grid: uniform [20, 300], 141 grid points.
- The 141-point H0 grid has spacing 2 km s^-1 Mpc^-1. The quadrature ladder tests inclusion-coefficient integration stability only; this pilot does not independently establish H0-grid resolution convergence. Posterior quantiles and coverage therefore retain a grid-discretization limitation.
- Inclusion-coefficient quadrature: orders 3, 5, 7, 9; convergence required before a cell is counted SUCCESS.

The inference coefficient prior is independent Normal with mean zero and standard deviations (5.0, 2.5) for intercept and slope. The fixed generating coefficients are not drawn from this prior (at gamma=0 the generating slope is -3.9). Consequently, the fixed-truth coverage arm is a conditional repeated-sampling diagnostic for this modular procedure, not prior-predictive or joint SBC.

The analysis-parent redshift realization and GW events are held fixed across gamma within replicate. The catalogue retention realization changes with gamma, as required by the DGP. Within each gamma/replicate, the catalogue, GW events, validation frame, and validation permutation are held fixed across validation sizes; only the prefix length changes.

## Primary summaries

For each gamma × validation-size cell, report:
- posterior-median bias relative to H0=70;
- RMSE of posterior medians;
- equal-tail 68%, 90%, and 95% interval coverage;
- 68% interval width;
- minimum catalogue Kish ESS and maximum normalized IPW weight;
- evaluability and failure counts.

Coverage is computed among SUCCESS cells. The failure ledger is reported alongside coverage; non-evaluable and failed cells must not silently disappear. Ten replicates yield wide Monte Carlo uncertainty, so coverage estimates are descriptive pilot quantities only.

## Validation and provenance

Each replicate is passed through the independent raw-artifact auditor, which reconstructs the quadratic-logit probabilities, catalogue retention, validation sample, GW generation trace, selected likelihood, detection efficiency, posterior normalization, and prior sensitivity. The aggregate auditor checks:
- exact planned-cell/ledger reconciliation;
- child-audit attestations and configuration;
- aggregate manifest hashes;
- same parent redshift realization and observed GW events across gamma;
- same analysis catalogue/events and validation frame across nested sizes;
- literal nested validation permutation prefixes;
- summary counts against the ledger.

The raw auditor regenerates inclusion-coefficient quadrature with the production fitter; its independent reconstruction focuses on the event-mixture likelihood and one-sided detection-efficiency calculation. This is not a fully independent second implementation of the coefficient posterior.

## Interpretation boundaries

- This is fixed-truth repeated-sampling coverage, not joint SBC.
- Gamma=0 is the correctly specified degree-one logistic reference; nonzero gamma values test this one quadratic omission pattern, not all misspecification.
- Ten replicates are insufficient for precise coverage claims or a confirmatory comparison.
- A converged quadrature ladder establishes numerical order stability only.
- The experiment does not establish estimator superiority, a universal IPW result, or the value of any method on real survey data.
- No LVK/GWTC data, real galaxy survey, or observational H0 measurement is used.

## Go/no-go

Use this 90-cell run to check runtime, failure rates, posterior-grid behavior, and whether the paired design behaves as intended. Do not scale to a larger confirmatory experiment until the artifacts pass the aggregate audit and the observed Monte Carlo uncertainty and failure rates justify it.
