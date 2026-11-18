# Paired quadratic-logit H0 calibration pilot — design note (2026-11-19)

## Question

How does finite validation size interact with quadratic misspecification in the catalogue-inclusion model when the resulting inverse-probability correction is propagated through the selected-event H0 likelihood?

This pilot is a software-validated, fixed-truth repeated-sampling experiment. It does not estimate the real Universe's H0.

## Design

- Generating H0: 70 km s^-1 Mpc^-1, fixed across all cells.
- Parent population: 1,200 low-redshift synthetic galaxies on [0.01, 0.40].
- GW process: four detected events per replicate; generated independently of catalogue retention. The design conditions on this fixed detected-event count; it does not model count-rate information or uncertainty in the number of detections.
- Catalogue inclusion: parent-measure-calibrated quadratic-logit law, beta=-3.9, gamma in {-1, 0, +1}.
- Inference: degree-one logistic inclusion model in all cells. Gamma=0 is the correctly specified reference; gamma=-1 and +1 are deliberate functional-form misspecification cases.
- Validation frame: 500 objects; nested SRS sizes 100, 300, 500.
- Replication: 10 independent replicate streams, shared across gamma and validation sizes for paired contrasts.
- H0 prior/grid: uniform [20, 300], 281 grid points.
- The 281-point H0 grid has spacing 1 km s^-1 Mpc^-1. The quadrature ladder tests inclusion-coefficient integration stability only; it does not establish H0-grid resolution convergence. A finer-grid posterior reconstruction remains a separate required check.
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
- Results are conditional on four detected events per replicate and do not include information from the GW detection count/rate.
- The 281-point H0 grid spacing is 1 km s^-1 Mpc^-1; coefficient-quadrature convergence is not an H0-grid refinement test.
- A converged quadrature ladder establishes numerical order stability only.
- The experiment does not establish estimator superiority, a universal IPW result, or the value of any method on real survey data.
- No LVK/GWTC data, real galaxy survey, or observational H0 measurement is used.



## Exact-head execution and audit (2026-11-19)

Workflow #44410603590 completed successfully on commit `c88ba40f8c23365e7ba32380c81dd86b335e9a68`. Artifact `14477511669` contains the complete run.

- 90/90 planned cells are SUCCESS; 90/90 child audits PASS.
- All 2,114 aggregate-manifest entries were independently SHA-256 checked against the downloaded artifact: no missing files and no hash mismatches.
- The aggregate auditor reports PASS, verifies 30 paired nested-validation blocks, and independently recomputes summary metrics.
- All 90 coefficient-quadrature ladders converged: 72 at order 7 and 18 at order 9.
- All cells condition on four detected events and use the 281-point H0 grid [20, 300] (1 km s^-1 Mpc^-1 spacing).

Descriptive pilot readout (10 replicates per cell; not confirmatory):

| gamma | validation n | mean posterior-median bias | RMSE | 68% coverage | median minimum catalogue Kish ESS |
|---:|---:|---:|---:|---:|---:|
| -1 | 100 | +3.95 | 10.29 | 7/10 | 14.49 |
| -1 | 300 | +3.95 | 10.21 | 7/10 | 37.57 |
| -1 | 500 | +3.46 | 10.23 | 7/10 | 63.83 |
| 0 | 100 | +8.21 | 11.51 | 7/10 | 10.04 |
| 0 | 300 | +8.22 | 11.72 | 8/10 | 20.43 |
| 0 | 500 | +7.50 | 11.38 | 9/10 | 30.49 |
| +1 | 100 | +13.30 | 16.40 | 5/10 | 9.67 |
| +1 | 300 | +13.63 | 16.73 | 5/10 | 16.28 |
| +1 | 500 | +12.87 | 16.07 | 5/10 | 24.61 |

The gamma=0 reference itself has substantial positive bias. Therefore the increase in bias at gamma=+1 cannot be interpreted as the full effect of quadratic misspecification: baseline modular-prior mismatch, four-event information content, and other features of the toy likelihood may contribute. Validation size improves weight ESS markedly but does not materially remove the observed H0 bias in this pilot. With only 10 replicates, coverage estimates have wide Wilson intervals and are unstable.

## Go/no-go

The artifact-integrity and execution gate is PASS. The scientific gate remains exploratory only. Before a confirmatory interpretation, the project still needs (i) a finer H0-grid reconstruction, (ii) a larger event-count sensitivity, (iii) a decomposition of the gamma=0 baseline bias including coefficient-prior sensitivity, and (iv) the matched five-arm A–E repeated-sampling comparison. This single-arm pilot cannot establish whether estimated inclusion correction improves H0 inference relative to catalogue-only, oracle, or known-probability controls.

