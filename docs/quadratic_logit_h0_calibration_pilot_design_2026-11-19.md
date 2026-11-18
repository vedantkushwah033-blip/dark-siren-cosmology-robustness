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
- The 281-point H0 grid has spacing 1 km s^-1 Mpc^-1. The quadrature ladder tests inclusion-coefficient integration stability only; it does not establish H0-grid resolution convergence. A separate 54-child grid-refinement diagnostic subsequently found 281-to-561 posterior-density and quantile tolerances passed in all 27 comparisons from the first three paired replicates; the full grid-refinement workflow is still running, so this is a partial numerical check rather than a full 90-cell grid audit.
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
- The 281-point H0 grid spacing is 1 km s^-1 Mpc^-1; coefficient-quadrature convergence is not an H0-grid refinement test. A separate first-three-replicate check found 281-to-561 stability in 27/27 comparisons, but full grid-refinement CI remains pending.
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



### Coefficient-prior scale sensitivity (same pilot artifact)

The saved within-cell sensitivity analysis changes the slope-prior standard deviation by multipliers 0.5 and 2.0 while holding the observed catalogue, events, selection rule, H0 grid, and base weights fixed. Across the 10 paired replicates in the gamma=0 reference, the mean posterior-median shift relative to multiplier 1.0 was -3.06 km s^-1 Mpc^-1 at validation n=100 and -0.72 at n=500 for multiplier 0.5; the corresponding shifts for multiplier 2.0 were +1.11 and +0.20. Thus this prior-scale perturbation changes the posterior, particularly at the smallest validation size, but does not by itself erase the +7.5 to +8.2 baseline mean median bias. This is a limited sensitivity check, not a complete prior robustness analysis: the prior mean, intercept scale, event count, and likelihood model were not varied, and the 10-replicate uncertainty remains substantial.



### Same-realization oracle decomposition on the 281-point grid

I independently recomputed the known-inclusion-probability and full-parent oracle posteriors directly from the current 281-point artifact, using the saved observable event distances, catalogue, parent truth, true inclusion probabilities, and the same selected-event likelihood with the one-sided detection-efficiency normalization. The calculation does not call the production shared-inclusion likelihood or the PR #89 post-hoc script. As a cross-check, the 141-point predecessor artifact and this 281-point artifact have byte-identical per-cell observable and truth CSVs (630 files compared, zero mismatches). The oracle posterior was recomputed on the current 281-point H0 grid, so the following table is on the same numerical grid as the fitted-IPW pilot.

Mean posterior-median bias in km s^-1 Mpc^-1; 10 paired replicates per cell:

| gamma | validation n | estimated IPW | known-pi IPW | full-parent oracle | estimated minus known-pi |
|---:|---:|---:|---:|---:|---:|
| -1 | 100 | +3.95 | +6.45 | +7.58 | -2.50 |
| -1 | 300 | +3.95 | +6.45 | +7.58 | -2.50 |
| -1 | 500 | +3.46 | +6.45 | +7.58 | -2.99 |
| 0 | 100 | +8.21 | +7.43 | +7.58 | +0.79 |
| 0 | 300 | +8.22 | +7.43 | +7.58 | +0.79 |
| 0 | 500 | +7.50 | +7.43 | +7.58 | +0.07 |
| +1 | 100 | +13.30 | +8.84 | +7.58 | +4.46 |
| +1 | 300 | +13.63 | +8.84 | +7.58 | +4.79 |
| +1 | 500 | +12.87 | +8.84 | +7.58 | +4.03 |

The full-parent oracle has mean median bias +7.58, Monte Carlo SE 3.15, and paired-t 95% interval [+0.46,+14.71] across the 10 shared event realizations; its 68% interval covers truth in 6/10 replicates. This directly shows that the large positive baseline bias is already present before catalogue thinning or estimated inclusion correction. At gamma=0, known-pi IPW bias is +7.43 and estimated-IPW bias ranges from +8.21 to +7.50; estimated-minus-known-pi shifts are small and their paired 95% intervals include zero. Thus the baseline bias should not be attributed mainly to the fitted inclusion prior or inclusion-model estimation.

Known-pi IPW minus full-parent median contrasts are -1.13 (paired 95% t interval [-5.30,+3.03]) at gamma=-1, -0.16 ([-2.98,+2.66]) at gamma=0, and +1.25 ([-0.16,+2.67]) at gamma=+1. These are uncertain finite-catalogue contrasts. Estimated-minus-known-pi contrasts show a clear curvature-dependent pattern in this small pilot: at gamma=-1 they are -2.50 at n=100 and -2.99 at n=500 (paired intervals [-4.27,-0.72] and [-4.74,-1.24]); at gamma=0 they are +0.79 and +0.07 (both intervals include zero); at gamma=+1 they are +4.46 and +4.03 (intervals [+2.38,+6.53] and [+2.66,+5.39]). The direction is consistent with curvature-dependent degree-one misspecification, but 10 replicates are not enough for a definitive conclusion.

This changes the next scientific priority: investigate the finite four-event realization and the selected-H0 likelihood/selection convention, then quantify how the estimated inclusion model adds to that baseline. The 100-replicate follow-up improves precision for the fitted-IPW arm but does not itself include oracle arms; an independent full-parent and known-pi decomposition on its 281-point raw artifact is still needed. This remains a synthetic fixed-truth diagnostic, not SBC or observational inference.

## Go/no-go

The artifact-integrity and execution gate is PASS. The scientific gate remains exploratory only. Before a confirmatory interpretation, the project still needs (i) completion of the full 281-to-561 H0-grid audit beyond the first three paired replicates, (ii) a larger event-count sensitivity, (iii) a larger-replicate full-parent/known-inclusion oracle reconstruction to quantify uncertainty in the gamma=0 baseline decomposition, and (iv) the matched five-arm A–E repeated-sampling comparison. This single-arm pilot cannot establish whether estimated inclusion correction improves H0 inference relative to catalogue-only, oracle, or known-probability controls.

