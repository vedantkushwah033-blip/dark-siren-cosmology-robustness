# Complete 90-cell quadratic-logit H0 calibration pilot: independent artifact readout

Date: 2026-11-19

## Executive result

A complete 90-cell artifact from the earlier pilot workflow has been recovered and independently checked. All 90 planned cells are present and marked SUCCESS; all 90 child audits report PASS; all 90 coefficient-quadrature runs report CONVERGED.

The original GitHub Actions run ended with **failure at the aggregate orchestration audit**, not during the 90 child simulations. The failure was caused by an invalid audit comparison: it required the validation-frame `selected_for_validation` indicator to be identical for validation sizes 100, 300, and 500. That indicator must change with sample size. The underlying validation frame itself is shared. The corrected audit now compares fixed frame IDs/redshifts and separately checks nested validation-ID prefixes.

This is a complete computational pilot, but it is still only a 10-replicate, four-event, fixed-truth synthetic calibration study. It is not confirmatory evidence, joint SBC, an observational H0 measurement, or a real GW-catalogue analysis.

## Independent corrected rerun

After fixing the orchestration auditor, the full 90-cell design was rerun on workflow #44408479569 at commit `3509ddcffec9a7807bc3867fe2553f6df90a5e7a`. The rerun completed SUCCESS, including the corrected independent orchestration audit and the complete-pilot assertion.

- Rerun artifact ID: 14476042847
- Rerun archive SHA-256: `cadf1a8ad7c9d2052a73972338fef493cbbcb2b85b2da89527c7e2fe5e0805ba`
- Rerun artifact: 90/90 SUCCESS, 90 child audits PASS, 90 quadrature runs CONVERGED.
- Aggregate audit report: PASS; 2,110 manifest entries verified; 30 nested validation blocks and 30 cross-curvature validation-pairing blocks checked; all summary metrics recomputed.
- Independent comparison of the original and rerun ledgers found identical cell IDs and all key numerical ledger metrics (H0 posterior medians/quantiles, interval widths, ESS, maximum normalized weights and median biases) within absolute tolerance (10^{-12}). Some raw CSV values differ at machine-roundoff scale (roughly (10^{-13})); this does not change the reported ledger results.

The first complete artifact below is retained as the original result; the successful rerun independently confirms that the corrected audit accepts the same 90-cell design and that the numerical summaries are reproducible.

## Artifact provenance

- Source workflow run: 44402834830
- Source branch: `research/quadratic-logit-h0-calibration-pilot-20261119`
- Source commit: `1ab60e2e7b7d446eb15f9b59efb3003f25965087`
- Artifact ID: 14473894295
- Artifact name: `quadratic-logit-h0-calibration-pilot-py3.12`
- Artifact archive SHA-256: `e22ae91eff8e0d1aa02f4d03699543df096518cf7d21cbe751ecee9ea66c5426`
- The downloaded archive SHA-256 was independently recomputed and matched the GitHub artifact digest.
- Archive contents: 2,201 files, 101,161,251 uncompressed bytes.
- Aggregate manifest: 2,110 entries; all paths and SHA-256 hashes verified.
- Child manifests: 90; all 1,656 listed child-file paths and SHA-256 hashes verified.
- Child audit reports: 90/90 PASS.
- Coefficient quadrature: 90/90 CONVERGED (72 at order 7, 18 at order 9).
- Planned cells: 90; ledger rows: 90; statuses: 90 SUCCESS, zero NON_EVALUABLE, zero FAILED, zero NOT_RUN.

The original aggregate auditor did not pass. Its failure was reproduced from the workflow log and traced to the size-dependent validation-selection indicator comparison. A separate replay using an intermediate auditor revision also exposed that the historical aggregate summary omitted an explicit `not_run` column. The artifact itself is complete: all 90 planned rows are SUCCESS, so NOT_RUN is zero. The current pilot runner writes the `not_run` count explicitly; the corrected auditor also distinguishes frame identity from sample-size selection.

## Design

| Component | Setting |
|---|---|
| Parent galaxies | 1,200 |
| Detected GW events | 4, fixed per replicate |
| Curvature values | -1, 0, +1 |
| Validation frame | 500 |
| Nested SRS validation sizes | 100, 300, 500 |
| Replicates | 10 per cell |
| Total cells | 90 |
| H0 truth | 70 km s⁻¹ Mpc⁻¹ |
| H0 prior/grid | Uniform [20, 300], 141 points; spacing 2 |
| Inclusion DGP | Quadratic-logit |
| Fitted inclusion model | Degree-one logistic |
| Prior coefficient SDs | [5, 2.5] |
| Curvature coefficient | beta=-3.9; gamma in {-1,0,+1} |
| Parent-average inclusion target | 0.5101919962291073 |
| Fixed event noise | Gaussian sigma=140 Mpc; detection threshold=250 Mpc |

At gamma=0 the degree-one fitted inclusion model has the correct functional form. At gamma=-1 and +1 the fitted degree-one form is deliberately misspecified. The number of detected events is fixed at four, so this experiment conditions on event count and does not assess GW event-rate/count calibration.

## Independent verification

### Artifact integrity
All aggregate and child manifest entries were verified against the downloaded artifact bytes. The archive digest matched GitHub's artifact digest exactly. All 90 child audit reports passed, and all 90 quadrature convergence reports were CONVERGED.

### Configuration checks
Across all 90 cells, the stored configuration agrees with the planned and ledger metadata:
- parent size 1,200 and H0 truth 70;
- validation frame size 500 and the cell's validation sample size;
- gamma matches the planned cell and beta is -3.9;
- degree-one inference and 141-point H0 grid.

### Paired-design checks
Independent cross-cell reconstruction across all 10 replicates found zero discrepancies:
1. Analysis-parent object IDs and redshifts are invariant across gamma and validation size.
2. Observed GW-event tables are invariant across gamma and validation size.
3. Validation-frame IDs and redshifts are invariant.
4. At fixed replicate and validation size, sampled validation IDs and redshifts are identical across gamma.
5. Within each replicate and gamma, validation samples are nested prefixes: 100 ⊂ 300 ⊂ 500.

Retention outcomes are expected to differ with gamma; they are not part of the shared-sampling-design equality check.

### Aggregate summary reconstruction
Every reported summary row was recomputed directly from the 90-row ledger. There were zero discrepancies for:
- planned and success counts;
- NON_EVALUABLE and FAILED counts;
- mean posterior-median bias;
- RMSE of posterior-median error;
- empirical 68/90/95% interval coverage;
- mean central 68% interval width;
- median minimum catalogue Kish ESS;
- median maximum normalized catalogue weight.

The historical summary file did not contain a `not_run` column. Because the full ledger contains all 90 planned cells and every cell is SUCCESS, the implied NOT_RUN count is exactly zero. This historical omission is an audit-schema defect, not a missing simulation cell.

## H0 calibration results

The outcome `median_bias` is posterior median minus the generating H0 value. RMSE is computed across the 10 replicate-level posterior-median errors. Coverage is the empirical fraction of posterior intervals containing H0 truth. Wilson 95% intervals are reported for 68% coverage to show the large uncertainty at n=10.

| γ | Validation n | N | Mean median bias | 95% t interval for mean bias | RMSE | 68% coverage (Wilson 95% interval) | 90% coverage | 95% coverage | Mean 68% width |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -1 | 100 | 10 | +3.95 | [-3.21, +11.11] | 10.29 | 7/10 = 0.70 [0.397, 0.892] | 0.90 | 1.00 | 24.78 |
| -1 | 300 | 10 | +3.95 | [-3.15, +11.05] | 10.21 | 7/10 = 0.70 [0.397, 0.892] | 0.90 | 1.00 | 24.33 |
| -1 | 500 | 10 | +3.46 | [-3.80, +10.72] | 10.23 | 7/10 = 0.70 [0.397, 0.892] | 0.90 | 1.00 | 24.21 |
| 0 | 100 | 10 | +8.21 | [+2.13, +14.29] | 11.51 | 7/10 = 0.70 [0.397, 0.892] | 1.00 | 1.00 | 25.22 |
| 0 | 300 | 10 | +8.22 | [+1.92, +14.52] | 11.72 | 8/10 = 0.80 [0.490, 0.943] | 0.90 | 1.00 | 24.69 |
| 0 | 500 | 10 | +7.50 | [+1.05, +13.95] | 11.38 | 9/10 = 0.90 [0.596, 0.982] | 0.90 | 1.00 | 24.74 |
| +1 | 100 | 10 | +13.29 | [+6.04, +20.54] | 16.40 | 5/10 = 0.50 [0.237, 0.763] | 0.90 | 0.90 | 25.48 |
| +1 | 300 | 10 | +13.62 | [+6.31, +20.94] | 16.72 | 5/10 = 0.50 [0.237, 0.763] | 0.80 | 1.00 | 24.76 |
| +1 | 500 | 10 | +12.86 | [+5.60, +20.12] | 16.07 | 5/10 = 0.50 [0.237, 0.763] | 0.80 | 1.00 | 24.97 |

All bias, RMSE and interval-width values are in km s⁻¹ Mpc⁻¹.

### Paired curvature contrasts

For each validation size, the following are paired differences in replicate-level posterior-median bias relative to gamma=0. All 10 replicate seeds are paired across curvature conditions. Intervals are paired t intervals with 9 degrees of freedom.

| Validation n | Contrast | Mean paired bias difference | 95% paired t interval |
|---:|---|---:|---:|
| 100 | γ=-1 minus γ=0 | -4.26 | [-7.20, -1.32] |
| 300 | γ=-1 minus γ=0 | -4.27 | [-6.67, -1.87] |
| 500 | γ=-1 minus γ=0 | -4.04 | [-6.33, -1.75] |
| 100 | γ=+1 minus γ=0 | +5.08 | [+2.32, +7.84] |
| 300 | γ=+1 minus γ=0 | +5.40 | [+2.39, +8.42] |
| 500 | γ=+1 minus γ=0 | +5.36 | [+2.38, +8.34] |

The direction of the median-bias shift is consistent across validation sizes: negative curvature lowers the median bias relative to gamma=0, while positive curvature raises it. These paired contrasts are exploratory and should not be read as a universal effect size.

For gamma=+1, the paired increase in absolute median error relative to gamma=0 was approximately +3.98, +4.13, and +3.97 at validation sizes 100, 300, and 500. The corresponding 95% paired t intervals were [0.49, 7.47], [0.14, 8.13], and [-0.06, 8.00]. Thus the absolute-error evidence is weaker than the signed-bias contrast, especially at n=500.

### Degree-one coefficient projection at validation n=500

As a diagnostic, I also summarized the degree-one inclusion-coefficient posterior means (prior multiplier 1) at the largest validation size. The generator's alpha values are gamma-specific; beta is -3.9 in every condition.

| γ | Generator α | Mean posterior α | Mean α difference from generator | Mean posterior β | Mean β difference from -3.9 |
|---:|---:|---:|---:|---:|---:|
| -1 | -0.125 | -0.052 | +0.073 [-0.023, +0.169] | -3.894 | +0.006 [-0.164, +0.177] |
| 0 | 0.000 | -0.048 | -0.048 [-0.155, +0.059] | -3.946 | -0.046 [-0.229, +0.136] |
| +1 | +0.125 | -0.050 | -0.176 [-0.293, -0.058] | -3.944 | -0.044 [-0.199, +0.112] |

Brackets are 95% t intervals across the 10 replicate-level posterior-mean coefficient differences. The slope posterior means remain near -3.9, while at gamma=+1 the degree-one model's posterior intercept remains below the generator alpha. This is consistent with the fitted linear form failing to represent the positive quadratic curvature. For gamma ±1, however, the degree-one model is misspecified, so its coefficients target a pseudo-true projection rather than the generator's alpha/beta; these differences are diagnostic only and must not be described as ordinary parameter-estimation bias.

### Prior-scale sensitivity

The child artifacts also evaluate slope-prior SD multipliers 0.5, 1, and 2, with multiplier 1 as the reference. The table gives the mean posterior-median shift from that reference across the 10 replicates in each cell (km s⁻¹ Mpc⁻¹):

| γ | Validation n | Multiplier 0.5 shift | Multiplier 2 shift |
|---:|---:|---:|---:|
| -1 | 100 | -2.62 | +0.96 |
| -1 | 300 | -1.10 | +0.32 |
| -1 | 500 | -0.63 | +0.17 |
| 0 | 100 | -3.06 | +1.11 |
| 0 | 300 | -1.26 | +0.37 |
| 0 | 500 | -0.72 | +0.20 |
| +1 | 100 | -3.04 | +1.00 |
| +1 | 300 | -1.26 | +0.36 |
| +1 | 500 | -0.75 | +0.20 |

The prior-scale shifts are smaller than the mean curvature-associated signed-bias differences, and they diminish as the validation sample grows. They do not explain the entire gamma=0 baseline offset or the gamma=+1 contrast. This is a limited sensitivity analysis over one coefficient-prior scale family, not a proof of prior robustness.

### Validation-size contrast

The paired change from validation n=100 to n=500 in signed median bias was:
- gamma=-1: -0.49 [-1.24, +0.26];
- gamma=0: -0.71 [-2.07, +0.65];
- gamma=+1: -0.43 [-1.36, +0.50].

All three intervals include zero. Mean 68% interval width changed by -0.57, -0.48, and -0.51 respectively; each change is small relative to the width itself. Meanwhile, median minimum catalogue Kish ESS increased from 14.49 to 63.83 at gamma=-1, 10.04 to 30.49 at gamma=0, and 9.67 to 24.61 at gamma=+1. Median maximum normalized weights fell with validation size in all three curvature conditions.

This suggests that a larger validation sample improves the stability diagnostics of estimated inclusion weights, but does not by itself remove the observed H0 median-bias shift in this four-event setup.

## H0 grid-refinement result

A separate fixed-observation experiment completed three child runs at 141, 281, and 561 H0 grid points. All three child likelihood/artifact audits passed, all used coefficient quadrature order 7, and the observable-input hashes match. The aggregate workflow initially failed because NumPy 2.5 removed `np.trapz`; the grid-audit script and regression fixture have been updated to `np.trapezoid`, and the corrected CI rerun is pending.

An independent local reconstruction of the archived grid artifact gives:

| H0 grid points | Spacing | Max absolute density difference vs 561 | Total variation vs 561 | Max absolute quantile difference vs 561 |
|---:|---:|---:|---:|---:|
| 141 | 2.0 | 2.23e-5 | 3.23e-4 | 0.0272 |
| 281 | 1.0 | 5.58e-6 | 6.47e-5 | 0.0055 |
| 561 | 0.5 | 0 | 0 | 0 |

The model and generator configurations match across grids after excluding only `h0_grid_points`. The finest-grid posterior median is approximately 130.036, with 95% interval approximately [86.389, 181.062], in this one synthetic realization. The difference between the 141-point and 561-point quantiles is negligible compared with the large posterior uncertainty and the roughly 60 km s⁻¹ Mpc⁻¹ median offset from truth.

This establishes numerical H0-grid stability for one fixed realization only. It does not establish repeated-sampling calibration, SBC, coverage, or estimator superiority.

## Interpretation

Three findings deserve further investigation:

1. **Curvature-direction sensitivity:** gamma=+1 produces about +5.1 to +5.4 km s⁻¹ Mpc⁻¹ additional median bias relative to gamma=0, while gamma=-1 produces about -4.0 to -4.3. The sign is stable across the three validation sizes in this pilot.
2. **A baseline offset remains at gamma=0:** even with the correctly specified degree-one inclusion functional form, mean posterior-median bias is around +7.5 to +8.2. Therefore the entire offset cannot be attributed to omitted quadratic curvature. The event-count conditioning, small number of events, event noise, selection normalization, finite catalogue and posterior model all remain relevant.
3. **More validation does not materially shift H0 in this design:** ESS improves and maximum normalized weights decline, but median bias and interval width change only modestly from n=100 to n=500.

The 68% coverage estimates are imprecise at 10 replicates. In particular, 5/10 coverage at gamma=+1 has a Wilson 95% interval [0.237, 0.763]. These results are not a calibrated estimate of long-run coverage.

## Limitations and next gates

- The H0 grid refinement is numerically stable for one fixed realization; the CI rerun of the aggregate comparison is pending, and repeated-sampling grid stability has not been established.
- Coefficient quadrature convergence (all 90 runs passed) is distinct from H0-grid convergence.
- Four detected events are fixed in every replicate; event-count/rate calibration is not tested.
- Ten replicates per cell are inadequate for confirmatory coverage or robust tail diagnostics.
- This is a low-redshift Euclidean synthetic benchmark, not LVK/GWTC/public-data analysis, real galaxy-survey validation, or an observational H0 measurement.
- Fixed-truth repeated-sampling coverage is not joint simulation-based calibration (SBC). SBC/PIT should be a separate subsequent experiment.
- The gamma=0 baseline bias needs diagnosis before assigning the curvature contrast solely to model misspecification.

### Next actions

1. Obtain passing CI replays of the corrected independent orchestration audit and the corrected grid-refinement auditor.
2. Diagnose the positive gamma=0 baseline bias using event-level likelihood and selection-normalization checks.
3. Increase replicates before any calibration claim; report Wilson/binomial uncertainty and preserve failed/non-evaluable cells in planned denominators.
4. Keep fixed event-count coverage separate from any future rate-informed analysis.

**Bottom line:** the full 90-cell simulation completed successfully and its artifacts, child audits, quadrature convergence, paired sampling design, and aggregate metrics have been independently checked. The preliminary result is an asymmetric curvature-related shift in H0 median bias on top of a nonzero gamma=0 baseline offset. This is a strong reason to proceed to grid refinement and larger repeated-sampling diagnostics—not yet a confirmatory scientific conclusion.
