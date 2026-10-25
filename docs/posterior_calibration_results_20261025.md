# Posterior-level H₀ calibration: audited results (2026-10-25)

**Status: independent audit PASS.** This report summarizes the reproducible synthetic two-host experiment. It is not a real gravitational-wave measurement.

- **Workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40172704350
- **Audited artifact:** https://api.github.com/repos/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/artifacts/12590881838/zip
- **Protocol:** [posterior calibration protocol](posterior_calibration_protocol.md)
- **Data:** 2,400 raw observations; 96 truth × noise × support × prior-measure × weighting summaries.
- **Injected truths:** 55, 65, 75, 85 km s⁻¹ Mpc⁻¹.
- **Noise fractions:** 3%, 10%, 20%.
- **Replicates:** 200 per truth/noise cell.
- **Audit:** independently reconstructed posterior medians and 68%/95% interval endpoints for every raw observation and inference arm; independently recomputed calibration metrics; verified hashes for both CSV inputs.

## Pooled posterior calibration across the prespecified design

The following figures pool the 2,400 synthetic observations equally across the 12 truth/noise cells. Each inference arm uses exactly the same synthetic observations. Bias is mean posterior-median error (median minus injected H₀); RMSE and MAE are in km s⁻¹ Mpc⁻¹. Coverage is the fraction of intervals containing the injected H₀.

| Prior support | Prior measure | Catalogue weights | Bias | MAE | RMSE | Mean 68% width | 68% coverage | Mean 95% width | 95% coverage |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| [50, 90] | Uniform H₀ | Correct | +0.60 | 6.22 | 8.31 | 16.68 | 68.9% | 26.87 | 98.0% |
| [50, 90] | Uniform H₀ | Reversed | +2.27 | 6.52 | 9.01 | 16.35 | 70.9% | 26.74 | 97.8% |
| [50, 90] | Log-uniform H₀ | Correct | −0.71 | 6.12 | 8.29 | 16.63 | 71.8% | 26.90 | 97.8% |
| [50, 90] | Log-uniform H₀ | Reversed | +1.00 | 6.31 | 8.50 | 16.68 | 69.8% | 26.92 | 97.8% |
| [20, 220] | Uniform H₀ | Correct | +13.69 | 29.51 | 36.94 | 74.62 | 60.0% | 109.77 | 93.5% |
| [20, 220] | Uniform H₀ | Reversed | +48.36 | 49.96 | 60.19 | 68.16 | 30.3% | 110.09 | 86.9% |
| [20, 220] | Log-uniform H₀ | Correct | +1.01 | 21.95 | 28.21 | 67.25 | 66.5% | 104.12 | 95.2% |
| [20, 220] | Log-uniform H₀ | Reversed | +39.83 | 42.33 | 51.55 | 73.40 | 52.5% | 108.96 | 93.8% |

## Stratified-bootstrap uncertainty

The workflow also performed 4,000 stratified bootstrap resamples for each inference arm and metric. Resampling was done within each of the 12 truth/noise cells, retaining the fixed design allocation. The intervals quantify Monte Carlo uncertainty conditional on this synthetic design.

**Workflow with uncertainty output:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40173199370  
**Artifact containing `uncertainty.csv` and `uncertainty_report.md`:** https://api.github.com/repos/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/artifacts/12591236798/zip

| Prior support | Prior measure | Weights | Bias (95% CI) | RMSE (95% CI) | 68% coverage (95% CI) | 95% coverage (95% CI) |
|---|---|---|---:|---:|---:|---:|
| [50,90] | Uniform H₀ | Correct | +0.60 [0.40, 0.80] | 8.31 [8.10, 8.52] | 68.9% [67.3%, 70.5%] | 98.0% [97.4%, 98.5%] |
| [50,90] | Uniform H₀ | Reversed | +2.27 [2.07, 2.46] | 9.01 [8.81, 9.22] | 70.9% [69.4%, 72.5%] | 97.8% [97.2%, 98.3%] |
| [50,90] | Log-uniform H₀ | Correct | −0.71 [−0.91, −0.50] | 8.29 [8.07, 8.52] | 71.8% [70.1%, 73.5%] | 97.8% [97.2%, 98.3%] |
| [50,90] | Log-uniform H₀ | Reversed | +1.00 [0.80, 1.19] | 8.50 [8.30, 8.70] | 69.8% [68.1%, 71.4%] | 97.8% [97.2%, 98.4%] |
| [20,220] | Uniform H₀ | Correct | +13.69 [12.44, 14.91] | 36.94 [36.02, 37.84] | 60.0% [58.1%, 62.0%] | 93.5% [92.5%, 94.4%] |
| [20,220] | Uniform H₀ | Reversed | +48.36 [46.98, 49.72] | 60.19 [59.13, 61.20] | 30.3% [28.4%, 32.1%] | 86.9% [85.5%, 88.2%] |
| [20,220] | Log-uniform H₀ | Correct | +1.01 [−0.04, 2.13] | 28.21 [27.49, 28.93] | 66.5% [64.7%, 68.2%] | 95.2% [94.3%, 96.0%] |
| [20,220] | Log-uniform H₀ | Reversed | +39.83 [38.57, 41.12] | 51.55 [50.61, 52.46] | 52.5% [50.5%, 54.5%] | 93.8% [92.8%, 94.7%] |

These are stratified percentile-bootstrap intervals over replicates within each truth/noise cell; they are not Bayesian credible intervals for H₀ and do not capture model misspecification beyond the simulation design.

## Interpretation

### 1. Prior support changes calibration dramatically in this toy model

Under the broad uniform prior, even the correctly weighted inference arm has a +13.69 km s⁻¹ Mpc⁻¹ pooled median bias and RMSE 36.94. With reversed weights, the pooled bias rises to +48.36 and RMSE to 60.19. The narrow [50,90] support gives much smaller pooled errors in this selected truth grid, but this does not establish that narrow priors are generally preferable: the prior is much more informative and restricts the range of allowed H₀ values.

### 2. The prior measure matters, especially under broad support

With correct weights and broad support, changing from uniform to log-uniform H₀ reduces pooled bias from +13.69 to +1.01 and RMSE from 36.94 to 28.21. With reversed weights, it reduces pooled bias from +48.36 to +39.83 and RMSE from 60.19 to 51.55, but inference remains badly biased. These are design-specific effects of this mixture and prior choices, not evidence that log-uniform priors are universally superior.

### 3. Incorrect catalogue weights can cause severe posterior miscalibration

Under the broad support, reversed weights reduce 68% coverage to 30.3% for a uniform prior and 52.5% for a log-uniform prior. Corresponding pooled 95% coverage is 86.9% and 93.8%. The experiment therefore shows that host-weight misspecification can substantially affect posterior calibration in this model. The effect is partly hidden by the narrow support, where all injected truths lie inside [50,90] and the restricted posterior cannot explore much of the high-H₀ tail.

### 4. Coverage estimates have Monte Carlo uncertainty

There are 2,400 observations per inference arm, but the design intentionally stratifies observations across 12 truth/noise cells rather than sampling one homogeneous population. The table is a descriptive pooled estimate; it should not replace truth- and noise-specific coverage tables. With 200 replicates in each cell, cell-specific coverage has appreciable Monte Carlo uncertainty, especially near nominal 95% coverage. The current workflow now reports stratified-bootstrap intervals for pooled bias, RMSE, interval widths, and coverage. A higher-replicate follow-up should narrow Monte Carlo uncertainty and report truth/noise-specific coverage intervals as well.

## Scope and limitations

This is a synthetic low-redshift, two-host Gaussian-distance model. It does not include realistic gravitational-wave event likelihoods, galaxy-catalogue selection or incompleteness, redshift errors, detector selection, or hierarchical multi-event population inference. The distance-noise scale is defined relative to the z₂ reference distance in each truth/noise cell. Results are valid only for the stated simulation design.

The appropriate next step is a higher-replicate, truth-stratified calibration run with uncertainty intervals for coverage and a more realistic event/catalogue selection model. Do not present these numbers as an observational H₀ result.

## Truth/noise-specific uncertainty extension

The workflow now additionally writes `cell_uncertainty.csv` and `cell_uncertainty_report.md`. These report bias, RMSE, 68%/95% coverage, and mean 68%/95% interval widths separately for each of the 12 injected-truth × noise cells and each of the eight inference arms. Each interval is calculated by resampling replicates within that one cell, rather than pooling across cells. The design yields 96 arm-by-cell rows.

This cell-level view is important because pooled calibration metrics can conceal opposite behavior at different injected H₀ values or noise levels. The cell-specific percentile intervals quantify Monte Carlo uncertainty under the selected simulation design only; they do not quantify model uncertainty or real-event uncertainty.

The extension is gated on the existing independent audit and SHA-256 checks. The automated test suite includes checks for complete cell coverage, report generation, interval ordering, and rejection of artifacts that have not passed the audit.


## Artifact verification and cell-level findings (workflow run 40173820045)

I downloaded and inspected the workflow artifact rather than relying only on the green workflow badge. The independent audit reports `PASS`, with 2,400 raw observations, 96 summary rows, reconstructed posterior summaries and calibration metrics, and verified SHA-256 hashes. The archive contains both `cell_uncertainty.csv` and `cell_uncertainty_report.md`.

The cell-level CSV has exactly 96 rows: 12 truth/noise cells × eight inference arms. The report has 96 data rows and two Markdown table header/separator rows. Each cell has 200 replicates and 4,000 within-cell bootstrap resamples for bias, RMSE, 68%/95% coverage, and mean 68%/95% interval widths. A consistency check found all 576 reported metric intervals ordered correctly and all point estimates within their corresponding percentile intervals. This is a post-generation artifact inspection, separate from the independent posterior reconstruction audit.

### Cell-level contrast: broad uniform prior, correct versus reversed host weights

Across the 12 truth/noise cells, reversing the host weights increased the *mean of cell-specific bias estimates* by 34.67 km s⁻¹ Mpc⁻¹ and the *mean of cell-specific RMSE estimates* by 24.80 km s⁻¹ Mpc⁻¹. Mean 68% coverage fell by 29.8 percentage points, and mean 95% coverage fell by 6.6 percentage points. These are unweighted averages of the 12 cell-specific metrics, not pooled-replicate estimates and not uncertainty intervals for paired contrasts.

The most adverse single cell for broad-uniform reversed weights was truth H₀ = 85 km s⁻¹ Mpc⁻¹ with 3% distance noise: bias +60.93 km s⁻¹ Mpc⁻¹ (cell bootstrap 95% interval [+55.41, +66.19]), RMSE 72.35 km s⁻¹ Mpc⁻¹ (95% interval [68.79, 75.62]), and 68% coverage 23.5% (95% interval 18.0%–29.5%). This illustrates why the pooled results alone are insufficient.

**Reproducibility note:** these verified artifact files were generated by successful workflow run [40173820045](https://github.com/vedantkushwaha033-blip/dark-siren-cosmology-robustness/actions/runs/40173820045), at commit `afe9a427e7bace84be7fedabe6e18b1ce69134db`. The PR branch currently points to a different commit, so this historical run should not be described as the latest-head check. A fresh run against the current branch head is required to establish that the current branch reproduces the same outputs.

These findings remain specific to the synthetic two-host model. The per-cell bootstrap intervals quantify replicate Monte Carlo uncertainty conditional on this design; they do not address realistic catalogue selection, gravitational-wave likelihoods, or model misspecification.


## Paired host-weight contrasts: exploratory result and reproducible workflow

A paired contrast script, `scripts/analyze_paired_host_weight_contrasts.py`, has been added. It compares correct and reversed weights on the **same simulated observation**, within each truth/noise cell. This pairing avoids treating the two inference arms as independent samples. The script requires the independent audit and SHA-256 checks to pass, then writes:
- `paired_host_weight_cell_contrasts.csv` (48 rows: 12 cells × four prior configurations);
- `paired_host_weight_contrasts.csv` (24 rows: four prior configurations × six metrics);
- `paired_host_weight_contrasts.md` (aggregate contrast table).

The aggregate is an equal-weight mean across the 12 design cells. Bootstrap intervals resample paired observations within each cell and then combine cell-level bootstrap draws. Coverage differences are proportions (e.g., −0.298 is −29.8 percentage points). This target is distinct from a pooled-observation difference.

### Exploratory paired estimates from the previously audited artifact

For the broad uniform-H₀ prior, reversed minus correct host weights produced:
- bias difference: +34.67 km s⁻¹ Mpc⁻¹ (stratified paired-bootstrap 95% interval [+34.33, +35.02]);
- RMSE difference: +24.80 km s⁻¹ Mpc⁻¹ (95% interval [+23.94, +25.66]);
- 68% coverage difference: −0.298 (95% interval [−0.320, −0.276]);
- 95% coverage difference: −0.0658 (95% interval [−0.0771, −0.0554]).

For broad log-uniform H₀, the corresponding differences were +38.82 in bias (95% interval [+38.41, +39.21]), +23.68 in RMSE (95% interval [+22.69, +24.70]), −0.140 in 68% coverage (95% interval [−0.160, −0.119]), and −0.0138 in 95% coverage (95% interval [−0.0217, −0.0054]).

Under the narrow support, the effect is smaller and metric-dependent. For narrow uniform H₀, reversed weights increased bias by +1.67 (95% interval [+1.55, +1.79]) and RMSE by +0.277 (95% interval [+0.140, +0.405]); 68% coverage increased by 0.0204 (95% interval [+0.0113, +0.0296]). For narrow log-uniform H₀, bias increased by +1.71 (95% interval [+1.59, +1.83]), while the RMSE contrast interval included zero. Thus, this simulation does not support the simplistic claim that reversed weights reduce coverage under every prior: the narrow uniform configuration is a counterexample for 68% coverage.

These intervals are conditional Monte Carlo intervals under this synthetic design. They do not include uncertainty about model structure, realistic GW likelihoods, catalogue selection, or astrophysical populations. The figures above were exploratory values from the previously audited artifact. A fresh run against the current branch head has now completed successfully; its artifact was downloaded and the generated paired-contrast files inspected. Current-head results are recorded below.


## Latest completed paired-contrast artifact verification

**Verified run:** [Posterior-level H₀ calibration, run 40185764559](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40185764559)  
**Commit:** `70dad55ea35e6f40c77bd32b4fc693b45bf7e090`  
**Artifact ID:** `12596401738` (GitHub-reported archive digest: `sha256:aec34e3c4d2d0543f0ae1e82d1a5fcecfd401f58a6c072342f4a5a512bf9932e`).

This run completed successfully. I downloaded and inspected the artifact archive. The independent audit reports PASS for 2,400 raw observations and 96 summary rows, confirms posterior summaries and calibration metrics were reconstructed from raw observations, and confirms hashes were verified. The artifact contains 48 paired cell-level rows and 24 aggregate rows. The run uses 200 paired observations per truth/noise cell and 4,000 bootstrap resamples. The associated test run at this commit also passed.

This is the **latest completed successful calibration artifact**, but it is not the latest repository head: subsequent commits vectorized the paired bootstrap implementation and added a regression test comparing vectorized results with scalar calculations. Those later commits have fresh workflows running, so their final statuses and artifact equivalence must still be confirmed.

### Paired contrasts: reversed minus correct host weights

| Prior support | Prior measure | Bias difference (km s⁻¹ Mpc⁻¹) | RMSE difference (km s⁻¹ Mpc⁻¹) | 68% coverage difference | 95% coverage difference |
|---|---|---:|---:|---:|---:|
| [20,220] | Uniform H₀ | +34.669 [34.320, 35.026] | +24.800 [23.935, 25.679] | −0.29792 [−0.31917, −0.27792] | −0.06583 [−0.07667, −0.05500] |
| [20,220] | Log-uniform H₀ | +38.818 [38.425, 39.210] | +23.684 [22.689, 24.680] | −0.13958 [−0.16083, −0.11833] | −0.01375 [−0.02168, −0.00542] |
| [50,90] | Uniform H₀ | +1.668 [1.548, 1.786] | +0.277 [0.147, 0.398] | +0.02042 [+0.01125, +0.02958] | −0.00208 [−0.00458, 0.00000] |
| [50,90] | Log-uniform H₀ | +1.709 [1.584, 1.837] | +0.083 [−0.061, +0.221] | −0.02000 [−0.02958, −0.01042] | +0.00042 [−0.00208, +0.00292] |

### Interval-width contrasts

| Prior support | Prior measure | Mean 68% width difference | Mean 95% width difference |
|---|---|---:|---:|
| [20,220] | Uniform H₀ | −6.464 [−6.588, −6.345] | +0.320 [+0.158, +0.478] |
| [20,220] | Log-uniform H₀ | +6.152 [+5.976, +6.328] | +4.839 [+4.646, +5.029] |
| [50,90] | Uniform H₀ | −0.329 [−0.459, −0.197] | −0.126 [−0.228, −0.016] |
| [50,90] | Log-uniform H₀ | +0.054 [−0.078, +0.187] | +0.020 [−0.075, +0.121] |

Bracketed values are stratified paired-bootstrap 95% intervals. Coverage differences are proportions: multiply by 100 to obtain percentage points. Width differences use the same H₀ units as the prior support. Negative values mean the reversed-weight posterior intervals are narrower; that is not automatically an improvement, because coverage and estimation error must be considered jointly.

The broad-support results show a pronounced degradation in point-estimation error and coverage under reversed weights for this synthetic model. The width results also show that the effect is not a uniform broadening: broad uniform weighting produces narrower 68% intervals but much worse coverage, whereas broad log-uniform weighting increases both interval widths. This is a useful warning against judging uncertainty quality from interval width alone.

The narrow-support effects are much smaller and depend on the prior measure and metric. For example, narrow uniform H₀ shows a positive 68% coverage contrast, while the 95% coverage interval reaches zero. These are descriptive, conditional Monte Carlo results—not universal claims about host weighting and not an observational H₀ measurement.

The successful run and artifact establish reproducibility for the stated synthetic design and commit. They do not address realistic gravitational-wave likelihoods, catalogue selection, peculiar velocities, detector selection, or model misspecification. The next milestone is to verify the newer vectorized implementation on the current branch head and then prioritize a separately prespecified extension that addresses a specific omitted physical or selection effect.

## Vectorized bootstrap reproduction check

**Successful calibration run:** [40185998682](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40185998682)  
**Code commit tested:** `a4cb98307adcf78eb760aadf7ff255263a60feec`  
**Artifact ID:** `12595443879`; GitHub-reported archive digest: `sha256:efd813266502c8daf25a70028b690bbb489355919e171542f3c430ed3253633d`.

This run exercised the vectorized paired-bootstrap implementation and completed successfully. Its test workflow also passed. I downloaded and independently inspected the artifact: the audit again reports PASS, the 2,400 raw observations and 96 summaries are present, and the 24 aggregate paired contrasts were generated.

I compared the paired-contrast CSV from this run with the preceding successful run at commit `70dad55ea35e6f40c77bd32b4fc693b45bf7e090`. All 24 aggregate rows matched in point estimates and both interval endpoints to a maximum absolute floating-point difference of **1.78 × 10⁻¹⁵**. The CSV bytes differ because of insignificant floating-point formatting/rounding, but the numeric results are equivalent at machine precision. This confirms that vectorization improved the calculation path without materially changing the statistical results for the fixed seed and design.

The currently newer branch head includes documentation updates after this verified code commit. The scientific implementation itself was tested successfully at the code commit above; any subsequent code change should trigger another validation run.

## Heterogeneity across the 12 truth/noise cells

As a descriptive check against a misleading macro-average, I summarized the 12 cell-specific paired point estimates from the verified vectorized-run artifact (run 40185998682). These are ranges of cell estimates, not confidence intervals, and no multiplicity-adjusted cell-wise hypothesis tests are implied.

| Prior support / measure | Bias difference range | RMSE difference range | 68% coverage difference range | Direction across cells |
|---|---:|---:|---:|---|
| Broad uniform H₀ | +15.78 to +65.19 | +12.35 to +48.94 | −0.415 to −0.235 | Bias/RMSE positive in 12/12; 68% coverage negative in 12/12 |
| Broad log-uniform H₀ | +20.96 to +66.86 | +11.49 to +46.58 | −0.215 to −0.070 | Bias/RMSE positive in 12/12; 68% coverage negative in 12/12 |

Bias and RMSE differences are in km s⁻¹ Mpc⁻¹; coverage differences are proportions. For broad uniform H₀, the 95% coverage contrast was also negative in all 12 cells (range −0.105 to −0.025). For broad log-uniform H₀, it was negative in 9 cells and positive in 3 (range −0.045 to +0.020). These cell-level signs strengthen the descriptive claim that the broad-support bias, RMSE, and 68% coverage effects are not driven by only one truth/noise setting in this design. They do not establish generality beyond this simulator or account for uncertainty in the choice of model.

## Latest completed validation after regression-test strengthening

**Successful calibration run:** [40186131023](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40186131023)  
**Code and test commit:** `1df6f450ec3bfaae72044951b4cd4e30069dbee1`  
**Artifact ID:** `12596377387`; GitHub-reported archive digest: `sha256:4637eafdbe0cf09ff4af3a4bcad7ba8b63fa20d243243c0bf932feda4daa3240`.

This later full calibration run also passed, as did the three-version test matrix. The tests now explicitly check that vectorized bootstrap results match scalar calculations, repeated analyses with the same seed are deterministic, and missing or duplicated replicate blocks are rejected.

I downloaded this run's artifact and compared its `paired_host_weight_contrasts.csv` with the artifact from run 40185998682. The files were **byte-for-byte identical**, and all 24 aggregate estimates and confidence-interval endpoints matched exactly. The audit again passed for 2,400 raw records and 96 summaries. Thus the final bootstrap implementation and strengthened tests are supported by a successful end-to-end calibration run at commit `1df6f450ec3bfaae72044951b4cd4e30069dbee1`.

The repository head subsequently advanced through documentation-only commits that record these checks, the per-cell heterogeneity summary, and a consistent displayed coverage percentage. No analysis code changed after the successful calibration commit described here; the latest documentation-only head's test workflow is queued separately.



## Posterior-grid numerical convergence checks

**Workflow:** [numerical diagnostics, run 40186240390](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40186240390)  
**Code revision:** `2072d254f51069d89cc75bc7cdd45053b0917ca1)  
**Environments:** Python 3.10, 3.11, and 3.12; all three numerical-diagnostic jobs completed successfully and uploaded result artifacts.

The diagnostic separates numerical-grid effects from changes to the prior itself:

- **Refining the grid on the same [50,90] prior:** increasing from 1,601 to 3,201 points changed the posterior median by about (1.5\times10^{-6}) and the 68%/95% interval endpoints by at most about (5.5\times10^{-6}) H₀ units in the reported baseline case.
- **Widening only the computational grid while keeping the prior fixed at [50,90]:** posterior summaries agreed to floating-point precision (maximum displayed differences around (1.4\times10^{-14})). This is evidence that the computational grid bounds outside the prior support are not truncating the posterior in this case.
- **Changing the prior support from [50,90] to [30,110]:** this is not a numerical-grid refinement. It materially changed the posterior: the reported median shifted from 70.42 to 68.20, posterior SD increased from 10.24 to 20.94, and the 95% interval widened from approximately [51.87, 88.76] to [31.72, 105.08]. Prior sensitivity is therefore much larger than the numerical-grid error in this example.
- **Targeted upper-tail stress case:** for an injected H₀ of 52, reversed weights, and 3% noise, the absolute change in the 95% upper endpoint between consecutive grids fell from 0.106 (1,001 to 2,001 points) to 0.0243, 0.00601, 0.00148, 0.000373, and 0.0000929 as the grid was successively doubled up to 64,001 points. Simpson-versus-trapezoid posterior-area discrepancy also decreased from (3.38\times10^{-6}) to (8.25\times10^{-10}).

These diagnostics support numerical stability for the tested examples, including a deliberately challenging upper-tail case. They do **not** prove uniform convergence for every possible observation, prior, or noise regime. The stress case also shows why reporting only a typical posterior is insufficient: upper-tail quantiles can converge more slowly than the median. Prior-support sensitivity is a separate modeling choice and must not be described as numerical error.

