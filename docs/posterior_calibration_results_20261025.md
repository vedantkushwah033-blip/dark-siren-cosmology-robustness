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
| [20,220] | Uniform H₀ | Reversed | +48.36 [46.98, 49.72] | 60.19 [59.13, 61.20] | 30.2% [28.4%, 32.1%] | 86.9% [85.5%, 88.2%] |
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

These intervals are conditional Monte Carlo intervals under this synthetic design. They do not include uncertainty about model structure, realistic GW likelihoods, catalogue selection, or astrophysical populations. The values above were calculated from the already audited artifact; the new automated workflow and tests are being run against the updated branch so the same outputs can be reproduced from the checked-in script.
