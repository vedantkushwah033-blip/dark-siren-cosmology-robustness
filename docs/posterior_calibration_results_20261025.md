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

## Interpretation

### 1. Prior support changes calibration dramatically in this toy model

Under the broad uniform prior, even the correctly weighted inference arm has a +13.69 km s⁻¹ Mpc⁻¹ pooled median bias and RMSE 36.94. With reversed weights, the pooled bias rises to +48.36 and RMSE to 60.19. The narrow [50,90] support gives much smaller pooled errors in this selected truth grid, but this does not establish that narrow priors are generally preferable: the prior is much more informative and restricts the range of allowed H₀ values.

### 2. The prior measure matters, especially under broad support

With correct weights and broad support, changing from uniform to log-uniform H₀ reduces pooled bias from +13.69 to +1.01 and RMSE from 36.94 to 28.21. With reversed weights, it reduces pooled bias from +48.36 to +39.83 and RMSE from 60.19 to 51.55, but inference remains badly biased. These are design-specific effects of this mixture and prior choices, not evidence that log-uniform priors are universally superior.

### 3. Incorrect catalogue weights can cause severe posterior miscalibration

Under the broad support, reversed weights reduce 68% coverage to 30.3% for a uniform prior and 52.5% for a log-uniform prior. Corresponding pooled 95% coverage is 86.9% and 93.8%. The experiment therefore shows that host-weight misspecification can substantially affect posterior calibration in this model. The effect is partly hidden by the narrow support, where all injected truths lie inside [50,90] and the restricted posterior cannot explore much of the high-H₀ tail.

### 4. Coverage estimates have Monte Carlo uncertainty

There are 2,400 observations per inference arm, but the design intentionally stratifies observations across 12 truth/noise cells rather than sampling one homogeneous population. The table is a descriptive pooled estimate; it should not replace truth- and noise-specific coverage tables. With 200 replicates in each cell, cell-specific coverage has appreciable Monte Carlo uncertainty, especially near nominal 95% coverage. A higher-replicate follow-up should add stratified bootstrap or Wilson intervals to coverage estimates.

## Scope and limitations

This is a synthetic low-redshift, two-host Gaussian-distance model. It does not include realistic gravitational-wave event likelihoods, galaxy-catalogue selection or incompleteness, redshift errors, detector selection, or hierarchical multi-event population inference. The distance-noise scale is defined relative to the z₂ reference distance in each truth/noise cell. Results are valid only for the stated simulation design.

The appropriate next step is a higher-replicate, truth-stratified calibration run with uncertainty intervals for coverage and a more realistic event/catalogue selection model. Do not present these numbers as an observational H₀ result.
