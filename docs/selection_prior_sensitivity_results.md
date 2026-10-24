# Selection-normalization prior-sensitivity audit

**Status:** Completed by GitHub Actions on 2026-10-24.  
**Workflow:** [Toy selection prior sensitivity, run 40071245383](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40071245383)  
**Artifact:** `selection-prior-sensitivity-results` (contains `per_run.csv` and `summary.csv`).  
**Seed:** 20261029.

## Question

Does dividing the toy dark-siren likelihood by the selection probability improve inference consistently when the prior support changes?

This is a controlled synthetic experiment, not an analysis of LVK detections or a measurement of the real Hubble constant.

## Design

- True values: H0 = 55, 70, 85 km s^-1 Mpc^-1.
- 600 detected mock observations per true value (1,800 observations total).
- Each observation is analysed with all three prior supports: [50, 90], [40, 100], [30, 120].
- Each observation and prior is analysed both without and with the selection normalization.
- Host redshifts: 0.02 and 0.04; nominal mixture weights 0.70 and 0.30.
- Gaussian distance uncertainty: 20 Mpc.
- Detection rule: 0 < observed distance < 180 Mpc.
- Posterior grid: 4,001 points spanning H0 = 25 to 125.
- Outputs: 5,400 per-run rows (1,800 observations x 3 priors) and 27 summary rows.

The paired metric is corrected absolute error minus uncorrected absolute error. Negative values favour the correction for that cell; positive values favour the uncorrected method. The reported MCSE is the standard error of the paired mean difference, not a posterior uncertainty.

## Results by truth and prior

Bias, RMSE, MAE and coverage are estimated over 600 simulated observations in each cell. Coverage values are empirical proportions, not exact guarantees.

| True H0 | Prior | Method | Bias | RMSE | MAE | 68% interval coverage | 95% interval coverage |
|---:|:---:|:---|---:|---:|---:|---:|---:|
| 55 | [50, 90] | Uncorrected | +11.288 | 12.527 | 11.288 | 57.7% | 95.5% |
| 55 | [50, 90] | Corrected | +9.404 | 10.867 | 9.404 | 64.2% | 97.0% |
| 55 | [40, 100] | Uncorrected | +9.392 | 11.819 | 9.392 | 72.2% | 96.5% |
| 55 | [40, 100] | Corrected | +6.696 | 10.143 | 7.370 | 77.3% | 97.7% |
| 55 | [30, 120] | Uncorrected | +13.385 | 16.185 | 13.851 | 67.0% | 95.5% |
| 55 | [30, 120] | Corrected | +9.076 | 13.918 | 10.938 | 72.3% | 96.5% |
| 70 | [50, 90] | Uncorrected | +1.510 | 5.974 | 5.104 | 94.8% | 100.0% |
| 70 | [50, 90] | Corrected | +0.105 | 6.142 | 5.185 | 96.5% | 100.0% |
| 70 | [40, 100] | Uncorrected | +1.852 | 8.673 | 6.712 | 89.5% | 99.8% |
| 70 | [40, 100] | Corrected | −0.718 | 9.684 | 7.989 | 91.5% | 100.0% |
| 70 | [30, 120] | Uncorrected | +2.508 | 16.439 | 13.807 | 81.3% | 98.7% |
| 70 | [30, 120] | Corrected | −1.333 | 18.829 | 15.387 | 79.3% | 98.8% |
| 85 | [50, 90] | Uncorrected | −11.096 | 12.836 | 11.096 | 53.2% | 100.0% |
| 85 | [50, 90] | Corrected | −12.491 | 14.357 | 12.491 | 41.3% | 100.0% |
| 85 | [40, 100] | Uncorrected | −10.133 | 14.766 | 11.662 | 91.0% | 100.0% |
| 85 | [40, 100] | Corrected | −12.529 | 17.608 | 13.912 | 77.2% | 100.0% |
| 85 | [30, 120] | Uncorrected | −6.020 | 18.564 | 15.252 | 85.8% | 99.8% |
| 85 | [30, 120] | Corrected | −8.975 | 21.929 | 17.650 | 82.8% | 99.3% |

H0 is in km s^-1 Mpc^-1; RMSE and interval widths use the same units.


## Monte Carlo uncertainty on empirical coverage

Each coverage estimate is based on n=600 independent replicates within a cell. The following are 95% Wilson intervals for the *simulation-estimated coverage probability*; they are not Bayesian credible intervals for H0.

| True H0 | Prior | Method | 68% coverage (95% Wilson interval) | 95% coverage (95% Wilson interval) |
|---:|:---:|:---|:---:|:---:|
| 55 | [50, 90] | Uncorrected | 0.577 [0.537, 0.616] | 0.955 [0.935, 0.969] |
| 55 | [50, 90] | Corrected | 0.642 [0.603, 0.679] | 0.970 [0.953, 0.981] |
| 55 | [40, 100] | Uncorrected | 0.722 [0.684, 0.756] | 0.965 [0.947, 0.977] |
| 55 | [40, 100] | Corrected | 0.773 [0.738, 0.805] | 0.977 [0.961, 0.986] |
| 55 | [30, 120] | Uncorrected | 0.670 [0.631, 0.706] | 0.955 [0.935, 0.969] |
| 55 | [30, 120] | Corrected | 0.723 [0.686, 0.758] | 0.965 [0.947, 0.977] |
| 70 | [50, 90] | Uncorrected | 0.948 [0.928, 0.963] | 1.000 [0.994, 1.000] |
| 70 | [50, 90] | Corrected | 0.965 [0.947, 0.977] | 1.000 [0.994, 1.000] |
| 70 | [40, 100] | Uncorrected | 0.895 [0.868, 0.917] | 0.998 [0.991, 1.000] |
| 70 | [40, 100] | Corrected | 0.915 [0.890, 0.935] | 1.000 [0.994, 1.000] |
| 70 | [30, 120] | Uncorrected | 0.813 [0.780, 0.842] | 0.987 [0.974, 0.993] |
| 70 | [30, 120] | Corrected | 0.793 [0.759, 0.824] | 0.988 [0.976, 0.994] |
| 85 | [50, 90] | Uncorrected | 0.532 [0.492, 0.571] | 1.000 [0.994, 1.000] |
| 85 | [50, 90] | Corrected | 0.413 [0.375, 0.453] | 1.000 [0.994, 1.000] |
| 85 | [40, 100] | Uncorrected | 0.910 [0.884, 0.930] | 1.000 [0.994, 1.000] |
| 85 | [40, 100] | Corrected | 0.772 [0.736, 0.803] | 1.000 [0.994, 1.000] |
| 85 | [30, 120] | Uncorrected | 0.858 [0.828, 0.884] | 0.998 [0.991, 1.000] |
| 85 | [30, 120] | Corrected | 0.828 [0.796, 0.856] | 0.993 [0.983, 0.997] |

These intervals reinforce that several 68% coverage estimates are far from the nominal 0.68, beyond what is plausibly explained by Monte Carlo noise alone. However, the table does not by itself identify the cause; prior boundaries, the finite grid, the simplified likelihood, and the selection-conditioning setup must be assessed separately.

## Paired change in absolute error

The metric below is mean(|error_corrected| - |error_uncorrected|). Its Monte Carlo standard error (MCSE) quantifies simulation uncertainty in that mean for the fixed design.

| True H0 | Prior | Mean paired change | MCSE | Interpretation within this toy design |
|---:|:---:|---:|---:|---|
| 55 | [50, 90] | −1.884 | 0.042 | Correction lowers mean absolute error |
| 55 | [40, 100] | −2.022 | 0.076 | Correction lowers mean absolute error |
| 55 | [30, 120] | −2.913 | 0.126 | Correction lowers mean absolute error |
| 70 | [50, 90] | +0.081 | 0.061 | Small positive difference; uncertain relative to MCSE |
| 70 | [40, 100] | +1.277 | 0.135 | Correction raises mean absolute error |
| 70 | [30, 120] | +1.580 | 0.173 | Correction raises mean absolute error |
| 85 | [50, 90] | +1.395 | 0.048 | Correction raises mean absolute error |
| 85 | [40, 100] | +2.250 | 0.122 | Correction raises mean absolute error |
| 85 | [30, 120] | +2.398 | 0.123 | Correction raises mean absolute error |

## Interpretation

1. **There is no uniform improvement.** In this simulation, selection normalization improves mean absolute error at H0=55 for all tested priors, but worsens it at H0=85 for all tested priors. At H0=70, results depend on prior support.
2. **Prior width changes uncertainty materially.** Wider support yields substantially wider intervals and, in several cells, larger RMSE. This experiment does not isolate a universal causal effect of prior width; it demonstrates sensitivity in this particular setup.
3. **Nominal 68% coverage is often far from 68%.** Some cells substantially under-cover, while others over-cover. The model and simulation/analysis assumptions need further calibration before any coverage claim is made.
4. **95% coverage near 100% is not proof of excellent calibration.** Intervals can be too wide, and the prior-boundary behavior plus the discrete finite grid need further study.
5. The H0=55, 70, and 85 cells must remain separate. A pooled score could hide the change in direction across truths.

## Limitations and next validation

- This is a deliberately simplified two-host, low-redshift Gaussian-distance model.
- The mock observations are selected by truncating observed distance, while inference uses a simplified selection normalization. Its assumptions must be documented and mathematically verified against the actual generative process.
- The study does not include a realistic gravitational-wave distance likelihood, sky localization, a real galaxy catalogue, redshift uncertainty, peculiar velocities, or an LVK event population.
- The nominal host weights and detection rule are fixed; no generalization to real catalogues is justified.
- 600 replicates per cell provide finite Monte Carlo precision. MCSEs above describe paired mean absolute-error differences only; coverage has its own binomial uncertainty.

**Recommended next step:** audit the generative and inferential selection model algebraically, then add a unit test for a case where the detected host mixture can be computed analytically. After that, run a small selection-threshold/noise sweep and check calibration before scaling the experiment. Do not claim that selection correction is universally better based on these results.
