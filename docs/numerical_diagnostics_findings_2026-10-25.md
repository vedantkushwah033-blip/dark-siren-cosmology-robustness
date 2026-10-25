# Numerical diagnostics findings — 2026-10-25

## Provenance

The values below were read from the uploaded `numerical-diagnostics-py3.12` artifact for successful workflow run [40217025570](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40217025570), commit `3c3823f9e5756df58a7cdb4f7f2e5264c28d344e`. The workflow completed successfully on Python 3.10, 3.11, and 3.12. The run contains posterior-grid stability and targeted convergence CSVs, plus the 1,000-replicate selected-likelihood PIT raw data, summary, and independent audit report.

## 1. Grid refinement: targeted upper-tail case

The targeted convergence diagnostic holds the synthetic observation and model setup fixed and increases the posterior grid from 1,001 to 64,001 points. The posterior normalization difference between trapezoidal and Simpson integration decreases by roughly four orders of magnitude, and the 95% upper endpoint stabilizes.

| Grid points | Posterior area-rule absolute difference | Posterior median | 95% upper endpoint | Change in upper endpoint from previous grid |
|---:|---:|---:|---:|---:|
| 1,001 | 3.38e-6 | 51.713915 | 60.848373 | — |
| 2,001 | 8.45e-7 | 51.713591 | 60.742249 | 0.106124 |
| 4,001 | 2.11e-7 | 51.713486 | 60.717920 | 0.024330 |
| 8,001 | 5.28e-8 | 51.713451 | 60.711905 | 0.006015 |
| 16,001 | 1.32e-8 | 51.713446 | 60.710428 | 0.001477 |
| 32,001 | 3.30e-9 | 51.713445 | 60.710055 | 0.000373 |
| 64,001 | 8.25e-10 | 51.713445 | 60.709962 | 0.000093 |

**Interpretation:** In this specific synthetic case, the posterior median is stable to the shown precision well before the upper 95% endpoint is equally stable. The coarse 1,001-point grid changes that endpoint by about 0.138 relative to the 64,001-point result; at 8,001 points the remaining difference is about 0.00194. This supports reporting grid sensitivity for tail endpoints rather than relying on posterior normalization alone.

The tested case uses a low-redshift host mixture, a 3% noise fraction, a fixed synthetic observed distance of about 119.47 Mpc, and reversed inference weights (0.3, 0.7). It is a numerical stress case, not an observed gravitational-wave event.

## 2. Grid support versus prior support

A separate stability diagnostic compares a 1,601-point baseline grid on [50, 90] with a 3,201-point finer grid on the same support, a wider numerical grid [30, 110] with the **same** prior support [50, 90], and a wider numerical grid with the **prior support itself widened** to [30, 110].

| Case | Numerical grid | Prior support | Median | Posterior SD | 68% interval | 95% interval |
|---|---|---|---:|---:|---:|---:|
| Baseline | 1,601 points on [50, 90] | [50, 90] | 70.415478 | 10.239961 | [59.723130, 82.573675] | [51.869119, 88.756257] |
| Finer baseline | 3,201 points on [50, 90] | [50, 90] | 70.415477 | 10.239958 | [59.723136, 82.573676] | [51.869120, 88.756256] |
| Wider grid, same prior | 6,401 points on [30, 110] | [50, 90] | 70.415477 | 10.239958 | [59.723136, 82.573676] | [51.869120, 88.756256] |
| Wider grid and wider prior | 6,401 points on [30, 110] | [30, 110] | 68.200082 | 20.939685 | [40.935987, 89.985520] | [31.720089, 105.078786] |

With prior support held fixed, the baseline-to-finer median difference is about (1.5\times10^{-6}), and widening only the numerical grid leaves the posterior summaries unchanged to floating-point precision. In contrast, widening the **prior support** changes the median by about 2.22 and doubles the posterior standard deviation from about 10.24 to 20.94.

**Interpretation:** The tested numerical integral is stable under grid refinement for this case, but inference is materially sensitive to prior support. Numerical convergence does not imply prior robustness. Prior-support sensitivity should be treated as a substantive modeling choice and reported separately from grid-convergence error.

## 3. Calibration artifact cross-check

The same successful run's selected-likelihood experiment used 1,000 replicates with seed `20261025`. Its raw records, summary, and independent audit report were archived. The audit passed. Across Python 3.10, 3.11, and 3.12, the recorded mean PIT was 0.4965598 (Monte Carlo SE 0.0088742), KS statistic 0.0248506, and descriptive KS p-value 0.5587675; all three environments produced identical summary values for this seed.

This is evidence of reproducible synthetic calculations, not proof of exact calibration or astrophysical validity.

## 4. Limitations and next tests

- These are controlled finite-mixture diagnostics; they do not validate real detector selection functions or galaxy-catalogue completeness.
- The targeted convergence case is one deliberately selected stress case. It should be extended across a prespecified matrix of noise levels, host configurations, true values, and prior supports.
- The prior-support comparison is a sensitivity result, not a reason by itself to prefer one support. A defensible prior needs scientific justification.
- A broader study should predefine numerical tolerances for posterior medians, credible endpoints, and integration-rule differences, then evaluate every prespecified case.
