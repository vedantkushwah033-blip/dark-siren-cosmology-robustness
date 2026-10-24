# Factorial sensitivity experiment: audited results (2026-10-24)

## Status and provenance

The factorial-sensitivity GitHub Actions run completed successfully:

- Workflow run: [40141336573](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40141336573)
- Workflow commit: `97a63a5b8308ce09666fc48798e0bb67b59894df`
- Pinned simulator commit: `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`
- Artifact: `factorial-sensitivity`, ID `12579151793`
- Artifact SHA-256: `a14f99b12267128ad1838b19d8ca074bd93c6fa96a7038caccaca150262185b5`
- Artifact size: 979,620 bytes; expiry shown by GitHub: 2027-01-22
- Main experiment: 9 scenarios × 500 mock universes = 4,500 raw rows.
- Missing-host controls: 6 × 500 = 3,000 rows.
- Full-population oracle: 2 × 500 = 1,000 rows.

The job passed the main matrix, structural validation, independent main-factorial audit, contrast analysis, report generation, missing-host control, and full-population oracle checks. The artifact was uploaded successfully.

This is a controlled synthetic diagnostic, not a real gravitational-wave analysis or an observational measurement of (H_0).

## Main factorial results

All incomplete-catalogue runs retained exactly 350 of 500 galaxies (70%). The table reports the project's defined `median_bias` (median across runs of each run's posterior-median error), RMSE of posterior medians, and empirical interval coverage.

| Scenario | 68% coverage | 95% coverage | Median bias (km s⁻¹ Mpc⁻¹) | Median RMSE | Mean 68% width | Mean 95% width |
|---|---:|---:|---:|---:|---:|---:|
| Complete, uniform host / uniform inference | 0.956 | 1.000 | 3.811 | 6.095 | 22.405 | 35.888 |
| Complete, luminosity host / luminosity inference | 0.952 | 1.000 | 3.907 | 6.328 | 22.243 | 35.690 |
| Complete, luminosity host / uniform inference | 0.932 | 0.996 | 4.124 | 6.241 | 22.389 | 35.792 |
| Random 70%, uniform / uniform | 0.932 | 0.996 | 3.480 | 6.308 | 22.181 | 35.641 |
| Random 70%, luminosity / luminosity | 0.936 | 0.994 | 3.068 | 6.163 | 22.154 | 35.582 |
| Random 70%, luminosity / uniform | 0.940 | 0.994 | 3.712 | 6.246 | 22.279 | 35.721 |
| Redshift-dependent 70%, uniform / uniform | 0.940 | 0.996 | 2.295 | 5.899 | 22.525 | 35.770 |
| Redshift-dependent 70%, luminosity / luminosity | 0.938 | 0.992 | 2.566 | 5.999 | 22.461 | 35.731 |
| Redshift-dependent 70%, luminosity / uniform | 0.942 | 0.998 | 2.882 | 5.878 | 22.674 | 35.906 |

Widths are in km s⁻¹ Mpc⁻¹. RMSE is in the same units. The 68% coverage estimates have approximate Monte Carlo standard errors around 1.0–1.1 percentage points; with 500 repetitions, small differences are not precise.

### Interpretation

1. **The intervals substantially over-cover in this simulation.** Nominal 68% intervals cover the truth in 93.2–95.6% of runs, and nominal 95% intervals cover in 99.2–100%. The immediate result is not evidence of calibrated, informative intervals; it indicates that the posterior intervals are very conservative under this design.
2. **The posterior-median estimator is positively shifted.** The defined median bias is positive in all nine scenarios (about 2.29–4.12 km s⁻¹ Mpc⁻¹). The average error of posterior medians is smaller than the median error in several cases, suggesting skewness in the across-run error distribution.
3. **There is no strong evidence that matched luminosity weighting is universally superior.** Most matched-versus-uniform and matched-versus-mismatched bootstrap intervals for coverage, width, and squared error include zero. Avoid selecting one favorable metric as a general winner.
4. **The redshift-dependent scenarios look somewhat different from random missingness in this particular synthetic model.** For matched luminosity weights, the median-error contrast (redshift-dependent minus random) is about -0.743 km s⁻¹ Mpc⁻¹, with a bootstrap 95% interval of [-1.485, -0.015]. The squared-error contrast interval includes zero. This is an exploratory result among multiple unadjusted contrasts, not confirmatory evidence that redshift-dependent incompleteness improves inference.
5. **The very broad intervals need investigation.** Mean 68% widths are about 22.15–22.67 and mean 95% widths about 35.58–35.91 km s⁻¹ Mpc⁻¹. Prior range, low-redshift distance approximation, mixture structure, and the simulator's noise model may contribute. This design does not isolate the cause.

## Missing-host control

The six missing-host scenarios compare observed-catalogue-only inference, a missing-host term weighted using inference weights, and an oracle missing-host term weighted using the hidden simulated population.

| Missingness | Method | 68% coverage | 95% coverage | Mean median error | Median of median errors | Median RMSE |
|---|---|---:|---:|---:|---:|---:|
| Random | Observed catalogue only | 0.942 | 1.000 | 2.174 | 3.860 | 6.070 |
| Random | Missing term, inference weights | 0.948 | 1.000 | 2.005 | 4.167 | 6.312 |
| Random | Missing term, oracle weights | 0.934 | 0.996 | 2.139 | 4.142 | 6.227 |
| Redshift-dependent | Observed catalogue only | 0.936 | 0.996 | 0.961 | 2.661 | 5.959 |
| Redshift-dependent | Missing term, inference weights | 0.946 | 0.996 | 1.622 | 3.271 | 6.069 |
| Redshift-dependent | Missing term, oracle weights | 0.938 | 0.994 | 1.736 | 3.899 | 6.228 |

The added missing-host term did not uniformly improve bias, RMSE, or coverage. Most unadjusted contrast intervals include zero. The oracle version uses simulated hidden-population information that a real incomplete survey does not have; it is a validation control, not a deployable estimator.

## Full-population oracle benchmark

Both oracle scenarios use all 500 galaxies and the true luminosity host weights, regardless of the simulated catalogue mask. The independently inspected raw artifact contains 1,000 rows. The posterior summaries for paired random and redshift-dependent masks were **exactly identical for every replicate** (zero difference in medians, means, posterior standard deviations, and interval endpoints), while host-observed diagnostics differed for 218 of 500 replicates.

That equality is the intended invariant: if inference ignores the observed mask and uses the same full population, changing only the mask must not change the posterior. This is a successful implementation control, not evidence that incomplete catalogues can be corrected perfectly in practice.

## Independent artifact checks performed after download

I downloaded and extracted the completed GitHub Actions artifact and independently recalculated the key summary statistics from the raw CSVs.

- Main factorial: 4,500 rows, 9 summaries, 500 rows per scenario; all numeric values finite.
- Missing-host control: 3,000 rows, 6 summaries, 500 rows per scenario; all posterior fields checked were finite.
- Full-population oracle: 1,000 rows, 2 summaries, 500 rows per scenario; posterior summaries were exactly equal across paired mask mechanisms.
- Recomputed core summaries (bias definitions as specified by each runner, RMSE, interval widths, coverage, host-observed fraction) matched the stored summaries to floating-point precision; the largest discrepancy was below (4\times10^{-15}).
- Main factorial realized completeness was exactly 1.0 for complete scenarios and exactly 0.7 for all incomplete scenarios.
- Oracle-control realized completeness was exactly 0.7 in all six missing-host scenarios.

The contrast bootstrap calculations in the artifact were generated by the committed analysis scripts with 5,000 resamples and fixed seeds. This post-download audit independently recomputed the raw-row summary metrics and oracle invariants; it did not independently reproduce every bootstrap interval draw-for-draw.

## Metric-definition consistency fix

The main factorial defines `median_bias` as the median across runs of `posterior median - true H0`. The oracle scripts had been inheriting `summarize_runs`'s *mean* error under the same field name. This inconsistency has now been corrected in the oracle runner code and regression tests have been added. The historical artifact above remains an immutable record of the original run; its oracle summary `median_bias` fields should be interpreted as mean errors until the updated workflow is rerun.

The updated oracle runners now expose a small helper with the consistent definition. New tests verify the difference between mean and median error and reject empty input. These changes require a fresh CI run and fresh oracle artifacts before updated summary values can be claimed.

## Limitations and next scientific steps

- This is a simplified synthetic model with a low-redshift distance relation, stylized luminosity proxies, a finite (H_0) grid/prior, and a simplified Gaussian distance likelihood.
- The redshift-dependent mask is a stress test, not a calibrated observational selection function.
- All contrasts are exploratory, and intervals are unadjusted for multiple comparisons.
- Coverage well above nominal levels is a calibration warning. Next steps should diagnose interval construction and posterior information, compare against a deliberately known/tractable generative case, and only then expand to more realistic selection and cosmological models.
- Do not present these results as a real (H_0) measurement, proof that a weighting scheme is generally better, or evidence of publication-level validation.

## Source links

- [Completed workflow and logs](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40141336573)
- [Factorial workflow definition](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/.github/workflows/factorial-sensitivity.yml)
- [Main factorial runner](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/scripts/run_factorial_sensitivity.py)
- [Main independent auditor](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/scripts/audit_factorial_sensitivity_outputs.py)
- [Contrast analysis](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/scripts/analyze_factorial_contrasts.py)
- [Missing-host control runner](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/scripts/run_oracle_missing_host_control.py)
- [Full-population oracle runner](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/scripts/run_full_population_oracle.py)
