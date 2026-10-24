# Factorial sensitivity experiment: audited results (2026-10-24)

## Status and provenance

The factorial-sensitivity GitHub Actions run completed successfully:

- Workflow run: [40142824335](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40142824335)
- Workflow commit: `8863e309a558b88cc54a8af27ec36ccdd9a38aa0`
- Pinned simulator commit: `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`
- Artifact: `factorial-sensitivity`, ID `12579731455`
- Artifact SHA-256: `ea14c964859531758a8f6da4ed1ba2db2bb3cd498bf022880f1d3bf4547af56c`
- Artifact size: 1,080,480 bytes; expiry shown by GitHub: 2027-01-22
- Main experiment: 9 scenarios × 500 mock universes = 4,500 raw rows.
- Missing-host controls: 6 × 500 = 3,000 rows.
- Full-population oracle: 2 × 500 = 1,000 rows.
- Known-true-host oracle: 1 × 500 = 500 rows, with paired calibration contrasts against the full-population oracle.

The job passed the main matrix, structural validation, independent main-factorial audit, contrast analysis, report generation, missing-host control, and full-population oracle checks. The artifact was uploaded successfully.

This is a controlled synthetic diagnostic, not a real gravitational-wave analysis or an observational measurement of (H_0).

## Main factorial results

All incomplete-catalogue runs retained exactly 350 of 500 galaxies (70%). The table reports the project's defined `median_bias` (median across runs of each run's posterior-median error), RMSE of posterior medians, and empirical interval coverage.

| Scenario | 68% coverage | 95% coverage | Median of median errors | Mean median error | Median RMSE | Mean 68% width | Mean 95% width |
|---|---:|---:|---:|---:|---:|---:|
| Complete, uniform host / uniform inference | 0.956 | 1.000 | 3.811 | 1.969 | 6.095 | 22.405 | 35.888 |
| Complete, luminosity host / luminosity inference | 0.952 | 1.000 | 3.907 | 1.995 | 6.328 | 22.243 | 35.690 |
| Complete, luminosity host / uniform inference | 0.932 | 0.996 | 4.124 | 2.144 | 6.241 | 22.389 | 35.792 |
| Random 70%, uniform / uniform | 0.932 | 0.996 | 3.480 | 1.394 | 6.308 | 22.181 | 35.641 |
| Random 70%, luminosity / luminosity | 0.936 | 0.994 | 3.068 | 1.646 | 6.163 | 22.154 | 35.582 |
| Random 70%, luminosity / uniform | 0.940 | 0.994 | 3.712 | 1.741 | 6.246 | 22.279 | 35.721 |
| Redshift-dependent 70%, uniform / uniform | 0.940 | 0.996 | 2.295 | 0.796 | 5.899 | 22.525 | 35.770 |
| Redshift-dependent 70%, luminosity / luminosity | 0.938 | 0.992 | 2.566 | 0.902 | 5.999 | 22.461 | 35.731 |
| Redshift-dependent 70%, luminosity / uniform | 0.942 | 0.998 | 2.882 | 1.152 | 5.878 | 22.674 | 35.906 |

Widths, errors, and RMSE are in km s⁻¹ Mpc⁻¹. “Median of median errors” is the median across runs of (posterior median − true H0); “mean median error” is the arithmetic mean of that same per-run error. The separate `mean_bias` field in the simulator is the mean error of the posterior mean, not the mean error of the posterior median. The 68% coverage estimates have approximate Monte Carlo standard errors around 1.0–1.1 percentage points; Wilson 95% intervals for individual scenario coverage are roughly 90.6–97.1% for nominal 68% coverage and 98.0–100% for nominal 95% coverage. With 500 repetitions, small differences are not precise.

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

## Known-true-host oracle calibration and paired contrast

The updated workflow adds a deliberately idealized control: inference is given the actual simulated host redshift for each event. The raw artifact contains 500 known-host runs. Its summary is:

| Diagnostic | Known-host oracle |
|---|---:|
| Nominal 68% interval coverage | 0.750 (375/500) |
| Wilson 95% interval for 68% coverage | [0.710, 0.786] |
| Nominal 95% interval coverage | 0.986 (493/500) |
| Wilson 95% interval for 95% coverage | [0.971, 0.993] |
| Median of posterior-median errors | 2.562 km s⁻¹ Mpc⁻¹ |
| Mean posterior-median error | 2.019 km s⁻¹ Mpc⁻¹ |
| Median RMSE | 7.257 km s⁻¹ Mpc⁻¹ |
| Mean 68% interval width | 17.025 km s⁻¹ Mpc⁻¹ |
| Mean 95% interval width | 29.862 km s⁻¹ Mpc⁻¹ |

The known-host intervals are narrower than the full-population oracle intervals, but they still over-cover relative to their nominal levels. This indicates that candidate-host mixture ambiguity is not the only source of overcoverage; the simplified likelihood, posterior construction, and estimator calibration also need investigation.

### Paired known-host versus full-population comparison

The two oracle runners use the same deterministic seed schedule, mock catalogues, host choices, observed distances, and distance uncertainties. The updated artifact stores those shared inputs, and the analysis script verifies exact event-level matching before resampling paired differences. Each contrast below is **known-host minus full-population oracle**; intervals use 50,000 paired bootstrap resamples.

| Metric | Full-population value | Paired difference | Paired bootstrap 95% interval |
|---|---:|---:|---:|
| 68% coverage | 0.944 | -0.194 | [-0.234, -0.154] |
| 95% coverage | 1.000 | -0.014 | [-0.026, -0.004] |
| Mean 68% width | 22.405 | -5.380 | [-5.637, -5.126] |
| Mean 95% width | 35.899 | -6.036 | [-6.399, -5.681] |
| Mean median error | 2.144 | -0.125 | [-0.664, 0.404] |
| Mean absolute median error | 5.342 | +0.804 | [+0.449, +1.162] |
| Mean squared median error | 37.237 | +15.425 | [+10.630, +20.304] |

The full-population mixture produces substantially wider intervals and higher coverage, but it does not reduce absolute or squared posterior-median error in this paired sample. The mean signed median-error difference is inconclusive. This is a simulation-only comparison: both oracle methods use hidden information and neither is directly implementable for a real dark siren with an unknown host.

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

The main factorial and missing-host contrast intervals were generated by the committed scripts with 5,000 resamples and fixed seeds; this audit did not independently reproduce every one of those intervals draw-for-draw. For the new paired known-host comparison, I independently recalculated all seven mean differences and all 50,000-resample bootstrap interval endpoints from the raw rows and matched the stored JSON to floating-point precision (largest endpoint discrepancy below 4×10⁻¹⁵).

## Metric-definition consistency fix

The main factorial defines `median_bias` as the median across runs of `posterior median - true H0`. The oracle scripts had previously been inheriting `summarize_runs`'s *mean* error under the same field name. The latest workflow run uses the corrected definition consistently across the missing-host, full-population, and known-host oracle summaries. Regression tests distinguish the median of errors from the mean and reject empty inputs.

## Limitations and next scientific steps

- This is a simplified synthetic model with a low-redshift distance relation, stylized luminosity proxies, a finite (H_0) grid/prior, and a simplified Gaussian distance likelihood.
- The redshift-dependent mask is a stress test, not a calibrated observational selection function.
- All contrasts are exploratory, and intervals are unadjusted for multiple comparisons.
- Coverage well above nominal levels is a calibration warning. The newly added known-true-host oracle is intended to diagnose whether overcoverage persists when candidate-host ambiguity is removed; its results remain pending until the updated workflow completes. Only then should the study expand to more realistic selection and cosmological models.
- Do not present these results as a real (H_0) measurement, proof that a weighting scheme is generally better, or evidence of publication-level validation.

## Source links

- [Completed workflow and logs](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40142824335)
- [Factorial workflow definition](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/.github/workflows/factorial-sensitivity.yml)
- [Main factorial runner](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/scripts/run_factorial_sensitivity.py)
- [Main independent auditor](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/scripts/audit_factorial_sensitivity_outputs.py)
- [Contrast analysis](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/scripts/analyze_factorial_contrasts.py)
- [Missing-host control runner](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/scripts/run_oracle_missing_host_control.py)
- [Full-population oracle runner](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/scripts/run_full_population_oracle.py)
- [Known-host oracle runner](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/scripts/run_known_host_oracle.py)
- [Paired oracle calibration analysis](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/main/scripts/analyze_known_host_oracle_contrasts.py)
