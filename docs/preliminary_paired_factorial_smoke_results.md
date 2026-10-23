# Preliminary paired factorial smoke-test results

## Status and reproducibility

The 30-replicate smoke experiment completed successfully on 2026-10-23.

- Workflow: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39955059874
- Workflow commit: `1738e09d242e6bfdec3beaedfbd465cdc518983e`
- Pinned simulator commit: `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`
- Design: 9 scenarios, 30 simulated universes per scenario, 270 per-run rows, 49 paired contrast rows.
- True value: $H_0=70$ km s$^{-1}$ Mpc$^{-1}$; 500 synthetic galaxies per parent catalogue.
- Outputs: `scenario_summary.csv`, `per_run_results.csv`, `paired_contrasts.csv`, and `scenario_summary.json`.

This is a smoke test to verify the experiment pipeline and inspect preliminary behavior. Thirty replicates are not enough for stable coverage estimates or definitive scientific conclusions.

## Descriptive scenario summaries

| Scenario | Median absolute error/bias summary (km s$^{-1}$ Mpc$^{-1}$) | RMSE of posterior median | 68% coverage | Host observed |
|---|---:|---:|---:|---:|
| Complete, uniform host generation and inference | 3.08 | 5.87 | 93.3% | 100% |
| Random 70%, uniform generation and inference | 3.15 | 5.92 | 93.3% | 73.3% |
| Redshift-dependent 70%, uniform generation and inference | 2.33 | 5.60 | 93.3% | 76.7% |
| Complete, luminosity generation and matched luminosity inference | 5.00 | 6.47 | 96.7% | 100% |
| Complete, luminosity generation and uniform inference | 6.04 | 6.43 | 93.3% | 100% |
| Random 70%, luminosity generation and matched luminosity inference | 5.20 | 6.47 | 96.7% | 83.3% |
| Random 70%, luminosity generation and uniform inference | 6.03 | 6.41 | 93.3% | 83.3% |
| Redshift-dependent 70%, luminosity generation and matched luminosity inference | 4.51 | 6.00 | 96.7% | 76.7% |
| Redshift-dependent 70%, luminosity generation and uniform inference | 4.69 | 6.00 | 93.3% | 76.7% |

The bias-summary column is the simulator's `median_bias` metric; the separately stored `mean_median_error` is the mean signed posterior-median error. Do not conflate these metrics.

## Paired contrasts: preliminary observations

- Under uniform host generation and uniform inference, random 70% incompleteness changed the mean signed posterior-median error by about +0.08 km s$^{-1}$ Mpc$^{-1}$ relative to the complete catalogue; the bootstrap 95% interval was approximately [-0.01, +0.17]. This small pilot does not establish a meaningful effect.
- Under uniform host generation, redshift-dependent 70% incompleteness changed the mean signed error by about -0.61 km s$^{-1}$ Mpc$^{-1}$ relative to the complete catalogue; the bootstrap 95% interval was approximately [-0.72, -0.49]. This is a result for the specific toy selection rule, not a general prediction for real catalogues.
- With luminosity-generated hosts, uniform inference weighting produced a higher 68% coverage indicator rate than matched luminosity inference in this small sample (93.3% versus 96.7% for matched weighting, in the paired scenarios); the per-run differences are discrete and uncertain. The smoke test does not establish which weighting strategy is better.
- The observed-host fractions differ from the nominal 70% catalogue completeness because the host is a single random draw and masks are stochastic transformations of each parent catalogue.

## Important limitations

1. These are synthetic data and one fixed true $H_0$, not a measurement or reanalysis of LVK events.
2. The simulator uses simplified low-redshift distance inference and synthetic luminosity proxies.
3. Coverage estimates from 30 replicates are noisy. For example, a 93.3% rate means 28 of 30 intervals covered the truth; it should not be interpreted as a precise calibration estimate.
4. The bootstrap intervals describe paired Monte Carlo contrasts under this setup. They do not account for uncertainty about whether the toy model resembles real astrophysics.
5. This pilot is useful for debugging and planning, not for making the paper's final claims.

## Next step

Run the prespecified 300-replicate version, inspect its artifacts and validation logs, then assess whether the main contrasts are stable. If conclusions change with replication, report that sensitivity rather than selecting a favorable run.
