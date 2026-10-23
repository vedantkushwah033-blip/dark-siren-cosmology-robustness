# Prior-boundary sensitivity audit: results from 12,600 posterior fits

**Workflow:** [Prior-boundary sensitivity audit, run 39997110013](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39997110013)  
**Artifact:** `prior-boundary-sensitivity` (per-run and summary CSVs)  
**Pinned simulator revision:** `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`  
**Seed:** 20261024; 300 paired replicates per true H0 and scenario  
**Grid:** 1,801 points on H0 = 30–120 km/s/Mpc  
**Status:** exploratory diagnostic, not a real-data result.

## Design

For each of seven injected values (H_0\in\{55,60,65,70,75,80,85\}), the audit generated 300 simulated catalogues/observations. Each same simulated observation was analysed under three prior ranges and two host scenarios:

- baseline prior: 50–90 km/s/Mpc;
- wider prior: 40–100 km/s/Mpc;
- very wide prior: 30–120 km/s/Mpc;
- known-host inference;
- unknown-host inference using equal weights for all catalogue galaxies.

This produces 7 × 300 × 3 × 2 = **12,600 posterior fits** and 42 summary rows. The successful workflow validated the row counts and uploaded both CSV files.

## Key unknown-host results

All errors and interval widths are in km/s/Mpc; coverage is a proportion.

| True H0 | Prior | Mean signed error | RMSE of median | 68% coverage | Mean 68% width |
|---:|---|---:|---:|---:|---:|
| 55 | 50–90 | +8.870 | 12.383 | 48.3% | 14.988 |
| 55 | 40–100 | +7.989 | 14.290 | 62.0% | 20.348 |
| 55 | 30–120 | +8.795 | 18.284 | 64.3% | 24.476 |
| 70 | 50–90 | +2.882 | 5.976 | 94.3% | 21.473 |
| 70 | 40–100 | +5.311 | 10.147 | 93.3% | 27.227 |
| 70 | 30–120 | +9.885 | 17.599 | 66.0% | 32.039 |
| 85 | 50–90 | −8.936 | 9.304 | 78.7% | 23.110 |
| 85 | 40–100 | −4.276 | 6.777 | 87.7% | 29.665 |
| 85 | 30–120 | +3.632 | 11.901 | 89.0% | 35.836 |

Coverage is estimated from 300 replicates in each cell. Its Monte Carlo standard error is approximately 2.5–2.9 percentage points when coverage is between 0.5 and 0.7, and smaller near 0 or 1. These estimates therefore have meaningful simulation noise, but the large shifts across prior configurations warrant further investigation.

## Paired prior comparison for unknown-host inference

The same simulated observations were analysed under each prior. Below are paired changes relative to the baseline 50–90 prior. “Coverage change” is the change in 68% coverage, in percentage points; “width change” is the change in mean 68% interval width.

| True H0 | Wider prior: error change | Wider prior: RMSE change | Wider prior: coverage change | Very-wide prior: error change | Very-wide prior: RMSE change | Very-wide prior: coverage change |
|---:|---:|---:|---:|---:|---:|---:|
| 55 | −0.88 | +1.91 | +13.7 pp | −0.08 | +5.90 | +16.0 pp |
| 60 | +0.20 | +2.84 | +4.3 pp | +2.00 | +7.75 | +4.0 pp |
| 65 | +1.12 | +3.45 | −16.0 pp | +4.10 | +9.66 | −28.7 pp |
| 70 | +2.43 | +4.17 | −1.0 pp | +7.00 | +11.62 | −28.3 pp |
| 75 | +2.80 | +3.29 | −0.3 pp | +7.89 | +10.64 | −17.3 pp |
| 80 | +3.85 | +0.65 | −0.3 pp | +10.39 | +7.74 | −2.3 pp |
| 85 | +4.66 | −2.53 | +9.0 pp | +12.57 | +2.60 | +10.3 pp |

Error change is the change in mean posterior-median error; RMSE change is calculated from the paired set of posterior-median errors. These are descriptive contrasts. No multiplicity adjustment was made across the 7 truth values and multiple metrics.

## Interpretation

1. **The unknown-host results are sensitive to prior boundaries.** At true H0=70, widening the prior from 50–90 to 30–120 changes mean signed error from +2.882 to +9.885 and RMSE from 5.976 to 17.599. At true H0=85, the mean signed error changes sign.
2. **Wider intervals do not guarantee better calibration.** At H0=70, the very-wide prior produces a larger mean 68% interval width but lower empirical 68% coverage than the baseline prior.
3. **Known-host results also change with the prior.** This means prior sensitivity is not exclusive to unknown-host inference, although its pattern and magnitude differ.
4. These results are conditional on this simulator, catalogue construction, distance-noise model, equal-weight unknown-host scenario, grid, and prior choices. They do not establish the behaviour of real gravitational-wave events or real galaxy catalogues.
5. Because the posterior likelihood and simulation assumptions have not yet been fully validated against a real selection function, we should not conclude that prior choice is the sole cause. The next step is to check the prior handling and posterior density numerically, and to compare known-host/unknown-host models under controlled likelihood cases.

## Recommended next checks

- Verify that `posterior_h0` applies the prior support exactly as intended and that posterior density is zero outside the support.
- For representative single-observation cases, compare the numerical posterior with an independently computed likelihood on the same grid.
- Record the fraction of posterior mass close to each prior boundary, rather than relying only on posterior medians and coverage.
- Add paired bootstrap intervals for prior contrasts if the same findings remain after implementation checks.
- Keep this audit separate from any claims about real LVK observations until the astrophysical selection model is validated.

## Reproducibility

The workflow artifact contains `per_run.csv` and `summary.csv`. The experiment is automated at [`.github/workflows/prior-boundary-sensitivity.yml`](../.github/workflows/prior-boundary-sensitivity.yml), with the driver in [`scripts/run_prior_boundary_sensitivity.py`](../scripts/run_prior_boundary_sensitivity.py).


## Paired-bootstrap uncertainty for prior contrasts

The automated paired-contrast workflow completed successfully in [run 39998157948](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39998157948). It resampled the 300 paired replicates 5,000 times for each truth/scenario/prior comparison and wrote `paired_contrasts.csv` to the workflow artifact.

Selected unknown-host contrasts (alternative prior minus baseline 50–90 prior):

| True H0 | Alternative prior | Change in mean median-error (95% bootstrap CI) | Change in RMSE (95% bootstrap CI) | Change in 68% coverage (95% bootstrap CI) |
|---:|---|---:|---:|---:|
| 55 | 40–100 | −0.881 [−1.258, −0.474] | +1.907 [+1.497, +2.300] | +13.7 pp [+8.3, +19.0] |
| 55 | 30–120 | −0.075 [−0.975, +0.893] | +5.901 [+4.668, +7.120] | +16.0 pp [+10.0, +22.0] |
| 65 | 30–120 | +4.102 [+2.995, +5.245] | +9.660 [+8.551, +10.739] | −28.7 pp [−34.0, −23.3] |
| 70 | 30–120 | +7.003 [+5.871, +8.151] | +11.624 [+10.584, +12.633] | −28.3 pp [−33.7, −23.3] |
| 85 | 40–100 | +4.660 [+4.329, +4.990] | −2.527 [−2.931, −2.146] | +9.0 pp [+6.0, +12.3] |
| 85 | 30–120 | +12.568 [+11.492, +13.649] | +2.597 [+2.089, +3.100] | +10.3 pp [+7.0, +14.0] |

These intervals reflect Monte Carlo variability under paired resampling of the chosen simulator and random seeds. They do not include uncertainty from model misspecification, and the many truth/metric comparisons are exploratory rather than a pre-registered family of confirmatory tests.

The effect is not a uniform improvement or deterioration: widening the prior changes the posterior error, interval coverage, and RMSE differently at different injected values. This is evidence that the reported conditional performance is prior-sensitive in this simplified setup, not evidence that a particular prior range is correct for real dark-siren analyses.
