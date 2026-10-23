# Noise-scale consistency audit: results

**Run:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39974976415  
**Commit:** `62a238c92454c2c4a9708e50c283c925eb6f2b8f`  
**Status:** workflow and test suites passed; artifact contains 14 summary rows and 14,000 posterior fits.  
**Scope:** isolated one-host synthetic experiment, not a realistic gravitational-wave or galaxy-catalogue analysis.

## Question

The toy simulator generates a distance uncertainty proportional to the true distance, `sigma_true = 0.15 * d_true`. This audit compares two likelihood specifications using the same observed distance for each replicate:

- **Fixed reported uncertainty:** treats the generated distance uncertainty as fixed while evaluating candidate (H_0) values.
- **H0-dependent uncertainty:** sets (sigma(H_0)=0.15,cz/H_0) and includes the Gaussian (-\log \sigma(H_0)) normalization term.

There are 7 injected (H_0) values, 1,000 replicates per value, and two modes. The results below use a flat prior over 30–120 km/s/Mpc. Coverage Monte Carlo standard errors are included in the artifact's `summary.csv`.

## Results

All errors and interval widths are in km/s/Mpc. Coverage values are proportions.

| True H0 | Likelihood mode | Mean signed median error | Median RMSE | 68% coverage | 95% coverage | Mean 68% width | Mean 95% width |
|---:|---|---:|---:|---:|---:|---:|---:|
| 55 | Fixed reported sigma | 4.454 | 10.972 | 65.7% | 92.9% | 20.277 | 41.612 |
| 55 | H0-dependent sigma | 2.739 | 9.552 | 67.9% | 94.6% | 16.602 | 32.556 |
| 60 | Fixed reported sigma | 4.741 | 11.933 | 63.6% | 92.5% | 21.577 | 43.388 |
| 60 | H0-dependent sigma | 3.052 | 10.791 | 65.9% | 92.7% | 18.118 | 35.577 |
| 65 | Fixed reported sigma | 4.246 | 11.643 | 68.4% | 94.0% | 22.368 | 44.452 |
| 65 | H0-dependent sigma | 2.678 | 10.888 | 69.4% | 93.9% | 19.411 | 38.064 |
| 70 | Fixed reported sigma | 4.525 | 11.954 | 67.2% | 94.9% | 23.417 | 45.502 |
| 70 | H0-dependent sigma | 3.069 | 11.440 | 68.0% | 94.9% | 20.841 | 40.760 |
| 75 | Fixed reported sigma | 3.947 | 11.916 | 66.2% | 95.0% | 23.749 | 45.312 |
| 75 | H0-dependent sigma | 2.743 | 11.808 | 66.9% | 95.1% | 21.931 | 42.659 |
| 80 | Fixed reported sigma | 3.717 | 11.933 | 66.9% | 95.1% | 23.978 | 44.754 |
| 80 | H0-dependent sigma | 2.771 | 12.100 | 67.1% | 94.6% | 22.899 | 44.214 |
| 85 | Fixed reported sigma | 2.995 | 11.465 | 68.9% | 95.5% | 23.843 | 43.687 |
| 85 | H0-dependent sigma | 2.301 | 11.814 | 68.4% | 96.1% | 23.594 | 45.090 |

## Interpretation

1. **The H0-dependent uncertainty specification reduced positive mean signed median error in all seven truth conditions.** The reduction ranged from about 0.69 to 1.72 km/s/Mpc. Bias remained positive in every condition, so this does not resolve the estimator's bias.
2. **68% interval coverage is broadly near nominal under the H0-dependent specification** (65.9–69.4%). At (H_0=70), it is 68.0% from 1,000 replicates. The 95% coverage ranges from 92.7–96.1%, with some deviation from nominal 95%.
3. **The H0-dependent intervals are generally narrower**, especially at lower and middle truth values. This is consistent with a likelihood that accounts for the changing scale of measurement uncertainty, but narrower intervals are not automatically better.
4. **RMSE is not uniformly improved.** It decreases for truths 55–75 but is slightly worse at 80 and 85. At (H_0=70), RMSE decreases from 11.954 to 11.440 km/s/Mpc.
5. This comparison is a controlled diagnostic of a one-host likelihood. It does not validate the unknown-host galaxy-mixture likelihood, the catalogue-generation process, or the original project's broader conclusions. The synthetic distance errors are large (15% fractional standard deviation), and the prior range remains bounded.
6. These summaries are descriptive; no multiplicity adjustment or formal paired uncertainty interval for differences is reported here. Each likelihood mode uses the same generated observation per replicate, so future comparisons should explicitly compute paired replicate-level differences and their uncertainty.

## Conclusion

The noise-scale specification is a material contributor to the behavior of this simplified single-host experiment, but it is **not a complete explanation** for the earlier unknown-host catalogue results. The positive median bias persists even when the likelihood uncertainty depends on (H_0). The next step is to investigate the estimator's transformation and likelihood behavior and then build a matched host-mixture benchmark with explicit paired comparisons.

Do not describe these results as an observational (H_0) measurement or as validation of an LVK analysis.
