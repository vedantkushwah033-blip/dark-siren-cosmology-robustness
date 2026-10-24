# Prior-support sensitivity audit: results

**Workflow:** [Prior-support sensitivity audit](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40073983427)  
**Commit evaluated:** `c269818aad6081686b742ddc375f7eab2a69e0c5`  
**Status:** experiment, focused tests, and CSV artifact upload passed on Python 3.10, 3.11, and 3.12. The run used 300 replicates per fixed truth and seed 20261024.

This is a synthetic two-host model sensitivity experiment, not an observational gravitational-wave result.

## Design

The same 300 detected observations at each truth were reused across all prior intervals. The truths were (H_0=52,70,88); prior supports were [50,90], [40,100], and [30,110]. Both uncorrected and selection-normalized posteriors were evaluated on a 4,001-point grid. The design uses host redshifts 0.02 and 0.04, weights 0.7 and 0.3, Gaussian distance noise 20 Mpc, and detection interval (0<d_{obs}<180) Mpc.

Coverage intervals are Wilson score intervals for binomial proportions. They quantify Monte Carlo uncertainty in empirical coverage, not uncertainty about which prior is scientifically appropriate.

## Empirical 68% coverage

Each entry is empirical coverage; the Wilson 95% interval is in parentheses.

| True (H_0) | Prior | Uncorrected | Selection-corrected |
|---:|---|---:|---:|
| 52 | [50,90] | 0.000 (0.000–0.013) | 0.000 (0.000–0.013) |
| 52 | [40,100] | 0.707 (0.653–0.755) | 0.727 (0.674–0.774) |
| 52 | [30,110] | 0.667 (0.612–0.718) | 0.707 (0.653–0.755) |
| 70 | [50,90] | 0.947 (0.915–0.967) | 0.950 (0.919–0.969) |
| 70 | [40,100] | 0.887 (0.846–0.918) | 0.910 (0.872–0.937) |
| 70 | [30,110] | 0.823 (0.776–0.862) | 0.823 (0.776–0.862) |
| 88 | [50,90] | 0.043 (0.025–0.073) | 0.043 (0.025–0.073) |
| 88 | [40,100] | 0.810 (0.762–0.850) | 0.587 (0.530–0.641) |
| 88 | [30,110] | 0.923 (0.888–0.948) | 0.873 (0.831–0.906) |

## RMSE and bias

| True (H_0) | Prior | Method | Bias | RMSE |
|---:|---|---|---:|---:|
| 52 | [50,90] | Uncorrected | 13.988 | 14.999 |
| 52 | [50,90] | Corrected | 11.917 | 13.071 |
| 52 | [40,100] | Uncorrected | 10.881 | 12.668 |
| 52 | [40,100] | Corrected | 7.895 | 10.499 |
| 52 | [30,110] | Uncorrected | 11.563 | 14.024 |
| 52 | [30,110] | Corrected | 7.962 | 12.027 |
| 70 | [50,90] | Uncorrected | 1.602 | 5.899 |
| 70 | [50,90] | Corrected | 0.147 | 6.108 |
| 70 | [40,100] | Uncorrected | 2.155 | 8.743 |
| 70 | [40,100] | Corrected | -0.285 | 9.583 |
| 70 | [30,110] | Uncorrected | 0.312 | 14.837 |
| 70 | [30,110] | Corrected | -2.997 | 17.419 |
| 88 | [50,90] | Uncorrected | -13.968 | 15.513 |
| 88 | [50,90] | Corrected | -15.351 | 17.017 |
| 88 | [40,100] | Uncorrected | -12.993 | 17.247 |
| 88 | [40,100] | Corrected | -15.216 | 19.961 |
| 88 | [30,110] | Uncorrected | -11.491 | 19.508 |
| 88 | [30,110] | Corrected | -13.874 | 22.532 |

## Findings

1. **The original [50,90] prior strongly constrains edge truths.** At (H_0=52), neither method's 68% interval covered the truth in these 300 draws. At (H_0=88), coverage was only 4.3%. The upper Wilson bounds remain below 0.08, so ordinary Monte Carlo fluctuation alone is not a plausible explanation for coverage near the nominal 68% level.
2. **Wider support substantially changes edge-truth coverage.** At (H_0=52), the wider priors yielded approximately 67–73% coverage. At (H_0=88), coverage increased to 81–92% for the uncorrected method and 59–87% for the corrected method, depending on support.
3. **Wider support is not uniformly better.** At (H_0=70), 68% coverage fell as the prior widened, while RMSE increased, especially for [30,110]. At (H_0=88), selection correction had worse RMSE than the uncorrected method for all three prior intervals in this experiment.
4. **Selection normalization is not uniformly beneficial for fixed-truth metrics.** It reduces RMSE for (H_0=52) in these runs, but increases RMSE at (H_0=70) and (H_0=88). This is consistent with the earlier paired experiments: the correction should be assessed across parameter space, not promoted as an automatic improvement for every truth and loss function.

## Interpretation limits

This experiment demonstrates prior-support sensitivity in the specified toy model; it does not prove that any one prior is scientifically correct. The observations are synthetic, the host mixture is simplified, and the model does not represent a calibrated gravitational-wave detector or galaxy survey. Coverage targets are frequentist properties under repeated sampling, whereas a prior is part of the Bayesian model. A wider prior may reduce truncation near one edge while increasing posterior uncertainty and RMSE elsewhere.

The later code revision adds a paired posterior-median-difference CSV to the artifact output. The numerical summaries above come from the completed run identified at the top; the new paired-contrast artifact must be verified from the newer workflow before being cited as validated.

## Reproduction

```bash
python scripts/run_prior_support_sensitivity.py --replicates 300 --seed 20261024 --grid-size 4001 --output-dir results/prior-support-sensitivity
python -m pytest tests/test_prior_support_sensitivity.py -q
```


## Paired posterior-median contrasts

The successful artifact set from workflow [40074060836](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40074060836) includes `paired_prior_differences.csv`. It compares each observation's posterior median under a wider prior with its median under [50,90]. Negative values mean the wider-prior median is lower; positive values mean it is higher. These are paired per-observation contrasts, not independent-sample differences.

| True H0 | Method | Comparison prior | Mean signed difference | Mean absolute difference | 90th percentile absolute difference |
|---:|---|---|---:|---:|---:|
| 52 | Uncorrected | [40,100] | -3.106 | 4.133 | 11.781 |
| 52 | Corrected | [40,100] | -4.022 | 4.791 | 14.320 |
| 52 | Uncorrected | [30,110] | -2.425 | 5.964 | 20.741 |
| 52 | Corrected | [30,110] | -3.955 | 6.203 | 22.768 |
| 70 | Uncorrected | [40,100] | +0.553 | 3.640 | 6.917 |
| 70 | Corrected | [40,100] | -0.432 | 4.275 | 7.765 |
| 70 | Uncorrected | [30,110] | -1.291 | 10.358 | 26.267 |
| 70 | Corrected | [30,110] | -3.144 | 11.224 | 33.404 |
| 88 | Uncorrected | [40,100] | +0.975 | 5.567 | 9.553 |
| 88 | Corrected | [40,100] | +0.135 | 6.247 | 13.398 |
| 88 | Uncorrected | [30,110] | +2.477 | 10.523 | 22.430 |
| 88 | Corrected | [30,110] | +1.476 | 11.038 | 25.720 |

### How to read these contrasts

- For truths near the lower prior edge ((H_0=52)), widening support tends on average to move the median downward, but the distribution of paired changes is broad and can include positive changes for individual observations.
- At (H_0=70), the wider [30,110] prior yields much larger median movement than [40,100], with 90th-percentile absolute changes above 26 in both methods.
- At (H_0=88), widening support tends on average to move medians upward, away from the upper edge of the original prior, but individual shifts remain heterogeneous.
- The selection correction does not consistently reduce the magnitude of prior sensitivity.

These results directly quantify prior sensitivity in posterior medians. They do not isolate the prior's causal effect on all metrics or establish a best prior; the posterior interval widths and full posterior shapes should also be considered.
