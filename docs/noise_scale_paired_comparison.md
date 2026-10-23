# Paired comparison of noise-scale likelihoods

**Input artifact:** [Noise-scale consistency audit, workflow run 39974976415](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39974976415)  
**Analysis:** paired replicate bootstrap, 10,000 resamples per injected (H_0), fixed seed 20261023.  
**Sign convention:** every delta is *H0-dependent uncertainty minus fixed reported uncertainty*. Negative values mean the H0-dependent mode is lower for that metric.

## Selected paired results

All errors and widths are in km/s/Mpc; coverage differences are percentage points. The intervals are percentile 95% bootstrap intervals based on paired replicate-level differences.

| True H0 | Mean signed error delta (95% CI) | RMSE delta (95% CI) | 68% coverage delta (95% CI) | 68% width delta (95% CI) |
|---:|---:|---:|---:|---:|
| 55 | −1.716 [−1.781, −1.650] | −1.420 [−1.550, −1.280] | +2.2 pp [+0.6, +3.9] | −3.674 [−3.902, −3.445] |
| 60 | −1.689 [−1.752, −1.627] | −1.142 [−1.287, −0.982] | +2.3 pp [+0.7, +3.8] | −3.460 [−3.666, −3.256] |
| 65 | −1.568 [−1.625, −1.511] | −0.754 [−0.882, −0.622] | +1.0 pp [−0.6, +2.6] | −2.957 [−3.128, −2.784] |
| 70 | −1.456 [−1.511, −1.399] | −0.514 [−0.633, −0.390] | +0.8 pp [−0.8, +2.4] | −2.576 [−2.722, −2.425] |
| 75 | −1.204 [−1.260, −1.147] | −0.108 [−0.206, −0.008] | +0.7 pp [−0.9, +2.3] | −1.818 [−1.934, −1.699] |
| 80 | −0.946 [−1.003, −0.888] | +0.167 [+0.091, +0.240] | +0.2 pp [−1.4, +1.7] | −1.079 [−1.172, −0.984] |
| 85 | −0.694 [−0.747, −0.640] | +0.349 [+0.303, +0.395] | −0.5 pp [−1.8, +0.8] | −0.249 [−0.323, −0.176] |

## Interpretation

- The H0-dependent specification reduces mean signed error at all seven truths, but the mean error remains positive in every cell. This is a partial improvement, not removal of bias.
- Paired RMSE decreases at truths 55–75, though only slightly at 75, and increases at 80 and 85. The effect is not uniform across the parameter range.
- The 68% coverage difference is small and its paired interval includes zero at truths 65–85. At 55 and 60 the observed increase is more clearly positive under this bootstrap comparison. Coverage under the new mode is nevertheless not a universal proof of calibration.
- The H0-dependent mode narrows 68% intervals at every tested truth. This can be expected from its different likelihood shape and should be considered together with coverage and error, not in isolation.

## Method and caveats

Each replicate's same observed distance is used in both likelihood modes; bootstrap resampling is paired by replicate. The bootstrap intervals quantify Monte Carlo uncertainty conditional on this toy design. They are not a substitute for model checking, and they do not establish validity for the unknown-host galaxy-mixture likelihood or real gravitational-wave observations.

The paired analysis is now automated by `scripts/analyze_noise_scale_paired_differences.py` and is wired into the noise-scale audit workflow so that future runs will publish `paired_bootstrap.csv` alongside the raw and summary outputs.
