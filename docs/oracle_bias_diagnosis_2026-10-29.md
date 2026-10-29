# Exact-density oracle bias diagnosis

**Scope:** synthetic low-redshift model; no finite parent sample or KDE in the oracle estimator. The numerical summaries below come from the local diagnostic artifacts included in the uploaded patch. A GitHub Actions rerun of the same fixed-seed diagnostic has been added; these summary numbers should be considered provisional until that workflow completes and matches them.

## Reported calibration summaries

A prior-drawn (H_0) run tests prior-predictive calibration averaged over the declared uniform prior. Fixed-truth runs instead condition on one value of (H_0). The two questions are distinct.

| Truth design | Replicates | Mean posterior-mean bias ± Monte Carlo SE | 68% coverage | 90% coverage | PIT/rank chi-square p-value |
|---|---:|---:|---:|---:|---:|
| (H_0) drawn from prior [50, 90] | 2,000 | −0.060 ± 0.066 | 0.675 | 0.900 | 0.451 |
| Fixed (H_0=70) | 2,000 | +0.574 ± 0.072 | 0.683 | 0.895 | (3.87\times10^{-5}) |
| Fixed (H_0=60) | 2,000 | +0.439 ± 0.058 | 0.686 | 0.897 | 0.00322 |
| Fixed (H_0=80) | 2,000 | +0.308 ± 0.073 | 0.688 | 0.918 | 0.000352 |

For a fine-grid fixed-(H_0=70) diagnostic (1,000 replicates, 1,601 redshift and 321 (H_0) points), reported bias was +0.655 ± 0.100, with 68%/90% coverage 0.691/0.901 and PIT p-value 0.00282. This is consistent with the coarse-grid bias within Monte Carlo uncertainty; it does not make the fixed-truth bias vanish.

## Interpretation

Bayesian prior-averaged calibration means that when the parameter is generated from the same prior used in inference, posterior expectations are calibrated on average over repeated prior-predictive datasets. It does **not** imply that the posterior mean is unbiased at every fixed true parameter value. With a finite interval prior and a skewed likelihood, conditional posterior-mean bias at a fixed truth can occur even when the model is coherent.

The fact that the generator and likelihood use the exact same continuous toy host density rules out the catalogue KDE as the source of this particular bias. The result is not a claim that every modeling component is correct: prior-averaged calibration alone cannot establish detector realism or all likelihood assumptions. The independent tests include a negative control that drops the (H_0)-dependent selection normalization and checks that the resulting miscalibration is detected.

The historical KDE-based oracle estimate of approximately +0.88 ± 0.35 km s(^{-1}) Mpc(^{-1}) was not precise enough to isolate an additional KDE contribution from this exact-density fixed-truth effect. A new diagnosis of the finite-parent/KDE contribution needs a paired exact-density versus oracle-KDE comparison with its own Monte Carlo interval.

## Provenance and limits

The local summary JSONs had SHA-256 checksums verified in the supplied patch; the local fine-grid result uses a distinct sample count and should not be interpreted as exact agreement from identical datasets. The current action reruns 2,000 replicates for prior-drawn and fixed truths with a declared seed. It will be the reproducible artifact of record once it finishes and the saved JSON is inspected.

This study is still a controlled toy likelihood with fixed distance-noise scale and threshold selection. It is not a real gravitational-wave likelihood or a cosmological measurement.
