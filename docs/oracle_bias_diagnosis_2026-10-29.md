# Exact-density oracle bias diagnosis

**Scope:** synthetic low-redshift model; no finite parent sample or KDE in the exact-density oracle estimator. This note distinguishes a 2,000-replicate local diagnostic snapshot (supplied with the patch) from a later, lower-resolution 300-replicate GitHub Actions screening rerun. They use different replicate budgets and grid settings, so their numbers are not expected to match exactly. Neither is a real-data result.

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


## CI-executed exact-density screening rerun (2026-10-29)

The dedicated workflow [#40797592794](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40797592794) completed successfully at code commit `5bbe5a434878af396e8c58e3ddfe7af0356e29ed`. Artifact `12865285154` contains `oracle_bias_exact_density.json`; the ZIP digest is `sha256:bf2d5f99faf1fe5cf8e6a028202463ee4546d51ec95dcf5c7681fbcf18d9a884`. It contains 300 replicate summaries per truth regime on a coarser 401-point redshift / 81-point H0 grid. The workflow saves aggregate summaries only; raw per-replicate data and seeds are absent, so this is a CI-executed screening result, **not** a row-level independently replayed calibration artifact.

| Truth design | Replicates | Mean posterior-mean bias ± MC SE | RMSE | 68% coverage | 90% coverage | Ten-bin PIT chi-square p-value |
|---|---:|---:|---:|---:|---:|---:|
| H0 drawn from prior [50, 90] | 300 | −0.221 ± 0.168 | 2.909 | 70.0% | 91.0% | 0.105 |
| Fixed H0=60 | 300 | +0.562 ± 0.146 | 2.584 | 68.3% | 91.0% | 0.276 |
| Fixed H0=70 | 300 | +0.724 ± 0.181 | 3.219 | 69.7% | 90.7% | 0.185 |
| Fixed H0=80 | 300 | +0.543 ± 0.189 | 3.306 | 69.7% | 91.3% | 0.324 |

The prior-drawn result is broadly compatible with calibration on this screening run. Fixed-truth posterior-mean bias is positive at all three tested values but is small relative to RMSE. The fixed-truth PIT chi-square p-values must not be interpreted using the prior-drawn SBC uniformity criterion: ranks at fixed truths do not have the same prior-predictive uniformity guarantee. This new run supports the distinction between prior-averaged calibration and conditional estimator bias, but the aggregate-only artifact is not sufficient for independent seed-level auditing. The 2,000-replicate local table is a different dataset and was not reproduced exactly by this coarser screen.


## Latest independently audited exact-density artifact

**Run:** [#40802415000](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40802415000), source commit `6d5f999eacab3dce3dc9eac5e444af22c25dc3c2`; the exact-density oracle job completed successfully. **Artifact:** ID `12867070639`, ZIP digest `sha256:4853794c2840980a8efc8bd6d82e1b82fc492aea3af4b291818db24464bd4134`.

The raw CSV contains 1,200 rows across four truth regimes (300 each), 1,129,316 bytes, SHA-256 `a38c64a2b9be9afce46c8f7b88659fa9f34f08a146d68376142b2184ba14cbee`. The summary JSON has 3,881 bytes, SHA-256 `a26c3b62d2612f04f34d9c2fb416ca1ae6024a69845ef7622666f1fe2efe3d72`. The independent report is **PASS**: all 1,200 seeds and accepted-distance vectors were replayed, all posteriors independently reconstructed, summary aggregates rederived, and file hashes/byte counts verified. Maximum errors were ≤1.85×10⁻¹³ for posterior summaries and ≤1.85×10⁻¹⁴ for ranks; true parameters and detected distances replayed exactly.

| Truth design | Replicates | Mean posterior-mean bias ± MC SE | RMSE | 68% coverage | 90% coverage | Mean posterior-CDF rank |
|---|---:|---:|---:|---:|---:|---:|
| H₀ drawn from prior [50,90] | 300 | −0.221 ± 0.168 | 2.909 | 70.0% | 91.0% | 0.524 |
| Fixed H₀=60 | 300 | +0.375 ± 0.151 | 2.639 | 72.0% | 89.3% | 0.480 |
| Fixed H₀=70 | 300 | +0.154 ± 0.189 | 3.271 | 68.7% | 88.0% | 0.505 |
| Fixed H₀=80 | 300 | +0.305 ± 0.191 | 3.314 | 67.0% | 92.7% | 0.474 |

The prior-drawn rank mean and ten-bin p-value (0.105) are broadly compatible with prior-predictive calibration at this moderate Monte Carlo resolution. Fixed-truth rows describe conditional estimator behavior and should not be judged using the same rank-uniformity criterion as prior-drawn SBC. This continuous-target control isolates finite-parent/KDE error from the remaining behavior of this toy likelihood; it is not a realistic gravitational-wave analysis.
