# Independent reproduction: catalogue-selection confirmatory study

**Date:** 29 October 2026  
**Scope:** independent numerical audit of the three 250-replicate synthetic selection-function conditions. This report concerns controlled toy simulations only. It is not a real LVK analysis, a real H0 measurement, or external scientific validation.

## Integrity checks

The three saved condition artifacts were checked independently of their aggregate summaries. Across the three conditions, 4,500 saved posterior grids (1,500 per condition) each contained 81 finite, nonnegative density values and integrated to one within floating-point tolerance. Recorded artifact manifest checksums matched the saved files, and the independent reanalysis reproduced the stored summary values within the reported numerical precision.

The three condition jobs used overlapping seed ranges. This was disclosed in the manuscript, protocol and condition-results report. The within-condition paired contrasts remain interpretable; the three conditions must not be treated as statistically independent replications.

## Primary paired contrast

Define D as the catalogue-only absolute posterior-mean error minus the correct-weighting absolute posterior-mean error. Positive D favors correct weighting under the specified synthetic truth law. The following estimates use an independently chosen bootstrap seed (20261029); small differences from the original stored bootstrap intervals are expected because the resampling seed differs.

| True inclusion condition | Mean paired D (km s^-1 Mpc^-1) | Independent percentile-bootstrap 95% interval |
|---|---:|---:|
| Low inclusion floor (0.05) | 7.4358 | [6.8050, 8.0406] |
| High inclusion floor (0.25) | 2.2462 | [1.7217, 2.7476] |
| Broad transition (width 0.08) | 2.3537 | [1.7767, 2.9286] |

These results show that the size of the catalogue-only error penalty varies materially across the tested toy inclusion conditions. They do not show that weighting is universally robust: the correct inclusion law is known by construction, and real catalogues have additional uncertainties.

## Arm-level calibration context

| Condition | Arm | Bias | RMSE | Nominal 68% interval coverage |
|---|---|---:|---:|---:|
| Low floor | Catalogue only | -10.6166 | 11.0401 | 4.8% |
| Low floor | Correct weighting | +1.1850 | 4.0952 | 65.6% |
| High floor | Catalogue only | -4.7169 | 5.8978 | 38.0% |
| High floor | Correct weighting | +0.9377 | 4.0556 | 69.6% |
| Broad transition | Catalogue only | -4.9943 | 6.1527 | 33.6% |
| Broad transition | Correct weighting | +1.3042 | 4.1541 | 70.8% |

All numbers are in the simulation's H0 units (km s^-1 Mpc^-1). Coverage is based on 250 replicates per condition, so it remains uncertain; it should not be presented as precise evidence of nominal calibration.

## Bias-compensation follow-up

A separate four-cell follow-up compared correct weighting with a deliberately shallow assumed transition at true H0 = 60 and 80. Let Delta_RMSE = RMSE(correct) - RMSE(shallow); positive values favor the shallow arm on RMSE in that cell.

| Condition | True H0 | Delta_RMSE | Paired-bootstrap 95% interval | Correct bias | Shallow bias |
|---|---:|---:|---:|---:|---:|
| High floor | 60 | +0.2275 | [0.1329, 0.3165] | +0.985 | +0.193 |
| Broad transition | 60 | +0.3270 | [0.1520, 0.4912] | +1.066 | -0.391 |
| High floor | 80 | -0.1172 | [-0.2338, +0.0043] | +0.513 | -0.550 |
| Broad transition | 80 | -0.7301 | [-0.9056, -0.5471] | -0.226 | -2.181 |

At H0 = 60, the shallow arm has lower RMSE in both conditions, with intervals above zero. At H0 = 80, that advantage disappears in the high-floor condition (interval crosses zero) and reverses in the broad-transition condition, where correct weighting has lower RMSE. The results are consistent with truth-dependent bias compensation, not proof that it is the sole mechanism. The high-floor H0 = 80 cell is inconclusive and must not be described as evidence for either arm.

The variance/bias decomposition reinforces the need for this qualification. In the broad-transition H0 = 80 cell, squared bias is approximately 0.051 for correct weighting and 4.757 for the shallow arm; the shallow arm's RMSE is worse despite a variance estimate that is not higher. A lower RMSE in one truth cell therefore cannot be promoted into a general recommendation for a misspecified inclusion curve.

## Reproducibility and limitations

The independent audit was run against the saved workflow artifacts, and the reanalysis script/test used for the follow-up passed locally (1 test). The current PR-head CI result is tracked separately in the GitHub Actions run linked from the pull request. The analysis is synthetic, model-conditional and limited by its generator, finite replicate counts and shared seed ranges across confirmatory conditions. No CHIMERA v2.2 execution is claimed: the required version-compatible incomplete-catalogue inputs remain unavailable/unverified. No real-data validation, domain-expert sign-off, submission, peer review or publication is claimed.
