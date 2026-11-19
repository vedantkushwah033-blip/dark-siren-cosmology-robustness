# Full-parent prior-predictive SBC: 1,000-replicate extension

**Status:** COMPLETE; independent raw-seed reconstruction PASS. This is a synthetic full-parent baseline only.

## Design
- 1,000 independent prior-predictive replicates.
- Draw each true H0 independently from Uniform(20, 300) km s^-1 Mpc^-1.
- Generate a fresh parent population of 1,200 galaxies and 30 detected events per replicate.
- Use the complete finite parent catalogue and the matching selected-event likelihood with its detection-efficiency normalization.
- Evaluate on 1,121 equally spaced H0 points.
- Use independent deterministic seed streams for truth draws and population/event generation.
- Independently regenerate all 1,000 prior draws, parent populations, event sequences, posterior ranks, quantiles, and summary statistics.

## Verified workflow and artifact
- GitHub Actions run: [44440102723](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44440102723) — SUCCESS.
- Artifact ID: `14493796508`, name `full-parent-h0-sbc-1000`.
- Artifact ZIP SHA-256: `60a5ac8e43b6e6863fc8bc0f88bc921137593fee14c39833269b48eef84f7cc6`.
- Independent audit: PASS; 1,000/1,000 truths, parent populations, detected-event sequences, posterior-rank and quantile rows regenerated; summary metrics independently recomputed.
- The workflow also checked the complete replicate ID range, finite ranks in [0,1], summary/coverage agreement, and the 1,000-row denominator.

## Results
| Diagnostic | Result |
|---|---:|
| Replicates | 1,000 |
| Mean posterior rank | 0.49731 |
| Rank SD | 0.28490 |
| 10-bin rank chi-square | 4.14 |
| Rank chi-square p-value | 0.90196 |
| Uniform KS statistic | 0.02081 |
| Uniform KS p-value | 0.77127 |
| Mean posterior-median bias (km s^-1 Mpc^-1) | +0.19654 |
| MCSE of mean median bias | 0.46881 |
| Posterior-median RMSE (km s^-1 Mpc^-1) | 14.81886 |
| 50% interval coverage | 51.0% (MCSE 1.58 pp) |
| 68% interval coverage | 69.9% (MCSE 1.45 pp) |
| 80% interval coverage | 81.1% (MCSE 1.24 pp) |
| 90% interval coverage | 90.0% (MCSE 0.95 pp) |
| 95% interval coverage | 94.4% (MCSE 0.73 pp) |

## Interpretation
The rank diagnostics and empirical interval coverages are broadly compatible with the nominal prior-predictive calibration expectations for this exact synthetic full-parent model. The mean posterior-median bias is small relative to its Monte Carlo standard error. These are finite-simulation diagnostics, not proof of exact calibration.

**Numerical-resolution caveat:** the 1,000-replicate run used 1,121 H0 grid points. A separate 200-replicate 561→1,121 grid audit found a small number of large tail/outlier discrepancies and failed its preregistered maximum-shift gates. Therefore, this result should not be used to declare 561 points adequate. A grid-refinement check for the 1,121-point 1,000-replicate run remains outstanding.

## Scope guardrails
- This tests only the full-parent baseline likelihood.
- It does not test catalogue-inclusion correction, finite-validation estimation, or the project's five-arm comparison.
- P-values are descriptive diagnostics, not binary proof of calibration.
- It remains a synthetic low-redshift Euclidean model, not LVK/GWTC data, real-survey validation, or an observational H0 measurement.
