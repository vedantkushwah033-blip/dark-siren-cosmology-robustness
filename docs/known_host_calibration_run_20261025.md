# Known-host calibration run: audited results

## Run identity

- GitHub Actions run: [`40167238554`](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40167238554)
- Workflow commit: `858ea0588de32615de0285cb41504db80fabf1be`
- Pinned simulator commit: `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`
- Base seed: `20261025`
- Replicates: 300 per truth/noise cell; 36 cells; 10,800 raw records.
- Prior: H₀ ∈ [50, 90] km s⁻¹ Mpc⁻¹.
- Raw calibration audit: **PASS**. The auditor independently reconstructed posterior summaries from raw observations, checked all 2,700 truth/replicate pairing groups, cross-checked the summary files, and recorded SHA-256 hashes.
- Paired-comparison audit: **PASS**. A separate audit script independently reconstructed all 135 paired estimates and 4,000-resample percentile-bootstrap intervals per contrast, checked raw-input and output hashes, and wrote `paired_comparison_audit.json`.
- Artifact: [known-host-calibration-surface](https://api.github.com/repos/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/artifacts/12588868361). GitHub reports expiry on 2027-01-23.

## Main results

The table below gives unweighted means across the seven interior truths (H₀ = 55–85), with each truth cell contributing equally. These are descriptive summaries, not pooled-event estimates.

| Distance error | Mean 68% coverage | Mean 95% coverage | Mean absolute cell median error | Mean cell RMSE of median | Mean 68% interval width | Mean 95% interval width |
|---:|---:|---:|---:|---:|---:|---:|
| 5% | 0.695 | 0.954 | 0.589 | 3.301 | 6.638 | 12.816 |
| 10% | 0.729 | 0.964 | 1.805 | 5.854 | 12.034 | 22.388 |
| 15% | 0.756 | 0.974 | 3.253 | 7.453 | 15.944 | 28.284 |
| 20% | 0.770 | 0.983 | 4.434 | 8.329 | 18.671 | 31.728 |

At the central truth H₀ = 70, the 68% coverage estimates are 0.703, 0.707, 0.787, and 0.860 for 5%, 10%, 15%, and 20% distance error, respectively. The corresponding 95% coverage estimates are 0.960, 0.973, 0.997, and 1.000. Higher noise broadens intervals and increases coverage in this cell; that is not evidence of uniformly correct calibration across the grid.

The two boundary truths (H₀ = 50 and 90) have zero empirical coverage for both central credible intervals at all four noise levels. This is the expected structural effect of trying to include an exact prior boundary in equal-tailed intervals from a continuous posterior supported on [50, 90]. These cells must not be mixed into interior-calibration claims.

## Paired comparison example

At H₀ = 70, comparing 20% distance error with the 5% reference, the paired mean absolute-median-error difference is **+3.167** (bootstrap 95% CI **[+2.972, +3.373]**). Mean interval-width differences are **+12.628** for the 68% interval (95% CI **[+12.342, +12.913]**) and **+19.177** for the 95% interval (95% CI **[+18.752, +19.597]**). The 68% truth-inclusion rate increases by 0.157 and the 95% rate by 0.040 in this cell. Contrasts are paired by replicate ID within each truth; they are descriptive and are not multiplicity-adjusted confirmatory tests.

## Scope and limitations

This is a **synthetic, single-known-host, matched low-redshift distance-model diagnostic**. It is not an observational gravitational-wave result, a measurement of the real Hubble constant, or a validation of catalogue incompleteness corrections, host-galaxy weighting under model mismatch, detector selection, or realistic LVK inference. The artifact contains raw per-run records, summary tables, the independent audit report, the cell-level report, and paired-comparison outputs. Reproduce claims from those files rather than from this compact summary alone.
