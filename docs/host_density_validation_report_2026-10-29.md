# Host-density validation: local paired sensitivity campaign

**Status:** Exploratory local result snapshot from the uploaded three-commit patch. The patch's CSV and summary SHA-256 checksums were checked, and a two-replicate subset passed the new independent reconstruction audit. The full 13-configuration campaign has now been wired into GitHub Actions for a clean rerun; this report is **not yet the CI-generated replacement artifact**. Do not conflate these local results with the separately locked Stage B study.

## Baseline Scott estimator

The local run used 100 paired replicates, seed 20261029, 6,000 parent galaxies, 40 detected events, true (H_0=70) km s(^{-1}) Mpc(^{-1}), distance-noise scale 140 Mpc, and a 250 Mpc detection threshold. The catalogue inclusion law had floor 0.1, midpoint 0.22, and width 0.05. The estimator used the corrected weighted-Scott bandwidth convention and the same fixed-noise likelihood as the host-density benchmark.

| Arm | Bias (km s(^{-1}) Mpc(^{-1})) | RMSE | Nominal 68% interval coverage |
|---|---:|---:|---:|
| Complete-parent oracle KDE | +0.88 | 3.64 | 0.61 |
| Incomplete catalogue without correction | −8.97 | 9.46 | 0.06 |
| Known inclusion probabilities | +0.93 | 3.66 | 0.59 |

The catalogue-only arm's mean absolute posterior-mean error exceeded oracle by 6.07 km s(^{-1}) Mpc(^{-1}), paired 95% interval [5.11, 6.90]. Known-completeness weighting differed from oracle by +0.01 [−0.09, +0.10] in that same absolute-error contrast. The severe catalogue-only degradation is specific to this constructed example.

## Bandwidth and grid sensitivity

The tested physical factors mean (h=f\,\times\,SD_w(\text{original sample})); the default Scott rule uses (h=SD_w n_{\rm eff}^{-1/5}). Factors 0.05, 0.1, and 0.2 gave results broadly similar to Scott's default in the local snapshot. Factors 0.5, 0.8, and 1.2 increasingly oversmoothed the density. At factor 1.2, catalogue-only bias was −1.93, RMSE 4.13, and 68% coverage 0.62. That does **not** establish a completeness correction: flattening the inferred host density can accidentally change posterior errors while leaving the selection problem unresolved.

Changing the posterior/redshift grids through the locally tested coarse, default, and fine settings did not change the reported summaries at the displayed precision. The current GitHub workflow will regenerate these configurations and produce the underlying tables/manifest in one downloadable artifact.

## Selection-law misspecification

The estimator's assumed inclusion law was deliberately changed while the data generator retained the true midpoint 0.22 and width 0.05.

| Assumed inclusion law | Bias | RMSE | 68% coverage | Paired absolute-error difference vs oracle (95% CI) |
|---|---:|---:|---:|---:|
| Midpoint 0.18 | little degradation in this sample | — | — | — |
| Midpoint 0.26 | little degradation in this sample | — | — | — |
| Width 0.10 (shallower transition) | −2.86 | 4.35 | 0.49 | +0.79 [ +0.17, +1.40 ] |
| Width 0.025 (steeper transition) | +2.31 | 4.28 | 0.51 | +0.43 [ +0.18, +0.68 ] |

These are finite-sample results from a single synthetic population and noise model. The midpoint variants produced little degradation in this sample; that is not evidence that midpoint misspecification is generally harmless. The two width variants clearly changed coverage, but the paired contrasts compare to oracle and should not be read as a universal penalty under all settings.

The inverse-inclusion weight diagnostics in the local data showed maximum-to-mean weight ratios around 2.8–5.0 and minimum effective-sample-size fractions around 0.46–0.88. These diagnostics are useful here because the inclusion-probability floor bounds extreme weights. They do not establish stable inverse weighting for selection laws with tiny probabilities.

## Independent audit scope

The accompanying independent auditor:
- resimulates selected replicate seeds from the same declared generator;
- recomputes reflected weighted KDEs with a separate kernel-sum implementation;
- computes selected posterior summaries with Simpson quadrature on finer grids;
- checks the saved manifest hashes when present;
- compares posterior means within 0.05 and medians/interval bounds within 0.10 km s(^{-1}) Mpc(^{-1}).

Only a small subset of the 100-replicate Scott run was re-audited in the original local work; the other configurations were checksum-checked but were not all independently reconstructed. The new campaign workflow is intended to close the reproducibility gap by rerunning all configurations and separately auditing a fresh smoke sample.

## Decision

**GO:** exploratory synthetic-methods reporting with explicit uncertainty and model limits.

**NO-GO:** claim of validated real LVK performance, a universal correction for incompleteness, or a measurement of real (H_0). This benchmark omits realistic detector response, uncertain redshifts, peculiar velocities, event-rate marginalization, and survey-level hierarchical inference.
