# Full 90-cell fixed-truth H0 pilot: independent artifact readout

Date: 2026-11-19

## Executive status

A complete 90-cell fixed-truth pilot artifact exists. The simulation step completed all 90 planned cells successfully. The original workflow was marked failed only because its aggregate auditor rejected a legitimate change in the size-specific `selected_for_validation` flags across nested validation sample sizes. A later replay exposed a second auditor compatibility defect: the historical summary has no explicit `not_run` column. Neither failure indicates a failed child simulation.

The archived raw artifact has now been independently checked outside that original workflow:

- 90/90 ledger rows are `SUCCESS`.
- 90/90 child audit reports are `PASS`.
- 90/90 coefficient-quadrature runs are `CONVERGED`.
- All 2,110 aggregate-manifest file hashes and paths verify.
- All 90 planned-cell identifiers match the ledger.
- 30 within-`(gamma, replicate)` nested-validation blocks pass.
- 30 cross-gamma paired-validation-design blocks pass.
- All 9 aggregate summary rows were independently recomputed from the raw ledger to floating-point precision.
- All 90 highest-order posterior densities integrate to 1 within approximately (7\times10^{-15}); independently recomputed quantiles agree with the ledger to below (9\times10^{-13}).

This is an integrity-checked exploratory fixed-truth pilot, not a confirmatory calibration study or SBC campaign.

## Provenance

- GitHub Actions run: `44402834830`
- Run commit: `1ab60e2e7b7d446eb15f9b59efb3003f25965087`
- Artifact ID: `14473894295`
- Artifact ZIP size: 37,660,291 bytes
- Artifact ZIP SHA-256: `e22ae91eff8e0d1aa02f4d03699543df096518cf7d21cbe751ecee9ea66c5426`
- Workflow simulation step: SUCCESS
- Original aggregate orchestration audit: FAILURE due to a false-positive comparison of size-specific selection flags.
- Corrected independent replay: subsequent compatibility failure because the auditor required a `not_run` summary field absent from this legacy complete-run summary. The auditor has since been patched to accept the legacy schema only when the ledger independently confirms zero NOT_RUN rows; CI replay is pending.

The 90-cell artifact is distinct from the earlier 70/90 partial artifact and the earlier 64-cell partial artifact. Do not combine their rows or denominators.

## Design

| Component | Setting |
|---|---:|
| Parent galaxies | 1,200 |
| Detected GW events | 4 per replicate |
| Curvature values | -1, 0, +1 |
| Validation frame | 500 |
| Nested SRS validation sizes | 100, 300, 500 |
| Replicates per cell | 10 |
| Planned / successful cells | 90 / 90 |
| Fixed truth | H0 = 70 km s⁻¹ Mpc⁻¹ |
| H0 prior | Uniform [20, 300] |
| H0 grid | 141 points; spacing 2 km s⁻¹ Mpc⁻¹ |
| Inclusion DGP | Quadratic logistic |
| Inclusion fit | Degree-one logistic model |
| Coefficient quadrature | Adaptive orders 3, 5, 7 (all converged) |

At γ=0 the degree-one logistic fit matches the inclusion model family; at γ=±1 the omitted quadratic term induces deliberate misspecification. The inference is conditional on four detected events and does not model event-count/rate calibration.

## Independent integrity and pairing audit

The audit reconstructed and checked:

1. Every aggregate-manifest path and SHA-256 digest: 2,110/2,110 passed.
2. Every child audit report and quadrature status: 90/90 passed.
3. Planned-cell / ledger identity: 90/90 matched.
4. Catalogue, parent, GW event observations, and underlying validation frame invariance across validation sizes within each gamma and replicate.
5. Nested validation SRS prefixes: 100 ⊂ 300 ⊂ 500 in all 30 gamma-by-replicate blocks.
6. Shared validation IDs and redshifts across gamma at fixed replicate and sample size: all 30 blocks passed. Retention outcomes are intentionally allowed to vary with gamma.
7. Aggregate summary recomputation from the replicate ledger: maximum absolute discrepancy was below (4\times10^{-15}) across reported numeric metrics.
8. Posterior density normalization and quantiles: all areas were between 0.9999999999999934 and 0.9999999999999973; maximum absolute quantile discrepancy was below (9\times10^{-13}).

The original auditor failure was a design-audit bug: it compared `selected_for_validation` as though it were invariant across sample sizes, although this flag must change as the SRS prefix grows. The corrected audit compares only invariant frame columns and separately verifies the selected IDs and nested-prefix property.

## Full 90-cell descriptive results

Bias is the mean replicate posterior-median error (posterior median minus 70). RMSE is computed across the 10 replicate median errors. Coverage denominators are 10 in every cell. ESS and maximum normalized weight are medians across the 10 replicates.

| γ | Validation n | Mean median bias | RMSE | 68% coverage | 90% coverage | 95% coverage | Mean 68% width | Median min Kish ESS | Median max normalized weight |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -1 | 100 | +3.951 | 10.289 | 0.70 | 0.90 | 1.00 | 24.781 | 14.493 | 0.2208 |
| -1 | 300 | +3.948 | 10.209 | 0.70 | 0.90 | 1.00 | 24.328 | 37.571 | 0.1266 |
| -1 | 500 | +3.461 | 10.231 | 0.70 | 0.90 | 1.00 | 24.210 | 63.834 | 0.0739 |
| 0 | 100 | +8.211 | 11.511 | 0.70 | 1.00 | 1.00 | 25.224 | 10.036 | 0.2550 |
| 0 | 300 | +8.218 | 11.719 | 0.80 | 0.90 | 1.00 | 24.693 | 20.433 | 0.1651 |
| 0 | 500 | +7.500 | 11.380 | 0.90 | 0.90 | 1.00 | 24.743 | 30.489 | 0.1275 |
| +1 | 100 | +13.290 | 16.401 | 0.50 | 0.90 | 0.90 | 25.481 | 9.672 | 0.2289 |
| +1 | 300 | +13.620 | 16.723 | 0.50 | 0.80 | 1.00 | 24.756 | 16.282 | 0.1517 |
| +1 | 500 | +12.860 | 16.066 | 0.50 | 0.80 | 1.00 | 24.972 | 24.608 | 0.1171 |

Bias, RMSE, and interval width are in km s⁻¹ Mpc⁻¹.

### Paired curvature contrast

At each validation size, the paired replicate-level difference in median bias (γ=+1 minus γ=-1) is:

| Validation n | Paired replicates | Mean paired difference | 95% paired t interval |
|---:|---:|---:|---:|
| 100 | 10 | +9.34 | [+5.10, +13.57] |
| 300 | 10 | +9.67 | [+5.95, +13.40] |
| 500 | 10 | +9.40 | [+5.81, +12.99] |

This is a consistent curvature-direction contrast in this particular synthetic design. With only 10 paired replicates, the intervals are exploratory and should not be treated as confirmatory evidence.

### Validation-size effect

Paired changes in median bias from n=100 to n=500 were:

- γ=-1: mean change -0.49; paired SD 1.05.
- γ=0: mean change -0.71; paired SD 1.90.
- γ=+1: mean change -0.43; paired SD 1.30.

By contrast, the median minimum Kish ESS rises and the maximum normalized weight falls as validation size increases in all three curvature settings. This indicates improved weight stability without a comparably large change in H0 median bias over this small pilot. It does not establish that validation size is unimportant.

The positive bias already present at γ=0 is important: not all bias can be attributed to curvature misspecification. Four detected events, finite catalogue sampling, distance noise, and the simplified selection/inference model can all contribute.

## H0 grid-refinement audit

The separate fixed-observation grid experiment completed all three child runs and their independent likelihood/artifact audits. Its aggregate comparison workflow initially failed because NumPy 2.5 removed `np.trapz`; the audit script and regression fixture have been updated to `np.trapezoid`. The rerun is pending.

An independent local reconstruction from the archived three-grid artifact found:

| H0 grid points | Grid spacing | Max absolute density difference vs 561 | Total variation vs 561 | Max absolute quantile difference vs 561 |
|---:|---:|---:|---:|---:|
| 141 | 2.0 | 2.23e-5 | 3.23e-4 | 0.0272 |
| 281 | 1.0 | 5.58e-6 | 6.47e-5 | 0.0055 |
| 561 | 0.5 | 0 | 0 | 0 |

The observable-input hashes match across all three grids, the quadrature order is 7 in each, and the model/generator configurations match after excluding only the H0 grid-point count. The posterior median is approximately 130.036 on the finest grid, far above the fixed truth 70 in this single realization.

This establishes numerical grid stability for one fixed synthetic realization only. It does not establish repeated-sampling calibration, coverage, SBC, or estimator superiority.

## Interpretation and limitations

The complete 90-cell pilot supports three cautious observations:

1. The curvature direction is associated with a substantial shift in H0 median bias in this synthetic experiment: γ=+1 is more positively biased than γ=-1 in paired replicates.
2. Increasing validation size improves IPW weight diagnostics (higher ESS, lower maximum normalized weight), but the average H0 median bias changes only modestly between n=100 and n=500.
3. There is also appreciable positive bias at γ=0, so the result cannot be framed as curvature misspecification being the sole cause.

Limitations remain decisive:

- Ten replicates per cell give coarse and uncertain coverage estimates.
- Four detected events per replicate make this a deliberately small conditional experiment.
- Fixed-truth coverage is not joint SBC; no prior-predictive parameter generation/PIT campaign is included.
- H0 grid stability was assessed for one realization, not across the repeated-sampling ensemble.
- The model is a low-redshift Euclidean synthetic toy, not LVK/GWTC or real galaxy-survey analysis.
- No observational H0 estimate, real-data validation, or general estimator superiority is claimed.

## Next gates

1. Obtain a passing CI replay of the corrected independent orchestration auditor against the archived 90-cell artifact.
2. Obtain a passing CI run of the corrected 141/281/561 grid-refinement auditor.
3. Run a larger pre-specified fixed-truth ensemble with adequate Monte Carlo precision and explicit failure accounting.
4. Design and execute joint SBC/PIT separately from fixed-truth coverage.
5. Only then decide whether the H0 propagation study is a meaningful central paper contribution or a supporting synthetic sensitivity analysis.
