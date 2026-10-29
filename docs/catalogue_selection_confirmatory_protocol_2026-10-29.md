# Confirmatory replication protocol: catalogue-selection stress conditions

**Protocol status:** prespecified follow-up to the 29 October 2026 exploratory 50-replicate stress screen. This is a larger Monte Carlo replication using the same toy model, not external validation or a new physical model.

## Question and primary estimand

Across three fixed catalogue-inclusion regimes, does inverse-inclusion weighting under the exactly known synthetic selection law reduce replicate-level absolute posterior-mean error relative to ignoring catalogue incompleteness? For each regime, the primary estimand is

`mean(|posterior_mean_catalogue_only - H0_true| - |posterior_mean_correct_weighting - H0_true|)`

in km s^-1 Mpc^-1. Positive values favor correct weighting under the specified synthetic model. Each replicate shares the same generated parent population, catalogue, and GW-like distances across arms; comparisons are paired by replicate.

## Conditions and fixed settings

All conditions use true H0 = 70 km s^-1 Mpc^-1; 4,000 parent galaxies per replicate; 30 detected toy events; 81 H0 grid points on [50, 90]; 401 redshift grid points on [0.01, 0.40]; distance-noise scale 140 Mpc; and detection threshold 250 Mpc. The six analysis arms remain analytic truth (diagnostic only), finite-parent oracle KDE, catalogue-only, correct weighting, shifted transition, and shallow slope.

| Condition | True inclusion floor | Midpoint | Width | Seed start | Replicates |
|---|---:|---:|---:|---:|---:|
| Low floor | 0.05 | 0.22 | 0.05 | 20271000 | 250 |
| High floor | 0.25 | 0.22 | 0.05 | 20271100 | 250 |
| Broad transition | 0.10 | 0.22 | 0.08 | 20271200 | 250 |

Seed ranges are disjoint from the pilot and the first stress screen. The three conditions will be run independently, with raw replicate rows, summary, configuration manifest, checksums, and the independent posterior-grid/aggregate auditor. No condition may be dropped based on its outcome.

## Analysis and reporting

For each condition, report the primary paired mean contrast and the runner's 10,000-resample percentile bootstrap 95% interval; arm-level bias, RMSE, nominal 68% coverage with Wilson interval, mean interval width, and weight effective sample size; and sensitivity-arm comparisons. Preserve all six arms. Do not pool conditions as if they were exchangeable. Report the 50-replicate screen and this larger replication separately.

The confirmatory interval quantifies Monte Carlo variation conditional on the simulator and analysis choices; it does not include model-form uncertainty. Because the generator and inference likelihood are deliberately simplified, a successful replication only strengthens internal evidence for the toy model. It does not establish validity for real GW catalogues, realistic detector selection, host-population systematics, or a real H0 measurement.

## Audit and stop rules

A condition counts as completed only if its workflow succeeds, the independent auditor returns PASS, artifact upload succeeds, and raw rows/summary/manifest are preserved. A failure must be diagnosed and rerun under a documented new seed/configuration; do not selectively exclude replicates. If a result disagrees materially with the exploratory screen, report the discrepancy rather than tuning the model post hoc.
