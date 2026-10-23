# Paired-contrast interpretation: 300-replicate factorial experiment

## Source and status

This interpretation was calculated from `paired_contrasts.csv` in artifact `paired-factorial-main-300-replicates`, produced by successful workflow run:

https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39959387009

The artifact contains 49 paired contrast-metric rows. Each contrast uses 300 paired synthetic universes. Bootstrap intervals are percentile intervals from 5,000 resamples. They quantify Monte Carlo uncertainty conditional on this simulation design, not uncertainty about model realism.

## Selected posterior-error contrasts

All differences are scenario B minus scenario A, in km s$^{-1}$ Mpc$^{-1}$. The metrics `median_error` and `mean_error` are signed errors relative to the fixed simulated truth, not absolute errors.

| Prespecified comparison | Metric | Mean paired difference | Bootstrap 95% interval |
|---|---|---:|---:|
| Random 70% minus complete; uniform hosts/inference | Posterior-median error | +0.057 | [+0.017, +0.097] |
| Random 70% minus complete; uniform hosts/inference | Posterior-mean error | +0.047 | [+0.012, +0.081] |
| Redshift-dependent 70% minus complete; uniform hosts/inference | Posterior-median error | −0.622 | [−0.660, −0.584] |
| Redshift-dependent 70% minus complete; uniform hosts/inference | Posterior-mean error | −0.494 | [−0.524, −0.464] |
| Redshift-dependent 70% minus random 70%; uniform hosts/inference | Posterior-median error | −0.679 | [−0.737, −0.623] |
| Redshift-dependent 70% minus random 70%; uniform hosts/inference | Posterior-mean error | −0.541 | [−0.587, −0.495] |
| Uniform inference minus matched luminosity inference; complete catalogue | Posterior-median error | +0.081 | [−0.002, +0.169] |
| Uniform inference minus matched luminosity inference; complete catalogue | Posterior-mean error | +0.034 | [−0.037, +0.102] |
| Uniform inference minus matched luminosity inference; random 70% catalogue | Posterior-median error | +0.063 | [−0.042, +0.171] |
| Uniform inference minus matched luminosity inference; redshift-dependent 70% catalogue | Posterior-median error | +0.147 | [+0.050, +0.243] |

## Host visibility and interval behavior

- Random 70% versus complete catalogue under uniform host generation: paired host-observation fraction changed by −0.293 (95% bootstrap interval [−0.347, −0.243]).
- Redshift-dependent 70% versus complete under uniform host generation: host-observation fraction changed by −0.290 (interval [−0.340, −0.240]).
- In the uniform-host block, redshift-dependent selection increased mean 68% interval width by +0.391 (interval [+0.349, +0.432]) and mean 95% interval width by +0.137 (interval [+0.117, +0.156]) relative to the complete catalogue.
- The observed difference in 68% coverage for random 70% versus complete was +0.0033, with a percentile bootstrap interval [0, 0.01]. This is a difference in empirical proportions, not evidence that random incompleteness improves calibration.
- The 95% coverage contrast is exactly zero in several comparisons because all or nearly all sampled intervals covered the fixed truth. A degenerate bootstrap interval for this observed paired difference does not establish equal true coverage probabilities.

## Interpretation guardrails

1. The random-incompleteness contrast in signed error is small in absolute magnitude even though its bootstrap interval excludes zero. Statistical separation from zero is not the same as scientific importance.
2. The redshift-dependent selection contrasts are larger in this simulator, but the selection rule is an artificial toy mechanism. Do not generalize the direction or magnitude to real galaxy catalogues.
3. There are 49 related contrast-metric intervals and no multiplicity adjustment. Treat them as prespecified exploratory diagnostics, not 49 independent confirmatory tests.
4. The posterior intervals are broad (mean 68% widths around 22 km s$^{-1}$ Mpc$^{-1}$) and empirical 68% coverage at fixed $H_0=70$ is about 94%–95%. Because the truth is fixed at the centre of the inference grid, this fixed-truth coverage result can be influenced by the grid/prior geometry and posterior behavior. The correct next step is a coverage audit over several fixed true values and an inspection of the posterior interval construction.
5. All results are from a synthetic, simplified, single-event model. They are not an LVK analysis or an estimate of the real Hubble constant.

## Recommended next analysis

- Re-run coverage diagnostics at multiple fixed true $H_0$ values for the same scenario definitions.
- Audit posterior quantile calculations and the relationship between posterior standard deviation and central credible-interval width.
- Report bias/RMSE and interval widths beside coverage; do not select a scenario based on coverage alone.
- Consider reporting the seven prespecified scientific comparison families separately from the seven metric-level diagnostics, with a clear distinction between exploratory and confirmatory analyses.
