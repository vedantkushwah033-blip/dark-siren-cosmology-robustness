# Host-mixture and detection-threshold sensitivity experiment protocol

**Status: implemented and completed; independent raw-to-summary audit passed on 2026-10-24.** See [the results report](host_mixture_threshold_sensitivity_results_2026-10-24.md), [factorial workflow run](https://github.com/vedantkushwaha033-blip/dark-siren-cosmology-robustness/actions/runs/40101614384), and [full test-suite run](https://github.com/vedantkushwaha033-blip/dark-siren-cosmology-robustness/actions/runs/40101614170). The report describes a synthetic diagnostic only, not an observational H0 measurement. This document prespecifies the next experiment to avoid choosing scenarios after seeing their outcomes. It is a synthetic statistical study, not a real gravitational-wave or galaxy-survey analysis.

## Motivation

The matched host-structure control found that the relative absolute-error performance of a two-host mixture and a collapsed effective-redshift model changes with the fixed true value of H0. Its detection interval and host mixture were held fixed. The next experiment should test whether this pattern is robust to host separation, intrinsic weights, and detection threshold rather than treating one configuration as universal.

## Factorial design

Use a fully crossed design with these factors:

| Factor | Levels |
|---|---|
| Fixed true H0 (km s^-1 Mpc^-1) | 60, 70, 80, 88 |
| Host redshifts | (0.02, 0.04), (0.02, 0.03), (0.02, 0.05) |
| Intrinsic host weights | (0.5, 0.5), (0.7, 0.3), (0.9, 0.1) |
| Upper detection threshold (Mpc) | 90, 120, 150 |
| Repetitions per cell | 300 detected observations |

This gives 4 × 3 × 3 × 3 = 108 design cells and 32,400 detected observations. Each generated observation is analysed by both the correctly specified two-host mixture and the collapsed single-host approximation at z_eff = sum_i w_i z_i. Use identical observations for both model arms. For each arm, report results with and without the corresponding toy selection-normalization term. The lower detection bound remains 0 Mpc, and accepted observations satisfy 0 < d < threshold.

## Controls held fixed

- Gaussian distance noise standard deviation: 20 Mpc.
- Uniform inference prior: H0 in [50, 90] km s^-1 Mpc^-1.
- Posterior evaluation grid: 4,001 points on [25, 125].
- Fixed seed recorded in the manifest; derive a distinct, deterministic RNG stream per design cell so changing cell ordering does not change other cells' results.
- At least 5,000 paired bootstrap resamples per cell for the primary paired contrast. Record bootstrap seed and quantile convention.
- Use the same detection threshold in the data-generation and analysis selection-normalization arms for the primary experiment. This isolates sensitivity to threshold value; it does not yet test threshold misspecification.

## Primary estimand

For each design cell and each correction status, define

**Delta = |error_collapsed| - |error_two-host|.**

Positive Delta favors the two-host model on absolute error; negative Delta favors the collapsed model. Report the mean paired Delta and percentile-bootstrap 95% interval. Do not infer a universal winner from the grand mean alone. Also report, by method and cell, bias, MAE, RMSE, 68% and 95% interval coverage with Wilson intervals, and mean interval widths.

Secondary outcomes: generated-host proportions among detected observations; analytic expected detected host proportions computed from the specified Gaussian model and truncation; and differences between empirical and expected proportions with Monte Carlo uncertainty. Intrinsic weights must never be mislabeled as detected-sample weights.

## Quality checks and acceptance criteria

1. Confirm exactly 32,400 rows in the primary per-run table and 108 unique design cells, each with 300 observations.
2. Verify every observed distance is finite and strictly inside its cell's detection interval.
3. Verify posterior medians, interval endpoints, errors, and all summary statistics are finite; posterior summaries must respect the numerical grid and stated prior support.
4. Recompute bias, MAE, RMSE, interval coverage, widths, and paired contrasts from the raw rows independently of the summary-writing code.
5. Check that paired model comparisons use the same observation and replicate identifier.
6. Compare detected-host fractions with analytic truncation probabilities. Large discrepancies must fail validation or be explicitly explained before interpretation.
7. Run focused tests and the full repository test suite; archive raw CSVs, summaries, manifest, code commit, workflow run, and artifact digest.
8. Report numerical-grid sensitivity on at least a small prespecified subset (e.g. grid sizes 2,001 versus 4,001) before attributing differences to model structure.

## Interpretation limits

- This is a synthetic stress test with a Gaussian distance proxy and hard truncation, not a detector-level likelihood or real catalogue.
- The collapsed redshift is a weighted mean, not a generally equivalent replacement for a mixture likelihood.
- Selection normalization and catalogue incompleteness are distinct mechanisms. This design changes the detection threshold but does not simulate missing galaxies or a realistic detector selection function.
- The fixed truths and equal cell sizes define a design average, not a posterior or population distribution over H0.
- A confidence interval excluding zero in a cell supports a paired difference under this simulator only; it does not establish astrophysical significance.
- Results must not be written into the README as findings until the implementation, tests, workflow, artifact, and raw-to-summary audit all pass.
