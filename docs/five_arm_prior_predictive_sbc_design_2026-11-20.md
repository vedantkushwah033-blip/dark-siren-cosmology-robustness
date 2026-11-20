# Five-arm prior-predictive SBC: design contract and smoke status

## Purpose

The existing five-arm curvature experiments used a fixed generating value H0=70. Their empirical coverage estimates are fixed-truth repeated-sampling coverage, not prior-predictive simulation-based calibration (SBC).

This follow-up changes the truth-generation layer only: each replicate draws H0 independently from the same Uniform(30,110) prior used by every inference arm. The same drawn truth, parent population, detected events, validation frame, and validation IDs are shared across the three curvature settings within a replicate. Each inference arm uses the identical prior support.

## Frozen design for the first engineering smoke

- Truth prior: Uniform(30,110), sampled once per replicate.
- Independent truth seed stream: SeedSequence([master_seed, replicate, 2]).
- Parent/event stream: stream 0; validation stream: stream 1.
- Master seed: 2026112020.
- Five replicates; gamma = -1, 0, +1.
- N_parent=1200; 30 detected events.
- Validation frame=500; nested validation sizes 100, 300, 500.
- Five inference arms A–E as in the fixed-truth experiments.
- All arms use Uniform(30,110), H0 grid of 281 nodes for smoke only.
- Each child run is independently audited before inclusion in summaries.
- Aggregate artifact has a fail-closed ledger and an independent top-level audit.

The five-replicate run is an engineering smoke, not sufficient for coverage or rank-uniformity claims. It must pass before any larger prior-predictive run is launched.

## Required SBC outputs for a later confirmatory-scale run

1. Prior-predictive posterior rank for each arm and replicate.
2. Rank histogram with predeclared bins and a uniformity diagnostic, interpreted with Monte Carlo uncertainty rather than a single p-value.
3. Coverage at central 50%, 68%, 80%, 90%, and 95% levels.
4. Bias, RMSE, interval width, and paired arm contrasts.
5. Separate summaries by curvature and validation size.
6. Failures/non-evaluable runs retained in the planned denominator with conservative bounds.
7. Rank/quantile and posterior-grid convergence checks on all replicates.
8. Independent reconstruction of truth seeds, raw parent/event/validation inputs, likelihoods, posterior summaries, and rank statistics.

## Interpretation boundaries

- A prior-predictive rank diagnostic is not a guarantee of frequentist calibration at every fixed H0 value.
- The oracle-parent arm is a diagnostic reference, not an implementable catalogue estimator.
- No result is evidence about LVK/GWTC data, real galaxy surveys, or an observational H0 measurement.
- No estimator is presumed superior; the aim is to identify conditions under which estimated inclusion weighting helps or fails.
- No IPW plus additional missing-host mixture double correction is introduced.

## Status

Pending CI smoke and independent artifact audit. No calibration conclusion is claimed until a larger run, rank diagnostics, and grid-convergence audits are completed.

## Rank and coverage implementation

The five-arm child runner now records the posterior CDF evaluated at the generating H0 for every arm and validation size. Its independent child auditor reconstructs that rank directly from the archived likelihood curve, separately from the production posterior-normalization helper. The same audit reconstructs central quantiles at 2.5%, 5%, 10%, 16%, 25%, 75%, 84%, 90%, 95%, and 97.5%.

A separate summarizer writes ten-bin rank histograms and coverage summaries for central 50%, 68%, 80%, 90%, and 95% intervals, with Wilson intervals for coverage. Rank-test p-values are explicitly descriptive and are not treated as confirmatory when expected bin counts are small. The rank/coverage summaries use successful runs only while reporting planned and failed counts alongside them; the independent lifecycle auditor retains the complete denominator.

The smoke workflow now runs focused truth-stream, posterior-rank, and summary tests; executes the five-replicate prior-predictive design; independently audits child and aggregate artifacts; and emits rank/coverage tables. These small-run summaries are an implementation check only.
