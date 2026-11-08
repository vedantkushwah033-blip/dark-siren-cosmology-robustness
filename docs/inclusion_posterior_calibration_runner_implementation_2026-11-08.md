# Multi-replicate inclusion-posterior calibration runner

**Status:** implementation added; no experiment has been executed or validated by CI at this commit.
**Scope:** synthetic low-redshift Euclidean benchmark only.

## Entry point

`scripts/run_inclusion_posterior_calibration.py` orchestrates the existing
single-replicate H0 runner and its independent artifact auditor. It keeps the
calibration experiments separate:

- `H0_PIT_FIXED_INCLUSION_TRUTH`: draw each replicate's generating H0 from the
  same declared uniform H0 prior used by inference, while holding the catalogue
  inclusion-truth parameters fixed, then calculate the posterior CDF at that
  truth. This is an H0 prior-predictive PIT conditional on fixed inclusion truth,
  **not joint SBC**, because inclusion coefficients are not drawn from their
  inference prior.
- `FIXED_TRUTH`: hold generating H0 at the declared fixed value (default 70)
  and calculate repeated-sampling equal-tail credible-interval coverage.

Example pilot command (20 replicates per experiment × floor × validation size):

```bash
python scripts/run_inclusion_posterior_calibration.py \
  --replicates 20 \
  --experiments H0_PIT_FIXED_INCLUSION_TRUTH,FIXED_TRUTH \
  --floors 0.0,0.10,0.25 \
  --validation-sizes 100,300,1000 \
  --n-validation-frame 4000 \
  --n-parent 2400 --n-detected 8 \
  --orders 7,11,15,19 --h0-points 81 \
  --output-dir artifacts/inclusion_posterior_calibration
```

This default design plans 360 independently audited replicate-size runs. It is
intentionally a long-running pilot; increase replication only after reviewing
runtime, numerical evaluability and the audit trail. The output directory must
be new or empty; choose a fresh path for every run so prior artifacts are never
overwritten.

## Pairing and seeds

For each experiment/floor/replicate, the runner derives deterministic,
separate analysis, validation and H0-prior-truth seed streams using
`numpy.random.SeedSequence`. The same analysis and validation seeds are reused
across the requested validation sizes. The single-replicate generator creates
one independent validation frame and selects prefixes of one random
permutation, so smaller validation samples are strict subsets of larger ones.
The analysis parent, observed catalogue and GW events are regenerated from
the same analysis seed and configuration at each size, so they are held fixed
within the paired validation-size contrast.

The validation frame remains independent of the analysis population. It is
not a disjoint subsample of the analysis parent.

## Status semantics

Every attempted cell is retained in `replicate_ledger.csv`:

- `SUCCESS`: smoke runner completed, quadrature converged under its declared
  rule, and the independent artifact audit passed.
- `NON_EVALUABLE`: artifacts passed audit but quadrature convergence was not
  achieved; such a posterior is excluded from primary PIT/coverage summaries.
- `FAILED`: runner or independent audit failed. Error stage and stderr tail
  are retained.
- `NOT_RUN`: reserved in the ledger vocabulary for planned work not attempted;
  the normal runner continues after per-cell failures rather than dropping
  failed attempts.

Each replicate directory contains the single-replicate observable/truth
artifacts, quadrature outputs, checksums and audit report. Runner/auditor
stdout and stderr are retained beside those directories.

## Outcomes

Conditional H0 prior-predictive PIT summaries include the empirical PIT mean, variance, fixed-bin counts,
expected equal-bin count and Monte Carlo standard error. They also report
prior-predictive inclusion rates for 68%, 90% and 95% equal-tail intervals,
with Wilson uncertainty intervals. These are descriptive diagnostics, not a
standalone proof of calibration.

Fixed-truth summaries include 68%, 90% and 95% equal-tail coverage, Wilson
95% intervals, binomial Monte Carlo standard errors, median interval widths,
posterior-median bias, and median coefficient-quadrature/IPW stability
diagnostics (Kish ESS and maximum normalized host weight). The 95% interval uses posterior q2.5% and q97.5%, now archived by the
single-replicate smoke runner.

All summary denominators use successful, numerically evaluable, independently
audited replicates; attempted, failed and non-evaluable counts are reported
alongside them. No failed replicate is silently removed from the ledger.

## Scientific limits

This runner does not establish calibration or coverage until it is executed
with adequate replication and independently reviewed. A 20-replicate pilot
is for orchestration and failure handling, not confirmatory inference.
Positive catalogue floors are controlled misspecification stress cases for
the fitted degree-one logistic model; they are not a universal proxy for
real-survey incompleteness. The benchmark is not LVK/GWTC analysis, real
galaxy-survey validation, detector selection validation or an observational
H0 measurement.


## Failure-inclusive bounds

For each binary interval-inclusion/coverage indicator, the summary now reports
an additional all-attempt partial-identification interval. Successful,
independently audited, numerically evaluable replicates contribute their
observed 0/1 outcome. Failed and non-evaluable attempts are treated as unknown:
the lower endpoint assigns all unknown outcomes 0, and the upper endpoint
assigns all unknown outcomes 1. NOT_RUN cells are excluded because they were
not attempted. These bounds make clear how much the successful-only coverage
or interval-inclusion estimate could move if numerical failures are
informative. They do not replace the evaluable-only estimate or its Wilson
interval, and a wide bound is a warning that calibration conclusions are
sensitive to evaluability.

The per-cell `attempted` count and `evaluable_fraction` denominator exclude
planned `NOT_RUN` rows; the separate `not_run` count preserves those planned
but unexecuted cells. This avoids understating evaluability when an experiment
is interrupted before all planned cells are reached.
