# Confirmatory simulation precision plan

**Date:** 30 October 2026  
**Status:** design recommendation; no confirmatory run has been executed.

## Why 200 replicates per cell are not enough for precise coverage claims

For an empirical coverage estimate based on (n) independent simulation replicates, the approximate standard error is

[
mathrm{SE}(hat p) simeq sqrt{p(1-p)/n}.
]

At nominal 68% coverage:

| Replicates per cell | Approx. SE | Approx. 95% half-width |
|---:|---:|---:|
| 200 | 0.0330 | 0.0647 |
| 500 | 0.0172 | 0.0337 |
| 1,000 | 0.0149 | 0.0292 |
| 2,000 | 0.0105 | 0.0206 |

At nominal 95% coverage:

| Replicates per cell | Approx. SE | Approx. 95% half-width |
|---:|---:|---:|
| 200 | 0.0154 | 0.0302 |
| 500 | 0.0097 | 0.0191 |
| 1,000 | 0.0069 | 0.0135 |
| 2,000 | 0.0049 | 0.0096 |

These are planning approximations, not exact confidence intervals. The analysis must report Wilson intervals (or a justified alternative) using actual successes and attempted replicates, including failed fits in the failure accounting. Coverage denominators are successful fits only, so differential fit failure can create selection into the coverage calculation and must be reported alongside coverage.

## Recommended staged design

### Stage A — engineering smoke test

- 20 replicates per cell.
- Confirm deterministic seeds, artifact schema, independent reconstruction, and all failure paths.
- No scientific efficacy or calibration claim.

### Stage B — variance/power pilot

- 100–200 replicates per cell.
- Estimate paired-contrast variability and failure rates.
- Use this stage only to plan the confirmatory campaign; do not reuse it as an unacknowledged confirmatory sample.
- Freeze the final factor grid, primary estimand, stopping rules, and analysis code before Stage C.

### Stage C — confirmatory campaign

For a paper whose central claim includes 68% coverage calibration, target at least 1,000 attempted replicates per primary cell if computationally feasible. This gives an approximate 95% coverage half-width near 0.03 when true coverage is around 0.68, before accounting for fit failures. If the paper needs substantially tighter calibration statements, increase the replicate count according to a prespecified precision target.

For each cell, report:
- attempted, successful, and failed replicates;
- fit-failure and boundary-limited rates;
- 68% and 95% coverage with Wilson intervals;
- posterior bias, RMSE, interval width, and paired uncertainty;
- positivity summaries, ESS, maximum and tail inverse weights;
- density error and inclusion-probability error (simulation scoring only);
- sensitivity to posterior-grid resolution.

Do not select the replicate count after seeing which count makes the result statistically significant. Use a precision target, not a significance target.

## Pairing and estimands

All arms within a replicate must share the same parent population, catalogue realization, GW observations, and random-number provenance wherever the estimand permits. The primary contrast should be declared before execution, for example:

[
Delta_r =
|\hat H_{0,r}^{\mathrm{catalogue-only}}-H_{0,*}|
-
|\hat H_{0,r}^{\mathrm{estimated}}-H_{0,*}|.
]

Positive values favor estimated weighting on absolute posterior-mean error. Report the mean paired contrast and a paired bootstrap interval over replicate IDs. Also report differential fit failure rates; a paired complete-case contrast is conditional on both arms succeeding and is not a substitute for the failure analysis.

Coverage is a separate estimand. Do not infer calibration from lower RMSE or narrower intervals.

## Bootstrap interpretation

The current validation bootstrap resamples paired validation records with replacement and holds the observed catalogue and GW observations fixed. It assumes the logistic selection family and does not apply a finite-population correction. It is a conditional sensitivity analysis, not a full posterior uncertainty analysis. The confirmatory study must report this scope explicitly.

A design-based or finite-population correction should be implemented only after its sampling assumptions and estimator are derived and independently tested. Do not apply a generic correction formula without checking the validation sampling design.

## Decision gate

Proceed from Stage B to Stage C only if:
1. all artifact and independent-audit checks pass;
2. the smoke test shows no unresolved numerical or schema failures;
3. positivity and weight-tail behavior are understood across the full design;
4. the primary contrast and coverage estimands are frozen;
5. the compute budget supports the stated precision target.

If any condition fails, revise the design and rerun the smoke stage. Do not quietly omit difficult cells.
