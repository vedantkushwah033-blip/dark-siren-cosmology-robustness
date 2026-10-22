# Repeated-sampling performance metrics

`summarize_interval_performance` computes several descriptive metrics from
one row per simulation replicate:

- **Bias:** mean estimate minus true value.
- **RMSE:** square root of mean squared estimation error.
- **Mean interval width:** average upper minus lower bound.
- **Empirical coverage:** fraction of intervals containing the true value.
- **Coverage Monte Carlo standard error:** approximate simulation uncertainty
  in the observed coverage fraction, treating replicates as independent.

For an observed coverage \hat p from (n) independent replicates, the
reported standard error is \sqrt{\hat p(1-\hat p)/n}. It is a
quick uncertainty summary, not an exact interval for coverage.

## Interpretation

For nominal 68% intervals, coverage near 94% is not nominal calibration. It
suggests conservatism or a mismatch between the interval procedure and the
repeated-sampling design, but the metric alone does not identify the cause.
Check the data-generating process, prior, interval construction, grid
resolution, boundaries, and whether the true parameter is fixed or sampled
from the prior.

These metrics do not replace simulation-based calibration. SBC asks a
different question: under a specified prior-generative model, are posterior
rank statistics behaving as expected? Nor do these metrics establish that a
simulation is realistic or that selection effects are correctly modelled.
