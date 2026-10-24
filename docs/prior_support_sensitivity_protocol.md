# Prior-support sensitivity audit

**Status:** implementation and CI workflow added; numerical results should not be reported until the workflow artifact has been inspected.

## Question

How sensitive are fixed-truth bias, RMSE, and credible-interval coverage to the hard prior support in the synthetic two-host model? This isolates prior support as an inferential design choice while holding each generated observation fixed across prior choices.

## Design

- Fixed truths: (H_0=52,70,88).
- Prior intervals: ([50,90]), ([40,100]), and ([30,110]).
- 300 detected observations per truth, generated once per truth and reused across all three priors.
- Random seed: 20261024.
- Posterior grid: 4,001 points on [25,125].
- Host redshifts: 0.02 and 0.04; pre-detection weights 0.70 and 0.30.
- Gaussian distance noise: 20 Mpc.
- Detection condition: (0<d_{obs}<180) Mpc.
- Both uncorrected and selection-normalized inference are evaluated for every observation and prior.

The same observed data are reused across priors to make prior comparisons paired rather than confounded by different random draws. The data-generating process is unchanged. The wider prior choices are sensitivity scenarios, not a claim that a wider prior is automatically preferable.

## Metrics and uncertainty

For each truth, prior, and method, the summary reports median-estimator bias, RMSE, empirical 68% and 95% credible-interval coverage, and Wilson score intervals for the coverage proportions. CSV files are read back and checked for schema and row-count consistency.

## Reproduction

```bash
python scripts/run_prior_support_sensitivity.py --replicates 300 --seed 20261024 --grid-size 4001 --output-dir results/prior-support-sensitivity
python -m pytest tests/test_prior_support_sensitivity.py -q
```

The CI workflow runs the experiment under Python 3.10, 3.11, and 3.12 and uploads the per-run and summary CSV files.

## Interpretation limits

This is a toy-model sensitivity analysis, not an observational gravitational-wave result. Changing the prior changes the inferential model but does not change the data-generating process. If coverage improves with wider support, that alone does not prove the wider prior is objectively correct; it only shows that conclusions under the original support were sensitive to that choice. Conversely, wider support can increase uncertainty or alter estimates. Compare all methods and truths, inspect Monte Carlo intervals, and avoid choosing a prior solely because it produces more favorable coverage.
