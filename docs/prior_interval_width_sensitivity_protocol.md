# Prior interval-width sensitivity protocol

## Purpose

The prior-support audit found that credible-interval coverage can change substantially when the hard prior support changes. This companion experiment quantifies interval width directly, so improved coverage is not mistaken for better precision when it may simply reflect wider intervals.

## Design

The script reuses the same paired synthetic observations and posterior summaries produced by `run_prior_support_sensitivity.py`. For each truth (52, 70, 88), prior ([50,90], [40,100], [30,110]), and method (uncorrected and selection-corrected), it reports mean, median, and 90th-percentile widths for 68% and 95% credible intervals. It also emits per-observation paired differences in width relative to the baseline prior [50,90].

## Artifacts

- `interval_width_summary.csv`: 18 rows (truth × prior × method).
- `paired_interval_width_differences.csv`: 2,400 rows for 300 observations per truth, two wider priors, two methods.
- The script checks row counts, CSV schemas, nonempty outputs, and basic interval-width ordering in tests.

## Reproduction

```bash
python scripts/run_prior_interval_width_sensitivity.py --replicates 300 --seed 20261024 --grid-size 4001 --output-dir results/prior-interval-width
python -m pytest tests/test_prior_interval_width_sensitivity.py -q
```

The workflow runs Python 3.10, 3.11, and 3.12 and uploads the CSV outputs.

## Interpretation

A wider interval can increase coverage without improving the quality of the likelihood or removing bias. Compare width, coverage, bias, and RMSE jointly. The prior is not to be selected solely because it makes coverage closer to a nominal target. This is a simplified synthetic model, not a real detector or galaxy-catalogue analysis.
