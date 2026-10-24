# Coverage-transition map protocol

## Question

How do fixed-truth coverage, bias, RMSE, and credible-interval widths change as the true parameter moves from the lower to the upper edge of a bounded prior? Does widening the support change those patterns?

## Design

- True values: 52, 55, 60, 65, 70, 75, 80, 85, 88.
- Prior supports: [50,90] and [40,100].
- Methods: uncorrected and selection-corrected likelihood.
- Repetitions: 300 per truth value.
- Random seed: 20261025.
- Posterior grid: 4,001 points.
- Data generator and likelihood implementations are reused unchanged from the established fixed-truth/prior-support experiment.
- For each true value, a single set of generated detected observations is reused across both prior supports and both methods.

The expected output contains 5,400 per-run rows (9 truths × 2 priors × 300 replicates) and 36 summary rows (9 × 2 × 2 methods).

## Outcomes

For each truth/prior/method cell, report bias of the posterior median, RMSE, empirical 68% and 95% equal-tailed interval coverage, Wilson confidence intervals for the coverage fractions, and interval-width summaries inherited from the shared fixed-truth generator.

## Reproducibility

```bash
python scripts/run_coverage_transition_map.py --replicates 300 --seed 20261025 --grid-size 4001 --output-dir results/coverage-transition
python -m pytest tests/test_coverage_transition_map.py -q
```

CI runs the experiment on Python 3.10, 3.11, and 3.12 and uploads per-run and summary CSV artifacts.

## Interpretation and limitations

This is a synthetic low-redshift toy model. The denser grid helps locate where prior-boundary effects change but cannot establish behavior for arbitrary truth values, priors, catalogues, or real gravitational-wave analyses. Wilson intervals quantify Monte Carlo uncertainty in estimated coverage; they do not account for model misspecification. Compare coverage alongside bias, RMSE, and interval width. Do not choose a prior solely to maximize coverage in this finite experiment.
