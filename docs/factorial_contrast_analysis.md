# Factorial contrast analysis

After the factorial sensitivity workflow has produced
`results/factorial-sensitivity/per_run_results.csv`, run:

```bash
python scripts/analyze_factorial_contrasts.py \
  --input results/factorial-sensitivity/per_run_results.csv \
  --output results/factorial-sensitivity/contrast_analysis.json \
  --bootstrap 5000 \
  --seed 20261023
```

The script reports scenario-level means for coverage, interval widths, errors,
and median pulls. It also estimates bootstrap 95% intervals for prespecified
scenario contrasts. Each difference is scenario B minus scenario A.

## Interpretation

- Coverage is computed per run by checking whether the fixed truth,
  (H_0=70) km/s/Mpc, lies within the corresponding credible interval.
- Widths are the posterior interval widths in km/s/Mpc.
- Median error is posterior median minus the fixed truth.
- RMSE is computed across posterior-median errors.
- Bootstrap intervals represent Monte Carlo variability conditional on the
  simulation design. They do not represent uncertainty about real-world H0.
- These are marginal, independent-scenario bootstrap comparisons. The current
  factorial script assigns separate seed blocks to scenarios, so do not
  interpret the contrasts as paired-run effects.
- Eight prespecified contrasts and seven metrics produce 56 intervals. Treat
  them as exploratory; intervals are not adjusted for multiple comparisons.
- Interpret coverage alongside bias, RMSE, and interval widths. A coverage
  increase can occur because intervals are wider.

This analysis is a diagnostic of the pinned synthetic simulator, not an
observational LVK analysis. It should not be used to claim a real H0 shift.
