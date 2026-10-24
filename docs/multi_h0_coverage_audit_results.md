# Multi-$H_0$ interval coverage audit

**Status:** completed; GitHub Actions passed 25,200 posterior fits across 84 design cells.  
**Workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40066479000  
**Scope:** synthetic low-redshift two-host model only. These are not observational constraints on $H_0$.

## Design

The injected Hubble parameter values were $H_0 \in \{52,55,60,70,80,85,88\}$ km s$^{-1}$ Mpc$^{-1}$ under a uniform inference prior $[50,90]$. The two candidate redshifts were $z_1=0.02$ and $z_2 \in \{0.025,0.04,0.06\}$. Two fixed noise scales (3% and 10% of the reference distance scale) and two inference weightings (correct 70/30 and reversed 30/70) were tested. Data were generated using 70/30 latent-host probabilities, with noise scale fixed within each injected-$H_0$ design and shared by generation and inference.

Each of the 7 × 3 × 2 × 2 design cells used 300 replicates: 25,200 posterior fits in total. The report gives average performance across the six redshift/noise combinations for each injected $H_0$ and inference weighting. The cell-level CSV files preserve the full detail.

## Aggregate results

Coverage and RMSE below are averaged across the six design combinations at each injected $H_0$. RMSE is in km s$^{-1}$ Mpc$^{-1}$.

| Injected $H_0$ | Weighting | Mean 68% coverage | Mean 95% coverage | Mean median-estimate RMSE |
|---:|---|---:|---:|---:|
| 52 | Correct | 55.2% | 94.2% | 6.40 |
| 52 | Reversed | 45.6% | 89.4% | 9.55 |
| 55 | Correct | 68.6% | 95.8% | 5.56 |
| 55 | Reversed | 57.2% | 93.6% | 8.40 |
| 60 | Correct | 71.6% | 96.3% | 5.19 |
| 60 | Reversed | 62.8% | 95.3% | 7.50 |
| 70 | Correct | 71.8% | 97.0% | 5.54 |
| 70 | Reversed | 69.3% | 97.2% | 6.20 |
| 80 | Correct | 77.6% | 96.8% | 5.51 |
| 80 | Reversed | 81.4% | 97.0% | 3.97 |
| 85 | Correct | 74.3% | 98.7% | 5.83 |
| 85 | Reversed | 78.7% | 98.8% | 4.19 |
| 88 | Correct | 52.2% | 95.7% | 7.17 |
| 88 | Reversed | 59.2% | 96.8% | 5.29 |

These are averages across six heterogeneous design cells, not pooled coverage from one homogeneous data-generating process. The nominal targets are 68% and 95%.

## Main findings

1. **Prior-boundary effects are visible.** Near the lower edge ($H_0=52$) and upper edge ($H_0=88$), the 68% interval coverage is well below the nominal 68% target in many settings. At these points the prior boundary truncates posterior support and affects interval construction and bias.
2. **Incorrect weights can worsen results, especially toward the lower boundary.** At $H_0=52$, reversed weights raise average RMSE from 6.40 to 9.55 and reduce mean 68% coverage from 55.2% to 45.6%. At $H_0=55$ and 60, the reversed-weight results are also worse on average.
3. **The effect is not uniform across the prior.** At $H_0=80$ and 85, reversed weights yield lower average RMSE in this particular design. This is a reminder that a misspecified model can sometimes look better on a chosen metric in a finite controlled experiment; it does not mean the wrong weights are preferable generally.
4. **95% coverage is often above target near the upper boundary.** Some of this pattern may reflect wider intervals and posterior shapes under the bounded prior. It needs further diagnosis rather than a blanket claim that intervals are calibrated.
5. **Monte Carlo uncertainty matters.** Each individual design cell has 300 replicates; for coverage near 68%, the binomial Monte Carlo standard error is about 2.7 percentage points, and near 95% it is about 1.3 points. Because aggregate values average heterogeneous cells, they should be interpreted descriptively.

## Limitations

- This uses the low-redshift approximation $d_L \approx cz/H_0$, two fixed candidate redshifts, Gaussian distance noise, a uniform bounded prior, and no realistic detector selection.
- The audit is a controlled stress test, not a complete simulation-based calibration test and not validation against LVK data or a real galaxy catalogue.
- Grid-based posterior summaries may introduce numerical error. A grid-refinement audit is still needed, particularly near boundaries and for multimodal posteriors.
- These aggregate averages must not replace inspection of cell-level results when writing the paper.

## Reproducibility

- Script: `scripts/run_multi_h0_coverage_audit.py`
- Tests: `tests/test_multi_h0_coverage_audit.py`
- Workflow: `.github/workflows/multi-h0-coverage-audit.yml`
- Artifact: `multi-h0-coverage-audit-results` (contains per-run and cell-level summary CSVs; retention is time-limited).
