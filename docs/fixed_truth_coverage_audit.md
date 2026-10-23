# Fixed-truth coverage audit — preliminary results

**Status:** completed computational audit; exploratory and not yet a calibration certification.  
**Workflow:** [Fixed-truth coverage audit, run 39961223917](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39961223917)  
**Artifact:** `fixed-truth-coverage-audit` (contains `summary.csv` and `per_run.csv`).  
**Pinned simulator revision:** `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`.

## Design

The audit runs 100 independently seeded replicates at each fixed input value H0 = 55, 60, 65, 70, 75, 80, and 85 km s⁻¹ Mpc⁻¹. It compares two simplified simulator conditions: `known_host` and `unknown_host_complete_catalogue`. That is 14 condition/truth cells and 1,400 posterior fits. It records posterior-median error, RMSE, central 68% and 95% posterior intervals, empirical coverage, and mean interval widths.

## Summary

All values below are fractions for coverage and km s⁻¹ Mpc⁻¹ for errors and widths.

| True H0 | Scenario | Mean signed median error | RMSE | 68% coverage | 95% coverage | Mean 68% width | Mean 95% width |
|---:|---|---:|---:|---:|---:|---:|---:|
| 55 | Known host | +5.97 | 8.29 | 0.68 | 0.98 | 15.99 | 29.37 |
| 55 | Unknown host, complete catalogue | +9.36 | 11.79 | 0.56 | 1.00 | 19.05 | 32.56 |
| 60 | Known host | +4.24 | 7.98 | 0.70 | 0.97 | 17.02 | 30.68 |
| 60 | Unknown host, complete catalogue | +8.01 | 11.09 | 0.60 | 1.00 | 20.56 | 33.88 |
| 65 | Known host | +2.10 | 7.19 | 0.78 | 1.00 | 17.38 | 31.06 |
| 65 | Unknown host, complete catalogue | +4.61 | 8.36 | 0.94 | 1.00 | 21.44 | 34.97 |
| 70 | Known host | +1.26 | 7.24 | 0.75 | 0.99 | 16.99 | 29.97 |
| 70 | Unknown host, complete catalogue | +1.87 | 6.75 | 0.92 | 0.98 | 22.20 | 35.51 |
| 75 | Known host | −1.39 | 7.15 | 0.74 | 0.98 | 16.15 | 28.39 |
| 75 | Unknown host, complete catalogue | −1.41 | 4.79 | 0.94 | 1.00 | 23.05 | 36.38 |
| 80 | Known host | −3.27 | 6.89 | 0.80 | 0.99 | 15.07 | 26.42 |
| 80 | Unknown host, complete catalogue | −6.06 | 7.64 | 0.85 | 1.00 | 22.83 | 36.16 |
| 85 | Known host | −4.80 | 6.76 | 0.78 | 0.97 | 13.31 | 23.19 |
| 85 | Unknown host, complete catalogue | −9.25 | 9.88 | 0.78 | 0.98 | 23.20 | 36.47 |

## Interpretation and important limitations

1. **The known-host results are not uniformly calibrated at nominal 68% coverage.** Observed coverage ranges from 68% to 80% across truth values. The nominal 68% interval has Monte Carlo standard errors of roughly 4–5 percentage points at n=100, so individual deviations are noisy; nevertheless, the truth-dependent signed errors and interval widths merit investigation.
2. **Unknown-host inference shows substantially broader intervals in these runs**, and often larger positive errors at lower truths and negative errors at higher truths. At H0=55, its 68% coverage is 56%, while at H0=70 it is 92%. This pattern is descriptive of this simulator and configuration only.
3. **The parameter grid is close to hard posterior bounds.** The saved per-run intervals visibly use a posterior support around 50–90 km s⁻¹ Mpc⁻¹. Truth values 55 and 85 are near those boundaries, so boundary truncation can affect posterior medians, interval widths, and coverage. This must be addressed before treating the coverage curve as an intrinsic property of the method.
4. The `unknown_host_complete_catalogue` condition is a simplified model, not a realistic representation of a full LVK dark-siren analysis. The audit does not use real GWTC/LVK events or a real galaxy survey.
5. Each cell has only 100 replicates. Coverage estimates should be accompanied by binomial Monte Carlo uncertainty; a follow-up with substantially more replicates is needed after the posterior-support and simulator assumptions are reviewed.
6. These results **do not establish performance on real observations** and should not be presented as a measurement of H0 or as evidence that real dark-siren inference has the same bias.

## Next audit steps

- Inspect and document the posterior support and prior in `posterior_h0`; test a wider support with the same simulations to quantify boundary effects.
- Run simulation-based calibration/rank diagnostics under a fully matched generative and inference model.
- Increase replicates only after fixing the support/prior question; report Monte Carlo uncertainty for coverage.
- Add tests for truths near and well inside the posterior support.
- Keep the synthetic study explicitly labelled as a controlled simulation until validated against a realistic selection function and public gravitational-wave/catalogue data.

## Reproducibility

The CI workflow validates 14 summary rows and 1,400 per-run rows and uploads the CSVs as an Actions artifact. The workflow and test suite completed successfully on 23 October 2026. The output artifact is the source of the numerical values in this report.
