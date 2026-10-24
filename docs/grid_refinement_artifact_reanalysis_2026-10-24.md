# Grid-refinement artifact reanalysis — 2026-10-24

## Source and reproducibility

This report reanalyses the archived `grid-refinement-audit-results` artifact from workflow run [40088129105](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40088129105), commit `ebfb7ad22d1aabf90d5269c39dcbff17dc12efc4`. The artifact contains `per_run.csv` (21,600 rows) and `summary.csv` (36 design-cell summaries). Its recorded SHA-256 digest is `df1e40bb421ab026fa49ef9f10aa8d3a0d6f46bb741f84097ab2995663092d31`.

The analysis below was computed from the raw per-run CSV, pairing each simulated observation and weighting method across grid sizes. There are 7,200 paired cases for each grid-size comparison (3 true H0 values × 3 second-host redshifts × 2 noise fractions × 200 replicates × 2 weighting schemes). Differences are absolute and in H0 units (km s^-1 Mpc^-1).

## Paired grid-size discrepancy summary

| Posterior summary | Comparison | Median absolute difference | 95th percentile | 99th percentile | Maximum | Cases > 0.001 | Cases > 0.01 |
|---|---|---:|---:|---:|---:|---:|---:|
| Median | 1,001 vs 4,001 | 0.000095 | 0.000589 | 0.000964 | 0.002502 | 56 / 7,200 (0.78%) | 0 |
| Median | 2,001 vs 4,001 | 0.000019 | 0.000118 | 0.000197 | 0.000572 | 0 / 7,200 | 0 |
| Lower 68% endpoint | 1,001 vs 4,001 | 0.000236 | 0.001272 | 0.001874 | 0.002487 | 562 / 7,200 (7.81%) | 0 |
| Lower 68% endpoint | 2,001 vs 4,001 | 0.000047 | 0.000250 | 0.000412 | 0.000586 | 0 / 7,200 | 0 |
| Upper 68% endpoint | 1,001 vs 4,001 | 0.000183 | 0.001261 | 0.002270 | 0.006801 | 660 / 7,200 (9.17%) | 0 |
| Upper 68% endpoint | 2,001 vs 4,001 | 0.000036 | 0.000261 | 0.000462 | 0.001508 | 2 / 7,200 (0.03%) | 0 |
| Lower 95% endpoint | 1,001 vs 4,001 | 0.000484 | 0.002339 | 0.003704 | 0.006135 | 1,986 / 7,200 (27.58%) | 0 |
| Lower 95% endpoint | 2,001 vs 4,001 | 0.000099 | 0.000460 | 0.000760 | 0.001213 | 16 / 7,200 (0.22%) | 0 |
| Upper 95% endpoint | 1,001 vs 4,001 | 0.000208 | 0.002159 | 0.004210 | 0.130453 | 1,463 / 7,200 (20.32%) | 1 / 7,200 (0.014%) |
| Upper 95% endpoint | 2,001 vs 4,001 | 0.000041 | 0.000426 | 0.000795 | 0.024330 | 43 / 7,200 (0.60%) | 1 / 7,200 (0.014%) |

The median discrepancies are small, and the 2,001-vs-4,001 differences are smaller than the 1,001-vs-4,001 differences for all listed summaries. However, the 95% endpoints are more sensitive than the posterior median. The maximum upper-95% discrepancy remains an important tail-specific exception. Reporting only the median discrepancy would conceal it.

The threshold counts are descriptive diagnostics, not universal acceptance criteria. A tolerance should be justified in advance by the precision required for the intended scientific interpretation.

## Exact largest upper-95% discrepancy

The largest upper-95% endpoint discrepancy occurs for:

- True synthetic H0: 52
- First-host redshift: 0.02
- Second-host redshift: 0.04
- Noise fraction: 0.03
- Replicate index: 138 (zero-based)
- Weighting: reversed, weights (0.3, 0.7)
- Observed synthetic distance: 119.47415214556634 Mpc

| Grid points | Posterior median | Lower 68% | Upper 68% | Lower 95% | Upper 95% |
|---:|---:|---:|---:|---:|---:|
| 1,001 | 51.713915 | 50.514052 | 53.692477 | 50.080249 | 60.848373 |
| 2,001 | 51.713591 | 50.513957 | 53.691894 | 50.080253 | 60.742249 |
| 4,001 | 51.713486 | 50.513928 | 53.691714 | 50.080256 | 60.717920 |

The median changes by only about 0.000429 between the coarsest and finest of these three grids, while the upper 95% endpoint changes by 0.130453. This is evidence that tail summaries need separate checks; stable posterior medians alone are insufficient to establish convergence of all credible-interval endpoints.

## Targeted finer-grid follow-up

The targeted artifact from numerical-diagnostics workflow run [40089986883](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40089986883) reports the following upper-95% sequence for this specific observation:

| Grid points | Upper 95% endpoint | Absolute change from previous grid |
|---:|---:|---:|
| 1,001 | 60.848373 | — |
| 2,001 | 60.742249 | 0.106124 |
| 4,001 | 60.717920 | 0.024330 |
| 8,001 | 60.711905 | 0.006015 |
| 16,001 | 60.710428 | 0.001477 |
| 32,001 | 60.710055 | 0.000373 |
| 64,001 | 60.709962 | 0.000093 |

For this case, the refinement sequence approaches approximately 60.710, with successively smaller changes. This supports grid discretization and quantile interpolation as the explanation for this particular discrepancy. It does not prove monotonic convergence for every observation or design cell.

## Interpretation and recommended decision

1. The 1,001-point grid should not be treated as uniformly adequate for tail summaries in this design.
2. A 4,001-point grid appears stable for typical cases but leaves a visible error in the identified tail case.
3. The targeted 64,001-point result is a useful diagnostic reference for this one case, not proof that every 4,001-point posterior is accurate to the same tolerance.
4. A prespecified tolerance should be evaluated for every posterior summary and every design cell, with explicit reporting of maxima and tail-specific failures.
5. If computational cost matters, consider a tested adaptive/local refinement strategy for tail quantiles, but benchmark it against fine-grid references over a broader set of difficult cases before adopting it.
6. The two integration rules applied to the same sampled posterior are numerical cross-checks, not fully independent posterior validation.

## Scope limitation

This is a numerical sensitivity audit of a simplified synthetic two-host model. It does not validate the model's astrophysical assumptions, realistic gravitational-wave selection effects, galaxy-catalogue completeness, or any observational estimate of H0.
