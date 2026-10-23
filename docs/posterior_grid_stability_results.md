# Posterior-grid stability diagnostic: verified output

**Status:** Executed in GitHub Actions on 2026-10-23. This is a deterministic numerical diagnostic for one illustrative synthetic posterior, not a result from the factorial simulation and not an observational measurement of $H_0$.

Workflow: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39898692225  
Commit: `75cccdd63587eb38ad966e25c0cfa8d7f520b59b`

## Results

| Case | Grid points | Prior | Posterior median (km s⁻¹ Mpc⁻¹) | 68% interval | 95% interval |
|---|---:|---|---:|---|---|
| Baseline | 1,601 | [50, 90] | 70.41548 | [59.72313, 82.57368] | [51.86912, 88.75626] |
| Finer baseline | 3,201 | [50, 90] | 70.41548 | [59.72314, 82.57368] | [51.86912, 88.75626] |
| Wider numerical grid, same prior | 6,401 | [50, 90] | 70.41605 | [59.72099, 82.57754] | [51.86303, 88.76206] |
| Wider numerical grid and wider prior | 6,401 | [30, 110] | 68.20008 | [40.93599, 89.98552] | [31.72009, 105.07879] |

All four numerical posterior integrals were 1 to displayed precision.

## Interpretation

1. **Baseline grid refinement is stable for this example.** Doubling the number of points from 1,601 to 3,201 changed the median by about $1.47\times10^{-6}$ km s⁻¹ Mpc⁻¹. The largest change among the reported 68% and 95% interval endpoints was about $5.46\times10^{-6}$ km s⁻¹ Mpc⁻¹.
2. **Extending the numerical grid while keeping the prior fixed has a small effect.** Compared with the finer baseline, the median changed by about 0.00057 km s⁻¹ Mpc⁻¹ and the largest reported interval-endpoint change was about 0.0061 km s⁻¹ Mpc⁻¹. This is consistent with small discretization differences for this specific setup.
3. **Changing the prior is not a grid-resolution test.** Widening the prior from [50, 90] to [30, 110] changes the median and substantially widens the credible intervals. This illustrates prior sensitivity for this example; it does not show that either prior is objectively correct.
4. **Boundary-band mass is diagnostic, not a pass/fail criterion.** The posterior mass in the lower and upper 2%-of-prior-width bands was about 1.08% and 1.60% for the baseline prior, versus about 2.31% and 0.72% for the wider prior. These figures should be interpreted alongside the posterior shape and the chosen prior.

## Scope and limitations

- This is one deterministic, illustrative posterior configuration; it is not a Monte Carlo coverage study.
- It does not establish robustness across all catalogue-incompleteness or host-weighting scenarios.
- It does not validate the low-redshift approximation against a full flat-$\Lambda$CDM inference.
- The prior-width comparison intentionally changes the inferential assumption, so it should be reported separately from numerical grid convergence.
- The output supports a limited statement: the posterior summaries for this example are numerically stable under baseline grid refinement, while the prior bounds materially affect the result.

## Reproducibility

Run from the repository root after installing the package and development dependencies:

```bash
python scripts/run_posterior_grid_stability.py --output artifacts/posterior_grid_stability.csv
```

The CI workflow runs this diagnostic on Python 3.10, 3.11, and 3.12 and uploads the CSV as the `posterior-grid-stability` artifact. The artifact is retained by GitHub Actions for a limited time; see the linked workflow run.
