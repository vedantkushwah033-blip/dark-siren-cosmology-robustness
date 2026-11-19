# Adaptive-CDF tail audit readout (2026-11-19)

**Workflow:** [44464048418](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44464048418)  
**Source artifacts:** immutable 200-replicate run 44440150116 and 1,000-replicate run 44446994455  
**Audit commit:** `eb516373afc37bdf3edec2bc3d423d6ec4a2f676`  
**Status:** targeted adaptive-CDF diagnostic completed; numerical-diagnostics CI on Python 3.10 and 3.11 passed, Python 3.12 was still in progress at the time of inspection.

## Scope

The audit independently regenerates selected parent populations, detected events, and prior-predictive H0 truth draws from the archived seed schedule. It evaluates the scalar full-parent selected-event log posterior and integrates its CDF using adaptive quadrature. It compares those adaptive references with a 4,481-point uniform-grid calculation.

This is a targeted diagnostic, not a full 200- or 1,000-replicate posterior-grid convergence audit. It does not alter the production likelihood, prior support, or tolerances, and it does not establish SBC calibration.

## 200-replicate source: four previously problematic draws

Selected replicates: 100, 174, 156, and 62. Source grid: 561 points. Adaptive/uniform comparison grid: 4,481 points.

| Diagnostic | Maximum across four replicates |
|---|---:|
| Archived 561-grid rank vs adaptive rank | 0.108485 |
| Archived 561-grid median vs adaptive median | 0.091211 km s^-1 Mpc^-1 |
| 4,481-grid rank vs adaptive rank | 0.000875 |
| 4,481-grid maximum quantile difference | 0.003161 km s^-1 Mpc^-1 |
| Adaptive reported relative quadrature error estimate | 1.16e-14 |

The original 561-point grid is materially inaccurate for these selected boundary cases. The 4,481-point grid agrees closely with the adaptive CDF for rank and reported quantiles in these four cases.

The adaptive posterior medians were 20.221, 21.483, 22.653, and 22.241 km s^-1 Mpc^-1, respectively. Their lower 2.5% quantiles were 20.009, 20.138, 21.030, and 20.177. These posteriors are concentrated close to the lower H0 prior boundary; this is a boundary-dominated prior-predictive stress case, not evidence of an observationally meaningful low H0 value.

## 1,000-replicate source: ten lowest-q975 draws

Selected replicates: 208, 17, 622, 538, 883, 89, 158, 618, 847, and 705. Source grid: 1,121 points. Adaptive/uniform comparison grid: 4,481 points.

| Diagnostic | Maximum across ten replicates |
|---|---:|
| Archived 1,121-grid rank vs adaptive rank | 0.010188 |
| Archived 1,121-grid median vs adaptive median | 0.006805 km s^-1 Mpc^-1 |
| 4,481-grid rank vs adaptive rank | 0.000610 |
| 4,481-grid maximum quantile difference | 0.003492 km s^-1 Mpc^-1 |
| Adaptive reported relative quadrature error estimate | 1.14e-14 |

The 4,481-point grid agrees closely with adaptive CDF ranks and quantiles in these ten selected draws. The 1,121-point source grid has larger rank discrepancies in some cases, even where median differences are small.

## Interpretation and limits

1. The adaptive reference strongly indicates that the earlier 561-point SBC grid is inadequate for the selected sharp, lower-boundary posteriors.
2. A 4,481-point grid resolves the rank and quantile discrepancies in the targeted four- and ten-replicate subsets under this scalar likelihood implementation.
3. This does **not** certify density-L1 convergence. The adaptive audit reports CDF ranks and quantiles, not a continuous-density L1 comparison against an independently represented adaptive density.
4. The very small quadrature error values are SciPy adaptive-integrator error estimates, not formal mathematical error bounds.
5. The targeted subsets were selected for problematic/boundary-tail behaviour. Their diagnostics must not be generalized to every replicate.
6. No full-sample SBC rank/coverage result is revised by this targeted audit. The original full-sample summaries remain as reported until a complete, independently audited recomputation is performed.
7. The prior boundary is scientifically important: many selected posterior draws place most mass within the first 1%–5% of the prior interval. This is a feature of the declared synthetic prior-predictive experiment and likelihood, and must be explicitly discussed.

## Next required numerical gate

Before claiming full SBC grid convergence:

- Run a full 200-replicate comparison with at least 4,481 points and preserve the complete lifecycle denominator.
- Extend the independent audit to a density-L1 comparison using a well-defined adaptive or independently refined density representation, not only CDF/quantile comparisons.
- Run the same targeted/full-sample checks for the 1,000-replicate experiment.
- Keep prior-boundary mass diagnostics and report boundary-dominated replicates rather than dropping them.
- Retain the original tolerances; any failure remains a failure.

**Conclusion:** targeted CDF/rank/quantile convergence is encouraging at 4,481 points, while full-sample density convergence and SBC recalculation remain open. This is not a final numerical-convergence pass.
