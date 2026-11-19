# Targeted 1,000-replicate SBC boundary-grid follow-up

**Status:** targeted numerical stress audit passed its final adjacent-grid criterion on the selected subset; not a full-source convergence certification and not a revised SBC result.

## Source and workflow evidence

The immutable source is the independently audited 1,000-replicate full-parent SBC run described in the earlier readout. Two targeted grid-audit runs completed:

- [Workflow 44467800575](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44467800575), artifact 14502044964: 1,121 / 2,241 / 4,481 / 8,961 point comparison.
- [Workflow 44468985931](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44468985931), artifact 14503542171: extended comparison through 17,921 points.
- [Workflow 44468985819](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44468985819): tests PASS.

Both grid-audit jobs independently verified the source run and reconstructed their selected replicates from the immutable source artifact. The extended run selected ten boundary-stressed replicates: 208, 622, 538, 17, 158, 883, 89, 847, 705, and 224.

## Adjacent-grid results

Maximum differences across the ten selected replicates:

| Grid comparison | Maximum density L1 | Maximum quantile shift (km s^-1 Mpc^-1) | Maximum rank shift |
|---|---:|---:|---:|
| 1,121 -> 2,241 | 0.0172854 | 0.045001 | 0.008240 |
| 2,241 -> 4,481 | 0.0045387 | 0.013749 | 0.001658 |
| 4,481 -> 8,961 | 0.0011363 | 0.002755 | 0.000437 |
| 8,961 -> 17,921 | 0.0002839 | 0.000905 | 0.000152 |

The previously used density-L1 threshold is 0.001. The 4,481 -> 8,961 comparison narrowly exceeds it (0.0011363); the 8,961 -> 17,921 comparison is below it (0.0002839). Quantile and rank shifts also decrease substantially. Thus 17,921 points are needed to clear the adjacent-grid density criterion for this selected stress subset. Do not infer that 4,481 points are sufficient for all full-parent SBC replicates.

## Adaptive-CDF cross-check

The separate adaptive-quadrature audit selected the ten lowest archived q97.5 cases (replicates 208, 17, 622, 538, 883, 89, 158, 618, 847, 705). Against its adaptive CDF reference, the 4,481-point calculation had:
- maximum absolute posterior-rank difference 0.0006099;
- maximum posterior-quantile difference 0.003492 km s^-1 Mpc^-1;
- maximum archived 1,121-grid versus adaptive rank difference 0.010188.

The adaptive integrator's reported relative quadrature error estimate was at most 1.14e-14. This is the numerical integrator's estimate, not a formal mathematical error bound. The adaptive subset and the grid-refinement subset are not identical (replicate 618 versus 224), so these summaries are complementary targeted diagnostics, not one matched ten-case table.

## Interpretation and remaining limitation

The selected lower-tail posteriors are sharply concentrated near the lower H0 prior boundary. Their narrow widths make density and rank calculations sensitive to grid resolution. The results justify using a substantially finer grid for boundary-stressed cases and motivate a full-source adaptive or high-resolution convergence audit.

This does **not** establish:
- density-L1 convergence for all 1,000 replicates;
- full-sample rank or SBC stability at 17,921 points;
- a new or revised SBC coverage estimate;
- convergence for the five-arm catalogue-inclusion experiment.

The next numerical gate is a full-sample audit or a predeclared, reproducibly selected stress set with a justified worst-case envelope, including density-L1, quantile, rank, and boundary-mass checks. Keep the existing tolerances; do not drop difficult replicates or relabel targeted evidence as full-sample validation.
