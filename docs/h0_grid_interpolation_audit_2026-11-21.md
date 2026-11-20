# H0 posterior-grid interpolation audit (2026-11-21)

## Scope and provenance

This diagnostic reuses the archived 30-event quadrature artifact from PR #141 (workflow 44806809662; artifact ID 14669447416). It processes all 3 curvature settings × 2 inclusion degrees × 3 quadrature orders × 2 validation sizes × 5 inference arms = **540 posterior likelihood curves**.

The audit starts from each archived 281-point H0 log-likelihood curve on [30,110]. It uses shape-preserving cubic interpolation (PCHIP) of the *log-likelihood curve* and re-integrates the resulting interpolant on 561-, 1,121-, and 2,241-point grids. The original 281-point trapezoidal normalization and quantiles are compared with the refined-grid integration of that same interpolant.

Reproduction:
```bash
python scripts/audit_h0_grid_interpolation.py degree-matched-30event-quadrature-sensitivity.zip --output-dir h0_grid_interpolation_audit
```

The summary table is in `results/h0_grid_interpolation_audit_summary_2026-11-21.csv`.

## Results

The maximum values below are taken across the 36 curves per arm at each refined resolution.

| Arm | Grid | Max posterior L1 change | Max abs median shift | Max abs 2.5% shift | Max abs 97.5% shift |
|---|---:|---:|---:|---:|---:|
| A catalogue-only | 561 | 1.75e-5 | 0.000604 | 0.001085 | 0.005507 |
| A catalogue-only | 1,121 | 2.19e-5 | 0.000699 | 0.001373 | 0.006844 |
| A catalogue-only | 2,241 | 2.30e-5 | 0.000734 | 0.001448 | 0.007181 |
| B parent oracle | 561 | 7.09e-7 | 0.000094 | 0.005493 | 0.003755 |
| B parent oracle | 1,121 | 7.52e-7 | 0.000123 | 0.007476 | 0.004778 |
| B parent oracle | 2,241 | 7.55e-7 | 0.000126 | 0.007693 | 0.005090 |
| C known-probability IPW | 561 | 8.45e-7 | 0.000149 | 0.007148 | 0.006899 |
| C known-probability IPW | 1,121 | 8.97e-7 | 0.000167 | 0.008227 | 0.008811 |
| C known-probability IPW | 2,241 | 9.00e-7 | 0.000175 | 0.008711 | 0.009082 |
| D MAP plug-in IPW | 561 | 1.00e-6 | 0.000143 | 0.008422 | 0.008425 |
| D MAP plug-in IPW | 1,121 | 1.06e-6 | 0.000169 | 0.009308 | 0.009648 |
| D MAP plug-in IPW | 2,241 | 1.07e-6 | 0.000181 | 0.009722 | 0.010098 |
| E integrated IPW | 561 | 7.79e-7 | 0.000706 | 0.007611 | 0.008938 |
| E integrated IPW | 1,121 | 8.27e-7 | 0.000782 | 0.008321 | 0.009712 |
| E integrated IPW | 2,241 | 8.29e-7 | 0.000804 | 0.008505 | 0.009898 |

Between 1,121 and 2,241 points, the largest absolute median change across all arms and curves is below 0.00006 km s^-1 Mpc^-1 for A and below 0.000031 for E. The largest 2.5th/97.5th percentile changes across all arms are below 0.00050/0.00055 km s^-1 Mpc^-1, respectively. The largest posterior L1 difference between the original 281-grid normalized curve and the 2,241-grid normalized interpolant is 2.31e-5 (arm A).

## Interpretation

For these archived curves, trapezoidal posterior integration on the 281-point grid is already very stable relative to finer integration of the PCHIP interpolant. This makes H0-grid discretization unlikely to explain the several-km s^-1 Mpc^-1 posterior-median offsets seen in the fixed-realization experiment.

This is **not** an exact likelihood recomputation at new H0 values. The likelihoods were not reevaluated from event-level data on the refined grids; only the archived 281-point curves were interpolated. Therefore this is a posterior-integration/interpolation audit, not the final exact H0-grid convergence gate. PCHIP can smooth between saved values and cannot reveal structure that the original grid failed to sample.

## Next gate

1. Add an exact-evaluation route that recomputes event likelihood and selection-efficiency terms at arbitrary H0 values from the saved event and host-density inputs.
2. Compare exact 281/561/1,121 (and, if needed, 2,241) point posteriors against the interpolation diagnostic.
3. Keep quadrature convergence, H0-grid convergence, and IPW weight stability as separate diagnostics.
4. Only after exact grid validation proceed to the independent prior-predictive calibration design.

No scientific superiority, calibration, or observational H0 claim is supported by this interpolation-only check.
