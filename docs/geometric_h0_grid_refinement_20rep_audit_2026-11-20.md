# Exact geometric H0-grid refinement audit: 20-replicate five-arm pilot

**Date:** 2026-11-20  
**Status:** PASS for the 1,121 → 2,241 point geometric-grid comparison.

## Frozen source

- Experiment workflow: [44731915717](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44731915717)
- Source artifact ID: 14628789845
- Source artifact SHA-256: `812d05e94f504e76ac80c9f2531ec4f99be389cf8367c2aa0a0485e6b56a314a`
- Design: 20 independent prior-predictive truths; 3 curvature strata (-1, 0, +1); validation sizes 100, 300, 500; five arms A–E; 30 detected events; Uniform(30,110) H0 prior.
- Archived H0 grid: 1,121 geometric points. Reference grid: 2,241 geometric points.
- Curves audited: 20 × 3 × 3 × 5 = **900**.

The source experiment itself passed its independent lifecycle and aggregate audit: all 900 planned cells succeeded, all 60 child runs passed, and the artifact contains 20 distinct simulated truths.

## Exact-grid method

For each frozen replicate/curvature/validation-size cell, the audit holds the archived KDE host-density estimates fixed, reconstructs the selected likelihood on the archived 1,121-point grid, and recomputes event numerators and selection-efficiency denominators on the 2,241-point geometric grid. Arm E is re-integrated over its 25 shared coefficient nodes. Posterior integration is with respect to **dH0**.

The batched implementation is algebraically equivalent to summing each event's log numerator and selection denominator separately. A regression test compares the batched result against an eventwise reference calculation.

This isolates H0-grid integration error conditional on archived KDE densities. It does **not** regenerate KDEs from raw galaxy weights, validate the KDE bandwidth, or establish posterior calibration.

## Results

Independent local execution of the exact-grid audit algorithm covered all 900 rows in ten disjoint two-replicate batches. Every batch passed; the maximum reconstructed archived selected-log-likelihood discrepancy was (1.14\times10^{-12}).

| Diagnostic: 1,121 → 2,241 geometric points | Maximum observed change | Tolerance |
|---|---:|---:|
| Absolute posterior-rank change | 0.00004596 | 0.002 |
| Posterior-median change (km/s/Mpc) | 0.0003490 | 0.05 |
| Maximum change among nine reported quantiles (km/s/Mpc) | 0.0009942 | 0.05 |
| Posterior-density L1 distance | 0.00007995 | 0.001 |

All 900 curves satisfy the declared refined-grid tolerances. The 1,121-point geometric grid is therefore adequate against this 2,241-point reference for the frozen KDE densities and this finite design.

### Descriptive rank pattern at the refined grid

Mean posterior ranks over the 20 independent truths, shown by curvature and arm:

| Curvature | Catalogue-only A | Parent oracle B | Known-π IPW C | MAP IPW D | Integrated IPW E |
|---:|---:|---:|---:|---:|---:|
| -1 | 0.9761 | 0.4393 | 0.4663 | 0.6515 | 0.6277 |
| 0 | 0.9692 | 0.4393 | 0.4413 | 0.4546 | 0.4527 |
| +1 | 0.9568 | 0.4393 | 0.4311 | 0.2942 | 0.3068 |

These are descriptive pilot averages, not formal evidence of rank uniformity. The repeated validation-size rows share the same simulated truth within each replicate/curvature block, so they must not be counted as independent truths. The parent-oracle arm is a simulation diagnostic only.

## Interpretation and limitations

1. The 1,121 → 2,241 geometric-grid comparison is numerically stable under the declared rank, quantile, and density-L1 tolerances.
2. The catalogue-only rank distribution remains concentrated near one, while the estimated degree-one IPW arms shift strongly with the true inclusion curvature.
3. This is consistent with the project's misspecification concern, but the 20-truth pilot is not large enough for confirmatory calibration claims.
4. The audit is conditional on archived KDE densities. It is not an independent re-generation of the KDE pipeline from raw galaxy weights.
5. No observational H0, LVK/GWTC, or real galaxy-catalogue validation claim is made.

## Reproducibility

The reusable audit is `scripts/audit_geometric_prior_predictive_h0_grid.py`. Its eventwise-equivalence regression is in `tests/test_audit_geometric_prior_predictive_h0_grid.py`. The manual audit workflow is `.github/workflows/exact-geometric-prior-predictive-h0-grid-audit.yml`.
