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

The batched implementation is algebraically equivalent to summing each event's log numerator and selection denominator separately, but avoids a Python loop over the 30 events. A regression test compares the batched result against an eventwise reference calculation.

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

### Prior-predictive ranks at validation size 500

The following means use the 20 independent prior truths within each curvature stratum. Parent-oracle B is the same across curvature because it uses the full parent population.

| Curvature | Catalogue-only A | Parent oracle B | Known-π IPW C | MAP IPW D | Integrated IPW E |
|---:|---:|---:|---:|---:|---:|
| -1 | 0.9761 | 0.4393 | 0.4663 | 0.6520 | 0.6392 |
| 0 | 0.9692 | 0.4393 | 0.4413 | 0.4614 | 0.4559 |
| +1 | 0.9568 | 0.4393 | 0.4311 | 0.2982 | 0.3012 |

Paired E-minus-C mean rank differences (20 paired truths) are +0.1729 at gamma -1 (95% t interval +0.1107 to +0.2351), +0.0145 at gamma 0 (−0.0212 to +0.0503), and −0.1299 at gamma +1 (−0.1779 to −0.0818). This is a curvature-dependent displacement of the fitted degree-one IPW ranks, not a universal ranking of estimators.

These intervals are exploratory and conditional on this synthetic design. The 20 truth draws are the independent units; the three validation sizes share the same truth and are not extra independent replicates. Mean rank near 0.5 alone is not a complete calibration test, and the small number of truths limits inference. The parent-oracle arm is a simulation diagnostic only.

## Interpretation and limitations

1. The 1,121 → 2,241 geometric-grid comparison is numerically stable under the declared rank, quantile, and density-L1 tolerances.
2. Catalogue-only posterior ranks are concentrated near one in this prior-predictive design.
3. The degree-one MAP and integrated IPW ranks shift with the true inclusion curvature; integrating over fitted coefficients does not remove the omitted-curvature effect in this pilot.
4. The 20-truth pilot is not large enough for confirmatory calibration claims.
5. The audit is conditional on archived KDE densities. It is not an independent re-generation of the KDE pipeline from raw galaxy weights.
6. No observational H0, LVK/GWTC, or real galaxy-catalogue validation claim is made.

## Reproducibility

The reusable audit is `scripts/audit_geometric_prior_predictive_h0_grid.py`. Its eventwise-equivalence regression is in `tests/test_audit_geometric_prior_predictive_h0_grid.py`. The automated audit workflow is `.github/workflows/exact-geometric-prior-predictive-h0-grid-audit.yml`.
