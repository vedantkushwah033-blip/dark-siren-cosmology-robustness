# Exact selected-likelihood H0-grid refinement audit — 2026-11-20

## Scope

This audit addresses a specific numerical concern: whether the 561-point uniform H0 grid used in the 50-replicate prior-predictive five-arm pilot resolves posterior ranks and credible intervals adequately.

The audit re-evaluates the selected likelihood at new H0 points, rather than interpolating the archived log-likelihood curves. It holds each archived KDE host-density estimate fixed, recomputes every detected-event numerator and selection-efficiency denominator, and re-integrates Arm E over its 25 shared coefficient nodes. The posterior prior remains Uniform(30,110), and posterior integration is with respect to **dH0**, not d(log H0).

This isolates H0-grid integration error conditional on the archived KDE density estimates. It does not independently regenerate the KDE densities from raw galaxy weights.

## Data and reproducibility

- Source run: [50-replicate workflow 44683807350](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44683807350)
- Source artifact ID: 14613021797
- Source artifact SHA-256: `1d8667e2f7391283237e27b0c6db26cac6e9dfd56d344649b6931ec047745c0c`
- Replicates: 50; curvature strata: -1, 0, +1; validation sizes: 100, 300, 500; arms: A–E.
- Total posterior curves audited: 2,250.
- Exact evaluation grids: 561-point original uniform grid; 1,121-point geometric grid; 2,241-point geometric reference grid.
- New reproducible script: `scripts/audit_prior_predictive_exact_h0_grid.py`.

Example:
```bash
python scripts/audit_prior_predictive_exact_h0_grid.py \
  --artifact-zip five-arm-prior-predictive-50rep.zip \
  --output-dir results/exact-h0-grid-audit
```

The script verifies that its exact 561-point reconstruction reproduces all archived A–E selected log-likelihood curves, then reports per-curve ranks, nine posterior quantiles, posterior-density L1 distances, and refinement summaries.

## Results from the independent numerical execution

All 2,250 curves were recomputed on both 1,121- and 2,241-point geometric grids. The 561-point archived curves were independently reconstructed on their original grid and compared with the 1,121-point results.

| Diagnostic | 561 uniform → 1,121 geometric | 1,121 geometric → 2,241 geometric |
|---|---:|---:|
| Maximum absolute posterior-rank change | 0.00066485 | 0.00006873 |
| Maximum posterior-median change (km/s/Mpc) | 0.003064 | 0.000542 |
| Maximum change among q2.5, q16, q84, q97.5 (km/s/Mpc) | 0.00810 | 0.00113 |
| Maximum posterior-density L1 difference | 0.002035 | 0.0001245 |

For the 1,121→2,241 comparison, the maximum quantile change above covers q2.5, q16, q84 and q97.5; the committed audit script computes and reports all nine quantiles.

### Interpretation

1. The original 561-point uniform grid has very small rank and quantile changes when refined, but its maximum posterior-density L1 difference from the 1,121-point grid is 0.002035, above the strict 0.001 diagnostic threshold. It should not be described as fully grid-converged under that density-shape criterion.
2. The 1,121-point geometric grid is strongly supported by comparison with the 2,241-point geometric reference: maximum rank change 0.00006873, median change 0.000542 km/s/Mpc, maximum change among the four reported quantiles 0.00113 km/s/Mpc, and maximum posterior-density L1 difference 0.0001245.
3. The geometric mesh is only an evaluation mesh. The prior is still uniform in H0, and all posterior normalization/CDF integration uses dH0. No log-H0 prior or d(log H0) integration is introduced.
4. This is an exact selected-likelihood recomputation conditional on the frozen KDE densities, not a complete end-to-end regeneration of all density estimates. The distinction is material and must remain in any paper or abstract.

## Decision

- **Legacy 561-point grid:** fails the strict posterior-density L1 threshold; rank and quantile summaries are nevertheless stable to the reported precision.
- **1,121-point geometric grid:** passes the reported comparison against the 2,241-point geometric reference on rank, median, selected quantiles and posterior-density L1 diagnostics.
- **Recommended analysis grid for the next calibrated five-arm runs:** at least 1,121 geometric H0 points over [30,110], retaining integration with respect to dH0.
- Do not retroactively label the 561-point calibration results fully converged. Re-run the final prior-predictive calibration on the refined grid before making confirmatory calibration claims.

The result remains a synthetic-methodology benchmark. It is not an observational H0 measurement or validation on real galaxy catalogues.
