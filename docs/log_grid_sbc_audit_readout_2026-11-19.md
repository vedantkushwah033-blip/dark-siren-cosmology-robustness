# Geometric H0-grid audit: full 200-replicate comparison (2026-11-19)

## Status

**Promising rank/quantile diagnostic; not a complete numerical-convergence pass.**

This readout records an independent local reconstruction of all 200 replicates from the immutable full-parent SBC source artifact, using a 2,241-node geometric grid in H0. It compares the resulting ranks and medians against the existing 4,481-node uniform-H0 audit. Four difficult lower-boundary replicates were also checked against the adaptive-CDF reference, and ten lower-tail replicates from the separate 1,000-replicate source were checked against its existing adaptive-CDF reference.

No prior support, likelihood, event-generation rule, selection normalization, or numerical tolerance was changed.

## Why test a geometric grid?

The prior is Uniform(H0_min, H0_max) with bounds [20, 300]. The observed distance uncertainty is fixed in Mpc, so the posterior can become very narrow in H0 near the lower prior boundary. A uniform H0 grid spends the same number of nodes per Mpc everywhere.

A geometric grid allocates more nodes near the lower boundary. Importantly, this implementation evaluates and integrates the posterior directly in H0:

- Nodes: `np.geomspace(20, 300, 2241)`.
- Integration measure: `dH0`, using nonuniform trapezoidal quadrature.
- Prior: unchanged Uniform(20, 300).

This is **not** a change of integration variable to t = log(H0). Therefore no log-coordinate Jacobian is applied. If the likelihood were instead integrated in t-space, the density would need the Jacobian factor exp(t).

## Full 200-replicate comparison

Source design:

- Source workflow: `44440150116`.
- Source artifact: `14490552961`.
- Source artifact SHA-256: `a26fa4cc1af2c568ea0d309a606eced37fc4288bb90648f4132bf8780e1ad1fc`.
- Master seed: `2026111907`.
- 200 prior-predictive replicates.
- 1,200 parent galaxies and 30 detected events per replicate.
- H0 prior Uniform(20, 300).
- Source grid: 561 points.

Comparison grid reference:

- Workflow: `44462233294`.
- Artifact: `14502499332`.
- Artifact SHA-256: `55429d28be1d6246f953f49a286810ee2e1f3d006632dff698c8f8759de6c78d`.
- Uniform grid: 4,481 points.
- This audit itself failed its original density-L1 and rank gates; it is a high-resolution comparator, not an exact reference.

### Results: 2,241 geometric nodes versus 4,481 uniform nodes

| Diagnostic | Maximum absolute difference | Median absolute difference |
|---|---:|---:|
| Posterior rank | 0.0005539 | 0.00000538 |
| Posterior median (km/s/Mpc) | 0.0014740 | 0.00006491 |

All 200 replicates were reconstructed. The four replicates with truth H0 below 25 had the largest discrepancies:

- Maximum rank difference: 0.0005539.
- Maximum median difference: 0.0014740.

For the other 196 replicates:

- Maximum rank difference: 0.0000938.
- Maximum median difference: 0.0009053.

The worst rank comparison was replicate 100 (truth H0 = 20.1273). Replicates 156 and 174 were the next most difficult.

The 561-point source grid versus the 2,241-point geometric grid had maximum rank difference 0.10816 and maximum median difference 0.09109. This confirms that the geometric grid materially changes the previously under-resolved lower-boundary results.

## Independent adaptive-CDF checks: four difficult replicates

The adaptive reference independently evaluates the scalar selected-event likelihood and integrates its CDF using adaptive quadrature.

| Replicate | Adaptive rank | Geometric-grid rank | Absolute rank difference | Maximum 11-quantile difference |
|---:|---:|---:|---:|---:|
| 100 | 0.3127145 | 0.3123934 | 0.0003211 | 0.0006596 |
| 174 | 0.1493825 | 0.1494075 | 0.0000250 | 0.0003853 |
| 156 | 0.2041507 | 0.2042046 | 0.0000538 | 0.0004607 |
| 62 | 0.8433806 | 0.8433413 | 0.0000393 | 0.0005056 |

For these four cases, the largest approximate density-L1 difference between the geometric-grid normalized density and the adaptive-normalized density, evaluated on the geometric grid and integrated by the trapezoid rule, was 0.0001789.

That last number is a **targeted grid-based approximation**, not a full continuous adaptive density-L1 certification.

## Independent adaptive-CDF checks: ten lower-tail replicates from the 1,000-replicate source

Source:

- Workflow: `44446994455`.
- Artifact: `14497517034`.
- Source ZIP SHA-256: `7b16872f339cea51f6c8e566f073aa1e2b9d6b4d6121f8820d6259c3cc1812ca`.

Adaptive diagnostic artifact:

- Artifact: `14502130062`.
- ZIP SHA-256: `8adb43b25731574378257c780fe746aa5138cab7d092395c8413d7ae4a92e8ac`.

The ten lowest archived q97.5% cases were reconstructed on the 2,241-node geometric grid and compared with the adaptive-CDF reference:

- Maximum rank difference: 0.0001333.
- Maximum median difference: 0.0000807 km/s/Mpc.
- Maximum difference over all 11 reported quantiles: 0.0007178 km/s/Mpc.

This is a selected ten-replicate tail diagnostic, not a full 1,000-replicate recalculation.

## Reproducibility and implementation

A reusable runner has been added:

`scripts/audit_full_parent_h0_sbc_log_grid.py`

It:

1. Verifies every file in the source artifact manifest.
2. Verifies the complete 200-replicate source design.
3. Reconstructs truth and simulation seeds from the declared SeedSequence schedule.
4. Regenerates parent populations and detected events.
5. Evaluates the selected-event full-parent likelihood on a geometric H0 grid.
6. Integrates with respect to H0, preserving the uniform-H0 prior.
7. Compares ranks and medians with the 4,481-point uniform-grid audit.
8. Writes a per-replicate CSV and summary JSON.

Example:

```bash
python scripts/audit_full_parent_h0_sbc_log_grid.py \
  --source-dir PATH_TO_EXTRACTED_200_REPLICATE_SOURCE \
  --baseline-csv PATH_TO_4481_GRID_AUDIT/sbc_grid_audit_replicates.csv \
  --output-dir /tmp/log-grid-audit \
  --grid-points 2241
```

Tests were added for:

- Posterior normalization on a geometric H0 grid.
- Quantile accuracy against a truncated Gaussian reference.
- Correct integration with respect to H0 rather than log(H0).

Local independent reconstruction environment: Python 3.13, NumPy 2.3.5, pandas 2.2.3, SciPy 1.17.0.

## What this does and does not establish

### Established by this diagnostic

- The original 561-point uniform grid materially under-resolves some low-H0 boundary posteriors.
- A 2,241-point geometric H0 grid agrees closely with the 4,481-point uniform grid for ranks and medians across all 200 source replicates.
- The geometric grid agrees closely with adaptive-CDF ranks and quantiles in the four difficult 200-replicate cases and ten selected 1,000-replicate tail cases.
- The original prior and likelihood can be retained while changing only the placement of quadrature nodes.

### Still unresolved

1. Full 200-replicate comparison of all 11 quantiles against independent adaptive CDF references.
2. Full-sample continuous density-L1 comparison against an independently represented adaptive density.
3. Full 1,000-replicate recalculation using the geometric grid.
4. Whether the new grid should replace the existing production SBC grid after CI and independent review.
5. The original 4,481-point uniform audit still fails density-L1 (maximum 0.006788 versus 0.001) and narrowly fails rank (maximum 0.002046 versus 0.002). These failures have not been waived.

**Do not label this a full numerical-convergence PASS.** It is a strong, reproducible diagnostic that justifies a complete geometric-grid audit with the existing tolerances preserved.

## Scope

This is a synthetic full-parent prior-predictive numerical audit. It is not catalogue-inclusion SBC, not a real galaxy-survey validation, and not an observational H0 measurement.
