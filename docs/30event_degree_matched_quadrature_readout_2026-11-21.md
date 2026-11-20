# 30-event degree-matched inclusion quadrature: audited readout

**Status:** completed numerical engineering pilot; not repeated-sampling calibration, SBC, estimator-superiority evidence, or observational inference.

**Repository experiment:** PR #141, branch `research/30event-degree-quadrature-convergence-20261121`  
**Workflow:** run `44806809662` (completed successfully)  
**Artifact:** `degree-matched-30event-quadrature-sensitivity` (artifact ID `14669447416`; SHA-256 `d0a8817b0716b2e6168ad4d69ad27b939d1359033fe49e6567b75c75077b77d0`)

## Design actually executed

- Synthetic parent population: 1,200; detected events: 30; fixed (H_0=70).
- H0 prior support [30, 110], 281-point grid.
- Curvature (gamma\in\{-1,0,+1\}); validation frame 500 and nested validation sizes 100 and 500.
- Inclusion-model degrees 1 and 2; Gauss-Hermite orders 5, 7, and 9.
- Five inference arms: catalogue-only, parent oracle, known-probability IPW, MAP plug-in IPW, and coefficient-posterior-integrated IPW.
- Eighteen child runs (3 curvatures × 2 degrees × 3 quadrature orders).

## Execution and artifact audit

The workflow completed successfully. The Python 3.12 job executed the 30-event experiment and its bundle verification; Python 3.10 and 3.11 completed their test jobs, while Python-3.12-only numerical experiment steps were correctly skipped in those two matrix jobs.

All 18 child independent audit reports say PASS. Each of the three curvature bundle manifests says PASS; each contains 114 listed files. I independently rechecked all 342 listed file byte counts and SHA-256 hashes after downloading the artifact, rechecked the 18 child audit statuses and degree/order configuration labels, and confirmed 35 paired-input checks per curvature bundle (105 total). This verifies artifact integrity and pairing, not the scientific adequacy of the model.

## Quadrature convergence

Across the 12 degree × curvature × validation-size comparisons, order 7 versus 9 gave:

- Maximum posterior-density L1 difference: **0.000654**.
- Maximum absolute posterior-median shift: **0.002071 km s⁻¹ Mpc⁻¹**.
- Maximum absolute 2.5th-percentile shift: **0.003339 km s⁻¹ Mpc⁻¹**.
- Maximum absolute 97.5th-percentile shift: **0.000830 km s⁻¹ Mpc⁻¹**.

Thus orders 7 and 9 are very close for this particular realization and grid. The 5-to-7 differences were larger in some n=100 cells (largest L1 difference 0.004187; largest median shift 0.014673), so order 5 is less reassuring at the smaller validation size. These are numerical convergence comparisons, not proof of convergence for all possible datasets or posterior geometries.

## Main H0 readout: integrated-IPW arm E

Posterior-median errors relative to the fixed truth (H_0=70), using order 9:

| Curvature | Validation n | Degree 1 error | Degree 2 error |
|---:|---:|---:|---:|
| -1 | 100 | +3.983 | +6.001 |
| -1 | 500 | +2.450 | +6.563 |
| 0 | 100 | +7.707 | +7.267 |
| 0 | 500 | +6.064 | +8.298 |
| +1 | 100 | +7.552 | +6.376 |
| +1 | 500 | +6.381 | +6.750 |

All values are km s⁻¹ Mpc⁻¹. Degree 2 does not uniformly improve the integrated-IPW posterior median. Its direction relative to degree 1 changes with curvature and validation size. In this same realization, the parent-oracle median error is +1.066 in every cell; known-probability IPW errors are +4.697, +5.333, and +3.281 for curvature -1, 0, and +1 respectively. The oracle and known-probability arms are paired controls, not independent replications.

## Important IPW stability warning

At order 9, degree-2 integrated-IPW arm E has 729 coefficient-quadrature nodes (9 × 9), but its posterior quadrature ESS is only about 41.6–42.6. More importantly, the minimum node-wise host-weight Kish ESS is:

- (gamma=-1): 2.84 (n=100), 4.50 (n=500).
- (gamma=0): 1.24 (n=100), 2.18 (n=500).
- (gamma=+1): 1.57 (n=100), 3.99 (n=500).

The maximum node-wise top-5%-weight mass ranges from 96.1% to 100.0% across these degree-2 cells. These are severe concentration diagnostics: a small number of catalogue objects dominate some node-specific IPW estimates. The quadrature ESS and the host-weight Kish ESS measure different things and must not be conflated.

## Interpretation and next gate

1. The numerical integration itself appears stable between orders 7 and 9 for this fixed dataset.
2. Stable quadrature does **not** imply stable IPW: degree-2 host weights remain extremely concentrated at some coefficient nodes.
3. Degree 2 does not show a consistent H0 advantage over degree 1 in this one paired realization.
4. No coverage, calibration, or estimator-ranking claim is justified from one fixed truth and one set of synthetic inputs.

**Next scientific gate:** retain the current result as a numerical pilot, then run independent prior-predictive truths with an explicit all-attempt lifecycle ledger. Report calibration and paired uncertainty separately from quadrature convergence, and retain node-wise IPW concentration diagnostics. Do not claim degree-2 superiority unless that broader experiment supports it.

---
