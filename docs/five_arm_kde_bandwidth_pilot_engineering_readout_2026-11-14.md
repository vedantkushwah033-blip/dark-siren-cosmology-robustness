# Five-arm KDE optimization and bandwidth pilot: engineering readout

**Date:** 2026-11-14  
**Source PR head:** `b4ec9a16065ae14ab3a8856299e00cbaa00535b9`  
**Workflow:** [run 43550061324](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/43550061324)  
**Artifact:** `five-arm-kde-bandwidth-pilot` (4,587,197 bytes; SHA-256 `d528c0fdb4cf904f8f8d6f33f8f5c633984b1fcc3760f3d0d35aaa63e6de20e6`)

## 1. Exact-head software validation

The exact PR head completed the three-version workflow successfully on Python 3.10, 3.11, and 3.12. The Python 3.12 path also completed the production-dimension kernel benchmark, the 100-replicate parent-density pilot, its raw-data reconstruction audit, and the separate training/holdout selection audit.

The new absolute-bandwidth estimator path is covered by the exact-head test matrix. This is software validation; it is not evidence of posterior calibration or estimator superiority.

## 2. Batched likelihood-kernel benchmark

Configuration: 321 H0 points × 801 redshift points × 30 events × 23 density rows (four fixed arms plus 19 coefficient nodes), three timed repetitions.

| Quantity | Result |
|---|---:|
| Scalar median | 6.07149 s |
| Batched median | 0.31343 s |
| Observed ratio (scalar / batched) | 19.3708× |
| Maximum absolute log-likelihood discrepancy | 1.42 × 10^-13 |
| Maximum absolute log-efficiency discrepancy | 2.66 × 10^-15 |
| Numerical status | PASS_NUMERICAL_EQUIVALENCE |

This benchmark isolates likelihood-kernel evaluation on precomputed densities; it excludes KDE construction. The timing is specific to the CI runner, BLAS, and machine load. It is not a universal speedup claim.

## 3. Parent-density Scott-multiplier pilot

The raw-data auditor independently reconstructed 1,200 replicate-by-candidate KDE/IAE calculations, covering 100 replicates and 400,000 parent objects. The maximum absolute IAE discrepancy was 7.49 × 10^-16. Raw audit status: **PASS**.

The train/holdout audit independently reconstructed the 70/30 split, all candidate scores and uncertainty summaries, source hashes, and selected value. Selection audit status: **PASS**.

However, the numerical result is not a frozen bandwidth choice:

- The training mean IAE decreased at every candidate from 0.03 through 0.40.
- The training minimum was at the largest candidate, 0.40 (training mean IAE 0.04626).
- The holdout mean IAE also decreased at every candidate, reaching 0.04727 at 0.40.
- Thus the minimum lies at the upper boundary of the searched range. The experiment did **not** bracket an interior optimum.

This is a scientifically important negative result about the candidate grid. The original selection artifact's value 0.40 must not be treated as an eligible frozen bandwidth. It is only the best value among the original candidates.

## 4. Corrective gate

The selection code now labels a candidate-grid endpoint minimum as `BOUNDARY_LIMITED_EXPAND_GRID`, sets `selection_eligible=false`, and requires the independent auditor to verify that status. The workflow candidate grid has been expanded to:

`0.03, 0.04, 0.05, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.80, 1.00, 1.20`.

The next exact-head workflow must determine whether this wider grid brackets an interior density-recovery minimum. If the minimum again occurs at an endpoint, the bandwidth remains ineligible and the range/design must be reconsidered before H0 propagation. The holdout remains descriptive and cannot be used to select.

## 5. Scope limits

This pilot uses a single synthetic parent DGP, N=4,000, a 801-point redshift grid, and a parent-sample KDE. It does not compare the five inclusion-correction arms, does not evaluate finite-validation estimation, and does not test H0 posterior calibration. The absolute physical bandwidth sensitivity protocol remains a separate, not-yet-executed stage. No observational or LVK/GWTC conclusion follows.
