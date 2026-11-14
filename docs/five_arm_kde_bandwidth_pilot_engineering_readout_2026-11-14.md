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


## 6. Relationship to the 200-replicate primary source experiment

The 100-replicate run in this PR is **not an independent replication** of the 200-replicate source experiment inherited from PR #44. Both use master seed `20261114`, the same parent-population DGP, and replicate IDs beginning at zero. The 100-run sample therefore overlaps the first 100 replicate IDs of the 200-run experiment. Its different candidate grid and train/holdout split make it a grid-sensitivity analysis on overlapping source replicates, not independent confirmation.

The 200-replicate experiment remains the primary source artifact for the original nine-candidate grid. Its independent raw-data audit reconstructed 1,800 replicate-by-candidate calculations with maximum absolute IAE discrepancy `4.16e-16`; the train/holdout audit passed with 140 training and 60 holdout replicates. Its training and holdout IAE both decreased monotonically through the largest candidate, `0.40` (training mean `0.04674`; holdout mean `0.04660`). This is also an unresolved upper-boundary minimum, not an interior optimum.

The expanded 100-replicate sensitivity workflow now searches a wider candidate grid through `1.20`. Because its replicate IDs overlap the primary 200-run source, it can diagnose whether the apparent boundary trend persists over a wider range, but it cannot be described as independent confirmation. Its expanded-grid output must be generated and audited at the current exact PR head before reporting any numerical result. If it still selects an endpoint, no bandwidth is eligible to freeze.

## 7. Resource and memory interpretation

The batched likelihood refactor bounds the **temporary event/H0/redshift kernel** by processing H0 chunks. It still retains full output arrays for log numerators and detection efficiencies with dimensions proportional to density count × event count × full H0-grid size. Therefore the implementation does not claim memory use independent of the final H0-grid resolution. The benchmark reports likelihood-kernel evaluation time only and excludes KDE construction.

## 8. Current-head validation boundary

The previously recorded successful run `43550061324` belongs to source head `b4ec9a16065ae14ab3a8856299e00cbaa00535b9`, not to later edits that expanded the candidate grid, tightened the selection audit, or added H0 chunk-boundary regression tests. Those later changes require a new exact-head workflow run. The old benchmark and audit results remain valid for the old artifact and configuration only; they must not be relabelled as validation of the current head.


## 6. Expanded candidate-grid result (workflow #43554565338)

The exact-head workflow for commit `2deee87dd7ccaad2b0c4e95c65c7825d0149a038` completed successfully on Python 3.10, 3.11 and 3.12. The expanded 17-candidate bandwidth source was generated and both independent auditors passed.

- Raw-data reconstruction: 1,700 replicate-by-candidate KDE/IAE calculations; 100 replicates; 400,000 parent objects; maximum absolute IAE discrepancy (7.49\times10^{-16}).
- Training/holdout split: 70/30, or 70 training and 30 holdout replicates.
- Selected training multiplier: **0.50**, now an interior candidate.
- Training mean IAE at 0.50: **0.0443212**.
- Holdout mean IAE at 0.50: **0.0439975** (MCSE **0.0016767**).
- Nearest competing candidate 0.60: training mean IAE **0.0453722**; holdout mean IAE **0.0441885**.
- Selection auditor: PASS; `boundary_limited=false`; `selection_eligible=true`.

The holdout result is consistent with the training-selected 0.50 candidate, but the holdout was not used to select it. This is a successful density-only selection pilot for the single parent DGP, N=4,000, grid=801, and this candidate set. It is not yet a common bandwidth selected across all five arms, catalogue-inclusion scenarios, or validation sizes. It must not be automatically substituted as a common Scott multiplier in the five-arm experiment.

The exact-head workflow also completed the production-like batched likelihood benchmark: observed ratio **19.531×**, scalar median **8.58972 s**, batched median **0.43980 s**, maximum log-likelihood discrepancy (1.42\times10^{-13}), and maximum log-efficiency discrepancy (2.66\times10^{-15}). The benchmark is runner-specific and excludes KDE construction.

## 7. Interpretation and next gate

The expanded Scott grid now has an interior density-recovery minimum. This closes the narrow *parent-only Scott-multiplier pilot* gate, not the broader five-arm bandwidth gate. Before any confirmatory five-arm H0 run, the project still needs:

1. common absolute-bandwidth density recovery across A–E, all E coefficient nodes, validation sizes, and preregistered inclusion DGP cells;
2. independent raw-input reconstruction and failure-inclusive accounting for that experiment;
3. a disjoint density-only training/holdout design for selecting a common physical bandwidth, with candidate-grid boundary checks;
4. downstream H0 propagation only after the density stage is frozen and audited.
