# H0 lower-prior-support sensitivity readout (2026-11-19)

**Workflow:** [44481137736](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44481137736) — SUCCESS  
**Tests:** support-sensitivity unit tests passed in the workflow.  
**Source:** immutable 1,000-replicate full-parent SBC artifact from [run 44446994455](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44446994455).  
**Adaptive reference:** 4,481-point coarse knots with adaptive quadrature; source-seed regeneration and four manifest entries verified.  
**Artifact:** `h0-prior-support-sensitivity-20261119`, ID `14509491021`, SHA-256 `5f949a2a636466deaa2d606193559d6fabca072bfa8b33414a6afd7ffbe9dafe`.

## Question

Do the ten lowest-`q975` posterior draws from the original synthetic SBC run owe their low-H0 posterior concentration to the original prior lower bound of 20 km s⁻¹ Mpc⁻¹?

## Design

For each of the ten deterministically selected replicates, the audit independently regenerated the archived H0 truth, 1,200-galaxy parent population, and 30 detected events from the immutable seed schedule. It then held the generated data and likelihood fixed and recomputed adaptive posteriors under three uniform supports:

- Original: [20, 300] km s⁻¹ Mpc⁻¹.
- Expanded: [15, 300] km s⁻¹ Mpc⁻¹.
- Expanded: [10, 300] km s⁻¹ Mpc⁻¹.

The expanded-support posterior is a conditional sensitivity analysis. The archived truths were generated from Uniform(20,300), so the expanded-support posterior is **not** the posterior under the original generating prior. Its ranks must not be mixed into the original SBC rank sample.

## Aggregate results

Across the 20 expanded-support comparisons (10 replicates × 2 alternative lower bounds):

| Diagnostic | Maximum absolute value |
|---|---:|
| Posterior median shift | 0.4687 km s⁻¹ Mpc⁻¹ |
| Posterior 2.5% quantile shift | 3.0060 km s⁻¹ Mpc⁻¹ |
| Posterior 97.5% quantile shift | 0.0785 km s⁻¹ Mpc⁻¹ |
| Posterior mass below original H0=20 bound | 31.19% |
| Adaptive integrator reported relative error estimate | 4.24 × 10⁻¹¹ |

## Interpretation

The original lower support can materially affect the lower posterior tail in at least some of these selected low-H0 cases: when support is expanded, as much as 31.2% posterior mass lies below the original lower bound, and the 2.5% quantile can move by about 3.0 km s⁻¹ Mpc⁻¹. The largest median shift is smaller (about 0.47 km s⁻¹ Mpc⁻¹), while the upper 97.5% quantile changes little.

This is evidence of **prior-support sensitivity in the lower tail**, not evidence that all ten cases are equally affected, nor a correction to the original SBC calibration. The summary maxima do not identify how many cases exceed any particular sensitivity threshold; per-replicate rows are preserved in the workflow artifact.

## Limitations and next gate

1. Only the ten lowest-`q975` replicates were studied; this is a targeted tail audit, not a full-sample sensitivity analysis.
2. The alternative supports change the prior and are not calibrated against the archived Uniform(20,300) truth-generating process.
3. The adaptive integrator's reported error is an estimate, not a formal mathematical bound.
4. The original 1,000-replicate SBC rank/coverage summary is not recomputed here.
5. A full-sample support-sensitivity study must use a newly matched prior-predictive design if calibration under the expanded support is to be assessed.

**Conclusion:** the lower H0 support is a material modelling choice for the selected low-tail posterior cases. Report the original SBC results as conditional on Uniform(20,300), and do not describe them as prior-support-robust.
