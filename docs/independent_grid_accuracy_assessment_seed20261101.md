# Independent grid-accuracy assessment: seed 20261101 (normalization-checked rerun)

## Executive result

The full 200-replicate grid-refinement audit completed successfully in GitHub Actions and the generated CSV artifacts were downloaded and independently inspected. The evaluator reports **PASS_FOR_GRID_DIFFERENCES_ONLY** for the proposed grid-difference criteria.

- Workflow run: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40091695141
- Seed: 20261101
- Replicates per design cell: 200
- Design cells: 36
- Matched case/weighting units: 7,200
- Posterior fits: 21,600 across 1,001-, 2,001-, and 4,001-point grids
- Cell-summary comparisons: 180 (36 cells × five posterior summaries)
- Failed cell-level 95th-percentile checks: 0
- Case-summary discrepancies above the 0.05 guardrail: 0

The normalization-checked rerun's artifacts were retained in GitHub Actions under grid-refinement-r200-seed20261101. The values below were calculated from the downloaded per_run.csv, case_accuracy_assessment.csv, and cell_accuracy_assessment.csv, not inferred from the workflow badge alone.

## Acceptance criteria applied

1. For each design cell and each of the five summaries, the 95th percentile of the absolute difference between the 1,001- and 4,001-point results must be ≤ 0.01 km s⁻¹ Mpc⁻¹.
2. Any individual difference > 0.05 km s⁻¹ Mpc⁻¹ must be flagged for targeted refinement.
3. Credible-interval endpoints must be finite and ordered.
4. The same simulated observation must be reused across grid resolutions and inference-weighting configurations.

All four checks passed. In addition, the runner now records an independently integrated posterior area on the actual prior support for every fit. The rerun checked 21,600 area values against a 1e-8 absolute tolerance: zero failures, with maximum absolute area error 4.440892098500626e-16.

## Quantitative results

### Largest cell-level 95th percentile

The table reports the largest cell-specific 95th percentile across all 36 design cells for each summary. The named cell is the cell attaining that maximum.

| Summary | Largest cell-level p95 absolute difference | Cell attaining maximum |
|---|---:|---|
| Posterior median | 0.00130534 | H0=52, z2=0.025, noise=0.03, correct weighting |
| Lower 68% endpoint | 0.00220844 | H0=52, z2=0.06, noise=0.03, correct weighting |
| Upper 68% endpoint | 0.00259805 | H0=52, z2=0.06, noise=0.03, correct weighting |
| Lower 95% endpoint | 0.00463673 | H0=52, z2=0.06, noise=0.03, correct weighting |
| Upper 95% endpoint | 0.00454303 | H0=52, z2=0.06, noise=0.03, correct weighting |

All values are below the proposed 0.01 km s⁻¹ Mpc⁻¹ p95 tolerance. The largest cell-level p95 across all 180 cell-summary checks is 0.00463673 km s⁻¹ Mpc⁻¹.

### Largest individual discrepancies

| Summary | Largest observed absolute difference | Case |
|---|---:|---|
| Posterior median | 0.00534972 | H0=88, z2=0.04, noise=0.10, replicate 31, reversed weighting |
| Lower 68% endpoint | 0.00664471 | H0=88, z2=0.04, noise=0.10, replicate 33, reversed weighting |
| Upper 68% endpoint | 0.00925907 | H0=52, z2=0.04, noise=0.03, replicate 160, reversed weighting |
| Lower 95% endpoint | 0.01073677 | H0=88, z2=0.04, noise=0.10, replicate 68, correct weighting |
| Upper 95% endpoint | 0.01164793 | H0=52, z2=0.04, noise=0.03, replicate 45, reversed weighting |

The maximum individual difference across all summaries was 0.01164793 km s⁻¹ Mpc⁻¹, below the 0.05 guardrail. Consequently, this run flagged no case for targeted refinement. This is an empirical result for this seed and design, not proof that no other seed or model configuration can produce a larger discrepancy.

## Interpretation

For this 200-replicate, 36-cell synthetic experiment, the 4,001-point grid meets the protocol's proposed cell-level p95 threshold for all five posterior summaries. The maximum individual discrepancy is also below the case-level flagging threshold. This is evidence of adequate grid resolution under the stated operational criteria for this run.

This result should be reported narrowly. It does **not** establish:
- independent normalization of every posterior density to tolerance 1e-8;
- repeated-sampling coverage or calibration;
- robustness to prior-bound changes or likelihood misspecification;
- realism of the synthetic two-host model relative to LVK observations;
- accuracy of an observational H0 measurement.

The result also does not contradict the earlier targeted case where the coarse-grid upper 95% endpoint shifted by about 0.13. The current audit uses a different seed (20261101) and a different random sample. The earlier outlier remains valid evidence that a different simulated case can exhibit a larger coarse-grid discrepancy; it should not be deleted or overwritten by this more favourable independent run.

## Reproducibility

The audit was run through the repository workflow with seed 20261101 and 200 replicates per cell. The workflow validates row counts, finite numerical values, interval ordering, complete grid triplets, matched observations, summary-statistic ordering, and posterior-area normalization. The retained artifact contains the per-run results, cell summaries, targeted-refinement CSV (header-only because no cases were flagged), and evaluator outputs.

Run the experiment and evaluator locally with:

    python scripts/run_grid_refinement_audit.py --replicates 200 --seed 20261101 --output-dir results/grid-refinement-audit-full
    python scripts/evaluate_grid_accuracy.py --per-run results/grid-refinement-audit-full/per_run.csv --output-dir results/grid-refinement-audit-full/accuracy-assessment

## Final assessment

**Status: PASS_FOR_GRID_DIFFERENCES_ONLY.** All 180 cell-summary p95 checks passed, all 21,600 posterior-area checks passed, and no case exceeded the individual refinement guardrail. The result supports using 4,001 grid points for the specified numerical sensitivity experiment under these criteria. Coverage, prior sensitivity, model misspecification, and astrophysical validity remain separate questions.
