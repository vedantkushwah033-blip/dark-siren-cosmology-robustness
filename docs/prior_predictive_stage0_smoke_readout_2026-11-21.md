# Prior-predictive Stage 0 smoke: independent readout (2026-11-21)

## Scope and claim limit

This is a ten-truth-per-design-cell prior-predictive engineering smoke for the synthetic low-redshift Euclidean dark-siren model. It verifies the five-stream simulation schedule, cross-degree pairing, child runs, and independent aggregate reconstruction. It is **not** a confirmatory calibration study, an SBC validation, an observational H0 measurement, an LVK/GWTC analysis, or evidence that degree two is generally superior.

## Frozen design and provenance

- Workflow: `44841068499` (completed successfully).
- Artifact: `prior-predictive-stage0-ten-replicate-smoke`, artifact ID `14684946325`.
- Master seed: `20261121`; ten replicate-specific H0 truths sampled from Uniform[30,110].
- Curvature: gamma in {-1, 0, +1}; validation sizes n={100,300,500}.
- Parent population: N=1200; 30 detected events; validation frame=500.
- H0 evaluation grid: 281 points on [30,110]; inclusion-coefficient quadrature order 5.
- Five independent streams: parent population (0), validation-frame generation (1), H0 truth (2), GW event generation (3), and within-frame validation sampling (4).
- Degree one and degree two used identical replicate-specific generated inputs and H0 truths. The separate cross-degree pairing audit completed with PASS.

## Independent audit

For each inclusion degree:

- 450/450 planned arm cells were SUCCESS (10 replicates × 3 curvatures × 3 validation sizes × 5 arms).
- 30/30 child runs received independent audit PASS.
- All 578 aggregate manifest entries were verified.
- All 45 aggregate metric rows and 720 paired-contrast rows were independently reconstructed.
- The independent auditor's status was PASS.

The cross-degree audit independently verified shared inputs and truth/seed pairing. These are lifecycle and reproducibility checks, not evidence of posterior calibration.

## Exploratory performance snapshot

The following values pool the 10 truths across the three curvature settings and three validation sizes (90 rows per arm). Because the same ten replicate truths are reused across gamma and validation-size cells, these 90 rows are **not 90 independent truths**. Use the cell-level table below for the design-specific view.

| Arm | Degree 1 RMSE | Degree 1 68% coverage | Degree 2 RMSE | Degree 2 68% coverage |
|---|---:|---:|---:|---:|
| Catalogue-only | 17.857 | 0.000 | 17.857 | 0.000 |
| Full-parent oracle | 2.835 | 0.800 | 2.835 | 0.800 |
| Known-probability IPW | 3.047 | 0.667 | 3.047 | 0.667 |
| Estimated MAP IPW | 4.356 | 0.589 | 3.600 | 0.678 |
| Estimated integrated IPW | 3.966 | 0.689 | 3.258 | 0.822 |

The pooled summaries suggest degree two can improve the estimated-IPW arms in this particular smoke. They do not establish a stable gain: only ten independent truths were drawn, and the cell-specific behavior varies with curvature and validation size.

### Integrated-IPW cell-level RMSE

| gamma | validation n | Degree 1 E RMSE | Degree 2 E RMSE | Degree 1 E 68% coverage | Degree 2 E 68% coverage |
|---:|---:|---:|---:|---:|---:|
| -1 | 100 | 5.34 | 4.39 | 0.60 | 0.80 |
| -1 | 300 | 5.03 | 3.71 | 0.80 | 0.80 |
| -1 | 500 | 4.52 | 3.02 | 0.80 | 0.70 |
| 0 | 100 | 3.71 | 4.01 | 0.70 | 0.80 |
| 0 | 300 | 3.21 | 2.90 | 0.70 | 0.90 |
| 0 | 500 | 2.94 | 2.41 | 0.60 | 0.90 |
| +1 | 100 | 2.96 | 3.26 | 0.70 | 0.70 |
| +1 | 300 | 3.52 | 2.78 | 0.70 | 0.90 |
| +1 | 500 | 3.67 | 2.18 | 0.60 | 0.90 |

## Critical stability warning

Degree-two integrated IPW has severe node-wise weight concentration, especially with n=100:

| gamma | validation n | Minimum node-wise Kish ESS | Maximum node-wise top-5% mass |
|---:|---:|---:|---:|
| -1 | 100 | 1.990 | 0.995 |
| -1 | 300 | 3.061 | 0.992 |
| -1 | 500 | 4.453 | 0.984 |
| 0 | 100 | 1.859 | 0.996 |
| 0 | 300 | 5.376 | 0.939 |
| 0 | 500 | 6.080 | 0.916 |
| +1 | 100 | 3.004 | 0.995 |
| +1 | 300 | 7.947 | 0.900 |
| +1 | 500 | 9.261 | 0.835 |

For degree one, the minimum node-wise ESS across these cells was 18.223 and the maximum top-5% mass was 0.847. The degree-two concentration is therefore a material warning: a lower RMSE in a small smoke does not establish a reliable correction when a few hosts dominate the weights. These are node-wise host-weight diagnostics, not coefficient-posterior quadrature ESS.

## Interpretation

1. The five-stream generator and auditor now pass a complete paired Stage 0 smoke across both inclusion degrees.
2. The degree-two integrated-IPW arm has better pooled RMSE/coverage than degree one in this particular ten-truth smoke, but cell-level gains are not uniform.
3. The degree-two IPW concentration warning is severe and must be treated as a central limitation, not hidden by pooled performance.
4. The full-parent oracle and known-probability IPW are useful synthetic reference arms; neither validates estimated inclusion probabilities in real catalogues.
5. The smoke's 10 independent H0 truths per curvature/validation-size cell are insufficient for confirmatory coverage, SBC, multiplicity-aware comparisons, or robust ranking.

## Next gate

- Keep the five-stream implementation and independent seed audit.
- Complete and validate the lifecycle ledger with explicit SUCCESS / NON_EVALUABLE / FAILED / NOT_RUN distinctions before confirmatory runs.
- Freeze the estimands, failure-inclusive summaries, and multiplicity plan.
- Run at least 100 independent prior-predictive truths per curvature × validation-size cell as an interim calibration stage, then expand toward 300.
- Increase H0-grid and quadrature/KDE convergence settings separately; this smoke's 281-point grid and order-5 quadrature are engineering settings, not final numerical-convergence evidence.
- Preserve weight-concentration diagnostics as first-class outputs. Do not declare a degree winner from this pilot.
