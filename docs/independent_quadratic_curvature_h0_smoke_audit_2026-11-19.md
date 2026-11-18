# Independent audit of the curvature-matched H0 smoke (2026-11-19)

## Scope and provenance

This note records an independent post-run audit of the synthetic artifact produced by the dedicated GitHub Actions workflow.

- Workflow run: [44401970010](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44401970010)
- Exact tested head: `3202a0ba8d119cbefac4ff5d1302d0a6ef44f70d`
- Artifact ID: `14473762897`
- Artifact name: `quadratic-curvature-h0-smoke-3202a0ba8d119cbefac4ff5d1302d0a6ef44f70d`
- Downloaded archive: `quadratic-curvature-h0-smoke.zip`
- Archive size: 64,160 bytes
- Archive SHA-256: `7edfd5470e980c55fc819bd919412ba956bcd78e1ee1f02a8994407a4de31ad8`

The exact-head workflow completed successfully. Its steps show that the complete pytest suite passed, the curvature-matched H0 smoke ran, and the observable/truth-separated artifact upload succeeded. This establishes execution and test success for this synthetic smoke; it does not establish calibration or estimator superiority.

## Independent checks performed

The downloaded ZIP was inspected locally rather than relying only on the workflow's success badge.

1. **Archive integrity:** all 8 files listed in `manifest.json` independently matched both recorded byte counts and SHA-256 hashes.
2. **Observable/truth separation:** catalogue observations contain only `object_id, redshift`; validation observations contain only `validation_id, redshift, retained`. Inclusion probabilities and the complete validation frame remain in truth-only files.
3. **Catalogue reconstruction:** all 589 observed catalogue rows exactly correspond to the retained subset of the 1,200-row parent truth table; IDs and redshifts match.
4. **Validation reconstruction:** all 300 observed validation rows match selected records in the 1,000-row validation truth frame, including retention labels and redshifts.
5. **DGP probability check:** independently evaluating the declared quadratic-logit model reproduced parent and validation inclusion probabilities to floating-point precision (maximum absolute discrepancies below (6\times10^{-16})).
6. **Seed replay:** independent replays of the population, validation, and event-generation random streams reproduced the retained labels, validation selection order, event host IDs, and observed event distances. Redshift arrays agreed to floating-point roundoff (maximum discrepancy below (10^{-16})).
7. **Event selection:** all four observed event distances satisfy the declared one-sided threshold of 250 Mpc. Their truth-only host IDs and redshifts match the parent population.
8. **Posterior normalization:** independent trapezoidal integration of all six stored H0 curves gave unit area to numerical precision. Recomputed 5%, 50%, and 95% quantiles agree with the summary values to the displayed precision.

## Run configuration

- Parent population: 1,200 galaxies.
- Independent validation frame: 1,000 galaxies; validation SRS size: 300.
- Quadratic inclusion curvature: (gamma=1.0).
- Logistic slope: (eta=-3.9).
- Parent-average inclusion target: 0.5101919962.
- Calibrated intercept: 0.1253055057.
- Detected events: 4, generated from the full parent population with the declared merger-weight and one-sided distance threshold.
- Truth: (H_0=70) km s(^{-1}) Mpc(^{-1}).
- Distance noise: 140 Mpc; threshold: 250 Mpc.
- H0 grid: 55–85 km s(^{-1}) Mpc(^{-1}), 61 points.
- Inclusion-coefficient quadrature orders: 5 and 7.

## Numerical results

| Arm | H0 5% | H0 median | H0 95% |
|---|---:|---:|---:|
| Catalogue only | 56.2105 | 67.2235 | 82.3342 |
| Degree-one estimated inclusion | 60.5424 | 77.0215 | 84.2846 |
| Degree-two estimated inclusion | 60.0396 | 76.5220 | 84.2286 |
| Oracle true inclusion (diagnostic only) | 59.0267 | 75.3051 | 84.0724 |

Downstream H0-posterior quadrature refinement from order 5 to 7 was small:

- Degree one: density L1 change 3.82e-5; maximum 5/50/95% quantile shift 5.49e-4.
- Degree two: density L1 change 3.24e-5; maximum 5/50/95% quantile shift 4.05e-4.

These are numerical convergence checks for this realization, not evidence of frequentist coverage or simulation-based calibration.

## Interpretation limits

- This is one synthetic realization with only four detected events. The difference between catalogue-only and corrected posterior medians is not an estimator-performance result.
- The degree-one and degree-two posteriors are conditional on their respective logistic models and priors. The modest difference between their medians is not a model-selection result.
- The oracle arm is available only because this is a simulation; it is a diagnostic reference, not an operational estimator.
- The toy likelihood uses common event noise and a one-sided threshold. It is not an LVK/GWTC analysis, a real galaxy survey analysis, or an observational H0 measurement.
- The artifact does not establish SBC, PIT calibration, bias, RMSE, or nominal interval coverage. Those require a predeclared multi-replicate experiment and independent raw-record audit.
- The artifact's truth-only outputs permit this audit but must not be passed into non-oracle inference.

## Decision

**PASS for artifact integrity, deterministic data-generation replay, observable/truth separation, posterior normalization, and exact-head smoke execution.**

**NOT EVALUATED for calibration, coverage, estimator superiority, or external/observational validity.**

The next scientific gate is a paired multi-replicate H0 experiment comparing degree-one, degree-two, catalogue-only, and oracle diagnostic arms under fixed seeds and declared validation sizes. It must record all failures/non-evaluable runs and independently audit raw records before summarizing bias, RMSE, PIT, and coverage.
