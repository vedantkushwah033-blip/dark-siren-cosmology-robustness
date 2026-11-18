# Independent-validation selected-H0 smoke rerun (2026-11-18)

## Status

**Workflow: PASS. Independent artifact audit: PASS. Quadrature-order stability: CONVERGED.**

This is a software/integration smoke for the existing modular inclusion-posterior selected-H0 prototype. It is not a calibration experiment and does not yet use the quadratic-logit curvature DGP from the seven-gamma density study.

- Workflow run: `44395542339`
- Workflow commit: `76be16935cf9fec71fc59b6512292240c24755a1`
- Artifact ID: `14470390389`
- Archive SHA-256: `c17feb3494f7feb57c7944cd7f8ee77940d83def9b807a58f94a00ee0f7a05a3`
- GitHub artifact digest: `sha256:c17feb3494f7feb57c7944cd7f8ee77940d83def9b807a58f94a00ee0f7a05a3`
- Independently downloaded archive SHA-256 matches the GitHub digest.

## Design actually run

- Analysis parent: 1,200; observed analysis catalogue: 633.
- Independent validation frame: 1,000; SRS validation sample: 150.
- Separate seeds: 2026111901 (analysis), 2026111902 (validation).
- Detected events: 4; generator attempts: 1,000; detected candidates: 936; rejected candidates: 64.
- H0 grid: 81 points over [50, 90], uniform H0 prior.
- Degree-one logistic inclusion model; proper independent Gaussian coefficient prior with SD [5, 2.5].
- Adaptive Gauss-Hermite quadrature orders requested: 5, 7, 9, 11. The convergence rule stopped after order 9 passed two consecutive transitions.
- Shared coefficient node across all events. The node-specific detection efficiency is computed and the selection correction is applied inside each node before marginalization.
- Prior slope-scale sensitivity multipliers: 0.5, 1, 2.

## Independent audit

The separate raw-artifact auditor reports:

- `status=PASS`
- 18 manifest entries verified.
- Analysis parent/catalogue counts: 1,200 / 633.
- Validation frame/SRS counts: 1,000 / 150.
- Generation attempt, detected-candidate, rejected-candidate and selected-event counts independently reconstructed.
- Selected-likelihood posterior normalization checks: 3 (orders 5, 7, 9).
- Maximum absolute posterior-density reconstruction differences: exactly 0 for orders 5, 7, 9.
- Prior-scale sensitivity posterior-density reconstruction differences: exactly 0 for multipliers 0.5, 1, 2.
- Quadrature status: `CONVERGED`.

The H0 posterior median at the final evaluated order 9 was 69.3213 km s^-1 Mpc^-1; the central 68% interval was approximately [58.4775, 80.6250]. This describes one deliberately small synthetic smoke realization only. It is not an estimate of the real Hubble constant and is not evidence of unbiasedness or calibration.

## Numerical convergence

- Order 5 → 7: max absolute density difference 3.0752e-6; max quantile difference 0.0011374.
- Order 7 → 9: max absolute density difference 1.3414e-7; max quantile difference 5.5722e-5.
- Both transitions passed the configured tolerances; two consecutive passes were achieved at order 9.

## Critical scientific limitation

The generator's catalogue-inclusion law is a sigmoid-with-floor function:
`floor + (1-floor)/(1+exp((z-midpoint)/width))`.
The floor was set to zero for this smoke. This is not the quadratic-logit inclusion DGP used by the corrected seven-gamma density grid. Thus the run validates the existing H0 propagation software path under its legacy synthetic generator, but it does **not** connect the current curvature-misspecification experiment to H0 posterior behavior.

The next scientific implementation gate is a matched quadratic-logit independent-validation H0 generator and auditor, using the same curvature settings, calibration targets, and paired seeds as the density grid. Then run multiple independent validation/GW replicates and separately assess fixed-truth coverage, PIT and SBC under clearly stated generative assumptions.

No LVK/GWTC data, real survey validation, observational H0 measurement, or universal estimator-superiority claim is supported.
