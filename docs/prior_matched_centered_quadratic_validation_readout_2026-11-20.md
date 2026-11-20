# Prior-matched centered quadratic validation readout — 2026-11-20

## Scope

Synthetic low-redshift Euclidean methodology only. This readout concerns a single degree-two H0 smoke realization and a three-replicate paired D1/D2 pilot. It is not calibration evidence, an observational H0 result, or evidence of universal degree-two superiority.

## CI and independent audits

On PR #133 head `851b9d7b81c3a4b4d01e1f3ed8edc4dc1b5d595c`, all four triggered workflows completed successfully:

- Python test matrix: run `44689830839`.
- Paired degree pilot: run `44689830956`.
- Quadratic H0 smoke: run `44689830842`.
- Degree-two quadrature sensitivity: run `44689831059`.

The smoke artifact's independent audit passed: 16 manifest entries verified; 1,000 generation attempts and 932 detected candidates reconstructed; 4 selected events; shared-node selected-likelihood reconstruction had zero density difference at quadrature orders 3 and 5; prior-scale sensitivity reconstructions also had zero density difference. The quadrature ladder reported CONVERGED.

The quadrature-sensitivity artifact independently audited all six degree/order children (D1 and D2, orders 5/7/9), with all six child audits passing. The independent audit reconstructed the posterior summaries and reported maximum log-likelihood discrepancy at approximately 1.4e-14. For integrated-IPW posterior density, order 7 to 9 L1 differences were 1.25e-5 (D1, n=100), 5.73e-7 (D1, n=300), 8.37e-6 (D2, n=100), and 7.02e-7 (D2, n=300). These are one realization's quadrature checks, not calibration.

The paired pilot outer audit passed: 90/90 lifecycle cells successful for each degree, 9 paired child blocks, 90 paired summary rows reconstructed, 354 outer and 574 nested manifest entries verified, and common A/B/C likelihood curves invariant to numerical tolerance.

## Paired D1/D2 result (only 3 paired replicates)

For uncertainty-integrated IPW arm E:

| Curvature | Validation n | D1 RMSE | D2 RMSE | D1 68% coverage | D2 68% coverage |
|---:|---:|---:|---:|---:|---:|
| -1 | 100 | 3.56 | 8.86 | 66.7% | 33.3% |
| -1 | 500 | 2.16 | 6.75 | 100% | 33.3% |
| 0 | 100 | 4.35 | 6.90 | 66.7% | 66.7% |
| 0 | 500 | 4.86 | 5.84 | 33.3% | 0% |
| +1 | 100 | 8.94 | 8.17 | 0% | 66.7% |
| +1 | 500 | 8.88 | 6.24 | 0% | 0% |

D2 is better in these +curvature examples, but worse in the negative- and zero-curvature examples. With only three paired replicates, the paired intervals are broad and do not support a model ranking. The result is a pilot, not a degree-two endorsement.

## Weight-instability warning

At quadrature order 9, D2 integrated-IPW minimum node-level Kish ESS was 3.17 for n=100 and 14.47 for n=300. The corresponding maximum node top-5%-weight mass was 92.3% and 62.5%. This is a major stability warning: quadrature convergence of the posterior integral does not imply stable IPW weights or reliable inference.

## Interpretation and next gate

The centered basis `[1, x, x^2-c_g]` now matches the D1 prior marginal on intercept and slope, removing the raw-versus-centered coefficient-prior confound. The remaining key question is whether this better-specified inclusion model can deliver calibrated, stable H0 inference over many independent populations.

Next: a larger paired prior-matched D1/D2 experiment with at least 50–100 replicates per curvature/validation-size cell, quadrature sensitivity, predeclared weight-stability reporting, and prior-predictive rank/coverage checks. Do not interpret the three-replicate pilot as confirmatory.
