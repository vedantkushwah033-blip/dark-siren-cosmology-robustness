# Five-arm prior-predictive SBC readout — 50 replicates — 2026-11-20

## Audit status and design

Workflow run: `44683807350`; artifact: `14613021797`; ZIP SHA-256: `1d8667e2f7391283237e27b0c6db26cac6e9dfd56d344649b6931ec047745c0c`.

The workflow completed successfully. The independent aggregate audit reports PASS, verifies 2,858 aggregate manifest entries, reconstructs 2,250 lifecycle cells (all SUCCESS), reconstructs 45 aggregate metric rows and 3,600 paired contrast rows, verifies 150 cross-curvature paired-validation blocks, and confirms 50 unique prior-predictive truths. All 150 child runs were audited.

The experiment uses 50 prior-predictive H0 truths, shared across curvature conditions within replicate, 30 detected events, curvature values -1/0/+1, validation sizes 100/300/500, and five arms. This is a finite 50-replicate synthetic SBC exercise, not a high-precision calibration study.

## Selected results at validation size 500

| Curvature | Arm | Mean rank | 68% coverage | 90% coverage | 95% coverage |
|---:|---|---:|---:|---:|---:|
| -1 | Parent oracle | 0.481 | 70% | 86% | 92% |
| -1 | Known-probability IPW | 0.525 | 42% | 74% | 80% |
| -1 | MAP plug-in IPW | 0.710 | 44% | 66% | 72% |
| -1 | Integrated IPW | 0.700 | 52% | 66% | 76% |
| 0 | Parent oracle | 0.481 | 70% | 86% | 92% |
| 0 | Known-probability IPW | 0.482 | 54% | 88% | 94% |
| 0 | MAP plug-in IPW | 0.507 | 50% | 78% | 92% |
| 0 | Integrated IPW | 0.501 | 52% | 82% | 96% |
| +1 | Parent oracle | 0.481 | 70% | 86% | 92% |
| +1 | Known-probability IPW | 0.465 | 62% | 84% | 92% |
| +1 | MAP plug-in IPW | 0.303 | 50% | 76% | 82% |
| +1 | Integrated IPW | 0.306 | 52% | 82% | 86% |

The catalogue-only arm is strongly miscalibrated in all three curvature conditions (mean ranks 0.98–1.00), as expected when catalogue incompleteness is not corrected. The parent-oracle arm is the closest reference in this table, although its 50-replicate coverage estimates still have material Monte Carlo uncertainty.

For the estimated degree-one arms, ranks are shifted high at gamma=-1 and low at gamma=+1. Integrating over the fitted coefficient posterior does not remove this curvature-dependent rank distortion. At gamma=0, integrated IPW is closer to uniform rank behavior, but coverage remains variable.

## Limits

- 50 replicates gives only coarse rank/coverage resolution; Wilson intervals are wide.
- This is a prior-predictive result under the implemented synthetic model, not fixed-truth calibration or real-data validation.
- The degree-one inclusion model is misspecified at nonzero quadratic curvature.
- These results do not establish that integrated IPW is superior to known-probability IPW.
- The independent audit verifies lifecycle, provenance, pairing and aggregate arithmetic; the calibration interpretation remains conditional on the simulation design being the intended generative model.

## Next gate

The key next experiment is a prior-matched centered degree-two prior-predictive comparison, with the same prior-predictive truths and events paired across D1/D2, larger replicate counts, quadrature checks, and explicit node-level weight stability diagnostics. The centered quadratic-basis correction is now available on the stacked PR #133 branch, but its three-replicate pilot is not enough to claim improved calibration.
