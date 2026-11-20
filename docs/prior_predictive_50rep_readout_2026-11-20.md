# Prior-predictive five-arm SBC 50-replicate pilot readout (2026-11-20)

## Execution and artifact integrity

- Workflow: [44683807350](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44683807350), completed successfully.
- Design: 50 independent prior-predictive truths, H0 ~ Uniform(30,110), 30 detected events, N_parent=1200, gammas {-1,0,+1}, validation sizes {100,300,500}, five paired inference arms.
- Lifecycle: 2,250/2,250 planned arm-cells SUCCESS; 2,250/2,250 audit_status PASS.
- Child audits: 150/150 successful.
- Independent aggregate auditor: PASS; 2,858 aggregate manifest entries verified, 45 aggregate metric rows reconstructed, 3,600 paired contrast rows reconstructed, 150 cross-curvature paired-validation blocks checked, 50 distinct truths reconstructed.
- Artifact ID: 14613021797; ZIP size 168,798,378 bytes.
- Downloaded ZIP SHA-256: `1d8667e2f7391283237e27b0c6db26cac6e9dfd56d344649b6931ec047745c0c`.

The archive was downloaded and inspected locally. The detailed readout independently checked unique keys and exact one-to-one coverage of all 2,250 SUCCESS/PASS ledger rows by posterior-summary and weight-diagnostic rows.

## Descriptive results

Across all 150 replicate-truths × curvature strata × validation-size blocks (the validation sizes are repeated within each replicate and are not independent):

| Arm | Mean posterior rank | Mean median bias (km/s/Mpc) | RMSE (km/s/Mpc) | 68% interval coverage |
|---|---:|---:|---:|---:|
| A Catalogue-only | 0.9909 | -15.81 | 17.38 | 0.7% |
| B Full-parent oracle | 0.4811 | +0.36 | 3.60 | 70.0% |
| C Known-π IPW | 0.4908 | +0.15 | 4.70 | 52.7% |
| D Degree-one MAP IPW | 0.5244 | -0.58 | 5.80 | 45.8% |
| E Degree-one integrated IPW | 0.5167 | -0.45 | 5.46 | 54.4% |

These are descriptive pooled summaries; the repeated validation sizes and common replicate truths must be respected in uncertainty calculations. The catalogue-only arm is strongly miscalibrated in this particular synthetic design. The full-parent oracle is the best-calibrated reference arm, as expected for a matched synthetic benchmark, not a deployable estimator.

At validation size 500:

| Curvature | Arm | Mean rank | Median bias | RMSE | 68% coverage |
|---:|---|---:|---:|---:|---:|
| -1 | B oracle | 0.481 | +0.36 | 3.60 | 70% |
| -1 | C known-π | 0.525 | -0.62 | 5.37 | 42% |
| -1 | D MAP | 0.710 | -3.53 | 6.30 | 44% |
| -1 | E integrated | 0.700 | -3.40 | 6.14 | 52% |
| 0 | B oracle | 0.481 | +0.36 | 3.60 | 70% |
| 0 | C known-π | 0.482 | +0.37 | 4.35 | 54% |
| 0 | D MAP | 0.507 | -0.09 | 4.80 | 50% |
| 0 | E integrated | 0.501 | -0.03 | 4.64 | 52% |
| +1 | B oracle | 0.481 | +0.36 | 3.60 | 70% |
| +1 | C known-π | 0.465 | +0.70 | 4.30 | 62% |
| +1 | D MAP | 0.303 | +3.05 | 5.49 | 50% |
| +1 | E integrated | 0.306 | +2.97 | 5.37 | 52% |

At n=500, E has lower RMSE than D in all three curvature strata, but E does not consistently beat C. Paired squared-error contrasts (E minus C; positive means E worse) are +8.86 [95% paired t interval -1.21, 18.92] at gamma=-1, +2.66 [-1.96, 7.27] at gamma=0, and +10.37 [3.43, 17.32] at gamma=+1. E is better than D on average at gamma=-1 and gamma=0; the E-D squared-error contrast is -2.06 [-3.37, -0.75] and -1.45 [-2.99, 0.09], respectively; at gamma=+1 it is -1.26 [-2.90, 0.38]. These intervals are exploratory and are not adjusted for the many strata/comparisons.

The E arm's minimum node-level Kish ESS at n=500 decreases with positive curvature: mean 94.2 (gamma=-1), 64.2 (gamma=0), 50.4 (gamma=+1). Mean maximum-node top-5%-weight mass rises from 44.0% to 51.4% to 58.1%. This is consistent with increasing weight concentration in the difficult positive-curvature condition; ESS associations are descriptive, not causal.

## Numerical grid sensitivity: preliminary only

I performed an interpolation-only sensitivity check on all 2,250 stored posterior curves: linearly interpolated each stored 561-point log-likelihood onto a 2,241-point geometric H0 grid over [30,110], then renormalized by trapezoidal integration with respect to dH0.

- Maximum absolute posterior-rank shift: 0.0008474.
- Maximum posterior-median shift: 0.00335 km/s/Mpc.
- Maximum q16/q84 shifts: 0.00418 / 0.00446 km/s/Mpc.

This is reassuring as a *representation/interpolation sensitivity check*, but it is **not an exact likelihood recomputation** on a refined grid. It cannot certify H0-grid convergence, because it interpolates the already sampled likelihood rather than re-evaluating the underlying selected likelihood and coefficient quadrature at new H0 points. The scientific calibration interpretation remains provisional until an exact refined-grid audit is implemented.

## Interpretation and limits

1. This is a synthetic prior-predictive benchmark only, not an observational H0 result or a real-catalogue validation.
2. The degree-one estimated inclusion model is misspecified when gamma is nonzero. Integrating its coefficients does not marginalize over omitted quadratic model structure.
3. Increasing validation size from 100 to 500 does not remove the positive-curvature D/E bias pattern in this pilot.
4. E improves on D's RMSE in these strata but is not uniformly superior to known-π IPW C; positive-curvature E remains worse than C in paired squared error.
5. Weight instability is material in E, particularly at positive curvature.
6. The 50-replicate results are a pilot. No universal estimator-superiority claim is supported.

**Status:** execution and artifact audit PASS; descriptive readout complete; exact H0-grid convergence remains OPEN. Keep this PR draft and unmerged.
