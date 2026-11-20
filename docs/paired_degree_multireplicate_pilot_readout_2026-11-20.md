# Prior-matched centered degree-two validation — 2026-11-20

## A. Completed quadrature sensitivity pilot

The initial workflow `44723868794` and latest-head rerun `44726785025` completed successfully on the prior-matched PR #133 branch. The latest artifact `14623863771` (degree-matched-quadrature-sensitivity) was downloaded and independently inspected; it reproduces the same summaries as the initial artifact `14623002861`.

- Latest artifact SHA-256: `1f2ffba18438161cea4c7eb44eb7037edc21b006df1039fd078c403a32a629cf`.
- Initial artifact SHA-256: `ae13f348a3bd600b1f1f67a998a75e7da3b0a243596545483111139e4c114f0a`.
- Outer manifest status PASS; all 114 listed files matched recorded byte counts and SHA-256 digests.
- Six nested child manifests (D1/D2 × quadrature orders 5/7/9) independently checked; all nested entries matched byte counts and SHA-256 digests.
- All six child audit reports PASS; each reconstructed 10 posterior summary rows.
- Maximum independently reconstructed log-likelihood discrepancy across children was at most approximately 1.1e-14.
- D1/D2 common arms A/B/C agreed within floating-point tolerance.
- The fitted D2 basis recorded in its child configuration was `[1, x, x^2 - c_g]`, with centre `c_g=0.2845005579191241`.

This artifact is **not** the 20-replicate paired experiment. It is a single-realization numerical sensitivity pilot with N_parent=400, four detected events, gamma=+1, validation frame 300, validation sizes 100/300, H0 truth 70, prior [30,110], 41-point H0 grid and quadrature orders 5/7/9.

## B. Quadrature convergence

For the integrated-IPW posterior density, order 7 to 9 L1 differences were:

| Degree | Validation n | Density L1 (7 vs 9) | Absolute median shift |
|---:|---:|---:|---:|
| 1 | 100 | 1.253e-5 | 0.000180 |
| 1 | 300 | 5.729e-7 | 0.000008 |
| 2 | 100 | 8.368e-6 | 0.000118 |
| 2 | 300 | 7.022e-7 | 0.000010 |

Thus integrated likelihood quadrature is numerically stable in this particular realization at orders 7/9.

## C. Weight instability is the limiting warning

Order-9 minimum node-level Kish ESS and maximum node top-5%-weight mass:

| Degree | Validation n | Minimum node Kish ESS | Maximum node top-5% mass |
|---:|---:|---:|---:|
| 1 | 100 | 9.60 | 79.6% |
| 1 | 300 | 21.87 | 57.4% |
| 2 | 100 | 3.17 | 92.3% |
| 2 | 300 | 14.47 | 62.5% |

The degree-two model is especially unstable at n=100: a small fraction of galaxies dominate the inverse-probability weighted estimate. Increasing quadrature order does not solve this; the minimum node ESS becomes lower as more posterior nodes reveal more extreme inclusion-probability configurations. D1 is also unstable, so this is not exclusively a degree-two problem.

**Key distinction:** numerical convergence of the coefficient-posterior integral is not evidence that the IPW estimator has acceptable finite-sample variance or calibrated H0 uncertainty.

## D. Posterior shifts are descriptive only

At order 9, the D2 integrated-IPW posterior median is about 0.981 km/s/Mpc below D1 for n=100 and 1.434 km/s/Mpc below D1 for n=300. Yet both are materially above the fixed truth 70 in this four-event realization; the parent-oracle median is also high. This one realization cannot establish that D2 improves inference. It mainly demonstrates that model degree changes the weighted likelihood while leaving common arms invariant.

## E. Prior-matching correction

The D2 basis is now centered as `[1,x,x^2-c_g]`. With independent Gaussian coefficient priors, this makes the D1 prior on (alpha,beta) exactly equal to the D2 marginal prior on (alpha,beta), while gamma is independent. The covariance audit is implemented in `scripts/audit_degree_prior_parameterization.py` and its regression test. This removes the earlier raw-versus-centered prior mismatch.

## F. Current 20-replicate paired experiment

The PR #133 workflow was expanded to 20 paired replicates with master seed `2026112020`, gamma {-1,0,+1}, N_parent=1200, 30 events, validation sizes {100,500}, fixed H0=70, prior [30,110], 281-point grid and quadrature order 5.

Workflow 44723868722 has started: the centered-prior covariance audit and paired-orchestrator regression tests passed, and the 20-replicate simulation step is still in progress. The Python 3.10/3.11/3.12 test matrix (run 44724793170) is also in progress. No 20-replicate result is claimed here until the simulation and independent outer audit both finish.

## Conclusion

The prior-matched parameterization is now in place and the one-realization quadrature audit passes. The remaining risk is weight instability, especially for degree-two integrated IPW at small validation size. The next experiment must report weight diagnostics alongside bias, RMSE and coverage. No estimator superiority or calibration claim is supported by this quadrature pilot.
