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

## F. Completed 20-replicate paired degree-one/degree-two batch

The prior-matched paired experiment completed successfully in GitHub Actions run **44773150331** on commit `489234ec26ef9f505047dccba5f7ab35f456bab4`. Artifact **14657642461** (`paired-degree-multireplicate-pilot`) has SHA-256 `e92cdbf589548a2113c6db8f9a764d577bdcbc410554722a36e1cffec43b76f1`.

The independent outer audit passed:
- 600/600 lifecycle cells SUCCESS for each degree;
- 60/60 replicate-by-curvature D1/D2 child pairs verified;
- 2,190 outer manifest entries and 3,736 nested manifest entries verified;
- 600 paired posterior-summary rows and 30 paired aggregate rows reconstructed;
- common-arm A/B/C likelihoods invariant across degrees.

The experiment used master seed `2026112020`, 20 paired parent/event replicates, gamma {-1,0,+1}, N_parent=1200, 30 detected events, validation sizes {100,500}, H0 truth 70, prior [30,110], 281 H0 points, and quadrature order 5. Both degrees used identical parent, event, catalogue and validation inputs. The earlier single-realization quadrature pilot is not counted as a replicate.

### Degree-specific results for fitted arms D and E

Bias is posterior median minus H0 truth; RMSE is across the 20 paired replicates. The paired squared-error difference is D2 minus D1; its 95% interval is a replicate-level paired t interval.

| gamma | n | Arm | D1 bias | D2 bias | D1 RMSE | D2 RMSE | D1 68% coverage | D2 68% coverage | Paired squared-error difference (95% CI) |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| -1 | 100 | D MAP | -2.644 | -0.325 | 4.771 | 4.191 | 0.60 | 0.60 | -5.20 (-16.63, 6.24) |
| -1 | 100 | E integrated | -1.949 | +0.029 | 3.989 | 4.202 | 0.65 | 0.75 | +1.75 (-7.12, 10.63) |
| -1 | 500 | D MAP | -1.365 | +1.650 | 3.523 | 5.104 | 0.80 | 0.60 | +13.65 (-1.39, 28.69) |
| -1 | 500 | E integrated | -1.219 | +1.505 | 3.376 | 4.639 | 0.80 | 0.60 | +10.12 (-0.35, 20.60) |
| 0 | 100 | D MAP | -0.007 | +0.879 | 4.196 | 4.139 | 0.60 | 0.55 | -0.48 (-10.89, 9.92) |
| 0 | 100 | E integrated | +0.373 | +0.936 | 3.852 | 4.355 | 0.65 | 0.75 | +4.13 (-7.55, 15.80) |
| 0 | 500 | D MAP | +1.343 | +1.957 | 3.786 | 4.588 | 0.70 | 0.55 | +6.71 (-4.42, 17.84) |
| 0 | 500 | E integrated | +1.395 | +1.928 | 3.809 | 4.434 | 0.65 | 0.55 | +5.15 (-3.74, 14.04) |
| +1 | 100 | D MAP | +3.054 | +2.238 | 5.156 | 5.167 | 0.40 | 0.35 | +0.12 (-19.14, 19.38) |
| +1 | 100 | E integrated | +3.027 | +1.916 | 4.950 | 4.831 | 0.45 | 0.70 | -1.16 (-13.70, 11.37) |
| +1 | 500 | D MAP | +4.212 | +2.224 | 5.605 | 4.976 | 0.30 | 0.45 | -6.66 (-18.35, 5.04) |
| +1 | 500 | E integrated | +4.103 | +2.202 | 5.530 | 4.810 | 0.40 | 0.45 | -7.45 (-17.46, 2.56) |

All 12 paired squared-error intervals include zero. At gamma=+1, D2 shifts D/E posterior medians downward by about 1.1–2.0 km/s/Mpc relative to D1, moving closer to the truth in this batch. At gamma=-1, the direction reverses at n=500: D2 has larger positive bias and RMSE. Gamma=0 shows small-to-moderate changes. Thus the degree effect is curvature- and validation-size-dependent; this 20-replicate batch does not establish an accuracy advantage for D2.

The common arms A/B/C are exactly invariant across degree, as intended. They are not additional independent comparisons.

### Node-wise weight concentration

Median across replicates of minimum node-wise Kish ESS and maximum node-wise top-5% weight mass for arm E:

| gamma | n | D1 minimum node ESS | D2 minimum node ESS | D1 maximum top-5% mass | D2 maximum top-5% mass |
|---:|---:|---:|---:|---:|---:|
| -1 | 100 | 46.86 | 2.86 | 0.569 | 0.964 |
| -1 | 500 | 73.53 | 8.64 | 0.453 | 0.784 |
| 0 | 100 | 38.24 | 2.93 | 0.641 | 0.976 |
| 0 | 500 | 65.73 | 12.16 | 0.520 | 0.772 |
| +1 | 100 | 27.88 | 4.67 | 0.698 | 0.960 |
| +1 | 500 | 44.04 | 20.72 | 0.566 | 0.683 |

This is the strongest warning from batch 1. Degree-two coefficient integration is numerically computable, but its posterior nodes can induce extreme inverse-probability weights, especially at validation n=100. Increasing quadrature order or refining the H0 grid does not repair this weight instability.

### Interpretation and next step

1. The H0-grid numerical gate passed for the frozen 30-event KDEs, and the D1/D2 20-replicate pairing/audit passed.
2. Degree two changes the fitted-arm posterior, but no paired squared-error interval excludes zero in this 20-replicate batch.
3. Degree-two arm E has severe node-wise weight concentration, so any apparent bias reduction must be evaluated alongside this instability.
4. Do not claim D2 superiority, calibration, or observational validity.
5. A second independent 20-replicate paired batch was launched on branch `research/paired-degree-20rep-batch2-20261121`, master seed `2026112121`, workflow run **44841133375**. At the time this readout was updated, the paired simulation step was in progress and no second-batch result was available.

The project remains a synthetic low-redshift Euclidean benchmark. No LVK/GWTC data, real galaxy survey, or observational H0 measurement is used.
