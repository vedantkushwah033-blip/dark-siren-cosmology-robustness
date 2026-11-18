# Seven-gamma curvature pilot execution record

This isolated branch exists to execute the pre-specified synthetic density-only quadratic-logit curvature pilot after the exact-head correction CI passed on PR #59.

## Frozen design

- Gamma: -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0.
- 20 paired replicates per gamma (140 planned total).
- Parent population: 4,000 galaxies per replicate.
- Independent validation frame: 1,500; nested validation sizes 100, 300, 1,000.
- Linear coefficient beta: -3.9; gamma-specific intercept calibrated to preserve parent-average inclusion.
- Quadrature order: 5.
- Density grid: 401 points on z in [0.01, 0.40].
- Absolute KDE bandwidth candidates: 0.006, 0.008, 0.010, 0.012, 0.016, 0.024, 0.032, 0.040, 0.050, 0.060, 0.080.
- Seed schedule: master seed 20261215, stride 10; validation seed is population seed + 1.
- Arms A-E and density metrics follow the source code and per-cell independent auditor.

## Trigger and isolation

The existing repository test workflow is conditionally configured on this branch to skip its broad three-version test matrix and run the seven gamma cells plus the cross-gamma auditor. The workflow is limited to this branch; no change is made to main. Raw gamma-cell artifacts are uploaded even when a cell fails. Cross-gamma summaries run only if every cell job succeeds.

## Claim limits

This experiment is synthetic density recovery only. It does not run the GW likelihood or H0 posterior, does not establish calibration, and does not validate an observational galaxy catalogue. A complete run requires all 140 replicate records, passing per-gamma raw-input audits, and a passing cross-gamma audit. Any failed or non-evaluable replicate must be accounted for; do not silently drop it.
