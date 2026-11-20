# Latest-head degree-matched quadrature sensitivity artifact — 2026-11-20

## Artifact provenance and independent verification

PR #133 workflows `44727139750` and latest-head rerun `44727858540` completed successfully. The latest artifact `14625221228` (`degree-matched-quadrature-sensitivity`) was downloaded and independently inspected; its configuration, weight diagnostics, and posterior summaries reproduce the prior artifact to floating-point precision.

- Latest ZIP SHA-256: `dea86991999b8b2d403f3d2d9826b88622572d1d89dc83f69fdec47c1f1058ce`.
- Prior artifact `14624950996` ZIP SHA-256: `1fc482483ec68bcbecbcb73f83b8edb455604747fb2d4365e9761aa8ea5bac4c`.
- Root manifest: 114/114 files matched declared byte counts and SHA-256 digests.
- Six child manifests: 13/13 files each verified, zero mismatches.
- All six independent child audit reports: PASS.
- Maximum child log-likelihood reconstruction discrepancy: `1.42e-14`.
- Each child audit reconstructed 10 posterior-summary rows.

This artifact is a one-realization engineering/numerical sensitivity pilot, not calibration or an estimator ranking. The exact config is N_parent=400, four detected events, gamma=+1, validation frame=300, validation sizes 100/300, H0 truth=70, prior [30,110], 41-point H0 grid, and quadrature orders 5/7/9.

Latest-head rerun checks: all 114 outer-manifest entries and all six 13-file child manifests verified with zero mismatches. All six child audit reports PASS; maximum reconstruction discrepancy was `1.42e-14`. Comparing the latest and prior artifacts, maximum numeric changes in posterior-summary CSVs were `4.27e-14`, while weight-diagnostic CSVs were byte-identical. The rerun is the same fixed-seed realization, not an additional independent sample.

## Posterior summaries at quadrature order 9

All values below are posterior medians for this one fixed realization. The reported median error is posterior median minus the fixed truth H0=70.

| Arm | D1, n=100 | D2, n=100 | D1, n=300 | D2, n=300 |
|---|---:|---:|---:|---:|
| Catalogue only A | 76.779 | 76.779 | 76.779 | 76.779 |
| Full-parent oracle B | 91.473 | 91.473 | 91.473 | 91.473 |
| Known-probability IPW C | 84.133 | 84.133 | 84.133 | 84.133 |
| Degree-specific MAP IPW D | 85.735 | 84.446 | 84.841 | 83.282 |
| Degree-specific integrated IPW E | 85.661 | 84.680 | 84.854 | 83.420 |

The degree-two fitted inclusion model is correctly specified for this quadratic-logit synthetic DGP, unlike degree one. As an additional paired-input invariance check, I compared the stored likelihood arrays directly: H0/redshift grids and all A/B/C likelihood and density arrays are bitwise identical between D1 and D2 at both validation sizes. Only the degree-specific D/E arms differ. In this one realization D2 moves the D/E posterior medians closer to truth, particularly at n=300. This is not evidence of average improvement: the sample contains only four events, one parent/catalogue realization, one curvature, and one truth.

The oracle arm is also far from truth in this realization (median 91.47). This reinforces that an oracle inclusion probability does not remove event-sampling and finite-catalogue variability. It should not be interpreted as an estimator failure on its own.

## Quadrature and weight stability

Integrated-IPW posterior summaries are effectively unchanged between quadrature orders 7 and 9 (maximum reported median change below 0.001 H0 units). The independent likelihood reconstruction discrepancies are at floating-point scale.

However, the integrated weights are highly concentrated:

| Degree | n | Minimum node Kish ESS | Maximum node top-5% weight mass |
|---|---:|---:|---:|
| D1 | 100 | 9.60 | 79.6% |
| D1 | 300 | 21.87 | 57.4% |
| D2 | 100 | 3.17 | 92.3% |
| D2 | 300 | 14.47 | 62.5% |

The D2 model's additional flexibility can make estimated inclusion probabilities extremely variable at small validation sizes. This is a key trade-off: improved functional specification does not automatically imply stable IPW.

## Interpretation and next gate

This artifact validates numerical implementation and provides one illustrative D1/D2 contrast. It does not establish that D2 is better, that either estimated-IPW arm is calibrated, or that the method works on observational catalogues.

The active next experiment is the paired 20-replicate D1/D2 study with common synthetic inputs across degrees, three curvatures, and validation sizes 100/500. Its simulation and independent outer audit must complete before any multi-replicate inference is made.
