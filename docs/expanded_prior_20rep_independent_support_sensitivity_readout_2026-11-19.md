# Expanded-support five-arm pilot: independent prior-truncation readout

**Status:** exploratory 20-replicate synthetic engineering pilot. This report is an independent post-run reconstruction from the immutable Actions artifact produced by workflow `44482164577` (artifact `14513112496`). It does not certify calibration or estimator superiority.

## Lifecycle and provenance

- Frozen support: [30,110] km s^-1 Mpc^-1; truth H0=70.
- Parent population: 1,200 galaxies; 30 detected events.
- Curvature: gamma = -1, 0, +1; validation sizes 100, 300, 500.
- Five arms A-E; 20 paired parent/event replicates per curvature.
- Planned cells: 900; observed statuses: 900 SUCCESS, 900 independent child-audit PASS.
- Independently verified child runs: 60; reconstructed posterior-summary rows: 900.
- The orchestration audit reports PASS, including 1,148 manifest entries verified and 45 aggregate metric rows reconstructed.
- The 1,121-point H0 grid has spacing 80/1120 = 0.0714286. The same-likelihood restriction [50,90] uses exact nodes from this grid.

The audit confirms lifecycle completeness and numerical reconstruction of the saved summaries. It does not establish posterior calibration.

## Same-likelihood prior-support sensitivity

For each of the 900 saved arm/replicate/validation-size likelihood curves, I independently normalized the same curve on [30,110] and again on [50,90]. All 900 expanded-support posterior medians reproduced the archived medians within 2e-7. This is a paired prior-support transformation, not a new simulation.

The effect is concentrated in catalogue-only arm A:

| gamma | validation n | Mean H0 median on [30,110] | Mean median on [50,90] | Mean q2.5% on [30,110] | Mean q2.5% on [50,90] | Mean posterior mass below 50 |
|---:|---:|---:|---:|---:|---:|---:|
| -1 | 100/300/500 | 52.872 | 55.120 | 44.112 | 50.525 | 0.277 |
| 0 | 100/300/500 | 53.296 | 55.900 | 44.091 | 50.676 | 0.260 |
| +1 | 100/300/500 | 54.227 | 57.135 | 44.257 | 50.980 | 0.239 |

The catalogue-only arm's mean posterior median shifts downward by about 2.25–2.91 when moving from [50,90] to [30,110], and its mean lower 2.5% quantile shifts downward by about 6.41–6.72. Between 50% and 70% of its 20 replicate curves per curvature place more than 10% posterior mass below 50. Thus the earlier [50,90] support materially truncated this arm's lower tail.

For comparison, across the same 20 replicates, the mean median shifts are small for the other arms (roughly 0.001–0.055 in absolute value), and their mean posterior mass below 50 is approximately zero. The prior truncation was therefore not a common-mode effect across arms.

Validation-size rows for arms A, B and C are identical in this pilot because those arms do not depend on the fitted validation sample size; D and E can vary with validation size. The repeated rows are not independent evidence.

## Descriptive five-arm metrics (not calibration claims)

The 20-replicate success-only summaries show substantial curvature dependence and broad uncertainty. Selected examples:

- At gamma=-1, n=500: D has median bias -2.697 and RMSE 5.008; E has bias -2.494 and RMSE 4.962.
- At gamma=0, n=500: D has bias +0.403 and RMSE 4.327; E has bias +0.552 and RMSE 4.363.
- At gamma=+1, n=500: D has bias +3.340 and RMSE 5.545; E has bias +3.367 and RMSE 5.508.
- Catalogue-only A remains strongly negatively biased in this design (bias from -15.77 to -17.13 across curvature settings), even after widening the support; the expanded support changes its posterior summaries but does not remove the catalogue-selection problem.

With only 20 parent/event replicates per curvature, coverage fractions are coarse and their Wilson intervals are wide. These results should not be used to rank estimators or claim coverage calibration.

## Limitations and next checks

1. The experiment is a 20-replicate pilot, not the confirmatory follow-up to the 200-replicate [50,90] experiment.
2. The wider [30,110] prior is a sensitivity setting, not a uniquely justified cosmological prior. A substantial lower-tail shift under A remains.
3. The 1,121-point grid has not yet passed a full difficult-case refinement audit against 2,241 or denser grids. The same-likelihood comparison inherits this discretization.
4. The 25-node tensor quadrature and fixed KDE bandwidth are unchanged; this pilot does not validate their approximation error.
5. Oracle arm B is diagnostic only. No LVK/GWTC events, public galaxy catalogues, or observational H0 measurements are used.

## Interpretation

The defensible result is narrow: in this synthetic setup, [50,90] materially truncates the catalogue-only posterior's lower tail, while the other four arms have negligible mass below 50 in these pilot realizations. Widening the common support is necessary for fairer display of the same likelihood curves, but it does not resolve catalogue-incompleteness bias, establish calibration, or demonstrate a universally superior estimator.
