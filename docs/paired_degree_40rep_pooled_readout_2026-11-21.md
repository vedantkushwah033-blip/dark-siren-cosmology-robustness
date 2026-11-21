# Paired degree-one / degree-two 40-replicate pooled readout (2026-11-21)

## Scope and claim limit

This is a fixed-truth, low-redshift, synthetic engineering pilot. It is not a prior-predictive calibration study, an observational H0 result, an LVK/GWTC analysis, or evidence that either inclusion model is universally superior. All contrasts below compare degree-two minus degree-one within matched replicate inputs.

## Design

- 40 paired analysis replicates, assembled from two independent 20-replicate batches.
- Curvature gamma in {-1, 0, +1}; validation sizes n={100, 500}.
- N_parent=1200; 30 detected events; H0 truth 70; H0 prior support [30,110].
- H0 grid: 281 points; inclusion-coefficient quadrature order 5; Scott multiplier 0.7.
- Degree one uses [1,x]. Degree two uses centered [1,x,x^2-c_g], with c_g=E_g[x^2], matched marginal priors on alpha and beta and an independent Gaussian prior on gamma.
- Both degree models share the same analysis and validation seeds within each replicate. Batch master seeds are 2026112020 and 2026112121.

## Artifact and lifecycle audit

Batch 2 workflow #44841133375 completed successfully. Its artifact is ID 14683823051, SHA-256 `10ad7f7a0c6528b82f8d4816a841c18ac164490e3731a8a5da6479ac612bb3d1`. The independent outer artifact audit reports PASS: 2,190 outer-manifest entries and 3,736 nested-manifest entries verified; 60 paired child blocks checked; 600 between-degree paired posterior-summary rows reconstructed; and 30 between-degree aggregate metric rows reconstructed. Separately, each degree-specific five-arm audit reports 600/600 planned cells SUCCESS, 60 successful child runs, 1,088 aggregate-manifest entries, 960 within-degree arm/reference contrast rows, and 30 within-degree aggregate metric rows. The 960 count is not the D2-minus-D1 contrast count.

The earlier batch-1 artifact SHA-256 is `e92cdbf589548a2113c6db8f9a764d577bdcbc410554722a36e1cffec43b76f1`. Across both batches, each degree has 1,200/1,200 SUCCESS summary cells with PASS audit status. Replicate seeds match across degrees within each batch; analysis and validation seeds have zero overlap between batches.

The cross-Python numerical-diagnostics workflow #44841133336 completed successfully on Python 3.10, 3.11, and 3.12. The long Python 3.12 numerical workflow steps, including paired-degree experiments and independent artifact audits, completed successfully.

## Pooled degree-two minus degree-one contrasts

The table reports the paired mean difference in squared median error, with a two-sided 95% Student-t interval across 40 matched replicates. Negative values favor degree two on this metric; positive values favor degree one. These are exploratory contrasts, not a pre-registered familywise inference.

| gamma | validation n | arm | mean Δ squared error | 95% CI | mean Δ absolute error |
|---:|---:|---|---:|---:|---:|
| -1 | 100 | D plug-in IPW | -4.652 | [-12.561, 3.256] | -0.230 |
| -1 | 100 | E integrated IPW | -1.092 | [-7.677, 5.493] | -0.099 |
| -1 | 500 | D plug-in IPW | 0.253 | [-10.022, 10.528] | 0.023 |
| -1 | 500 | E integrated IPW | -1.125 | [-9.810, 7.560] | -0.010 |
| 0 | 100 | D plug-in IPW | 2.038 | [-5.547, 9.623] | 0.086 |
| 0 | 100 | E integrated IPW | 2.979 | [-3.843, 9.801] | 0.378 |
| 0 | 500 | D plug-in IPW | 3.019 | [-2.771, 8.809] | 0.201 |
| 0 | 500 | E integrated IPW | 2.211 | [-2.543, 6.966] | 0.229 |
| +1 | 100 | D plug-in IPW | 1.591 | [-9.135, 12.318] | 0.124 |
| +1 | 100 | E integrated IPW | -2.731 | [-9.847, 4.384] | -0.325 |
| +1 | 500 | D plug-in IPW | -5.937 | [-12.043, 0.170] | -0.778 |
| +1 | 500 | E integrated IPW | -5.946 | [-11.278, -0.613] | -0.705 |

Only the gamma=+1, n=500 integrated-IPW contrast has a 95% interval excluding zero in this exploratory 12-contrast table. The corresponding plug-in interval narrowly includes zero. The remaining intervals include zero. Do not interpret this isolated result as a robust degree-two advantage: it is one of multiple comparisons, under one fixed truth and a narrow design.

## Pooled D/E performance snapshots

Across the 40 replicates per cell, selected degree-one / degree-two RMSE values (and 68% central interval coverage) were:

| gamma | n | arm | degree 1 RMSE / coverage | degree 2 RMSE / coverage |
|---:|---:|---|---:|---:|
| -1 | 100 | D | 5.135 / 0.575 | 4.660 / 0.525 |
| -1 | 100 | E | 4.347 / 0.625 | 4.219 / 0.750 |
| -1 | 500 | D | 4.203 / 0.650 | 4.233 / 0.675 |
| -1 | 500 | E | 4.041 / 0.675 | 3.899 / 0.750 |
| 0 | 100 | D | 4.406 / 0.525 | 4.631 / 0.550 |
| 0 | 100 | E | 3.904 / 0.675 | 4.268 / 0.750 |
| 0 | 500 | D | 3.550 / 0.750 | 3.952 / 0.650 |
| 0 | 500 | E | 3.517 / 0.725 | 3.818 / 0.700 |
| +1 | 100 | D | 4.774 / 0.450 | 4.938 / 0.425 |
| +1 | 100 | E | 4.542 / 0.550 | 4.231 / 0.775 |
| +1 | 500 | D | 4.780 / 0.425 | 4.112 / 0.600 |
| +1 | 500 | E | 4.719 / 0.475 | 4.040 / 0.625 |

The patterns are mixed. Degree two improves several positive-curvature cells but worsens some zero-curvature cells; degree-one vs degree-two differences are not uniformly directional.

## Persistent numerical warning: degree-two integrated-IPW concentration

For arm E under degree two, median across replicates of the minimum node-wise Kish ESS and maximum node-wise top-5%-weight mass were:

| gamma | validation n | minimum node Kish ESS | maximum node top-5% mass |
|---:|---:|---:|---:|
| -1 | 100 | 2.482 | 0.971 |
| -1 | 500 | 6.877 | 0.758 |
| 0 | 100 | 3.013 | 0.964 |
| 0 | 500 | 10.975 | 0.771 |
| +1 | 100 | 5.302 | 0.960 |
| +1 | 500 | 19.154 | 0.717 |

Thus, even when degree-two integrated IPW has lower H0 RMSE in some cells, its coefficient-node-specific inverse weights can be extremely concentrated, especially at n=100. This is a material stability limitation. Do not confuse these node-wise IPW Kish ESS values with the coefficient-posterior quadrature ESS.

## Interpretation and next work

1. The second independent batch reproduces the lifecycle and pairing structure and expands the fixed-design pilot from 20 to 40 replicates.
2. No general degree-two advantage is established. The strongest apparent gain is localized to gamma=+1, n=500, and is not sufficient to claim superiority.
3. Degree-two integrated IPW concentration remains a central warning, especially for small validation samples.
4. Next scientific step should be a broader, pre-specified stress design with more curvature and validation sizes, stronger attention to overlap/positivity, and explicit calibration diagnostics. Increase replicate count only after defining the target estimand, multiplicity plan, and failure-inclusive summary policy.
5. Before any scientific claim, separate fixed-truth frequentist coverage, prior-predictive SBC, and density-recovery metrics; none should be substituted for another.

