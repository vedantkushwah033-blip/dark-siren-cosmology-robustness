# Stage A selection-misspecification pilot: audited results (2026-10-29)

## Status

**Exploratory pilot completed; artifact-audit job passed.** This report records the 200-replicate Stage A workflow run for commit `18863fcda677ec5df755314b1044f73a344a1bb4`. It is a controlled synthetic study, not a gravitational-wave catalogue analysis and not a measurement of the real Hubble constant.

- Workflow: [selection-misspecification-pilot run #1](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40793766654)
- Three-version test matrix on the same commit: [tests run #1877](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40793766726) — Python 3.10, 3.11, and 3.12 all succeeded.
- Saved artifact: `selection-misspecification-stage-a-18863fcda677ec5df755314b1044f73a344a1bb4` (raw replicate table, summary JSON, manifest; 90-day retention).
- The independent posterior-grid/artifact audit is a successful workflow step. Raw artifact SHA-256 checksums were verified again during this report's analysis.

## Design

The workflow used 200 paired replicate seeds beginning at 20261030, 4,000 simulated parent galaxies per replicate, 30 detected events, fixed true $H_0=70$ km s$^{-1}$ Mpc$^{-1}$, an $H_0$ grid of 81 points over [50, 90], and 401 redshift grid points over [0.01, 0.4]. The toy detection threshold was 250 Mpc, with Gaussian distance noise scale 140 Mpc. The true catalogue inclusion curve had floor 0.1, midpoint 0.22, and width 0.05. KDE used the corrected Scott weighted-bandwidth default; an explicit bandwidth factor instead means the factor times the weighted standard deviation of the original sample.

The six paired arms were: analytic-truth density (inaccessible diagnostic control), oracle KDE (full simulated parent population), incomplete catalogue without correction, inverse-inclusion weighting under the declared correct selection law, a shifted transition midpoint (0.27), and a shallower transition width (0.075).

## Main results

Bias is the mean posterior mean minus true $H_0$. Coverage is for nominal 68% posterior intervals. Coverage intervals are Wilson 95% binomial intervals; the pilot is deliberately small, so these intervals remain wide.

| Arm | Bias (km s$^{-1}$ Mpc$^{-1}$) | Bias SE | RMSE | 68% coverage | Wilson 95% interval | Mean interval width |
|---|---:|---:|---:|---:|---:|---:|
| Analytic truth | +0.069 | 0.266 | 3.754 | 136/200 = 68.0% | [61.2%, 74.1%] | 7.223 |
| Oracle KDE | +0.290 | 0.267 | 3.771 | 137/200 = 68.5% | [61.8%, 74.5%] | 7.245 |
| Catalogue only | −9.468 | 0.224 | 9.984 | 15/200 = 7.5% | [4.6%, 12.0%] | 8.654 |
| Correct weighting | +0.351 | 0.267 | 3.785 | 142/200 = 71.0% | [64.4%, 76.8%] | 7.265 |
| Shifted transition | −0.409 | 0.265 | 3.756 | 139/200 = 69.5% | [62.8%, 75.5%] | 7.653 |
| Shallow slope | −1.824 | 0.260 | 4.091 | 127/200 = 63.5% | [56.6%, 69.9%] | 7.737 |

The density error measure also sharply separates ignoring catalogue selection from the other arms: mean density ISE was 0.8724 for catalogue-only, compared with 0.0383 for correct weighting, 0.0189 for oracle KDE, and zero by construction for analytic truth. These ISE values are specific to this toy density/grid and should not be read as a universal ranking.

## Paired contrasts

For each replicate, calculate $|\hat H_{0,\mathrm{arm}}-H_{0,\mathrm{true}}|-|\hat H_{0,\mathrm{reference}}-H_{0,\mathrm{true}}|$. Negative values favour the comparison arm. Intervals below are percentile-bootstrap 95% intervals resampling the 200 paired replicate rows (10,000 bootstrap draws; fixed analysis seeds).

| Contrast | Mean paired absolute-error difference | Bootstrap 95% interval |
|---|---:|---:|
| Analytic truth − catalogue only | −6.563 | [−7.128, −5.942] |
| Oracle KDE − catalogue only | −6.541 | [−7.110, −5.930] |
| Correct weighting − catalogue only | −6.536 | [−7.129, −5.929] |
| Shifted transition − correct weighting | +0.041 | [−0.070, +0.152] |
| Shallow slope − correct weighting | +0.313 | [+0.028, +0.594] |
| Oracle KDE − analytic truth | +0.022 | [−0.032, +0.077] |
| Correct weighting − analytic truth | +0.027 | [−0.074, +0.129] |

In this pilot, ignoring catalogue incompleteness caused a large error in the inferred host density and a substantial downward shift in inferred $H_0$. Correct weighting greatly improved performance under the stipulated, known inclusion law. A shallower-than-true transition produced a modest but detectable increase in paired absolute error relative to correct weighting in this particular simulation, whereas the shifted-midpoint comparison is inconclusive at this sample size.

The corrected-weighting and oracle/analytic-control arms have similar performance within paired uncertainty here. This does **not** show that weighting is generally unbiased or robust: all arms use the same simplified low-redshift toy likelihood, one true $H_0$, one event count, one incompleteness regime, and one latent host-population shape.

## Integrity checks performed

The downloaded artifact was checked for:
- exact SHA-256 match to both entries in `manifest.json`;
- 1,200 rows = 200 replicates × 6 arms, with one common seed and true $H_0$ per replicate;
- exactly one row per replicate/arm, all six arms present for every replicate;
- finite posterior summaries, unit posterior integrals, and consistency of interval widths with endpoints;
- posterior-mean reconstruction from the persisted posterior density grids (median absolute reconstruction difference below $10^{-6}$; maximum difference about 0.0093 in $H_0$ units due to discretized integration conventions).

The workflow's independent audit step also completed successfully. These checks support artifact integrity and calculation reproducibility; they do not validate realism of the scientific model.

## Limitations and decision

**Decision: retain exploratory status; no confirmatory or real-data claim.** The pilot is useful evidence that catalogue incompleteness can severely distort this particular synthetic inference problem, and that misspecifying the inclusion curve can degrade results. It does not establish performance for realistic LVK events.

Still outstanding before any confirmatory claim:
1. a prospective precision calculation from paired pilot variability and a new, non-overlapping seed stream;
2. bandwidth and grid-resolution sensitivity after the KDE fix, with raw output and independent audit for each configuration;
3. calibration over multiple fixed $H_0$ values, event counts, completeness strengths, and host-population shapes;
4. implementation and calibration of a dataset-level selection-conditioned likelihood with the intended $H_0$-dependent distance-noise model;
5. literature comparison and a final scientific review of model assumptions.

Do not combine this result with older KDE-based runs produced under the previous SciPy bandwidth semantics.
