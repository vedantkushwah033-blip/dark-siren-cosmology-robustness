# Locked Stage B selection-misspecification confirmation: results and audit

**Status: Stage B simulation completed; its independent artifact audit passed.** This is a confirmatory experiment for a predeclared endpoint *within the stated synthetic model*, not confirmatory evidence for a real gravitational-wave analysis.

## Provenance and integrity

- Locked simulation and analysis code commit: `12ed103e83bc398df37db8eba8584a9b2b9bc9ad`.
- [Stage B GitHub Actions run](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40794942601): simulation, independent artifact audit, and upload all completed successfully.
- [Stage A follow-up re-run](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40794942515): completed successfully, including independent audit.
- The same-commit Python 3.10–3.12 CI matrix succeeded for an earlier source commit; exact-head CI for the latest primary-estimand changes was still running when this report was initially written. Check the current [Actions page](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions) before using the latest software state.
- Stage B artifact: `selection-misspecification-stage-b-12ed103e83bc398df37db8eba8584a9b2b9bc9ad` (ID 12864010259, 90-day retention).
- Raw CSV SHA-256: `916172113502e11c58ead35ab205259bafbe7cfca05b9248deedf01eec7c5dce`.
- Summary JSON SHA-256: `1082ebbd23c4405aa9043d942b5180c177a046044668166e5a80942e210411d9`.
- 3,000 rows = 500 replicates × six arms; seed range 20262030–20262529; all six arms are present and paired in every replicate. The runner recorded a 90.475-second runtime under Python 3.12.15, NumPy 2.5.4, pandas 3.1.0, and SciPy 1.18.1.

The artifact hashes were independently rechecked after downloading the ZIP. Bias, RMSE, coverage, and interval-width summaries were recomputed from the raw CSV and matched the saved summary; the workflow's posterior-grid and primary-estimand reconstruction auditor passed.

## Locked design

The run used 4,000 parent galaxies per replicate, 30 detected events, fixed truth (H_0=70) km s(^{-1}) Mpc(^{-1}), an 81-point (H_0) grid over [50, 90], a 401-point redshift grid over [0.01, 0.40], and the corrected weighted Scott KDE convention. Detection threshold was 250 Mpc and the Gaussian distance-noise scale was 140 Mpc. The true catalogue inclusion curve used floor 0.1, midpoint 0.22, and width 0.05. The analyst arms were analytic truth, oracle KDE, incomplete catalogue without correction, correct inverse-inclusion weighting, shifted transition midpoint 0.27, and width 0.075 (shallower transition). Stage B seeds do not overlap Stage A's 20261030–20261229 range.

## Primary predeclared result

For replicate (r), the locked primary estimand is

[
D_r=|widehat H_{0,mathrm{catalogue only},r}-H_{0,mathrm{true}}|
-|widehat H_{0,mathrm{correct weighting},r}-H_{0,mathrm{true}}|.
]

Positive values mean that correct weighting reduced absolute posterior-mean error. Across 500 paired replicates:

| Statistic | Result |
|---|---:|
| Mean (D) | **+5.9519 km s(^{-1}) Mpc(^{-1})** |
| Paired SD | 4.6784 |
| Monte Carlo SE | 0.2092 |
| Normal-approximation 95% interval | [5.5418, 6.3620] |
| Percentile-bootstrap 95% interval, 10,000 resamples | **[5.5410, 6.3573]** |
| Analysis/bootstrap seed | 20366759 |

The primary interval excludes zero by a large margin. The catalogue-only arm's posterior mean was closer to the fixed truth in 58 of 500 replicates; correct weighting had lower absolute error in 442/500 = 88.4%, Wilson 95% interval [85.3%, 90.9%]. This is strong evidence for the *predeclared contrast under this data-generating process*, not a universal guarantee that inverse weighting improves real analyses.

## Per-arm performance

Bias is mean posterior mean minus true (H_0). Coverage intervals are Wilson 95% score intervals for the empirical coverage of nominal 68% posterior intervals.

| Arm | Bias ± SE | RMSE | 68% coverage | Wilson 95% interval | Mean 68% interval width | Mean density ISE |
|---|---:|---:|---:|---:|---:|---:|
| Analytic truth | +0.612 ± 0.161 | 3.647 | 360/500 = 72.0% | [67.9%, 75.8%] | 7.351 | 0.0000 |
| Oracle KDE | +0.918 ± 0.160 | 3.695 | 364/500 = 72.8% | [68.7%, 76.5%] | 7.376 | 0.0200 |
| Catalogue only | −8.847 ± 0.138 | 9.369 | 48/500 = 9.6% | [7.3%, 12.5%] | 8.745 | 0.8624 |
| Correct weighting | +1.017 ± 0.163 | 3.787 | 348/500 = 69.6% | [65.4%, 73.5%] | 7.382 | 0.0386 |
| Shifted transition | +0.253 ± 0.161 | 3.615 | 362/500 = 72.4% | [68.3%, 76.1%] | 7.761 | 0.1280 |
| Shallow slope | −1.177 ± 0.157 | 3.709 | 352/500 = 70.4% | [66.3%, 74.2%] | 7.843 | 0.0992 |

The dominant result is clear: ignoring the stipulated catalogue incompleteness drives severe downward bias in this toy case and nominal 68% intervals cover the fixed truth only 9.6% of the time. The corrected and oracle/analytic-control arms perform much better. Density ISE is evaluated only on this declared toy redshift grid; it is not a universal estimator-quality metric.

## Secondary paired contrasts

These supplementary post-run comparisons use ( |mathrm{error}_{m first}|-|mathrm{error}_{m second}|), so negative values favor the first-listed arm. Percentile-bootstrap intervals use 10,000 paired resamples; RNG seed start is 20485637. These comparisons were not the single locked primary endpoint and should be treated as secondary.

| Contrast | Mean absolute-error difference | Paired bootstrap 95% interval | Interpretation in this toy setting |
|---|---:|---:|---|
| Correct weighting − catalogue only | −5.9519 | [−6.3567, −5.5448] | Large improvement under known inclusion law |
| Oracle KDE − analytic truth | +0.0228 | [−0.0170, +0.0627] | Finite-parent/KDE penalty unresolved at this precision |
| Correct weighting − analytic truth | +0.0946 | [+0.0240, +0.1650] | Small absolute-error penalty versus inaccessible truth control |
| Shifted transition − correct weighting | −0.0756 | [−0.1468, −0.0062] | Mild apparent improvement, consistent with bias compensation rather than evidence that misspecification is generally safer |
| Shallow slope − correct weighting | +0.0781 | [−0.1003, +0.2514] | Inconclusive; no stable degradation established by this run |
| Shallow slope − analytic truth | +0.1727 | [+0.0183, +0.3310] | Modest degradation versus ideal diagnostic control |
| Shifted transition − analytic truth | +0.0190 | [−0.0670, +0.1043] | No clear difference at this precision |

An important caution is that the 200-replicate Stage A pilot suggested a larger shallow-slope penalty than the 500-replicate Stage B run supports. The Stage B paired interval includes zero for shallow-slope versus correct weighting. Therefore the pilot's apparent shallow-slope effect is **not a robust conclusion** and must not be advertised as one. A misspecified transition can accidentally compensate other finite-sample biases for one chosen estimand.

## Scientific interpretation and stop/go

**GO:** report the preregistered Stage B primary contrast as a reproducible, independently audited result within the fixed-H0 synthetic experiment.

**NO-GO:** any claim that this is a realistic LVK selection model, that inverse weighting universally fixes catalogue incompleteness, or that the study measures real (H_0). The study still fixes one true H0, one event count, one host population, one noise setup and a declared known inclusion curve. It omits realistic detector response, event-rate marginalization, peculiar velocities, uncertain redshifts and survey-level population inference.

Before a broader methods or publication claim, run predeclared calibration over multiple true H0 values, event counts, completeness levels and host-population shapes; compare grid/bandwidth sensitivity; and validate against a realistic external likelihood/selection benchmark. The new H0-dependent noise module is a separate synthetic prototype and has its own CI gate; it is not used in this locked Stage B run.

## Reproduction

Download and unzip artifact ID 12864010259, then run:

`python scripts/audit_selection_misspecification.py /path/to/unpacked/artifact`

`python scripts/analyze_selection_misspecification_pairs.py /path/to/unpacked/artifact --output /tmp/stage_b_pair_analysis.json`

The first command reconstructs posterior grids and verifies the saved primary result. The second independently recomputes aggregate performance, the locked primary estimand, and the secondary paired contrasts from the raw CSV and checks the manifest hashes. The supplementary bootstrap seed is fixed for reproducibility.
