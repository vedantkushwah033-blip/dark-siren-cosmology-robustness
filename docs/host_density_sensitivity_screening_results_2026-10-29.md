# Host-density sensitivity screening: verified results and limits

**Date:** 29 October 2026  
**Status:** completed exploratory screen; no real-data inference and no confirmatory claim.

## Provenance

- GitHub Actions run: [#40797592794](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40797592794) — both campaign and exact-density diagnostic jobs succeeded.
- Campaign source commit: `5bbe5a434878af396e8c58e3ddfe7af0356e29ed`.
- Campaign artifact ID: `12865051001`; archive digest `sha256:bfd9132783c354a929f0e633e0b0536361e8427cbb06ae1f94e451aae48470e8`.
- Oracle screening artifact ID: `12865285154`; archive digest `sha256:bf2d5f99faf1fe5cf8e6a028202463ee4546d51ec95dcf5c7681fbcf18d9a884`.
- Campaign design: 13 configurations × 40 paired replicates; seeds 20261029–20261068, 6,000 parent galaxies, 40 detected events, fixed H0=70 km s⁻¹ Mpc⁻¹, Gaussian distance-noise scale 140 Mpc, detection threshold 250 Mpc.
- Integrity review of earlier campaign artifact ID `12865052874`: all manifest-listed campaign files passed SHA-256 and byte-count checks. The new `scripts/audit_host_density_campaign.py` was then run separately for all 13 configurations with 40 replicates × three arms each. All 1,560 fit rows passed independent generator replay, reflected-KDE reconstruction, selected-posterior reconstruction, and raw-summary/aggregate checks. The unified clean reproduction and one-shot audit passed in [run #40802415000](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40802415000). Canonical artifact ID `12866753211`, ZIP digest `sha256:0ec51d1c96e9f86d7f1ef2a6dbcb9e57a0c67fa5be2459617fcf51ccbc5565fe`. Its audit JSON reports PASS: 13/13 configurations; 40 replicates per configuration; 1,560/1,560 rows independently reconstructed; all 14 manifest files verified; 403 summary aggregate fields checked; maximum aggregate absolute difference 0.0; and no failed rows. One borderline coverage flag differed in oracle/`bw_0.5`/replicate 25 because an endpoint was only `9.37e-5` from the truth. It is explicitly counted as boundary-sensitive and does not trigger failure because it is inside the prespecified endpoint tolerance.
- The older compact two-replicate smoke audit remains a fast cross-implementation test but no longer defines the full scope of local reconstruction. All results remain exploratory and dependent on the declared synthetic model.

## Primary fixed setting: weighted Scott bandwidth

Bias is posterior-mean minus truth; RMSE and bias are in km s⁻¹ Mpc⁻¹. The known-inclusion arm uses the true declared logistic inclusion curve.

| Arm | Bias | RMSE | 68% coverage | Density L1 error |
|---|---:|---:|---:|---:|
| Complete-parent oracle KDE | +1.123 | 3.936 | 52.5% | 0.058 |
| Incomplete catalogue, uncorrected | −8.634 | 9.165 | 10.0% | 0.537 |
| Known inclusion correction | +1.180 | 3.896 | 55.0% | 0.082 |

The paired absolute-error difference for known-inclusion correction minus oracle was −0.057, bootstrap 95% interval [−0.199, +0.083]. The interval includes zero, and 40 replicates provide low precision for coverage. The defensible result is that catalogue-only inference is strongly distorted in this toy and its density error is much larger; the corrected arm tracks oracle closely at the default bandwidth. The screen does not resolve small differences between those two arms.

## Bandwidth sensitivity

An explicit bandwidth factor (f) means (h=f,SD_w) for the original sample; default Scott is (h=SD_w n_{eff}^{-1/5}).

| Bandwidth | Catalogue-only bias | RMSE | 68% coverage | Density L1 | Known-inclusion bias | RMSE | 68% coverage | Density L1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Scott default | −8.634 | 9.165 | 10.0% | 0.537 | +1.180 | 3.896 | 55.0% | 0.082 |
| (f=0.05) | −8.786 | 9.297 | 7.5% | 0.537 | +0.951 | 3.787 | 57.5% | 0.075 |
| (f=0.10) | −8.757 | 9.271 | 7.5% | 0.536 | +0.990 | 3.807 | 55.0% | 0.055 |
| (f=0.20) | −8.625 | 9.157 | 10.0% | 0.537 | +1.153 | 3.883 | 55.0% | 0.076 |
| (f=0.50) | −7.400 | 8.122 | 20.0% | 0.531 | +2.137 | 4.440 | 55.0% | 0.231 |
| (f=0.80) | −4.983 | 6.202 | 40.0% | 0.488 | +2.954 | 4.990 | 55.0% | 0.297 |
| (f=1.20) | −1.736 | 4.238 | 62.5% | 0.413 | +3.446 | 5.330 | 50.0% | 0.313 |

Heavy smoothing makes the catalogue-only posterior mean look less biased at this fixed truth, but the inferred host density remains substantially wrong. It also increases H0 RMSE in the known-inclusion arm. This is an example of bias compensation: posterior error can improve by accident without the nuisance distribution becoming accurate. It is a reason to report density diagnostics and prespecify bandwidth rather than choosing it by apparent H0 performance at one truth value.

The paired known-inclusion-versus-oracle absolute-error differences were: Scott −0.057 [−0.199, +0.083]; (f=0.05) −0.092 [−0.241, +0.055]; (f=0.10) −0.085 [−0.231, +0.059]; (f=0.20) −0.073 [−0.209, +0.065]; (f=0.50) +0.014 [−0.090, +0.117]; (f=0.80) +0.015 [−0.056, +0.083]; and (f=1.20) +0.008 [−0.032, +0.045]. All these paired intervals contain zero at this sample size.

## Grid stability

Compared with the 161 × 801 default grid:

| Alternative grid | Maximum absolute posterior-mean difference | Maximum absolute 68% interval-width difference |
|---|---:|---:|
| Coarse, 81 × 401 | 0.000851 km s⁻¹ Mpc⁻¹ | 0.0300 km s⁻¹ Mpc⁻¹ |
| Fine, 201 × 1001 | 0.000103 km s⁻¹ Mpc⁻¹ | 0.00476 km s⁻¹ Mpc⁻¹ |

These paired values indicate that grid discretization is not driving the main effects in this particular setup. They are not proof of asymptotic convergence or a guarantee for sharper likelihoods, different priors, or other grids.

## Inclusion-law misspecification

| Assumed logistic curve | Bias | RMSE | 68% coverage | Density L1 | Paired absolute-error difference vs oracle (95% interval) |
|---|---:|---:|---:|---:|---:|
| Midpoint 0.18 (truth 0.22) | +0.750 | 3.801 | 50.0% | 0.118 | −0.118 [−0.297, +0.064] |
| Midpoint 0.26 (truth 0.22) | +0.709 | 3.719 | 55.0% | 0.155 | −0.252 [−0.475, −0.043] |
| Width 0.10 (truth 0.05) | −2.536 | 4.245 | 60.0% | 0.259 | +0.293 [−0.714, +1.232] |
| Width 0.025 (truth 0.05) | +2.578 | 4.669 | 42.5% | 0.227 | +0.613 [+0.201, +1.021] |

The midpoint 0.26 arm looks better for one fixed-truth error contrast; that does not imply the wrong midpoint is preferable. The width configurations yield noticeably different biases and coverages; with 40 replicates their uncertainty is broad, although the width-0.025 paired contrast interval excludes zero. That is a candidate for prospective replication, not a general result.

## Exact-density oracle screening

The same workflow ran 300 exact-density replicates for each of four truth regimes on a coarser 401-point redshift / 81-point H0 grid. The output artifact stores aggregate summaries only, not raw rows and seeds, so it is CI-executed but not independently replayed row-by-row.

| Truth design | Bias ± MC SE | RMSE | 68% coverage | 90% coverage | Ten-bin PIT chi-square p-value |
|---|---:|---:|---:|---:|---:|
| Prior-drawn H0 | −0.221 ± 0.168 | 2.909 | 70.0% | 91.0% | 0.105 |
| Fixed H0=60 | +0.562 ± 0.146 | 2.584 | 68.3% | 91.0% | 0.276 |
| Fixed H0=70 | +0.724 ± 0.181 | 3.219 | 69.7% | 90.7% | 0.185 |
| Fixed H0=80 | +0.543 ± 0.189 | 3.306 | 69.7% | 91.3% | 0.324 |

The prior-drawn result is broadly compatible with calibration on this screening run. Fixed-truth posterior-mean bias is positive at all three values but small relative to RMSE. Fixed-truth PIT p-values are descriptive, not tests under the same uniformity criterion as prior-predictive ranks. The earlier 2,000-replicate local table is a different dataset and was not reproduced exactly by this coarser screen.

## What is supported—and what is not

**Supported:** in the declared toy setting, ignoring the inclusion curve damages host-density and fixed-truth H0 summaries at default bandwidth; the true-inclusion control tracks the complete-parent oracle at the default and narrow bandwidths; tested grid refinement has a negligible effect; bandwidth and selection-law assumptions materially change results.

**Not supported:** universal advantages of inverse-inclusion weighting; choosing a bandwidth based on H0 performance at a single truth value; realistic GW/LVK selection validity; real-H0 inference; or superiority to GP reconstruction or operational pipelines. The study fixes one truth, one parent/host law, one event count, a noise scale, and a threshold. A stronger follow-up should prespecify multiple H0 truth values, event counts, inclusion floors/transition shapes, and a larger replicate count based on paired precision.
