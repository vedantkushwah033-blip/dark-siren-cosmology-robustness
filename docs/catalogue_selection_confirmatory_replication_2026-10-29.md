# Confirmatory catalogue-selection replication — 29 October 2026

## Status

**The three-condition replication completed successfully and all three artifact auditors passed.** This is a larger follow-up to the 50-replicate stress screen. It is still a synthetic, low-redshift Euclidean toy study, not an LVK analysis or an estimate of the real Hubble constant.

- Workflow: [catalogue-selection confirmatory replication, run #40819870231](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40819870231)
- Source commit used for all conditions: `9aadda0ffde113ebb2584aa5425424dbe154dd2b`
- Each condition: 250 paired replicates, 4,000 parent galaxies per replicate, 30 detected events, fixed true H0 = 70 km s⁻¹ Mpc⁻¹, 81-point H0 grid over [50,90], 401-point redshift grid over [0.01,0.40].
- Six arms per replicate: analytic truth, oracle KDE, catalogue only, correct weighting, shifted transition, shallow slope.
- Each auditor reconstructed all 1,500 raw rows, six arm aggregates, and 1,500 posterior grids. All three reported `PASS`.
- I additionally recomputed bias, RMSE, coverage counts, the paired primary contrast, and the 10,000-resample percentile bootstrap from each downloaded CSV. Recomputed point estimates matched the stored summaries to floating-point precision; raw CSV and summary JSON byte counts and SHA-256 digests matched the manifest.

## Primary result

The prespecified paired contrast is
`D_r = |posterior_mean_catalogue_only - H0_true| - |posterior_mean_correct_weighting - H0_true|`.
Positive values mean the correct-weighting arm has smaller absolute posterior-mean error. Confidence intervals below are percentile-bootstrap 95% intervals using 10,000 paired resamples.

| True catalogue inclusion curve | Catalogue-only bias | Catalogue-only RMSE | Catalogue-only 68% coverage | Correct-weighting bias | Correct-weighting RMSE | Correct-weighting 68% coverage | Mean paired D (km s⁻¹ Mpc⁻¹) | Bootstrap 95% interval |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Floor 0.05, midpoint 0.22, width 0.05 | −10.617 | 11.040 | 12/250 (4.8%) | +1.185 | 4.095 | 164/250 (65.6%) | +7.436 | [6.798, 8.053] |
| Floor 0.25, midpoint 0.22, width 0.05 | −4.717 | 5.898 | 95/250 (38.0%) | +0.938 | 4.056 | 174/250 (69.6%) | +2.246 | [1.711, 2.764] |
| Floor 0.10, midpoint 0.22, width 0.08 | −4.994 | 6.153 | 84/250 (33.6%) | +1.304 | 4.154 | 177/250 (70.8%) | +2.354 | [1.779, 2.923] |

The primary paired interval excludes zero in all three tested conditions. The absolute-error reduction varies materially by inclusion regime: the low-floor case has the largest catalogue-only bias and the largest weighting benefit, while the high-floor case has a smaller but still positive contrast. This supports the limited conclusion that, under these declared synthetic conditions, omitting the known inclusion correction can substantially distort inference and the effect depends on the inclusion curve.

## Context from the control arms

| Inclusion curve | Analytic-truth RMSE | Oracle-KDE RMSE | Correct-weighting RMSE | Shifted-transition RMSE | Shallow-slope RMSE |
|---|---:|---:|---:|---:|---:|
| Floor 0.05, width 0.05 | 3.938 | 3.970 | 4.095 | 3.846 | 4.201 |
| Floor 0.25, width 0.05 | 3.964 | 3.999 | 4.056 | 4.052 | 3.873 |
| Floor 0.10, width 0.08 | 4.043 | 4.110 | 4.154 | 3.931 | 3.864 |

Correct weighting does not dominate every alternative on every metric. In particular, a misspecified selection arm can have lower RMSE than correct weighting in a finite simulation because bias and variance can compensate. Analytic truth and oracle KDE are diagnostic controls that are not available as ordinary estimators in a real incomplete-catalogue problem. These comparisons must not be presented as evidence of universal superiority.

## Interpretation and boundaries

1. **Replicated synthetic pattern:** in all three 250-replicate conditions, the catalogue-only arm is more negatively biased and has larger RMSE than the correct-weighting arm. The paired absolute-error contrast is positive with a bootstrap interval excluding zero in each selected condition.
2. **Regime dependence:** the mean absolute-error reduction is about 7.44, 2.25, and 2.35 km s⁻¹ Mpc⁻¹ across the low-floor, high-floor, and broad-transition conditions. These are simulation-specific effects, not physical uncertainty estimates for H0.
3. **Coverage:** catalogue-only coverage is poor in these settings (4.8%, 38.0%, 33.6%), whereas the correctly weighted arm gives 65.6%, 69.6%, and 70.8%. Coverage is estimated over 250 replicates; it is not proof of exact calibration, and the posterior intervals omit many real-world uncertainties.
4. **Known-law assumption:** the correct-weighting arm uses the true inclusion curve. Real surveys generally require estimating or externally validating the selection function; this experiment does not solve that problem.
5. **Model limitations:** simplified host population, fixed H0 and event count, toy distance/noise model, and no realistic detector response, network selection, event-rate marginalization, peculiar velocities, redshift uncertainty, or full hierarchical likelihood.
6. **Reproducibility is not realism:** independent reconstruction and checksums show that the saved toy results are internally consistent with the implemented specification. They do not establish that the specification accurately represents real gravitational-wave data.

## Artifact provenance

All three artifacts are retained by GitHub Actions for 90 days from upload:

- Low-inclusion floor: [artifact 12873663060](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40819870231/artifacts/12873663060) — CSV SHA-256 `39e4fddc6c2fa6c49700fecfb9adce27a5f12ba8cf2917e04a0747b007a617c8`; summary SHA-256 `7a8aa44a8a0580f302b9ce681ae5c203be7815e381f719127224343dc1541271`.
- High-inclusion floor: [artifact 12873737617](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40819870231/artifacts/12873737617) — CSV SHA-256 `add705e823021dce6fb69aafe9b8ef68aa9b547802b28de88855923b9e2ba28f`; summary SHA-256 `a49d4a66b0f0606f1cc4c59b92e35e8e70933e80e55460e6acdd41996fa18a41`.
- Broad transition: [artifact 12873982241](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40819870231/artifacts/12873982241) — CSV SHA-256 `d7e2c8028780ed69ec58f1b4f785b727bb5f0f59753dda756c293d0976a94d65`; summary SHA-256 `eb5f52758159b9b72830e47bdd7a74275283b25bde57af65fc2db61c7cc4bfa5`.

## Release decision

**GO:** include this as a prespecified, clearly scoped synthetic replication in the internal manuscript and use it to motivate further selection-function work.

**NO-GO:** claims of a new H0 measurement, real LVK validation, a realistic end-to-end dark-siren pipeline, universal weighting superiority, or external-method superiority. The exact CHIMERA v2.2 benchmark remains blocked until compatible inputs are obtained, versioned, and run with provenance. Domain-expert review is still needed before external submission.
