# Literature and scientific-positioning audit

**Date:** 29 October 2026. This is a focused literature audit to establish a defensible novelty boundary, not a systematic review.

## Closely relevant research

| Source | Relevant contribution | Implication for this project |
|---|---|---|
| Barbieri, Kalomenopoulos, Tagliazucchi, Gair & Khochfar, “Not all those missing are lost: leveraging galaxy clustering in incomplete catalogs to unleash dark sirens cosmology,” [arXiv:2608.30599](https://arxiv.org/abs/2608.30599) (31 Aug 2026). | Reconstructs missing galaxies using a spatial correlation-function model and sky pixels; mock studies report 2–4× tighter posteriors in selected 1–10% completeness regimes, depending on completeness and GW localization/uncertainty. | Unlike this method, the project uses a one-dimensional low-redshift toy, a known logistic inclusion law, and inverse-probability-weighted KDE. It has no sky map or spatial galaxy-clustering reconstruction; the models/data regimes are not directly comparable. |
| Ghosh & More, “Recovering Large-Scale Clustering of Missing Galaxies for Galaxy–Gravitational-Wave Cross-correlations,” [arXiv:2609.07348](https://arxiv.org/abs/2609.07348) (posted 7 Sep 2026). | Develops a clustering-preserving completion framework for missing galaxies in flux-limited surveys and tests it on simulated catalogues, with motivation from GW–galaxy cross-correlations and dark-siren cosmology. | This project neither reconstructs angular/large-scale clustering nor builds a completed galaxy catalogue; it tests one-dimensional host-density estimation and inverse-inclusion weighting in a deliberately simple toy likelihood. The estimands and models do not support a direct ranking. | Do not claim to be the first to reconstruct missing galaxies or preserve clustering, or to outperform clustering-informed completion. |
| Tagliazucchi, Gair, Barbieri & Moresco, “Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs,” Physical Review D 114, 083523 (published 13 Oct 2026), [DOI 10.1103/vmzf-1wxq](https://doi.org/10.1103/vmzf-1wxq), [publisher page](https://journals.aps.org/prd/abstract/10.1103/vmzf-1wxq), [arXiv:2608.24992](https://arxiv.org/abs/2608.24992). | Reconstructs a clustered galaxy redshift distribution from an incomplete catalogue with a selection-aware hierarchical likelihood and a Gaussian process. Its mock studies report improved H0 precision relative to homogeneous completion. | Catalogue incompleteness and reconstructed host density are already active research topics. Do not claim this project is the first to study incomplete galaxy catalogues, nonparametric reconstruction, or selection correction. The defensible distinction is a transparent controlled benchmark contrasting corrected KDE, inverse-inclusion weighting, known/unknown selection-law assumptions and posterior calibration under a deliberately simple model. |
| “Implementing a robust test of galaxy catalogue completeness for dark siren measurements of the Hubble constant,” Monthly Notices of the Royal Astronomical Society 549(1), stag828 (published 1 May 2026), [publisher page](https://academic.oup.com/mnras/article/549/1/stag828/8666390). | Implements a robust magnitude-completeness-limit test in the pixelated gwcosmo pipeline and assesses its effect on Hubble-constant inference, including GWTC-1/GWTC-3 catalogues. | Completeness correction is already being addressed in real-analysis pipelines. This synthetic project should not claim that practical completeness estimation has been solved. Testing known inclusion probabilities under controlled conditions is a narrower question. |
| Gray, Williams & Papadopoulos, “Expanding the scope of dark siren cosmology: Inferring the population properties of gravitational wave-hosting galaxies,” [arXiv:2609.04991](https://arxiv.org/abs/2609.04991), submitted 4 Sep 2026. | Describes a dark-siren methodology jointly inferring host-galaxy population properties with cosmological/population parameters using a sampled redshift prior and gwcosmo. | Host-population misspecification is an active methodological issue. This project's host weighting is an intentionally simplified control, not a competing full hierarchical host-population model. |
| Palmese et al. (2023) and Dalang & Baker (2024), discussed as prior completeness approaches in the 2026 MNRAS article linked above. | Establish earlier approaches to catalogue incompleteness and dark-siren inference. | Read the underlying papers directly before submitting a manuscript; position the project relative to their assumptions rather than presenting it as a first application. |
| CHIMERA v2.2, [tagged release](https://github.com/CosmoStatGW/CHIMERA/releases/tag/v2.2), tag commit fd35f897e36850d865c760c13325c5a46ebe94db. | The public release includes a changed likelihood API, restructured population modules and complete/incomplete galaxy-catalogue tutorials. | Do not use an older CHIMERA interface or legacy notebook as if it were v2.2 compatible. An external benchmark needs a pinned commit, compatible notebook/API, reproducible environment and exact data provenance. If an expected input file is unavailable or mismatched, the comparison is not run; another version’s file is not a silent substitute. |

## Contribution that can be defended now

The project is best described as a **reproducible synthetic robustness benchmark**, with several narrow methodological/software outputs:

1. It exposes a KDE bandwidth-convention error: SciPy's gaussian_kde bw_method argument is not a direct physical-bandwidth multiplier, especially after reflection augmentation. The corrected code computes a reflected weighted Gaussian sum with a stated physical bandwidth and a direct-sum regression test.
2. It separates catalogue-only inference, inverse-inclusion weighting under known selection, oracle KDE and an analytic-density control, so selection misspecification and finite-parent/KDE error can be compared.
3. It records raw replicate rows, seed pairing, saved posterior grids, SHA-256 manifests and auditors that reconstruct calculations instead of trusting internally consistent summaries.
4. It tests a separate discrete-host likelihood where distance-noise scale and detection efficiency both depend on H0, with independent seed replay and posterior reconstruction.

These are useful reproducibility and robustness contributions. They do not by themselves establish a new astrophysical result or a superior method to GP reconstruction, gwcosmo or CHIMERA.

## Safe claims for an abstract

- “We built and audited a controlled synthetic experiment to quantify how catalogue incompleteness and selection-law misspecification affect posterior summaries in a simplified dark-siren model.”
- “We corrected and directly regression-tested the physical KDE bandwidth convention used by the simulation.”
- “In the declared toy setup, ignoring catalogue incompleteness strongly degraded fixed-H0 error and interval coverage, while inverse-inclusion weighting under the exactly known inclusion law substantially reduced this particular error.”
- “A separate H0-dependent-noise toy likelihood passed deterministic simulation replay and posterior-reconstruction checks; this is internal software/model validation, not detector realism.”

Avoid claims that the project solved the Hubble tension, measured a new H0, is the first to handle incomplete galaxy catalogues, universally fixes inverse-probability weighting, validates LVK/GWTC or CHIMERA, or has run the CHIMERA v2.2 benchmark unless the exact inputs, commands, versions and logs have been inspected.

## Remaining literature/external validation

- Compare the GP reconstruction, spatial-correlation reconstruction and clustering-preserving catalogue-completion work (including Ghosh & More, arXiv:2609.07348) assumption by assumption, including the inference target, selection-aware likelihood, sky/line-of-sight information and mock-data settings.
- Read the methods and supplement of the 2026 `gwcosmo` completeness paper before comparing the project's known-probability toy with its apparent-magnitude-limit estimator. Also position it beside the older CHIMERA clustering-based incompleteness work and 2025 completeness/host-weighting literature.
- Read the 2026 host-population preprint and identify which missing-host and host-weight assumptions differ.
- Complete the CHIMERA v2.2 input inventory. This audit does not certify the presence of the previously requested 100-event SNR-25 mock or pcat_incomplete.h5; those exact inputs remain unverified.
- Claim an advantage only after a prespecified same-data, same-model comparison with uncertainty intervals and independent reproduction.

## Bottom line

**GO** for a narrowly framed synthetic methods/reproducibility report if exact-head CI and the regenerated campaign finish cleanly. **NO-GO** for real-data, astrophysical-validation, superiority or Hubble-tension claims. Recent published work makes precise claim discipline essential.
