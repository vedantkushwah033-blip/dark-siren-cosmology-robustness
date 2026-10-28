# Preliminary literature and novelty-risk memo — 2026-10-28

## Status

This is a focused first-pass literature comparison, not a systematic review and not a novelty claim. It records why broad claims about solving galaxy-catalogue incompleteness or improving dark-siren H0 inference would be unsupported by the current synthetic model.

## Closest relevant work inspected

### 1. Gray, Messenger & Veitch (2022): direction-dependent catalogue completeness

**Reference:** Gray, R., Messenger, C. & Veitch, J. (2022), “A pixelated approach to galaxy catalogue incompleteness: improving the dark siren measurement of the Hubble constant,” *Monthly Notices of the Royal Astronomical Society* 512, 1127–1140. DOI: 10.1093/mnras/stac366.

Publisher page: https://academic.oup.com/mnras/article/512/1/1127/6526882

**Relevance:** This work models spatially varying catalogue completeness in a dark-siren analysis and demonstrates its effect on H0 inference. Catalogue incompleteness itself, and the idea that its modelling affects H0 inference, are not novel contributions available to claim for this project.

**Difference from current project:** The current repository is a simplified synthetic statistical diagnostic, not a sky-pixelated catalogue-completeness pipeline. That difference is a scope difference, not automatically a novel contribution.

### 2. Gair et al. (2023): methodological guide to galaxy-catalogue dark sirens

**Reference:** Gair, J. R. et al. (2023), “The Hitchhiker’s Guide to the Galaxy Catalog Approach for Dark Siren Gravitational-wave Cosmology,” *The Astronomical Journal* 166, 22. DOI: 10.3847/1538-3881/acca78.

Accessible version: https://eprints.gla.ac.uk/334400/

**Relevance:** This paper explains the galaxy-catalogue dark-siren method and highlights how inconsistent mock-data generation, Bayesian implementation, and selection-effect modelling can create spurious bias. The current project’s emphasis on consistency audits is relevant to this existing methodological concern, so the manuscript must explain what its tests add beyond rechecking known implementation safeguards.

**Difference from current project:** A small toy model may be useful for isolating a mechanism, but it cannot be described as validating the general dark-siren framework without a much stronger connection to the full model and a carefully defined new estimand.

### 3. Datrier & Hendry (2026): robust catalogue-completeness test

**Reference:** Datrier, L. E. H. & Hendry, M. A. (2026), “Implementing a robust test of galaxy catalogue completeness for dark siren measurements of the Hubble constant,” *Monthly Notices of the Royal Astronomical Society* 549(1), stag828. DOI: 10.1093/mnras/stag828.

Publisher page: https://academic.oup.com/mnras/article/549/1/stag828/8666390

Preprint: https://arxiv.org/abs/2502.14164

**Relevance:** This work applies a statistical test of galaxy-catalogue completeness in the gwcosmo pipeline. It is especially important to distinguish the current project’s proposed study of diagnostic distinguishability under controlled misspecification from a claim to develop or improve a catalogue-completeness test.

**Difference from current project:** The current project has not shown that its diagnostic improves completeness estimation or H0 inference in gwcosmo or another realistic pipeline.

### 4. Leyde, Baker & Enzi (2026): reconstruction of incomplete galaxy catalogues

**Reference:** Leyde, K., Baker, T. & Enzi, W. (2026), “Cosmic cartography. Part II. Completing galaxy catalogs for gravitational-wave cosmology,” *Journal of Cosmology and Astroparticle Physics* 01 (2026) 013. DOI: 10.1088/1475-7516/2026/01/013.

Publisher page: https://doi.org/10.1088/1475-7516/2026/01/013

**Relevance:** This work reconstructs the underlying galaxy field and magnitude distribution, with validation on simulated data. It reinforces that physically motivated catalogue completion and spatial structure are active, developed research areas.

**Difference from current project:** The current project does not reconstruct a galaxy field and does not model spatial correlations. It should not imply that it offers a competing catalogue-completion method.

### 5. Tagliazucchi et al. (accepted 2026): Gaussian-process reconstruction

**Reference:** Tagliazucchi, M., Gair, J., Barbieri, R. & Moresco, M. (2026), “Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs,” *Physical Review D*, accepted 29 September 2026. DOI: 10.1103/vmzf-1wxq.

Publisher page: https://journals.aps.org/prd/accepted/10.1103/vmzf-1wxq

**Relevance:** The paper describes a Gaussian-process approach to line-of-sight reconstruction of incomplete galaxy catalogues, using a hierarchical Bayesian likelihood and explicitly accounting for the survey selection function. The repository’s toy model is substantially narrower and should not claim to improve this kind of reconstruction.

**Difference from current project:** The current project is not a GP catalogue reconstruction or a selection-function-aware hierarchical analysis. Its possible contribution must be a distinct, explicitly defined statistical diagnostic or a reproducibility finding.

## Preliminary novelty-risk assessment

| Proposed claim | Risk | Reason | Defensible alternative |
|---|---|---|---|
| “We solve galaxy-catalogue incompleteness in dark-siren cosmology.” | Extremely high | Established catalogue-completeness and reconstruction methods already exist; the current model is a simplified mixture. | “We study specified consequences of host-mixture misspecification in a controlled synthetic model.” |
| “We improve real H0 inference.” | Extremely high | No realistic GW event likelihood, detector-selection normalization, or validated reference-pipeline result has been demonstrated. | Make no real-data performance claim until a reference-pipeline comparison exists. |
| “We introduce a new universal diagnostic for misspecification.” | High | The proposed diagnostic-comparison question is not yet demonstrated to be novel or identifiable. | Define a narrow diagnostic, test whether misspecification classes can actually be distinguished, and report conditions where they cannot. |
| “We provide an independently auditable synthetic study of interactions between host and distance-noise misspecification.” | Plausible but unproven | This could be a defensible narrow contribution if the factorial design, independent audit, uncertainty analysis, and literature comparison are completed. | Present as a model-conditional simulation study, with the exact estimand and failure regimes stated. |
| “We fixed a mixture-likelihood bug and audited affected artifacts.” | Evidence-based software contribution if fully documented | The repository contains a correction, regression tests, and a reconciliation report; current evidence must be tied to exact commits and artifacts. | Report the precise defect, tests, impacted outputs, and correction provenance without implying a new astrophysical method. |

## Recommended primary question

A defensible candidate for further investigation is:

> Under a prespecified synthetic host-mixture model, can a declared set of calibration and posterior-predictive diagnostics distinguish host-population misspecification from distance-noise misspecification at finite sample size, and in which conditions are the two error sources observationally indistinguishable to those diagnostics?

This question is only worth pursuing if the project:
1. Defines each diagnostic and its reference distribution mathematically.
2. Includes host-only, noise-only, and joint-misspecification conditions using common observations where pairing is intended.
3. Quantifies diagnostic classification or distinguishability performance with uncertainty, rather than relying only on visual plots.
4. Includes negative controls and cases where the diagnostics are expected to fail.
5. Prespecifies the primary metric and simulation budget.
6. Shows that the result is not a restatement of existing Bayesian calibration or posterior-predictive-checking literature.
7. Remains explicitly conditional on the simplified simulator unless an external reference-pipeline study is completed.

Do not add a classifier merely to make the project sound more advanced. If the diagnostics cannot identify the error source from the available summaries, that non-identifiability may itself be a valid result, provided it is defined and established rigorously.

## Required next literature work

- Read the full methods and results of each cited paper, not just the abstract.
- Add relevant work on posterior predictive checks, Bayesian calibration, simulation-based calibration, misspecified Bayesian models, and diagnostic identifiability.
- Record explicit overlaps and differences in a claim-by-claim matrix.
- Check for additional papers and software released in 2026.
- Do not use “first,” “novel,” or “state of the art” in a manuscript unless the completed review supports the exact wording.

## Conclusion

The current strongest defensible direction is a rigorous, reproducible, narrowly scoped statistical simulation study—not a claim to solve dark-siren catalogue incompleteness. A stronger contribution might emerge from a well-defined diagnostic-identifiability result, but its novelty and usefulness remain unproven and require further work.
