# Dark-siren literature positioning update (31 October 2026)

**Purpose:** update the project's novelty boundary and next-stage validation plan against recent peer-reviewed work.  
**Status:** literature and research-design note; no new numerical experiment or astrophysical result is claimed.

## Executive assessment

The field has moved quickly. Gaussian-process reconstruction of incomplete line-of-sight galaxy catalogues is now a published *Physical Review D* result (13 October 2026), and variance-completion methods have been applied to GLADE+ and gwcosmo, including an O3 analysis, in a 2026 JCAP paper. A manuscript that presents generic GP reconstruction, clustering-aware completion, or a simple comparison of homogeneous versus multiplicative completion as its central novelty would therefore be poorly positioned.

The defensible opportunity for this repository is narrower and more testable: **a controlled, independently auditable study of when dark-siren posterior calibration and robustness fail under interacting catalogue-selection, host-population, and weighting misspecification, with explicit propagation of uncertainty in the incompleteness correction.** This is a proposed research direction, not a demonstrated result. It must be benchmarked against current completion methods before any novelty claim is made.

The linear finite-population bootstrap branch is a separate statistical-methods pilot. It can strengthen reproducibility and design-based reasoning, but it is not itself evidence that a nonlinear dark-siren posterior is calibrated.

## Recent literature that changes the comparison set

### 1. Gaussian-process reconstruction of incomplete line-of-sight catalogues

Tagliazucchi, Gair, Barbieri and Moresco, “Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs,” *Physical Review D* **114**, 083523 (published 13 October 2026), DOI: https://doi.org/10.1103/vmzf-1wxq.

The paper models the line-of-sight galaxy redshift distribution with a Gaussian process, fitting incomplete catalogues with a hierarchical likelihood that accounts for survey selection. On mock catalogues to redshift 0.4, it reports average constraint-precision improvements over homogeneous completion of 23% at 24% catalogue completeness and 37% at 8% completeness, with a largest reported configuration-specific improvement of 66%.

**Implication:** GP-based clustering reconstruction is an existing method and must be included as a comparator if the repository makes claims about catalogue completion. The paper also identifies a useful unresolved issue: its conclusion notes that a fully Bayesian treatment jointly marginalizing over GP realizations and cosmological parameters remains to be developed. That is a possible future direction, but this repository must not claim to have solved it without implementing and validating such a model.

### 2. Variance completion and large-scale-structure-informed corrections

Dalang, Fiorini and Baker, “Large scale structure prior knowledge in the dark siren method,” *Journal of Cosmology and Astroparticle Physics* **01** (2026) 034, DOI: https://doi.org/10.1088/1475-7516/2026/01/034.

This work describes practical variance completion, compresses predictions for missing galaxies into a correction relative to homogeneous completion, incorporates the correction into line-of-sight calculations in gwcosmo, and reports applications involving GLADE+, GW190814, and O3 data.

**Implication:** a credible benchmark should compare the repository's catalogue-completion assumptions with homogeneous, multiplicative, and structure-informed/variance-completion baselines where reproducible inputs and compatible implementations are available. A synthetic oracle alone is not an external validation.

### 3. Clustering of invisible host galaxies

Dalang and Baker, “The clustering of dark sirens' invisible host galaxies,” *Journal of Cosmology and Astroparticle Physics* **02** (2024) 024, DOI: https://doi.org/10.1088/1475-7516/2024/02/024.

This paper studies how the clustering structure of missing host galaxies affects completion schemes and contrasts homogeneous and multiplicative approaches.

**Implication:** missingness mechanisms should be tied to spatial/redshift structure, not only to a scalar catalogue completeness fraction. At minimum, the simulation design should distinguish missing completely at random, redshift-dependent selection, and structure-dependent selection, while preserving a known truth for calibration.

### 4. Computational scaling and pipeline validation

Papadopoulos, Chapman-Bird, Gray, Messenger and Bertheas, “Scalable dark siren cosmology with gwcosmo: GPU acceleration, validation, and systematics,” *Physical Review D* **114**, 063528 (published 15 September 2026).

The paper describes GPU acceleration and validation/systematics work in gwcosmo, motivated by growing event catalogues and the computational demands of hierarchical inference.

**Implication:** runtime and numerical stability matter, but speed is not a substitute for statistical validity. If this project later integrates an external inference pipeline, it should record software versions, configuration, numerical tolerances, and independent reproduction checks.

### 5. Bright-galaxy subsets and weighting robustness

Naveed, Turski and Ghosh, “Dark standard siren cosmology with bright galaxy subsets,” *Monthly Notices of the Royal Astronomical Society* **548** (2026), stag605.

The paper includes a uniform-weighting robustness comparison against luminosity-based weighting across brightness-selected subsets.

**Implication:** host weighting should be treated as a model choice with sensitivity analysis, not as a universally correct prescription. A useful experiment compares correctly specified, misspecified, and uniform weights while separately reporting bias, interval coverage, interval width, and failure rates.

## Research-positioning matrix

| Topic | Existing literature baseline | What this repository may test | Claim boundary |
|---|---|---|---|
| Homogeneous completion | Established baseline | Reproduce as a reference arm | Not novel |
| Multiplicative completion | Established baseline | Matched comparison under controlled selection | Not novel |
| Clustering/variance completion | Published method and applications | Compare when compatible code/data can be reproduced | Do not claim priority |
| GP line-of-sight reconstruction | Published PRD 2026 method | Independent benchmark or uncertainty-propagation extension | Do not present generic GP completion as new |
| Host-galaxy weighting | Existing luminosity/stellar-mass and uniform-weighting studies | Calibration under weight misspecification and selection interactions | Novelty requires a precise, literature-checked estimand |
| Posterior calibration under joint misspecification | A potential cross-cutting evaluation target | Repeated-simulation coverage, bias, interval score, rank calibration, and failure accounting | Proposed until executed and audited |
| Finite-population bootstrap | Wang, Peng & Kim (2022) theory and current pilot implementation | Reproducible linear design-based benchmark | Does not validate nonlinear dark-siren posterior intervals |

## Recommended high-value research question

> Under what combinations of catalogue-selection mechanism, missing-host clustering, and host-weight misspecification do dark-siren posterior summaries become miscalibrated, and which completion/weighting strategy remains reliable when uncertainty in the reconstructed host distribution is propagated rather than plugged in?

This question is more defensible than claiming a new catalogue-completion algorithm before the necessary comparisons exist.

### Minimum confirmatory design

1. **Lock the estimands:** posterior-mean error, bias, RMSE, nominal 50/68/90% interval coverage, interval width, proper interval score, numerical failure probability, and calibration diagnostics.
2. **Separate data-generating factors:** catalogue completeness; random versus redshift-dependent versus structure-dependent missingness; host-rate weighting law; weight-model misspecification; GW distance uncertainty/localization; and number of events.
3. **Use paired simulations:** every method receives the same generated GW events and observed catalogue for a replicate. Preserve replicate-level seeds and truth.
4. **Compare meaningful baselines:** homogeneous completion, multiplicative completion, a structure-informed method, and GP reconstruction if an independently reproducible implementation is available. Include oracle/full-catalogue and deliberately misspecified controls.
5. **Propagate reconstruction uncertainty:** compare plug-in summaries with posterior/ensemble marginalization over plausible missing-host distributions. Do not label a plug-in GP median as full Bayesian marginalization.
6. **Use independent audits:** regenerate truth, selections, seeds, posterior grids, summaries, and intervals from raw replicate records. Serialize every attempted fit and every internal resampling failure.
7. **Plan Monte Carlo precision before running:** determine replicate counts from target uncertainty for coverage and paired contrasts; report binomial or paired Monte Carlo uncertainty, not only point estimates.
8. **Run a held-out validation stage:** freeze method settings on development simulations, then evaluate on prespecified stress regimes not used for tuning.

## Immediate repository actions and release gates

### Active PR #17: linear finite-population bootstrap pilot

- Keep the PR draft and unmerged until exact-head CI and review are complete.
- Preserve the v1 artifact unchanged as a historical pilot.
- Implement the documented v2 internal-attempt accounting before calling the linear benchmark fully audited.
- Add multiple independently generated finite frames per distribution regime before drawing distribution-level conclusions.
- Do not call the pilot confirmatory or claim a method winner from 100 repetitions per cell.

### Main nonlinear dark-siren programme

- Keep the synthetic nature of current results explicit.
- Add this literature update to the manuscript's related-work and limitations sections before submission.
- Do not claim novelty for GP reconstruction, variance completion, or clustering-aware catalogue completion.
- The next substantive scientific contribution should be a preregistered, matched-method calibration benchmark with uncertainty propagation and a complete independent audit.
- External-data validation remains a separate gate; no real LVK measurement is implied by the synthetic studies.

## References

1. Tagliazucchi, M., Gair, J., Barbieri, R. & Moresco, M. (2026). Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs. *Physical Review D*, 114, 083523. https://doi.org/10.1103/vmzf-1wxq
2. Dalang, C., Fiorini, B. & Baker, T. (2026). Large scale structure prior knowledge in the dark siren method. *JCAP*, 01, 034. https://doi.org/10.1088/1475-7516/2026/01/034
3. Dalang, C. & Baker, T. (2024). The clustering of dark sirens' invisible host galaxies. *JCAP*, 02, 024. https://doi.org/10.1088/1475-7516/2024/02/024
4. Papadopoulos, A., Chapman-Bird, C. E. A., Gray, R., Messenger, C. & Bertheas, T. (2026). Scalable dark siren cosmology with gwcosmo: GPU acceleration, validation, and systematics. *Physical Review D*, 114, 063528.
5. Naveed, K., Turski, C. & Ghosh, A. (2026). Dark standard siren cosmology with bright galaxy subsets. *MNRAS*, 548, stag605.
6. Wang, Z., Peng, L. & Kim, J. K. (2022). Bootstrap inference for the finite population mean under complex sampling designs. *JRSS Series B*, 84(4), 1150–1174. https://doi.org/10.1111/rssb.12506
