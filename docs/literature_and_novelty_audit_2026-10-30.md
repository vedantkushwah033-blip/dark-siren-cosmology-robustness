# Literature and novelty audit: estimated inclusion and incomplete dark-siren catalogues

**Audit date:** 30 October 2026  
**Scope:** galaxy-catalogue incompleteness, selection-function estimation, host weighting, large-scale-structure reconstruction, and uncertainty propagation in dark-siren cosmology.  
**Status:** focused literature audit, not a systematic review. Search coverage is not exhaustive.

## Executive verdict

The broad claim that catalogue incompleteness, host weighting, selection-function misspecification, or nonparametric reconstruction is a new dark-siren problem is not defensible. Several directly relevant 2025–2026 papers address these subjects, including GP-based line-of-sight reconstruction, spatial/magnitude catalogue completion, empirical completeness-limit tests, and joint inference of host-population properties.

The current repository's estimated-inclusion pilot is best described as a **controlled statistical-methods benchmark**: it studies what happens when a retention probability is estimated from a separate parent-frame validation sample, then used for inverse-probability host weighting, under finite validation size and deliberate selection-model misspecification. That precise combination may be a useful experimental niche, but this audit does **not** establish that it is novel. Inverse-probability weighting and validation-sample uncertainty are established statistical ideas outside this application.

The project should not claim a new dark-siren estimator, a new cosmological measurement, or superiority over GP/large-scale-structure completion methods.

## Closest and most consequential literature

### 1. GP reconstruction of incomplete line-of-sight catalogues

Tagliazucchi, Gair, Barbieri & Moresco, *Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs*, Physical Review D 114, 083523 (published 13 October 2026), DOI: 10.1103/vmzf-1wxq.

The paper models the line-of-sight galaxy redshift distribution with a Gaussian process and explicitly accounts for the survey selection function. It reports improved precision over homogeneous completion on mock catalogues at 24% and 8% completeness. It also identifies fully Bayesian joint marginalization over the reconstructed redshift distribution and H0 as a computationally challenging future direction.

**Implication:** a broad proposal to use flexible/nonparametric reconstruction to preserve clustering information is already occupied. The current project's one-dimensional weighted KDE is not a substitute for this method, and comparisons must not imply otherwise.

### 2. Spatial and magnitude-aware catalogue completion

Leyde, Baker & Enzi, *Cosmic cartography. Part II. Completing galaxy catalogs for gravitational-wave cosmology*, JCAP 01 (2026) 013, DOI: 10.1088/1475-7516/2026/01/013.

This work reconstructs the underlying true galaxy field while respecting spatial correlations and infers the true magnitude distribution alongside the spatial distribution. It is validated using simulated data based on the Millennium simulation and provides the LIGHT codebase.

**Implication:** any future benchmark needs to acknowledge spatial clustering and magnitude-dependent selection. A redshift-only simulator cannot support claims about full three-dimensional survey completion.

### 3. Empirical catalogue-completeness testing

Datrier & Hendry, *Implementing a robust test of galaxy catalogue completeness for dark siren measurements of the Hubble constant*, MNRAS 549 (2026), stag828, DOI: 10.1093/mnras/stag828.

This work applies a statistical test to estimate the apparent-magnitude completeness limit of magnitude-redshift samples and reports changes in H0 inference for GWTC-1/GLADE.

**Implication:** estimating catalogue selection from data is already an active line of work. The present pilot's synthetic external validation frame is a different design, but not a general novelty claim.

### 4. Joint host-population and cosmological inference

Gray, Williams & Papadopoulos, *Expanding the scope of dark siren cosmology: Inferring the population properties of gravitational wave-hosting galaxies*, arXiv:2609.04991 (4 September 2026).

The preprint describes jointly inferring host-galaxy weighting/population properties, GW population properties, and cosmology, and reports an application using gwcosmo and GWTC-5.0-related results.

**Implication:** host weighting should not be treated as a nuisance that is known or solved by a single plug-in regression. A publication-quality study must distinguish uncertainty in the catalogue-retention mechanism from uncertainty in the astrophysical host-rate model.

### 5. Incompleteness and host-weighting stress tests

Borghi et al., *Echoes from the Dark: Galaxy Catalog Incompleteness in Standard Siren Cosmology*, arXiv:2509.18243 (2025).

This study examines catalogue incompleteness and host-weighting effects with a self-consistent framework in CHIMERA and mock BBH populations.

Alfradique, Bom & Castro, *Systematic bias in dark siren statistical methods and its impact on Hubble constant measurement*, Physical Review D 112, 063561 (2025), DOI: 10.1103/vd36-3mys.

This work studies biases from catalogue incompleteness and incorrect host weighting.

VanWyngarden et al., *How Low Can You Go: Constraining the Effects of Catalog Incompleteness on Dark Siren Cosmology*, arXiv:2511.04786 (2025).

This work studies how clustering and host-halo assumptions affect the completeness needed for useful dark-siren inference.

**Implication:** incompleteness × host-weighting sensitivity is not a novelty gap by itself.

### 6. Other large-scale-structure and cross-correlation methods

Leyde, Baker & Enzi (JCAP 2026) and Dalang, Fiorini & Baker, *Large scale structure prior knowledge in the dark siren method*, JCAP 01 (2026) 034, DOI: 10.1088/1475-7516/2026/01/034, discuss approaches that use large-scale-structure information to improve incomplete-catalogue inference.

Cheng & Gair, *Unified harmonic framework for dark siren cosmology*, accepted by Physical Review D in October 2026, develops a harmonic-space/cross-correlation perspective and clarifies its relationship to galaxy-catalogue methods.

**Implication:** this project should state explicitly that it is not a spatial reconstruction, variance-completion, or angular cross-correlation method.

## What is actually distinct in the current repository?

The current controlled experiment has a particular data-generating and validation design:

1. A simulated complete parent galaxy frame is generated.
2. Catalogue retention is generated from a known probability law.
3. A simple random validation sample is drawn from the parent frame independently of catalogue retention; both redshift and retention indicator are observed for validation units.
4. A logistic retention model is fitted using only those validation records.
5. Estimated inverse-inclusion weights are applied to the retained catalogue to construct a one-dimensional host-redshift density.
6. The resulting density is passed into a simplified selection-conditioned GW likelihood.
7. Oracle probabilities, catalogue-only inference, and estimated-probability inference are compared on paired simulated populations and GW observations.
8. A bootstrap utility now measures sensitivity to validation-record resampling conditional on the observed catalogue and GW data.

This design is not the same as estimating a magnitude completeness relation from a real survey, reconstructing a 3D galaxy field, or jointly inferring the host astrophysical population. It is a deliberately simplified external-validation experiment.

## Claims that are currently supportable

Subject to exact-head CI and artifact audits, the project may support claims such as:

- The code implements a reproducible synthetic experiment with a separately sampled parent-frame validation set.
- The fitted inclusion model is kept separate from hidden simulation truth during inference.
- The study quantifies sensitivity to validation sample size and to selected forms of selection-model misspecification.
- The study reports positivity/weight-tail and effective-sample-size diagnostics alongside posterior outcomes.
- The bootstrap utility quantifies approximate iid validation-resampling sensitivity conditional on the fixed catalogue, GW data, and assumed logistic family.

These are implementation/design claims. Numerical claims require completed runs and independently audited artifacts.

## Claims that are not supportable

Do not claim that:

- inverse-probability weighting is a new method;
- catalogue incompleteness or host weighting is an unexplored problem;
- the project outperforms GP, LIGHT, variance completion, CHIMERA, or gwcosmo;
- the bootstrap is a full Bayesian or exact finite-population uncertainty treatment;
- the model is validated for realistic 3D galaxy surveys or LVK observations;
- the synthetic posterior is a new H0 measurement;
- the present results solve host-galaxy population uncertainty.

## Recommended paper framing

A defensible working title is:

**Finite-Validation Uncertainty and Selection-Model Misspecification in Inverse-Probability Host Weighting: A Controlled Dark-Siren Simulation Benchmark**

A narrower alternative:

**When Does Validation-Based Catalogue Reweighting Fail? Positivity, Model Misspecification, and Cosmological Posterior Calibration in a Synthetic Dark-Siren Model**

The central question should be framed as a statistical reliability question, not as a claim to invent a new astrophysical completion method:

> Under what combinations of validation sample size, selection-model misspecification, and weak positivity does validation-estimated inverse-probability weighting produce unstable host-density estimates or miscalibrated cosmological posteriors?

The primary outcomes should be prespecified and reported jointly:
- posterior bias and RMSE;
- fixed-truth interval coverage with uncertainty intervals;
- interval width;
- inclusion-probability estimation error (simulation scoring only);
- host-density error (simulation scoring only);
- effective sample size and inverse-weight tails;
- fit failure and boundary-limited rates;
- validation-bootstrap sensitivity of posterior summaries;
- differential failure rates across paired arms.

The experiment must include correctly specified and misspecified retention laws, multiple validation sizes, and positivity regimes. A zero-misspecification control is essential. All comparisons must use matched parent populations and GW realizations.

## Necessary next scientific upgrades

1. **Do not stop at the current linear-vs-quadratic pilot.** Add a nonlinear selection surface with enough independent replicates to estimate coverage and paired effects precisely.
2. **Separate model uncertainty from validation sampling uncertainty.** The current bootstrap holds the logistic family fixed; it cannot reveal bias from choosing the wrong functional form.
3. **Add an explicitly named flexible selection-model comparator** (for example, a prespecified spline logistic model) only after a protocol and independent tests exist. Do not tune it against H0 outcomes.
4. **Use a design-aware validation uncertainty analysis.** The current ordinary bootstrap resamples with replacement and does not apply a finite-population correction. For substantial sampling fractions, implement a justified finite-population/design-based method or clearly restrict the interpretation.
5. **Benchmark positivity failures rather than concealing them.** Report low-probability mass, ESS, maximum and tail weights, and failure rates. Any weight truncation must be a separately preregistered estimator arm.
6. **Avoid claiming full astrophysical relevance.** The current simulator is one-dimensional in redshift and uses a simplified GW distance/selection likelihood. A later external benchmark against a realistic public simulator and mock survey is necessary before astrophysical claims.
7. **Freeze a confirmatory design before scaling.** Use a smoke gate, then a sufficiently powered replicate count based on precision targets for the primary paired contrast and coverage estimates.

## References and navigation

- Tagliazucchi et al. (2026), PRD 114, 083523: https://doi.org/10.1103/vmzf-1wxq
- Leyde, Baker & Enzi (2026), JCAP 01 013: https://doi.org/10.1088/1475-7516/2026/01/013
- Datrier & Hendry (2026), MNRAS 549, stag828: https://doi.org/10.1093/mnras/stag828
- Gray, Williams & Papadopoulos (2026): https://arxiv.org/abs/2609.04991
- Borghi et al. (2025): https://arxiv.org/abs/2509.18243
- Alfradique, Bom & Castro (2025), PRD 112, 063561: https://doi.org/10.1103/vd36-3mys
- VanWyngarden et al. (2025): https://arxiv.org/abs/2511.04786
- Dalang, Fiorini & Baker (2026), JCAP 01 034: https://doi.org/10.1088/1475-7516/2026/01/034
- Cheng & Gair (2026): https://arxiv.org/abs/2603.13053

## Bottom line

The project can become a rigorous, useful computational-statistics study if it demonstrates where validation-estimated weighting works, where it fails, and how users can diagnose that failure. Its contribution should be judged by the quality of the experimental design, calibration analysis, independent reproducibility, and clarity of limitations—not by a claim that the underlying weighting idea is new.
