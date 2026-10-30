# Related-work and novelty-boundary audit

**Audit date:** 30 October 2026  
**Scope:** Focused candidate on finite-parent validation-sampling uncertainty in inverse-probability host weighting. This is a targeted positioning audit, not a systematic review and not proof that no related work exists.

## 1. Established neighbouring work

1. **Galaxy-catalogue incompleteness and host weighting are already active research topics.** Borghi et al., *Echoes from the Dark: Galaxy Incompleteness in Standard Siren Cosmology* (2025), develop a self-consistent incompleteness and host-weighting framework in CHIMERA and study mock catalogues and multiple completeness levels. https://arxiv.org/abs/2509.18243
2. **Nonparametric reconstruction of incomplete line-of-sight galaxy distributions is not new.** A recent Physical Review D paper develops Gaussian-process reconstruction of the galaxy redshift prior while accounting for catalogue selection. https://journals.aps.org/prd/abstract/10.1103/vmzf-1wxq
3. **Incorrect host weighting can bias dark-siren inference.** Hanselman et al. (2024) study galaxy-weighting systematics with mock catalogues and discuss the dependence on localization, distance uncertainty, and galaxy correlations. https://arxiv.org/abs/2405.14818
4. **Selection and weighting systematics have direct H0 consequences.** A 2025 Physical Review D study analyses systematic bias from catalogue incompleteness and incorrect host weighting. https://journals.aps.org/prd/abstract/10.1103/vd36-3mys
5. **Host-galaxy population inference continues to evolve.** Gray, Williams & Papadopoulos (2026) describe joint inference of cosmology and host-galaxy population properties in a dark-siren framework. https://arxiv.org/abs/2609.04991

## 2. What this means for the present project

The paper must not claim novelty for:
- catalogue incompleteness itself;
- inverse-probability or stellar-mass host weighting in general;
- selection-aware dark-siren likelihoods;
- Gaussian-process reconstruction of missing-galaxy structure;
- generic host-population misspecification or generic robustness testing.

Those topics have substantial prior literature and must be treated as context, not as the contribution.

## 3. Narrow candidate contribution

The defensible candidate is a **design-based uncertainty-calibration benchmark** for a particular workflow:

- a finite synthetic parent frame of paired galaxy redshift and catalogue-retention records;
- external validation samples drawn exactly without replacement (SRSWOR);
- a nonlinear estimator that fits inclusion probabilities, forms inverse-probability weights, estimates a reflected weighted KDE, and propagates it into a simplified selection-aware H0 posterior;
- direct comparison of ordinary paired-record bootstrap uncertainty with the empirical repeated-SRSWOR distribution, including a finite-population-correction candidate;
- complete replayable sampling indices, seeds, failure accounting, and independent reconstruction of numerical outputs.

The target is not to invent a new dark-siren likelihood. It is to test whether a commonly used resampling uncertainty estimate is calibrated for this finite-parent validation design, and to map where nonlinear fitting, positivity, and extreme weights affect that calibration.

## 4. Claims that remain unproven

A targeted search has not established that this exact finite-parent bootstrap-calibration question has never been studied in adjacent survey-sampling, missing-data, semiparametric, or astronomical contexts. Before submission, search those fields systematically and ask a statistician to review the estimand and variance design.

The ordinary bootstrap may perform adequately in some regimes; the FPC may not repair nonlinear-estimator behaviour. Both are empirical hypotheses. No result should be written into the abstract until the confirmatory experiment and independent artifact audit are complete.

## 5. Mandatory decision gates

1. Green exact-head CI on every supported Python version.
2. Independently audited finite-parent smoke artifact.
3. Deterministic tests for paired resampling, SRSWOR, FPC edge cases, and failure accounting.
4. Pilot to estimate runtime, optimizer failures, variance components, and bootstrap Monte Carlo error.
5. Freeze confirmatory design and primary endpoints before confirmatory outcomes are inspected.
6. Independent reconstruction of all successful outer and bootstrap fits and all reported aggregate statistics.
7. External statistical and domain review before describing the work as publication-ready.

## 6. Appropriate positioning

Use language such as **“a controlled synthetic benchmark of finite-validation uncertainty calibration in a simplified dark-siren host-weighting pipeline.”**

Do not describe it as a new cosmological measurement, an LVK-data result, a solution to the Hubble tension, a new general dark-siren method, or a proven first-of-its-kind contribution.


## 7. Additional survey-sampling bootstrap literature identified on 30 October 2026

The resampling question has a substantial adjacent literature outside astronomy. This materially narrows the novelty claim and changes the experimental comparator set.

1. **Bootstrap inference for finite populations under complex designs.** *Bootstrap Inference for the Finite Population Mean under Complex Sampling Designs*, Journal of the Royal Statistical Society: Series B (2022), 84(4), 1150 onward. The paper explicitly includes simple random sampling without replacement and describes a method that generates bootstrap finite populations, reapplies the sampling design, and uses studentization. DOI and article record: https://doi.org/10.1111/rssb.125ಿಂಗ್
2. **Direct bootstrap methods for complex finite-population designs.** The JASA article *A Direct Bootstrap Method for Complex Sampling Designs From a Finite Population* (published online 2012) develops resampling methods that combine with- and without-replacement selections and discusses classical sampling designs including SRSWOR. DOI: https://doi.org/10.1198/jasa.2011.tm09767
3. **Bootstrap variance estimation for model parameters from survey data.** *Bootstrap variance estimation with survey data when estimating model parameters* (Computational Statistics & Data Analysis, 2012) highlights that both finite-population model variability and sampling-design variability may matter; ignoring one component can materially affect variance estimates. DOI: https://doi.org/10.1016/j.csda.2012.03.011
4. **Two-phase and without-replacement design adjustments.** Beaumont & Émond (2022), *A Bootstrap Variance Estimation Method for Multistage Sampling and Two-Phase Sampling When Poisson Sampling Is Used at the Second Phase*, Stats 5(2), 339–357, discusses design-specific bootstrap weights and the role of finite-population corrections. https://doi.org/10.3390/stats5020019

### Consequence for the proposed study

A scalar FPC multiplier is only a simple baseline, not a substitute for a design-consistent survey bootstrap. The revised protocol therefore requires three resampling arms—ordinary paired bootstrap, ordinary bootstrap with scalar FPC, and a literature-grounded design-aware finite-population bootstrap—compared with the empirical repeated-SRSWOR reference. The exact design-aware method must be selected and implemented from its published specification, with its assumptions checked before use. It must not be claimed to transfer automatically to the nonlinear logistic/KDE/posterior pipeline.

This makes the paper's potential contribution more modest but more credible: an independently auditable empirical comparison of established resampling approaches when propagated through one explicitly defined nonlinear synthetic dark-siren host-weighting workflow. It does not establish novelty for finite-population bootstrap theory or for survey-sampling resampling methods generally.
