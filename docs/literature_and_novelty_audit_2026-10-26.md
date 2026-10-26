# Targeted literature and novelty audit (updated 2026-10-27)

## Executive conclusion

**The broad research question is already well established.** Galaxy-catalogue incompleteness, host-galaxy weighting, selection functions, clustering, and their effects on dark-siren H0 inference have been studied in peer-reviewed papers and recent simulation frameworks. This repository must not claim that modelling missing hosts or testing catalogue incompleteness is itself a novel method.

The work completed here is most defensibly described as a **reproducible, independently audited synthetic robustness suite** that studies specific interactions among discrete-host misspecification, fixed assumed hidden-population fraction/redshift, prior support, and repeated-sampling calibration. Whether that narrower combination is novel has **not** been established by this targeted search. It requires a more systematic review and a direct comparison against realistic reference pipelines.

This is a targeted literature audit, not a systematic review: it is intended to prevent obvious overclaims and guide the next decision, not to establish priority exhaustively.

## Closest relevant work

| Work | What it already covers | Implication for this repository |
|---|---|---|
| Gray, Messenger & Veitch (2022), [A pixelated approach to galaxy catalogue incompleteness](https://doi.org/10.1093/mnras/stac366) | Direction-dependent catalogue completeness and its impact on dark-siren H0 inference, implemented with a realistic catalogue pipeline. | Catalogue incompleteness and the need to model its direction dependence are not new topics. Our one-dimensional discrete-host toy model is much simpler. |
| Gair et al. (2023), [The Hitchhiker's Guide to the Galaxy Catalog Approach for Dark Siren Gravitational-wave Cosmology](https://doi.org/10.3847/1538-3881/acca78) | Methodological overview of galaxy-catalogue dark sirens, potential-host marginalization, and common sources of bias from inconsistent simulations or selection modelling. | Correctly matching the data generator and likelihood is a known core requirement; internal simulation agreement alone is not evidence of astrophysical validity. |
| Finke et al. (2021), [Cosmology with LIGO/Virgo dark sirens](https://doi.org/10.1088/1475-7516/2021/08/026) | Dark-siren inference methodology, galaxy-catalogue analysis, and a detailed study of systematic effects. | Selection and population assumptions must be grounded in the established hierarchical framework before this toy model can support a cosmological claim. |
| Datrier & Hendry (2026), [Implementing a robust test of galaxy catalogue completeness for dark siren measurements of H0](https://doi.org/10.1093/mnras/stag828) | A robust completeness-limit test integrated into the gwcosmo pipeline and applied to GLADE/GLADE+ and GWTC datasets. | Catalogue completeness is already being measured and improved in real-data pipelines. Our synthetic fraction parameter is not a replacement for a survey completeness model. |
| Borghi et al. (2025), [Echoes from the Dark: Galaxy Catalog Incompleteness in Standard Siren Cosmology](https://arxiv.org/abs/2509.18243) | A self-consistent catalogue-incompleteness and host-weighting framework in CHIMERA, with joint cosmological and astrophysical population constraints over mock event populations. | Joint incompleteness and host weighting are already being studied with substantially more realistic simulation and inference. |
| VanWyngarden et al. (2025), [How Low Can You Go: Constraining the Effects of Catalog Incompleteness on Dark Siren Cosmology](https://arxiv.org/abs/2511.04786) | Mock catalogues that vary incompleteness and clustering; shows the impact depends on host-halo selection and localization volume. | The effect of incompleteness is conditional on clustering and localization, neither of which is represented realistically in our current simulator. |
| Tagliazucchi et al. (2026), [Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs](https://doi.org/10.1103/vmzf-1wxq) | Gaussian-process reconstruction of the line-of-sight galaxy field from incomplete catalogues, with a hierarchical selection-aware likelihood and mock catalogues. | More flexible latent-field reconstruction already exists; a single fixed hidden-host redshift is not a competitive replacement. |
| [Clustering effects on the dark siren determination of H0: a simulation study](https://doi.org/10.1093/mnras/stag1662) (2026) | High-resolution cosmological simulations with clustering, catalogue completeness variation, and tests of modelling choices; reports robustness to prior-range choices in its setup. | Our observed prior-support sensitivity is a property of our simplified setup and must not be generalized. The difference is a reason to compare generators and likelihood assumptions, not to claim disagreement with realistic analyses. |
| [Scalable dark siren cosmology with gwcosmo: GPU acceleration, validation, and systematics](https://doi.org/10.1103/53c5-mtk5) (2026) | Describes a validated, accelerated gwcosmo pipeline, including injection-recovery on 2,000 simulated events and comparison against the earlier implementation on 141 GWTC-4.0 events. | A mature pipeline already has substantial simulation and real-event validation. Our own toy-model reconstruction is a lower-level unit test, not a substitute for pipeline-level validation. |
| [Dark standard siren cosmology with bright galaxy subsets](https://doi.org/10.1093/mnras/stag605) (2026) | Tests the use of bright galaxy subsets with GWTC and GLADE+ data. | Host-selection/weighting alternatives have already been tested with observational data; this repository needs a distinct, explicit estimand to add value. |

| Leyde, Baker & Enzi (2026), [Cosmic cartography. Part II. Completing galaxy catalogs for gravitational-wave cosmology](https://doi.org/10.1088/1475-7516/2026/01/013) | Reconstructs the latent true galaxy field and magnitude distribution while respecting large-scale spatial correlations; validates on Millennium-simulation mocks. | A latent missing-galaxy population should not be presented as novel by itself. A single fixed or finite hidden-redshift mixture is a much simpler control model and needs to be compared with clustered field reconstruction. |
| Gray, Williams & Papadopoulos (2026), [Expanding the scope of dark siren cosmology: Inferring the population properties of gravitational-wave-hosting galaxies](https://arxiv.org/abs/2609.04991) | Updates gwcosmo methodology with a sampled redshift prior and jointly infers host-galaxy weighting/population properties alongside cosmology, reproducing a subset of GWTC-5.0 results. | Jointly learning host-population properties is an active direction; this repository's fixed assumed hidden fraction/redshift is a controlled misspecification study, not an alternative operational inference method. |
| Cheng & Gair (2026), [Unified harmonic framework for dark siren cosmology](https://doi.org/10.1103/295h-s6ss) | Clarifies the relationship between galaxy-catalogue inference and two-point cross-correlation approaches, including marginalization over the unknown galaxy field. | The field has multiple established routes for handling unknown galaxy structure. A credible external benchmark should state which estimand and latent-field assumptions it is comparing against. |

## What our results do and do not establish

### Established within this repository

- Seeded synthetic data can be regenerated and posterior results reconstructed independently for the declared likelihood.
- Paired contrasts reduce noise from using different observations for different inference assumptions.
- In the declared toy model, omission of a simulated hidden component degrades RMSE and interval coverage as the true hidden fraction grows.
- Fixed assumed hidden fraction/redshift values can differ from the oracle arm; adding a spurious hidden component can be harmful when the true hidden fraction is zero.
- Prior support changes posterior error, interval width, and coverage in this finite-host model. The magnitude and direction depend on the truth, hidden fraction, and inference arm.

### Not established

- That these exact interaction effects are novel relative to all prior literature.
- That the toy model is an adequate approximation to a GW likelihood or a real galaxy survey.
- That the oracle arm is usable for real data; it is given the true hidden fraction and redshift by construction.
- That any prior support is scientifically preferred.
- That the method improves H0 inference on public gravitational-wave/catalogue data.
- That the project is publication-ready.

## Recommended contribution framing

Until a more complete literature review and benchmark exist, describe the contribution as:

> A reproducible simulation and audit study of how prior support and hidden-host fraction/redshift misspecification interact in a simplified discrete-host dark-siren likelihood, using common random numbers, paired bootstrap uncertainty, negative controls, and independent numerical reconstruction.

This is deliberately a **candidate framing**, not a novelty claim. The primary technical contribution should be assessed against existing completeness-aware hierarchical methods, not against a baseline that omits effects already known to matter.

## Decision gates before any novelty or publication claim

1. **Systematic literature matrix:** search ADS/arXiv and journal indexes for catalogue completeness, latent galaxy-field reconstruction, host-weight misspecification, prior sensitivity, selection normalization, and simulation-based calibration. Record exact estimands, data generators, likelihoods, public code, and limitations.
2. **Baseline reproduction:** reproduce at least one published, open-source completeness-aware method on its documented mock data before adding our misspecification factor.
3. **Model compatibility:** make event detection selection, catalogue inclusion, latent host membership, and H0 prior support explicit in one joint generative model; derive the likelihood and check limiting cases.
4. **Realistic validation:** compare with realistic clustered mock catalogues or an established public pipeline. The current fixed-redshift mixture is not sufficient.
5. **Independent review:** ask a researcher familiar with dark-siren hierarchical inference to challenge the novelty framing, conditioning assumptions, and interpretation before presenting the work as original research.

The next defensible milestone is therefore not another isolated toy-model result. It is a reproducible comparison against a published completeness-aware baseline, with our current toy model retained as a unit-test/control environment.
