# Novelty scope and publication gate

**Date:** 29 October 2026  
**Status:** conservative positioning memo; not a systematic literature review  
**Project:** synthetic dark-siren H0 robustness study

## Executive decision

**No broad novelty gap is established.** Catalogue incompleteness, hidden-host populations, host weighting, selection normalization, and population misspecification are established topics in dark-siren cosmology. The project must not claim to be the first study of these topics or imply superiority to established pipelines.

The defensible candidate contribution is narrower: **a reproducible synthetic calibration and identifiability study that maps when specific simplified host-density, hidden-host, and selection-model assumptions miscalibrate H0 posteriors, and measures whether prespecified diagnostics detect those failures.** This contribution is not yet established by a title or protocol alone. It requires completed, auditable results that demonstrate useful diagnostic power beyond simply restating expected model misspecification.

## Conservative contribution statement

> This project develops a reproducible synthetic Bayesian framework to examine when simplified dark-siren H0 inference—using KDE-based host-density estimation, hidden-host fractions, inverse-inclusion weighting, and selection-normalized likelihoods—remains calibrated under controlled misspecification. It evaluates dataset-level coverage, bias, interval width, RMSE, and diagnostic behaviour, while explicitly separating oracle information from realistically available information.

Scope limits:
- The generative model is controlled and simplified; it is not a realistic end-to-end LVK detector and galaxy-survey model.
- Results do not constitute a measurement of the real Hubble constant.
- Results do not validate CHIMERA v2.2, CHIMERA generally, or gwcosmo.
- Known selection probabilities and oracle host-density information are diagnostic controls, not realistic estimators available in every application.
- Inverse-inclusion weighting is not a universal correction: it depends on correct inclusion probabilities, support/positivity, appropriate weights, and the estimand/model.
- Independent artifact reconstruction is evidence of reproducibility and software quality, not automatically a novel scientific discovery.

## Title recommendation

**When Do Simplified Dark-Siren H0 Inference Methods Remain Calibrated? A Reproducible Synthetic Study of Hidden-Host and Selection-Function Misspecification**

Alternative:

**A Reproducible Synthetic Calibration Study of Galaxy-Catalogue Incompleteness and Host-Population Misspecification in Dark-Siren Cosmology**

Avoid titles that imply a real H0 measurement, a validated observational pipeline, or a first-of-its-kind method.

## Claim-to-evidence gates

| Proposed claim | Current position | Evidence required before strengthening |
|---|---|---|
| Novel study of catalogue incompleteness/host weighting | Unsupported as a broad novelty claim | Precise, sourced distinction from existing methods and a substantive new result |
| KDE bandwidth sensitivity is novel | Not established | Show how bandwidth interacts with calibration/identifiability and that diagnostics detect a consequential failure regime |
| Inverse-inclusion weighting corrects incompleteness | Too broad | Restrict claim to the declared inclusion mechanism and verify positivity, weights, support, and correct model assumptions |
| Selection normalization affects inference | Established principle | Demonstrate a controlled, reproducible effect size and a negative control that catches the omission |
| Synthetic calibration generalizes to LVK observations | Unsupported | Realistic injection/selection, localization, redshift-error, clustering, catalogue and population-model validation |
| CHIMERA v2.2 is validated | NO-GO | Exact pinned compatible inputs, checksums, environment, end-to-end execution and independent output reconstruction |
| Tests passing establishes scientific validity | False | Combine exact-head CI with numerical, statistical, model, provenance, and interpretation audits |
| Independent artifact reconstruction is a major discovery | Not by itself | Demonstrate substantive diagnostic or inferential value, beyond an engineering contribution |

## Experimental design and reporting requirements

Before confirmatory runs, freeze:
1. Generative model and exact data-selection mechanism, including whether selection acts on latent or noisy measured distance.
2. The pre-selection meaning of the hidden-host fraction and the detected-sample mixture implied by selection.
3. The exact conditional-on-detected-event likelihood versus any rate-based point-process likelihood.
4. True and assumed selection laws, baseline cutoff, and cutoff offsets in physical units.
5. KDE bandwidth-selection rule, support, boundary treatment, and oracle/non-oracle information split.
6. Primary estimands, strata, replicate seeds, failed-fit handling, multiplicity/contrast policy, and Monte Carlo precision.
7. Raw-artifact schema, manifest, independent auditor, and acceptance criteria.

Report signed bias, RMSE, nominal 68% and 95% coverage, covered counts, Wilson intervals, interval widths, failed-fit rates, boundary posterior mass, and stratified results for endpoint, near-boundary, and interior truths. Pooled SBC or coverage must not replace conditional diagnostics. A nominal coverage band is a prespecified decision rule, not proof of calibration.

The proposed 5 × 3 × 3 × 3 design has 135 cells; 500 replicates per cell means 67,500 datasets. Pilot runtime and Monte Carlo precision first. Do not change cells after inspecting favourable results without clearly labelling the change exploratory and preregistering a new confirmatory stage.

## Publication readiness

**Publication readiness: NO-GO at this checkpoint.** A scoped manuscript may be developed, but it must not be described as submitted, accepted, peer reviewed, or publication-ready without the corresponding facts and completed evidence gates.

A software-journal route such as JOSS should be considered only if the repository becomes a substantial, maintainable research-software package with clear research use and evidence of impact. A collection of notebooks or a one-off simulation report may not satisfy that bar. Venue selection should follow the demonstrated contribution, not prestige.

## External validation gates

### CHIMERA v2.2

**NO-GO.** Existing provenance audit records that the pinned notebook expects mock-specific files absent from the tagged source tree and that exact compatible inputs have not been verified. Do not silently substitute a v2.0 bundle and call it a v2.2 reproduction. Re-entry requires version-matched assets, provenance, checksums, licensing/access review, pinned environment, successful execution, and independent reconstruction.

### Real observational analysis

**NO-GO.** A defensible observational analysis requires a compatible GW event/posterior set, selection injections, galaxy catalogue and completeness model, consistent versioning, event-selection rule, population assumptions, redshift-error treatment, convergence checks, and independently auditable outputs. No real-data H0 claim is permitted until these are implemented and verified.

## Literature handling note

The repository already contains a focused literature/positioning audit dated 29 October 2026. This memo is a synthesis of that audit and the external review supplied for this project; it is **not** a fresh systematic search and does not independently verify every citation in the supplied bibliography. Before submission, verify all bibliographic metadata, publication status, version, DOI/arXiv record, and the exact claim supported by each source. Clearly distinguish peer-reviewed articles from preprints and official software documentation.

## Release checklist

- [ ] Exact-head tests pass on supported Python versions.
- [ ] Likelihood derivation and implementation agree, including conditioning and normalization.
- [ ] Confirmatory protocol frozen before confirmatory outputs are inspected.
- [ ] Raw data, configuration, seeds, software versions, and hashes are archived.
- [ ] Independent reconstruction reproduces posterior summaries and metrics.
- [ ] Boundary, near-boundary, and interior performance are reported separately.
- [ ] Negative controls detect known deliberately incorrect models.
- [ ] All null/adverse outcomes and failed fits are reported.
- [ ] Novelty claims are narrowed to a source-backed, demonstrable contribution.
- [ ] CHIMERA and observational claims remain explicitly NO-GO until their own gates pass.
