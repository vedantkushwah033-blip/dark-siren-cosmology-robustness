# Claim-by-claim literature overlap matrix — 2026-10-28

**Status:** targeted bibliographic triage based on publisher/preprint records and abstracts, not a systematic review. The matrix is for research planning and novelty-risk control; it is not evidence that the cited methods have been reproduced or fully compared. Full papers, appendices, code, and exact estimands must be checked before manuscript submission.

## Why this matrix exists

The repository's possible contribution must be framed narrowly. Catalogue incompleteness, host weighting, selection effects, spatially structured catalogue completion, and misspecification in dark-siren inference are active existing research areas. Merely simulating one of these effects is not a novelty claim.

## Literature overlap

| Work | Problem / approach visible in the source record | Overlap with this repository | Difference that may matter | What must be checked in full text |
|---|---|---|---|---|
| Gray, Messenger & Veitch (2022), [MNRAS 512, 1127](https://academic.oup.com/mnras/article/512/1/1127/6526882) | Pixelated, direction-dependent galaxy-catalogue completeness in dark-siren inference. | Catalogue incompleteness affects the host/redshift contribution and H0 posterior. | Repository is a simplified synthetic host-mixture / robustness study, not a sky-pixelated catalogue pipeline. Scope difference alone does not establish novelty. | Completeness parameterization, calibration tests, misspecification scenarios, and released code. |
| Gair et al. (2023), [AJ 166, 22](https://eprints.gla.ac.uk/334400/) | Methodological guide to galaxy-catalogue dark sirens; emphasizes consistent mock generation, Bayesian implementation, and selection effects. | Strong relevance to inference/generator consistency and audit design. | Repository may offer a more focused controlled diagnostic, but only if it isolates a distinct estimand and tests beyond known consistency safeguards. | Exact failure modes, recommended validation tests, and overlap with the repository's mixture/selection controls. |
| Datrier & Hendry (2025 preprint; journal record 2026), [arXiv:2502.14164](https://arxiv.org/abs/2502.14164) | Robust test for the limiting apparent magnitude / completeness of galaxy catalogues, implemented in `gwcosmo`. | Directly overlaps catalogue-completeness diagnostics. | The repository does not currently implement that survey-level test or validate it in `gwcosmo`. | The published version's assumptions, statistical guarantees, data requirements, and failure regimes. |
| Leyde, Baker & Enzi (2026), [JCAP 01 (2026) 013](https://doi.org/10.1088/1475-7516/2026/01/013) | Reconstructs the true galaxy field while respecting spatial correlations and infers the galaxy magnitude distribution, tested on simulated catalogues. | Overlaps the scientific motivation that incomplete catalogues alter host/redshift priors. | The repository's toy model does not reconstruct a spatially correlated galaxy field. | Generator details, validation metrics, and relation to other catalogue-completion methods. |
| Tagliazucchi et al. (accepted 2026), [Physical Review D record](https://journals.aps.org/prd/accepted/10.1103/vmzf-1wxq) | Gaussian-process line-of-sight redshift reconstruction with a hierarchical likelihood accounting for survey selection; tested on mock GW/catalogue data. | Overlaps GP-related exploration and incomplete-catalogue robustness. | A repository GP pilot is not equivalent to a selection-aware hierarchical reconstruction. | Full model, priors, calibration, simulation configurations, and software implementation. |
| Borghi et al. (2025 preprint), [arXiv:2509.18243](https://arxiv.org/abs/2509.18243) | CHIMERA-based study of incomplete spectroscopic catalogues, stellar-mass host weighting, and an O5-like 100-event simulated population. | Closely related to the project's planned external benchmark and host-weighting sensitivity. | The current repository's synthetic toy model is not a reproduction of the paper or CHIMERA v2.2. | Exact version, mock products, host-weighting definitions, selection treatment, and available replication artifacts. |
| VanWyngarden et al. (2025), [arXiv:2511.04786](https://arxiv.org/abs/2511.04786) | Studies when bright-galaxy subsets can produce informative H0 posteriors under specific clustering and localization assumptions. | Shows catalogue incompleteness effects depend on data-generating regime. | The repository's conclusions must be similarly conditional on its declared generator; no universal direction of bias should be claimed. | Assumptions about clustering, localization, host selection, and what constitutes unbiased/informative inference. |
| Gray, Williams & Papadopoulos (2026 preprint), [arXiv:2609.04991](https://arxiv.org/abs/2609.04991) | Jointly infers host-galaxy weighting, GW population properties, and cosmology using a sampled redshift prior; abstract reports an application to GWTC-5.0. | Direct overlap with host-population assumptions and cosmological inference. | The repository does not jointly infer a realistic host population and cosmological model. | Exact hierarchical structure, identifiability, prior sensitivity, selection normalization, and reproducibility materials. |

## Claim disposition

### Claims that should be rejected in the current manuscript

- “We solve galaxy-catalogue incompleteness.”
- “We introduce a new general host-weighting method.”
- “We improve real-data H0 inference.”
- “Our method is validated against CHIMERA v2.2.”
- “Our diagnostic is the first” or “state of the art.”

The current evidence does not support these claims.

### Narrow candidate worth testing

A possible question is:

> In a declared synthetic low-redshift host-mixture model, when do prespecified posterior-calibration and predictive diagnostics distinguish host-population misspecification from distance-noise misspecification, and when are the sources observationally indistinguishable at finite dataset sizes?

This remains a **candidate question, not an established novel contribution**. Before running a large experiment, the project must:
1. Resolve the event-level versus dataset-level likelihood design gate in `docs/diagnostic_identifiability_design_and_go_no_go_2026-10-28.md`.
2. Prespecify the primary diagnostic, operating threshold, false-positive rate, sensitivity, and uncertainty intervals.
3. Include matched controls, host-only, noise-only, joint mismatch, and negative controls.
4. Replicate across multiple generator configurations and random seeds.
5. Compare the resulting estimand against full-text work on Bayesian calibration, simulation-based calibration, posterior predictive checks, and misspecified Bayesian models.
6. Report non-identifiability or low power honestly if that is what the experiment shows.

## External-validation gate

The planned CHIMERA comparison remains blocked until the exact reference version and required input products are available and verified. Do not substitute another version or a simplified synthetic calculation and label it exact validation.

## Sources inspected and limits

The new 2026 records for Tagliazucchi et al. and Gray et al. were found through publisher/preprint records; the other records above were also checked at publisher or preprint sources. This triage uses abstracts and metadata, not a full-text systematic search. It can identify obvious overlap and prevent overclaiming, but it cannot prove novelty. The next literature stage must read the complete relevant papers and their supplementary material and record explicit comparison by estimand, model, diagnostic, and validation standard.
