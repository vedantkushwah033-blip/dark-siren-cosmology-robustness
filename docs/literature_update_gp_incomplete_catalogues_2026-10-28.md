# Literature update: GP reconstruction of incomplete galaxy catalogues (2026-10-28)

## Source

Matteo Tagliazucchi, Jonathan Gair, Riccardo Barbieri, and Michele Moresco, “Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs,” *Physical Review D* **114**, 083523 (published 13 October 2026), DOI: [10.1103/vmzf-1wxq](https://doi.org/10.1103/vmzf-1wxq). Official publisher page: https://journals.aps.org/prd/abstract/10.1103/vmzf-1wxq.

This is directly relevant literature for catalogue-incompleteness modelling. The paper proposes reconstructing the line-of-sight galaxy redshift distribution with a Gaussian process (GP), fitting an explicitly selection-aware galaxy likelihood. In mock catalogues, the authors report average H0 interval-precision improvements of 23% for a 24%-complete catalogue and 37% for an 8%-complete catalogue relative to homogeneous completion. These are results for their particular simulated configurations, not a universal guarantee.

## Why it matters to this repository

Our existing robustness experiments study controlled synthetic finite-host mixtures and how omitted hosts, host weights, distance uncertainty, and prior support affect H0 inference. The current dataset-level prototype is intentionally much narrower: independent events share H0, but it does not implement a survey-selection likelihood, event-rate inference, or joint inference of the galaxy population and H0.

The paper reinforces an important distinction: adding a scalar missing-host fraction is not equivalent to reconstructing a structured redshift distribution, and neither is automatically equivalent to a selection-aware hierarchical model. Catalogue incompleteness and GW detection selection are separate mechanisms and should be represented explicitly.

## Research implications and proposed sequence

1. **Keep the current synthetic results scoped.** Describe them as conditional sensitivity analyses under the stated finite-host assumptions; do not claim they validate a GP completion method or reproduce the paper.
2. **Add a prespecified GP-versus-homogeneous benchmark only after the generative model is explicit.** Generate a known latent galaxy redshift intensity with overdensities/underdensities, apply a magnitude/redshift-dependent catalogue detection function, and compare (a) complete-catalog oracle, (b) homogeneous-completion baseline, and (c) GP reconstruction. Evaluate reconstruction error and downstream H0 bias, RMSE, interval coverage, and interval width.
3. **Separate information and selection.** State which likelihood conditions on detected GW events and include the corresponding detection-efficiency normalization when modelling a detected-event population. Use the same population assumptions consistently in event terms and selection normalization.
4. **Avoid oracle leakage.** A method supplied the true missing-host fraction, latent redshift intensity, or GP hyperparameters is a diagnostic oracle, not a deployable estimator. Label such arms separately.
5. **Validate uncertainty, not just point estimates.** Use independent simulation replicates, fixed seeds, calibration checks, grid-convergence checks, and an independent reconstruction of reported summaries. Predefine the primary estimand and contrasts before running the full experiment.
6. **Keep external validation gated.** The separate CHIMERA v2.2 reproduction remains NO-GO until the exact notebook inputs and derived catalogue products have verifiable provenance and checksums. This paper does not supply those CHIMERA-specific files and is not a substitute for them.

## Limitations of this literature note

This is a source-informed research-planning note, not an implementation of the GP method, a replication of the paper, or an independent review of its numerical results. The paper's reported precision gains are summarised from the publisher's abstract and should be checked against the full methods and tables before being used in a quantitative comparison. No real-data H0 result is produced by this repository note.

## Decision

Treat the GP reconstruction paper as a priority reference for designing a more realistic, selection-aware catalogue-incompleteness experiment. Do not fold its published percentages into our project's results or present our current model as a reproduction.
