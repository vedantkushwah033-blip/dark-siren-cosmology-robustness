# ASA manuscript and poster blueprint

**Working status:** manuscript architecture only; not a submission-ready paper.  
**Evidence freeze:** not yet declared.  
**Scope:** synthetic low-redshift dark-siren statistical simulation; discrete-host and continuous-hidden-population model families must be reported separately.

## Working title

**Finite-Sample Reliability of Dark-Siren Hubble-Constant Inference under Hidden-Host Omission and Mixture Misspecification: A Controlled Simulation Study**

This title is intentionally narrower than a claim to develop a new galaxy-catalogue completion method.

## Research question

In controlled synthetic dark-siren models, how do omission and fixed misspecification of a latent host component affect H0 posterior point error and fixed-truth credible-interval coverage, and how do these effects depend on prior support and finite-sample estimation of hidden-mixture parameters?

## Candidate contributions (claims remain conditional on completed audits)

1. A paired simulation design separates hidden-component omission from fixed hidden-component misspecification, using common simulated observations across inference arms.
2. It quantifies both point-estimation error and fixed-truth posterior interval coverage, rather than treating posterior width as a sufficient measure of reliability.
3. It distinguishes prior-predictive simulation-based calibration under an assumed model from repeated-sampling coverage when the data-generating process and fitted model differ.
4. It provides a transparent continuous-mixture extension to examine finite-sample estimation of hidden-component parameters, explicitly as a simplified statistical diagnostic rather than a realistic galaxy-field reconstruction method.

Do not describe these as first-ever contributions until a complete paper-level literature review confirms the gap.

## Abstract (do not finalize until evidence freeze)

Dark-siren inference combines gravitational-wave distance information with redshift information from potential host galaxies. Incomplete catalogues and uncertain host-population assumptions can affect the resulting posterior for the Hubble constant. We study these effects in controlled synthetic models, separating a discrete positive-truncated Gaussian-distance experiment from a continuous hidden-redshift log-distance extension. The discrete experiment compares omission, oracle knowledge, and fixed misspecification of a hidden host component, and examines prior-support sensitivity. The continuous extension evaluates finite-sample recovery of hidden-mixture parameters and compares oracle, fixed-wrong, plug-in, and integrated inference. Performance is assessed using posterior point error, RMSE, interval width, fixed-truth coverage, and numerical diagnostics, with Monte Carlo uncertainty reported explicitly. The analysis is restricted to simplified synthetic models and does not use public gravitational-wave events or realistic detector-selection and galaxy-field models. Final numerical results and conclusions will be inserted only after independent reproduction and validation.

## Manuscript structure

### 1. Introduction
- Explain dark-siren redshift ambiguity and why galaxy catalogues enter.
- Summarize the recent catalogue-completion and host-population literature (2024–2026).
- State that this work is a statistical reliability study, not a new astrophysical completion framework.
- End with the research question and scoped contributions.

### 2. Related work
Organize by:
1. Catalogue incompleteness and host weighting.
2. Spatially structured/large-scale-structure completion.
3. Hierarchical inference of host-population properties.
4. Bayesian computation validation and posterior behaviour under misspecification.

Compare the estimand, data model, missing-population model, selection treatment, inference target, simulation realism, and evaluation metrics in each paper.

### 3. Discrete hidden-host model
- Listed redshifts and weights.
- Point-mass hidden redshift and hidden fraction.
- Positive-truncated Gaussian distance generator.
- Conditional component likelihood and unchanged mixture weights under the stated redraw-same-host generator.
- H0 posterior and prior support.
- Explicit distinction from GW detection selection and catalogue-inclusion selection.

### 4. Continuous hidden-population extension
- (u=\log d) observation model.
- Truncated-normal hidden-redshift distribution on fixed support.
- Gauss-Legendre marginalization.
- Parameter vector ((H_0,f,\mu_h,s_h)).
- Which parameters are fixed versus estimated in each experiment.
- Numerical integration and prior-grid contracts.
- Local information diagnostic and its limitations.

### 5. Simulation design and inference arms
- Truths, sample sizes, noise levels, random seeds, replication counts.
- Oracle, omitted, fixed-wrong, plug-in, fully integrated methods.
- Common random numbers / paired observations.
- Primary and exploratory cells declared before confirmatory runs.

### 6. Statistical evaluation
- Posterior median bias and RMSE.
- 68% and 95% fixed-truth coverage.
- Mean interval widths.
- Monte Carlo standard errors and Wilson intervals for coverage.
- Paired contrasts with bootstrap intervals where prespecified.
- SBC as a separate prior-predictive computation diagnostic.
- Grid and quadrature convergence.

### 7. Results
Use separate subsections for:
- Existing discrete hidden-host experiments.
- Continuous hidden-population recovery.
- Local identifiability diagnostics.
- Prior-support sensitivity.
- Numerical convergence.
No pooling of the older calibration experiment with the newer hidden-host experiments.

### 8. Discussion
- Interpret results only within the synthetic assumptions.
- Explain why oracle performance is a benchmark, not a deployable method.
- Discuss when fixed-wrong assumptions can be worse than omission.
- Explain posterior coverage versus interval width.
- Compare carefully with the recent literature.

### 9. Limitations
State clearly:
- no real LVK/GWTC event data;
- no realistic detector selection function;
- no sky-position or galaxy clustering model;
- no realistic luminosity/stellar-mass host-rate model;
- no peculiar velocities or redshift measurement errors;
- simplified low-redshift distance relation;
- known event-noise model;
- simplified hidden-host distribution;
- finite grid and Monte Carlo uncertainty;
- no claim of external astrophysical validity.

### 10. Conclusion
Summarize only validated results. Avoid claims of a universal correction, a real H0 measurement, or superiority over current dark-siren pipelines.

## Planned tables

1. Model assumptions and notation.
2. Simulation factorial and inference arms.
3. Discrete hidden-host omission results.
4. Fixed hidden-assumption misspecification contrasts.
5. Prior-support sensitivity.
6. Continuous hidden-parameter recovery and Monte Carlo uncertainty.
7. Local information/correlation diagnostics.
8. Numerical convergence audit.

## Planned figures

1. Model schematic: listed versus hidden host component.
2. RMSE and coverage against true hidden fraction.
3. Paired effect of wrong hidden fraction/redshift assumptions.
4. Prior-support coverage-width trade-off.
5. Continuous-model parameter-recovery intervals by sample size.
6. Local standardized information eigenvalues/correlation matrix.
7. Grid and quadrature convergence summary.

## Poster storyboard (single landscape poster)

- **Top band:** title, authors, affiliation only as actually authorized.
- **Left column:** motivation, research question, related-work gap.
- **Middle column:** two clearly separated model families and simulation design.
- **Right column:** three validated results, uncertainty bars, limitations, reproducibility QR/link.
- **Bottom strip:** synthetic-only scope and acknowledgements.

The poster must not use institutional or conference logos without permission. It must not imply endorsement or present unvalidated analogue results as project findings.

## Evidence needed before final manuscript

- CI green on the exact evidence-freeze commit.
- Reproduction of the continuous-model unit tests and parameter-recovery smoke.
- Confirmatory parameter-recovery runs with prespecified sample sizes and replication counts.
- Full-grid and quadrature convergence on primary cells.
- Independent reconstruction of raw simulation metrics.
- Literature comparison based on full papers, not abstracts alone.
- Written domain review of the simplified generative assumptions.
- Final figure/table scripts and source hashes.
