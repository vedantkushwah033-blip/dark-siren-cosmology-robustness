# Project use, intended audience, and limits

## One-sentence description

This project is a reproducible synthetic study of how prior boundaries, host-galaxy weighting, catalogue incompleteness assumptions, and simplified detection-selection corrections can affect Bayesian inference of the Hubble constant, H0, in a dark-siren-like setting.

## Who may find it useful?

### Students and early-career researchers
The repository can serve as a worked example of a simulation-based methods study: define a data-generating process, compare inference assumptions on matched synthetic observations, assess numerical grid error, test interval coverage, and preserve the full workflow in version control.

### Researchers in Bayesian inference and statistical methodology
The project can be used as a small testbed for studying prior sensitivity, selection normalization, repeated-sampling coverage, paired comparisons, and reproducibility. Any reuse should begin by reading the model assumptions and running the tests.

### Gravitational-wave and cosmology researchers
The code may be useful as a pedagogical prototype or a starting point for designing controlled robustness checks. It is not a replacement for a realistic gravitational-wave likelihood, detector selection function, galaxy-catalogue model, or collaboration-grade analysis.

### Scientific software reviewers
The repository provides testable numerical routines, simulation scripts, structured outputs, CI workflows, and artifact-integrity checks. Reviewers can inspect how assumptions map to code and whether outputs can be reproduced.

## Where could the work be used?

- A student research portfolio, with an honest description of the project as synthetic computational research.
- A methods-focused research note or preprint, if the completed results support a clear and nontrivial conclusion and the limitations are explicit.
- A computational astrophysics or Bayesian statistics presentation.
- Teaching material for prior sensitivity, Monte Carlo uncertainty, numerical integration, credible intervals, and selection effects.
- As a preliminary design study before building a more realistic analysis with observational data.

These are possible uses, not claims of acceptance, publication suitability, or endorsement.

## What the current project can legitimately investigate

The current toy setup can compare methods under its own stated assumptions, including:
- How changing prior support changes posterior summaries.
- How a simplified assumed detection threshold changes a selection-corrected posterior.
- How repeated-sampling bias, RMSE, interval coverage, and interval width vary across the designed cells.
- How grid resolution changes posterior summaries.
- Whether paired comparisons are computed and exported consistently.

Conclusions should be conditional on this specified data-generating model. A method that performs well under the toy model is not automatically reliable for real data.

## What it cannot yet establish

The repository does not, by itself, establish:
- A measured value of H0 from gravitational-wave observations.
- The accuracy of an LVK or other collaboration's published inference.
- Realistic performance under detector SNR thresholds, duty cycles, waveform-dependent detection, or the full selection function.
- Robustness to real galaxy-catalogue incompleteness, redshift errors, host-galaxy evolution, or correlated population uncertainties.
- Publication readiness, peer-reviewed validity, or a guarantee that an academic programme or researcher will adopt the work.

## What would make it more useful to practitioners?

1. Complete the planned factorial experiment and inspect the actual artifact files.
2. Report Monte Carlo uncertainty for coverage and paired contrasts, not only point estimates.
3. Add sensitivity checks for simulation assumptions and demonstrate that conclusions are not driven by one chosen configuration.
4. Compare the numerical implementation against an independent reference calculation and test posterior quantiles at finer resolution.
5. If moving toward real-data relevance, implement a documented, validated event likelihood and selection model using appropriate public data and carefully controlled catalogue assumptions.
6. Have an experienced researcher review the scientific framing and the gap between the toy model and real analyses.

## Recommended public description

“An open, reproducible simulation study examining prior-boundary and simplified selection-threshold sensitivity in a Bayesian dark-siren-like H0 inference model. The project focuses on numerical convergence, repeated-sampling calibration, paired robustness comparisons, and transparent limitations. It is a synthetic methods study, not an observational H0 measurement.”

