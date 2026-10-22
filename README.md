# Dark-Siren Cosmology Robustness

A reproducible Bayesian study of galaxy-catalogue incompleteness, host-population misspecification, and statistical calibration in dark-siren cosmology.

**Status: early research scaffold. No scientific results have been established by this repository yet.**

## Research question

Under which controlled conditions do galaxy-catalogue incompleteness and host-population misspecification affect the bias, uncertainty, and calibration of Bayesian inference for the Hubble constant, and when can explicitly modelling missing hosts improve inference?

## Why this project exists

Gravitational-wave sources without uniquely identified electromagnetic counterparts are often called *dark sirens*. Statistical information about possible host galaxies can help constrain cosmological parameters. This project will use controlled simulations to examine how incomplete galaxy catalogues and incorrect host-population assumptions affect inference.

The first stage is validation: the existing research implementation must be audited before results or scientific claims are carried into this extension.

## Scientific guardrails

- Synthetic results will not be described as measurements from real gravitational-wave observations.
- Galaxy-catalogue inclusion and gravitational-wave detection selection are distinct effects and must be modelled separately.
- Posterior coverage, simulation-based calibration, bias, and precision are different quantities.
- Working hypotheses are not findings; both negative and inconclusive results will be reported.
- No claim of novelty or publication readiness will be made without appropriate literature review and independent scrutiny.

## Repository layout

- `src/dark_siren_robustness/`: reusable Python utilities.
- `tests/`: unit and sanity tests.
- `docs/research_questions.md`: research questions and working hypotheses.
- `docs/statistical_model.md`: model components that must be specified.
- `docs/validation_plan.md`: validation and calibration plan.
- `.github/workflows/tests.yml`: automated test workflow.

## Setup

Python 3.10 or newer is required.

```bash
python -m venv .venv
```

On Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pytest
```

On macOS/Linux:

```bash
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pytest
```

## Current implementation

The initial scaffold contains:
- a low-redshift luminosity-distance approximation;
- a numerical flat-ΛCDM luminosity-distance utility;
- an illustrative distance-dependent catalogue-inclusion probability;
- basic unit tests for these utilities.

These utilities are not a complete dark-siren likelihood or inference pipeline. The toy completeness function is not calibrated to any real galaxy survey. The cosmology implementation requires independent numerical validation over the intended redshift range.

## Planned milestones

1. Audit the original dark-siren project and the reported nominal-68% interval coverage result.
2. Validate simple likelihoods and known-answer cases.
3. Specify and test simulation-based calibration and repeated-sampling coverage.
4. Run controlled experiments varying catalogue incompleteness and host-weight misspecification.
5. Evaluate any proposed correction against baselines using pre-specified metrics.
6. Publish reproducible results and a manuscript only after the analysis supports the conclusions.

## Citation

There is no release or research result to cite yet. A citation file and formal release can be added when authorship, versioning, and the intended citation are settled.
