# Dark-Siren Cosmology Robustness

A reproducible Bayesian simulation study of galaxy-catalogue incompleteness, host-population misspecification, and statistical calibration in dark-siren cosmology.

**Research status:** early-stage simulation and validation scaffold. The unit-test suite passes, but the main factorial experiment has not yet produced validated scientific results. Do not treat this repository as a real gravitational-wave measurement or a validated (H_0) result.

## Research question

Under which controlled conditions do galaxy-catalogue incompleteness and incorrect host-population assumptions affect the bias, uncertainty, and calibration of Bayesian inference for the Hubble constant? When does explicitly modelling missing hosts help in the controlled simulations?

## Scientific guardrails

- Synthetic results will not be described as measurements from real gravitational-wave observations.
- Galaxy-catalogue incompleteness and gravitational-wave detection selection are distinct effects and must not be conflated.
- Posterior coverage, simulation-based calibration, bias, and precision are different quantities.
- Working hypotheses are not findings; negative and inconclusive results will be reported.
- Oracle controls that use hidden simulated galaxies are simulation-only benchmarks, not directly available corrections for real surveys.
- No claim of novelty or publication readiness will be made without appropriate literature review and independent scrutiny.

## Repository layout

- `src/dark_siren_robustness/`: reusable numerical and statistical utilities.
- `tests/`: unit and numerical-validation tests.
- `scripts/`: reproducible simulation, artifact validation, and contrast-analysis scripts.
- `docs/`: model assumptions, design audits, experiment protocols, and interpretation limits.
- `.github/workflows/tests.yml`: automated test suite.
- `.github/workflows/factorial-sensitivity.yml`: manually triggered factorial experiment and control analyses.

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

The repository contains:
- a low-redshift luminosity-distance approximation, (d_L \approx cz/H_0);
- a numerical flat-ΛCDM luminosity-distance utility with independent quadrature tests;
- calibration and repeated-sampling coverage utilities;
- a nine-scenario catalogue-completeness/host-weighting experiment;
- missing-host and full-population oracle control experiments;
- bootstrap contrast analysis and artifact-integrity checks.

These pieces are not a complete observational dark-siren inference pipeline. The toy completeness and host-population mechanisms are synthetic and are not calibrated to a real galaxy survey or gravitational-wave detection pipeline.

## Running the main experiment

Open the [factorial-sensitivity workflow](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/workflows/factorial-sensitivity.yml), select **Run workflow**, choose branch `main`, and initially use:

- `replicates`: `500`
- `seed`: `2026`

The workflow runs nine primary scenarios, six missing-host control scenarios, and two full-population oracle scenarios. At 500 replicates per scenario, that is 8,500 per-run records across these experiments, before summary and contrast outputs. It may take several minutes.

When it finishes, download the `factorial-sensitivity` artifact and preserve the raw per-run CSV files, summaries, configuration, and contrast analyses. Check that the run succeeded and inspect the outputs before writing conclusions. The 70% catalogue scenarios use fixed-size masks retaining exactly 350 of 500 galaxies; this conditions on the retained count rather than modelling uncertainty in the count of missing galaxies.

See [the factorial experiment protocol](docs/factorial_sensitivity_protocol.md), [the design audit](docs/factorial_design_audit.md), and [the missing-host model audit](docs/missing_host_model_audit.md) for details and limitations.

## Cosmological distance-model extension

A flat-ΛCDM distance function and numerical tests have been added. The planned model-robustness experiment is still protocol-only: first run and validate the existing factorial experiment, then implement and run matched flat-ΛCDM and deliberate distance-model-mismatch scenarios. The initial protocol fixes Ωm = 0.3; it does not marginalize over matter-density uncertainty. See [the cosmological-model robustness protocol](docs/cosmological_model_robustness_protocol.md). The separate [selection-normalization protocol](docs/selection_normalization_protocol.md) defines a toy detection-selection extension. The analytic implementation and numerical tests are in `src/dark_siren_robustness/selection.py`. The current repository also contains [selection-conditioned SBC results](docs/selection_sbc_audit_results.md), [paired fixed-truth results](docs/selection_paired_comparison_results.md), [prior sensitivity results](docs/selection_prior_sensitivity_results.md), [nonzero lower-threshold results](docs/nonzero_lower_threshold_results.md), and [lower-threshold sensitivity results](docs/lower_threshold_sensitivity_results.md). A new fixed-truth coverage map is running to examine pointwise coverage across the prior support. These are controlled toy-model audits, not a complete observational dark-siren inference pipeline or a real detector model.

## Interpretation limits

The current simulations use a simplified low-redshift relation in the original inference experiment, a Gaussian distance likelihood, a finite (H_0) grid, and a generated galaxy population. Redshift-dependent missingness is a designed stress test, not a calibrated survey selection function. A passing test suite establishes only the software properties actually covered by tests; it does not establish scientific realism or the truth of any hypothesis.

## Citation

There is no formal release or validated research result to cite yet. A citation file and formal release can be added once the authorship, version, and intended citation are settled.
