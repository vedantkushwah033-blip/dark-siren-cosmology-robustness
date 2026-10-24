# Factorial sensitivity experiment

This experiment separates three dimensions using the original simulator pinned
to commit `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`:

- Catalogue: complete versus 70% target completeness.
- Missingness mechanism: random removal versus redshift-dependent removal.
- Host population versus inference weights: matched uniform, matched luminosity,
  or luminosity-weighted host generation analyzed with uniform inference weights.

The nine scenarios are defined in `scripts/run_factorial_sensitivity.py`.
Each scenario uses 500 mock universes by default and writes scenario-level
metrics plus per-run results. The optional manual workflow can raise the
replicate count.

## Why this design

A completeness change alone cannot tell us whether any observed shift is
caused by which galaxies are missing or by how host galaxies are weighted.
The matrix includes complete controls, matched-weight controls, and a
deliberate host-generation/inference-weight mismatch. The redshift-dependent
selection mechanism tests whether the pattern of missing galaxies matters,
not just the nominal completeness fraction.

## Run

Open:
https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/workflows/factorial-sensitivity.yml

Select **Run workflow**, initially keep replicates at `500` and base seed at
`2026`, then download the `factorial-sensitivity` artifact.

## Outputs

- `scenario_summary.csv`: one row per scenario, with coverage, interval
  widths, bias, RMSE, pulls and realized catalogue completeness. For these fixed-size masks, the realized
completeness should be exactly 0.70 in the 70% scenarios (350/500), so this
column is also a useful integrity check.
- `per_run_results.csv`: each simulated run for audit and plotting.
- `scenario_summary.json`: settings, pinned source commit and metrics.

## Interpretation cautions

These are synthetic experiments, not observations and not a measurement of
the real Hubble constant. The simulation's host population and catalogue are
generated from the same underlying mock universe, and galaxy luminosity/mass
proxies are stylized. All selection mechanisms in this matrix retain exactly 350 of 500 galaxies
at 70% target completeness; redshift-dependent selection randomizes *which*
galaxies are retained, favoring lower redshifts. The host can still be absent
from the observed catalogue. No
missing-host likelihood term is enabled in this first matrix; it is an
intentional diagnostic of the naive observed-catalogue inference. A follow-up
should compare this with the original model's oracle missing-host term, while
clearly labeling that term as a simulation-only validation control.

Do not write conclusions until the workflow has run and both summary and
per-run outputs have been checked.

## Independent output audit

The factorial workflow now runs `scripts/audit_factorial_sensitivity_outputs.py`
after producing the CSV/JSON files and before computing contrasts. The audit
checks the nine prespecified scenario labels, requested and realized
replicate counts, exact 70% completeness for incomplete-catalogue scenarios,
finite posterior summaries, credible-interval nesting, and consistency of
reported interval widths with their endpoints. Its unit tests include
deliberately corrupted outputs.

This audit is a structural/data-integrity check. It does not rerun the original
simulator independently, establish that its missingness mechanism is a
realistic model of galaxy catalogues, or prove frequentist coverage. Those
questions require separate simulation design and validation. A successful
workflow should therefore be described as internally checked synthetic
evidence, not as an observational result.
