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

A flat-ΛCDM distance function and numerical tests have been added. The planned model-robustness experiment is still protocol-only: first run and validate the existing factorial experiment, then implement and run matched flat-ΛCDM and deliberate distance-model-mismatch scenarios. The initial protocol fixes Ωm = 0.3; it does not marginalize over matter-density uncertainty. See [the cosmological-model robustness protocol](docs/cosmological_model_robustness_protocol.md). The separate [selection-normalization protocol](docs/selection_normalization_protocol.md) defines a toy detection-selection extension. The analytic implementation and numerical tests are in `src/dark_siren_robustness/selection.py`. The current repository also contains [selection-conditioned SBC results](docs/selection_sbc_audit_results.md), [paired fixed-truth results](docs/selection_paired_comparison_results.md), [prior sensitivity results](docs/selection_prior_sensitivity_results.md), [nonzero lower-threshold results](docs/nonzero_lower_threshold_results.md), and [lower-threshold sensitivity results](docs/lower_threshold_sensitivity_results.md), [fixed-truth coverage-map results](docs/fixed_truth_coverage_map_results.md), and [single-host SBC results](docs/single_host_selection_sbc_results.md). The fixed-truth coverage map and single-host selection-conditioned SBC control have now run; see [fixed-truth coverage-map results](docs/fixed_truth_coverage_map_results.md) and [single-host SBC results](docs/single_host_selection_sbc_results.md). The adaptive-quadrature normalization check is documented in [selection quadrature audit results](docs/selection_quadrature_audit_results.md). The multi-resolution comparison is implemented in [the grid-refinement script](scripts/run_selection_grid_refinement.py), with results emitted as a schema-checked CSV and tested in [the grid-refinement workflow](.github/workflows/selection-grid-refinement.yml). and automated by [the quadrature audit workflow](.github/workflows/selection-quadrature-audit.yml). These are controlled toy-model audits, not a complete observational dark-siren inference pipeline or a real detector model.


## Prior-support × detection-threshold factorial experiment

A separate paired experiment varies three prior intervals, four assumed detection thresholds, three fixed synthetic truths, and corrected versus uncorrected inference. Read the [prespecified protocol](docs/prior_threshold_factorial_protocol.md) before interpreting any output. The full protocol run uses 300 replicates per truth, a 4,001-point posterior grid, and seed 20261028; it is distinct from the older factorial-sensitivity workflow above.

To launch the full factorial run, open the [full prior-threshold factorial workflow](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/workflows/full-prior-threshold-factorial.yml), choose **Run workflow**, and use the protocol settings. The workflow is manually dispatched; editing or committing code does not launch the full experiment. A successful run should generate 21,600 raw rows, 72 summary rows, and 396 paired-contrast rows, plus `experiment_manifest.json`. These are expected counts, not claims that a full run has already completed. Inspect the run conclusion and download the artifacts before reporting results.

## Interpretation limits

The current simulations use a simplified low-redshift relation in the original inference experiment, a Gaussian distance likelihood, a finite (H_0) grid, and a generated galaxy population. Redshift-dependent missingness is a designed stress test, not a calibrated survey selection function. A passing test suite establishes only the software properties actually covered by tests; it does not establish scientific realism or the truth of any hypothesis.

The fixed-truth coverage map uses a bounded uniform prior (H_0\in[50,90]). When the true value is near either edge, posterior credible intervals are constrained by the prior support; pointwise frequentist coverage can therefore be poor even when the selection-conditioned likelihood is implemented correctly. This is a property of the bounded inferential setup that must be reported, not hidden by changing the data-generating process. The prior-predictive SBC audit and fixed-truth coverage map answer different questions and must not be treated as interchangeable validation tests. Coverage fractions are Monte Carlo estimates; consult their Wilson intervals before interpreting small differences between methods.

The Wilson interval helper explicitly returns exact boundary values for zero successes and all successes. Independent adaptive quadrature checks are being added to compare posterior normalization from grid integration against a separate numerical integration route. Their status must be verified from the corresponding GitHub Actions run before claiming that check passed.

## Citation

There is no formal release or validated research result to cite yet. A citation file and formal release can be added once the authorship, version, and intended citation are settled.


### Recent validation notes

- [Validation status as of 2026-10-24](docs/validation_status_2026-10-24.md): distinguishes verified software checks from outstanding experiment and scientific-validation work.
- [Project use and intended audience](docs/project_use_and_audience.md): explains appropriate use cases and explicitly states what the synthetic study cannot establish.
- [Raw grid-refinement artifact reanalysis](docs/grid_refinement_artifact_reanalysis_2026-10-24.md): inspects the archived 21,600-fit artifact, quantifies paired discrepancies across grid sizes, and documents the extreme upper-95% tail case and targeted finer-grid follow-up.
- [Factorial experiment interpretation plan](docs/factorial_experiment_interpretation_plan.md): predefines how bias, RMSE, interval coverage, widths, paired contrasts, and Monte Carlo uncertainty will be interpreted without overstating the synthetic results.

The raw-artifact reanalysis finds that typical posterior-median discrepancies are small, but tail credible endpoints can be more grid-sensitive. In the worst identified case, the upper 95% endpoint changes by 0.13045 between 1,001 and 4,001 grid points; a targeted 64,001-point follow-up approaches approximately 60.710 for that specific synthetic observation. This is not a universal error bound. Do not describe the grid as uniformly converged until prespecified accuracy criteria have been assessed across relevant cases.

### Posterior grid refinement

The completed deterministic grid-refinement audit and its interpretation are documented in [grid-refinement results](docs/selection_grid_refinement_results.md). The workflow compared 1,001-, 4,001-, and 8,001-point posterior grids for three fixed toy observations; it should not be interpreted as a universal error guarantee.


### Prior-support sensitivity

A paired sensitivity experiment reuses the same detected synthetic observations across three prior intervals. See the [prior-support sensitivity protocol](docs/prior_support_sensitivity_protocol.md); numerical conclusions will be added only after the CI artifacts are inspected.


### Prior-support sensitivity results

The completed paired prior-support audit found substantial sensitivity near the hard prior boundaries. See [prior-support sensitivity results](docs/prior_support_sensitivity_results.md) for coverage with Wilson intervals, bias, RMSE, and limitations. These are synthetic-model results and do not establish an optimal prior.


### Posterior interval-width sensitivity

The companion audit quantifies changes in 68% and 95% posterior interval widths across prior supports, paired on identical observations. See the [interval-width sensitivity protocol](docs/prior_interval_width_sensitivity_protocol.md). Results are provisional until the corresponding CI artifacts are inspected.


Validated interval-width results are documented in [prior interval-width sensitivity results](docs/prior_interval_width_sensitivity_results.md). The experiment and artifact checks passed on Python 3.10, 3.11, and 3.12; results remain specific to the synthetic model.


The [joint prior-performance scorecard protocol](docs/joint_prior_scorecard_protocol.md) combines coverage, bias, RMSE, and credible-interval width across prior supports. Its results remain specific to the synthetic model.


### Joint prior-performance scorecard results

The completed joint scorecard compares bias, RMSE, fixed-truth coverage with Wilson intervals, and posterior interval widths across the tested prior supports and both inference methods. See [validated scorecard results](docs/joint_prior_scorecard_results.md). The scorecard workflow passed on Python 3.10, 3.11, and 3.12; these results are synthetic-model diagnostics, not observational H0 measurements.


### Coverage-transition map

A denser paired fixed-truth experiment tests nine true values from 52 to 88 and compares the original [50,90] prior with [40,100], reusing each generated observation across prior supports. It reports bias, RMSE, 68%/95% coverage with Wilson intervals, and interval widths. See [coverage-transition protocol](docs/coverage_transition_map_protocol.md) and [validated results](docs/coverage_transition_map_results.md). This remains a synthetic diagnostic, not an observational H0 analysis.


### Paired method contrasts with bootstrap uncertainty

The coverage-transition experiment now reports corrected-minus-uncorrected paired differences for posterior medians, absolute/squared errors, 68%/95% interval widths, and coverage indicators, with percentile-bootstrap 95% intervals. See the [paired contrast results](docs/paired_method_contrast_results.md). The workflow uploads `paired_method_contrasts.csv` with the per-run and summary artifacts. All findings are conditional on the synthetic model and should not be generalized directly to real LVK data.


### Detection-threshold misspecification

A new paired experiment generates detected observations under a 120 Mpc toy threshold and compares selection-corrected inference using the correct threshold against deliberately incorrect assumed thresholds (150, 180, and 210 Mpc). The [protocol](docs/threshold_misspecification_protocol.md) documents the design and paired bootstrap contrasts; CI uploads per-run, summary, and paired-contrast CSV artifacts. This remains a synthetic model, not a real detector-selection function.


**Threshold misspecification results:** [recomputed performance summaries and paired-bootstrap contrasts](docs/threshold_misspecification_results.md). These are synthetic-model results, not real LVK measurements. The report is based on the successful per-run artifact; the workflow upload-path fix is being validated separately.


### Prior-support × detection-threshold factorial audit

The [factorial protocol](docs/prior_threshold_factorial_protocol.md) defines a paired experiment varying three prior intervals, four assumed detection thresholds, three true H0 values, and both corrected/uncorrected methods. It reports coverage Wilson intervals and paired bootstrap contrasts. The experiment is a synthetic toy-model sensitivity analysis, not a real detector model; results should be interpreted only after inspecting the generated artifacts and CI status.


### Prior-support × detection-threshold factorial artifact review

The compact CI artifacts have been downloaded and audited for row counts, expected configuration coverage, finite numeric values, and interpretation risks. The artifact review documents selected provisional results and limitations; the compact run is a smoke test, not the full 300-replicate analysis. See [artifact review and compact-run diagnostics](docs/prior_threshold_factorial_results.md). Do not present the compact CI results as final or observational H0 measurements.


A follow-up [paired-contrast interpretation audit](docs/prior_threshold_factorial_contrast_interpretation.md) explains the uncorrected-threshold invariance, selected paired differences, zero-width bootstrap intervals, multiplicity caveats, and limitations of the compact run.


### Paired-bootstrap runtime improvement

The factorial experiment's percentile bootstrap now resamples paired differences in bounded vectorized batches of 256 draws. This preserves the paired-resampling definition while avoiding a Python-level loop for every bootstrap draw and limiting temporary array memory. The random-number stream and exact numerical bootstrap endpoints can differ from earlier script revisions because vectorized sampling consumes random draws differently; comparisons must therefore be tied to the script commit and seed, not assumed bitwise identical across revisions. The change does not alter the experiment's synthetic-model scope or increase the effective sample size.


### Additional regression tests for the factorial audit

The focused tests now verify that all uncorrected posterior summaries (median and 68%/95% interval endpoints) remain invariant across assumed-threshold labels, and exercise paired-bootstrap determinism on a small complete matched fixture. These checks complement the end-to-end compact experiment test and guard against accidental changes to the paired design.


### Full-resolution factorial workflow

A manually dispatched workflow is available at `.github/workflows/full-prior-threshold-factorial.yml`. From GitHub's **Actions** tab, select **Full prior-threshold factorial experiment** and choose 300 replicates, a 4,001-point grid, and seed `20261028` for the protocol-default run. The workflow runs focused regression tests, executes the factorial script, independently recomputes summary metrics/Wilson intervals and paired mean contrasts from the raw rows using `scripts/validate_prior_threshold_factorial.py`, checks expected CSV row counts and finite numeric values, and uploads all three CSV artifacts for 90 days. It is not automatically dispatched by this documentation commit; a user must explicitly start the workflow in GitHub Actions. The workflow's successful completion and its artifact contents must be verified before the full run is described as complete.


### Full numerical grid-convergence workflow

The manually dispatched workflow `.github/workflows/grid-refinement-full.yml` runs the matched-observation grid audit on 1,001-, 2,001-, and 4,001-point posterior grids. From **Actions → Grid refinement numerical convergence audit**, choose 200 replicates and seed `20261030` for the protocol default. It runs the focused test, checks expected row counts and finite numeric values, and uploads per-run and cell-summary CSVs for 90 days. This workflow is independent of the prior-threshold factorial workflow; passing one does not establish convergence in the other analysis.


### Grid-refinement regression hardening

The grid-refinement tests additionally check the ordering of posterior credible-interval endpoints, finiteness of all five reported posterior summaries, reuse of the exact same simulated distance observation across grid resolutions and weighting choices, and consistency of median/95th-percentile/maximum discrepancy summaries. These are structural integrity checks; they do not assert that the grid error is scientifically negligible. Numerical adequacy must be judged from the actual full-run discrepancy distributions against a prespecified tolerance.


The completed 200-replicate grid-refinement run is documented in [the full convergence results and interpretation](docs/grid_refinement_full_results.md). Typical discrepancies were small, but one upper-95%-endpoint maximum was materially larger than the rest and requires case-level inspection before claiming uniform convergence.


### Targeted upper-tail convergence diagnostic

The worst upper-95% endpoint discrepancy from the archived grid audit can be reproduced independently with:

```bash
python scripts/run_targeted_grid_convergence.py --output results/targeted-grid-convergence.csv
```

This fixed synthetic observation is evaluated on 1,001 through 64,001 grid points, with posterior summaries written to a CSV. A focused regression test checks credible-interval ordering and the decrease in successive upper-endpoint changes. The test is a case-specific numerical diagnostic, not a universal convergence guarantee; the full design must still be evaluated against a tolerance specified before examining results. The script's output should be retained as an artifact before making numerical-accuracy claims.


### Prospective grid-accuracy criteria

The independent numerical acceptance protocol is documented in [`docs/prospective_grid_accuracy_protocol.md`](docs/prospective_grid_accuracy_protocol.md). It defines cell-level percentile limits, a per-case guardrail, targeted refinement requirements, and normalization checks. These criteria were written after the previous audit, so they must be applied to a new run and must not be described as preregistered for the earlier results.


### Independent full-design grid-accuracy result

The 200-replicate audit with seed 20261101 completed for 36 design cells and passed all 180 cell-level 95th-percentile checks for the five posterior summaries; no individual discrepancy exceeded the 0.05 refinement guardrail. See the artifact-backed [independent grid-accuracy assessment](docs/independent_grid_accuracy_assessment_seed20261101.md) and [GitHub Actions run](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40091695141). The normalization-checked rerun also recomputed posterior area on the prior support for all 21,600 fits using the same trapezoidal rule as the inference routine (maximum absolute area error 4.44e-16; zero failures at 1e-8 tolerance). This is a consistency check, not independent-quadrature validation. These checks do not establish repeated-sampling coverage, prior robustness, or scientific/astrophysical realism.
