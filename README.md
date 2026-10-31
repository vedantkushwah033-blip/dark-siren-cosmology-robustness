# Dark-Siren Cosmology Robustness

> **Research status (29 October 2026; controlled synthetic methodology only):** The locked 500-replicate Stage B experiment independently passed its posterior-grid and raw-table audits; its primary paired absolute-error reduction was 5.9519 km s⁻¹ Mpc⁻¹ (bootstrap 95% interval [5.5410, 6.3573]). A distinct 13-setting KDE bandwidth/grid/selection screen passed its full audit: 1,560/1,560 posterior fits were independently reconstructed, 403 aggregate fields were rederived with zero discrepancy, and 14 manifest-listed files were checked. One 68% coverage indicator differed only at an interval boundary 9.37×10⁻⁵ from truth and is explicitly disclosed as boundary-sensitive. A separate 1,200-dataset exact-density oracle audit also passed. The configurable catalogue-selection study adds a 200-replicate baseline pilot and three predeclared conditions; the larger follow-up completed **250 paired replicates per condition (750 total)**, with all six inference arms retained and each condition's independent posterior-grid/aggregate audit and artifact upload passing in [run #40819870231](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40819870231). The paired absolute-posterior-mean-error reductions for correct known-law weighting versus catalogue-only inference were 7.436 [6.798, 8.053], 2.246 [1.711, 2.764], and 2.354 [1.779, 2.923] km s⁻¹ Mpc⁻¹ for the low-inclusion-floor, high-inclusion-floor, and broad-transition regimes. These are results conditional on simplified synthetic models—not real LVK/GWTC data, a real H₀ measurement, or proof of detector realism. The literature audit now includes Gaussian-process and clustering-based catalogue reconstruction, the 2026 practical `gwcosmo` completeness test, host-population inference, and clustering-preserving galaxy completion. **The exact CHIMERA v2.2 external benchmark remains NO-GO** until compatible, versioned input HDF5 products are obtained and executed with provenance; no older-version data are substituted. The source-bearing commit `1832dd8bbe562cbc95442f6c97ff3f73d1b25b2c` passed the Python 3.10/3.11/3.12 numerical-diagnostics matrix in push run [#40865223580](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40865223580) and PR-triggered run [#40865232383](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40865232383). Subsequent documentation-only commits are not covered by that earlier green result; verify the exact PR head in Actions before release. PR #7 remains open and unmerged and currently targets `experiment/catalogue-density-estimators`, not the default `main` branch; review that base-branch choice before merging.

- **Dark-siren literature positioning update (31 October 2026):** [recent GP line-of-sight reconstruction, variance completion, clustering-aware methods, pipeline validation, weighting comparisons, novelty boundaries, and a matched-method calibration roadmap](docs/dark_siren_literature_positioning_update_2026-10-31.md). Literature/design note only; no new numerical result is claimed.
- **Linear-to-nonlinear validation bridge (31 October 2026):** [claim ladder and separate validation requirements for finite-population means, catalogue-validation uncertainty, and nonlinear dark-siren posterior functionals](docs/linear_to_nonlinear_validation_bridge_2026-10-31.md). The v1 linear pilot's outer artifact audit passed (3 frames, 1,800 method attempts, 18 summary cells); its inner-bootstrap attempt serialization gap and required v2 release gate are documented in the [v2 remediation protocol](docs/inner_bootstrap_attempt_audit_gap_and_v2_protocol_2026-10-31.md). Full pilot results and artifact fingerprints are in the [audited pilot report](docs/audited_linear_finite_population_bootstrap_pilot_results_2026-10-31.md).\n- **Estimated-inclusion nonlinear stress protocol (30 October 2026):** [prespecified curvature × validation-size design, diagnostics, estimands, and stop/go gates](docs/estimated_inclusion_nonlinear_stress_protocol_2026-10-30.md). Protocol only; no factorial results are claimed.
- **Linear finite-population bootstrap pilot manuscript:** [exploratory methods note, v1 audited results, v2 auditability changes, conditional versus unconditional coverage, and explicit limitations](manuscript/linear_finite_population_bootstrap_pilot_note.md). The v1 numerical pilot is reported; v2 remains pending exact-head CI, artifact generation, and independent audit. The multi-frame design remains prospective.
- **Confirmatory multi-frame preregistration protocol:** [hierarchical estimands, precision-driven frame/repetition planning, frame-aware uncertainty, multiplicity, stopping rules, and release gates](docs/confirmatory_multiframe_preregistration_protocol_2026-10-31.md). Prospective protocol only; no confirmatory outcomes are claimed.
- **Executable frame-precision planner:** `src/dark_siren_robustness/frame_precision_planning.py` computes a planning-only frame count from an externally justified pilot SD and target half-width, using an upper confidence bound for the SD. Unit tests run in CI; it does not estimate pilot variability or lock the confirmatory design.
- **Pre-posterior weight-stability diagnostic protocol (30 October 2026):** [leakage-safe predictors, held-out evaluation, baselines, calibration, and stop/go gates](docs/preposterior_weight_stability_diagnostic_protocol_2026-10-30.md). Protocol only; no predictive performance is claimed.
- **Integrated draft manuscript:** [controlled synthetic catalogue-selection benchmarks and H0-dependent-noise calibration, draft v0.3](docs/draft_manuscript_2026-10-29.md). This is an internal draft, not peer reviewed.
- **Research-readiness report:** [artifact-verified achievements, open CI/data gates, and safe claim language](docs/research_readiness_report_2026-10-29.md).
- **H0-dependent-noise SBC report:** [results from 1,200 seeded synthetic datasets, independent replay audit, and fixed-truth calibration](docs/h0_dependent_noise_sbc_results_2026-10-29.md). The run passed independent reconstruction at [GitHub Actions run #40796895426](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40796895426).
- **Literature and claim-positioning audit:** [GP and clustering-based catalogue reconstruction, practical completeness testing, host-population inference, and CHIMERA v2.2 status](docs/literature_and_positioning_audit_2026-10-29.md).
- **Host-density validation campaign:** [protocol and audit status](docs/host_density_validation_methods_2026-10-29.md); campaign code in `scripts/analyze_host_density_campaign.py`; campaign-wide independent simulator/KDE/posterior audit in `scripts/audit_host_density_campaign.py`. All 13 configurations passed local per-configuration reconstruction checks (1,560 fit rows); the clean CI reproduction and all-row audit passed at [#40802415000](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40802415000), independently reconstructing all 1,560 fits.
- **Host-density campaign results:** [canonical 13-setting screening report](docs/host_density_sensitivity_screening_results_2026-10-29.md) records fixed-seed results, hashes, audit scope, uncertainty, and claim limits. The old local-snapshot report is now explicitly marked historical/superseded.
- **Exact-density oracle calibration:** [prior-averaged versus fixed-truth bias](docs/oracle_bias_diagnosis_2026-10-29.md). The current run #40802415000 saves 1,200 raw rows and its independent oracle auditor replayed every seed/distance, reconstructed every posterior and aggregate, and verified hashes; that oracle job is **PASS**. A separate additive H0-dependent-noise prototype is documented in [the selection-conditioned model note](docs/h0_dependent_noise_selection_model.md); it is not the likelihood used by the locked Stage B run.
- **Stage A/Stage B selection-misspecification study:** [locked 500-replicate result and paired contrasts](docs/selection_misspecification_stage_b_report_2026-10-29.md). Stage B's independent posterior-grid and raw-table paired-analysis artifacts passed at run [#40795993174](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40795993174). The primary result applies only to the declared fixed-(H_0) synthetic model.
- **Configurable selection-function study and larger replication:** [complete pilot, exploratory screen, and 750-replicate audited replication results](docs/selection_function_stress_screen_results_2026-10-29.md), with the larger replication predeclared in [the protocol](docs/catalogue_selection_confirmatory_protocol_2026-10-29.md). The separate three-condition 250-replicate-per-condition confirmation is documented in [the confirmatory replication report](docs/catalogue_selection_confirmatory_replication_2026-10-29.md), including independently rederived paired contrasts, bootstrap intervals, checksums, and scope limits. All three 250-replicate conditions passed the independent artifact audit and upload in [workflow #40819870231](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40819870231). The paired absolute-error reduction from catalogue-only to correct weighting was 7.44 [6.80, 8.05], 2.25 [1.71, 2.76], and 2.35 [1.78, 2.92] km s⁻¹ Mpc⁻¹ across low-floor, high-floor, and broad-transition settings. These remain results conditional on a simplified synthetic model; the shallow-slope misspecified arm sometimes had lower RMSE than correct weighting.

- **Independent confirmatory reproduction and truth-dependent bias-compensation follow-up:** [recomputed results and limitations](docs/independent_reproduction_confirmatory_2026-10-29.md), [interpretation protocol](docs/bias_compensation_followup_protocol_2026-10-29.md), and paired bootstrap utility `scripts/reanalyze_bias_compensation_followup.py` with unit tests. The 750 confirmatory replicates have overlapping seed ranges across conditions (150 low/high, 150 high/broad, 50 low/broad); within-condition paired contrasts remain valid, but the three conditions are not independent replications.
- **Pinned CHIMERA v2.2 execution contract:** [notebook configuration, input inventory, and go/no-go gate for the incomplete-catalogue benchmark](docs/chimera_v2_2_notebook_execution_contract_2026-10-26.md).

- **Known-host calibration protocol:** [design, estimands, pairing, audit gates, and interpretation limits](docs/known_host_calibration_protocol.md).
- **Known-host calibration analysis plan:** [prespecified cell summaries, paired comparisons, Monte Carlo uncertainty, and interpretation rules](docs/known_host_calibration_analysis_plan.md).
- **Known-host calibration results:** [audited 10,800-run experiment, boundary-aware coverage results, and paired contrasts](docs/known_host_calibration_run_20261025.md).
- **Known-host calibration synthesis and next decision gate:** [cross-experiment interpretation, boundary/prior caveats, and prespecified next steps](docs/known_host_calibration_synthesis_2026-10-25.md).
- **Prospective joint catalogue-incompleteness × event-selection protocol:** [staged model-compatibility gate, factorial design, estimands, and independent audit requirements](docs/joint_catalogue_incompleteness_detection_selection_protocol.md).
- **Joint likelihood derivation checkpoint:** [separating event-selection normalization from the latent missing-host population, with compatibility tests and a go/no-go gate](docs/joint_selection_incompleteness_likelihood_derivation_2026-10-25.md).
- **Joint-model prototype:** `posterior_h0_selected_with_missing_hosts` combines a finite missing-host mixture with the same selection denominator used for the complete host mixture; limiting-case tests are in [test_joint_selection_missing_hosts.py](tests/test_joint_selection_missing_hosts.py). This is a synthetic finite-mixture prototype, not a validated real-survey correction.
- **Selected-likelihood calibration runner:** `scripts/run_selected_likelihood_sbc.py` produces 1,000 prior-predictive PIT records, decile Wilson intervals, mean-PIT Monte Carlo standard error, and a descriptive KS statistic/p-value. The KS p-value is not treated as proof of calibration; the raw records are archived by the numerical-diagnostics workflow. `scripts/audit_selected_likelihood_sbc.py` independently reconstructs each PIT from the raw truth/host/distance records, replays the seeded RNG protocol to verify truth/host/distance/rejection counts, validates finite and physically ordered numerical metadata, recomputes prior-draw and PIT summary statistics plus decile intervals, and rejects tampered/inconsistent records before CI artifacts are uploaded. The completed calibration workflow passed on Python 3.10, 3.11, and 3.12; the audit establishes internal synthetic consistency only, not astrophysical validity. Regression tests also challenge malformed prior support, integration metadata, host redshifts, noise, and detection bounds.
- **Paired host-weight sweep protocol:** [prespecified prior-support factor, paired estimands, and audit gates](docs/paired_weight_sweep_protocol.md).
- **Paired host-weight sweep results:** [18-cell independently audited results showing prior-support sensitivity](docs/paired_weight_sweep_results_20261025.md).
- **Prior-measure sensitivity protocol:** [paired uniform-H₀ versus log-uniform-H₀ host-component evidence experiment and independent audit gates](docs/prior_measure_sensitivity_protocol.md).
- **Prior-measure sensitivity results:** [9,000-observation experiment, all 36 paired contrasts, and independent audit](docs/prior_measure_sensitivity_results_20261025.md).
- **Posterior-level calibration protocol:** [multi-truth H₀ bias, RMSE, interval-width, and coverage design](docs/posterior_calibration_protocol.md).
- **Posterior-calibration raw provenance audit:** `scripts/audit_posterior_calibration_generation.py` independently replays the fixed-seed host and distance generation before posterior reconstruction; tampered raw observations are rejected even if their CSV hash is refreshed.
- **Paired host-weight contrast audit:** `scripts/audit_paired_host_weight_contrasts.py` independently reconstructs all 48 cell-level and 24 aggregate paired contrasts, including stratified-bootstrap intervals, from the audited raw posterior rows.
- **Posterior-level calibration results:** [2,400-observation experiment with independently audited bias, RMSE, interval width, and coverage](docs/posterior_calibration_results_20261025.md).
- **Distance-uncertainty likelihood sensitivity protocol and interpretation limits:** [three explicitly defined likelihood arms, oracle-like benchmark caveat, generator mismatch, prior-boundary effects, and scope](docs/distance_uncertainty_model_sensitivity_interpretation.md).
- **Distance-uncertainty full-run results:** [32,400 posterior records, independent reconstruction audit, truth-dependent coverage, stratified paired A-minus-B contrasts, and grid-resolution check](docs/distance_uncertainty_full_run_results.md). Reproduce the run with `scripts/run_distance_uncertainty_sensitivity.py`, independently audit it with `scripts/audit_distance_uncertainty_sensitivity.py`, summarize paired contrasts with `scripts/summarize_distance_uncertainty_paired_contrasts.py`, and recompute all posterior summaries on a 3,201-point grid with `scripts/check_distance_uncertainty_grid_stability.py`.
- **Host-population misspecification stress test:** [independent audit protocol](docs/host_misspecification_independent_audit_protocol_2026-10-26.md) and `scripts/run_host_misspecification_audit.py` / `scripts/audit_host_misspecification_audit.py`. The audit reconstructs all synthetic fits and summary metrics, and now reports 95% Wilson intervals for empirical coverage to expose finite-replicate Monte Carlo uncertainty. These intervals are scenario-wise uncertainty summaries, not paired scenario-comparison tests; all conclusions remain limited to the declared toy model.
- **Latest host-misspecification run:** [1,000 replicates per scenario, 4,000 independently reconstructed fits, results and interpretation limits](docs/host_misspecification_results_2026-10-26.md); [workflow run and raw artifact](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40278173997).
- **Prior-support × hidden-host interaction protocol:** [prespecified factors, paired estimands, negative controls, and audit gates](docs/host_prior_support_factorial_protocol.md).
- **Prior-support × hidden-host interaction results:** [86,400 posterior fits, independent reconstruction, equal-cell paired-bootstrap intervals, and the prior/coverage trade-off](docs/prior_support_hidden_host_interaction_results_2026-10-26.md); [successful workflow and artifact](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40283661686).
- **Hidden-host assumption-sensitivity protocol:** [paired design for incorrect hidden fraction/redshift, including a zero-hidden-population negative control](docs/hidden_host_assumption_sensitivity_protocol.md).
- **Hidden-host assumption-sensitivity results:** [105,600 fits, independent reconstruction, and paired-bootstrap effects of incorrect assumed fractions/redshifts](docs/hidden_host_assumption_sensitivity_results_2026-10-26.md); [successful workflow and raw artifact](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40284995493).
- **Literature and novelty audit:** [targeted review of established incompleteness/host-weighting methods, explicit claim limits, and gates before any originality claim](docs/literature_and_novelty_audit_2026-10-26.md).
- **Published-pipeline benchmark plan:** [CHIMERA/gwcosmo data prerequisites, reproducibility risks, and go/no-go criteria for a fair baseline comparison](docs/reference_pipeline_benchmark_plan_2026-10-26.md).
- **CHIMERA benchmark target/version audit:** [blocking mismatch between the pinned v2.0 source tree, current-main incomplete-catalogue example, and required data products](docs/chimera_benchmark_target_audit_2026-10-26.md).

A reproducible Bayesian simulation study of galaxy-catalogue incompleteness, host-population misspecification, and statistical calibration in dark-siren cosmology.

**Research status (updated 2026-10-27):** The direct missing-host mixture likelihood was corrected after a transient extra catalogue-mass factor was found in a log-space rewrite. The four additional directly identified hidden-host experiments have now been rerun from corrected revision `ecb260886716c54fdd2d5e47390701949f5a18d5`; their raw CSV artifacts were compared row-by-row with historical artifacts, with all categorical fields identical and numerical differences limited to floating-point precision. General numerical-diagnostics tests passed on Python 3.10, 3.11, and 3.12 at that revision. See the [corrected rerun provenance report](docs/corrected_hidden_host_rerun_provenance_2026-10-27.md) and [impact audit](docs/missing_host_mixture_likelihood_impact_audit_2026-10-27.md). These checks establish reproducibility for the declared synthetic implementation, not astrophysical validity or novelty. The literature audit does not establish novelty; the CHIMERA/reference-pipeline benchmark remains blocked on verified matching inputs and provenance. No real gravitational-wave/catalogue data have been validated here, no observational H₀ result is claimed, and publication readiness has not been established.  Every result must remain tied to its exact source revision. The proposed next-stage design is documented in the [prospective diagnostic-comparison protocol](docs/prospective_diagnostic_comparison_protocol_2026-10-27.md); it is a protocol, not completed evidence.

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

A flat-ΛCDM distance function and numerical tests have been added. The matched flat-ΛCDM versus low-redshift distance-model experiment has been implemented and run, and the joint distance-model × host-mixture experiment has completed with an independent artifact audit. See [cosmological-model robustness results](docs/cosmological_model_robustness_results_2026-10-24.md), [joint distance-model and host-mixture results](docs/joint_distance_model_host_mixture_results_2026-10-24.md), and [posterior-grid resolution audit](docs/joint_grid_resolution_audit_results_2026-10-24.md). These are synthetic diagnostics; they do not constitute an observational dark-siren analysis. The initial protocol fixes Ωm = 0.3; it does not marginalize over matter-density uncertainty. See [the cosmological-model robustness protocol](docs/cosmological_model_robustness_protocol.md). The separate [selection-normalization protocol](docs/selection_normalization_protocol.md) defines a toy detection-selection extension. The analytic implementation and numerical tests are in `src/dark_siren_robustness/selection.py`. The current repository also contains the completed paired host-weight sweep with narrow and broad prior-support arms; its independently audited results show that weight sensitivity can be hidden when alternate-host H₀ modes are excluded by a narrow prior. See [the sweep protocol](docs/paired_weight_sweep_protocol.md) and [audited sweep results](docs/paired_weight_sweep_results_20261025.md). The current repository also contains [selection-conditioned SBC results](docs/selection_sbc_audit_results.md), [paired fixed-truth results](docs/selection_paired_comparison_results.md), [prior sensitivity results](docs/selection_prior_sensitivity_results.md), [nonzero lower-threshold results](docs/nonzero_lower_threshold_results.md), and [lower-threshold sensitivity results](docs/lower_threshold_sensitivity_results.md), [fixed-truth coverage-map results](docs/fixed_truth_coverage_map_results.md), and [single-host SBC results](docs/single_host_selection_sbc_results.md). The fixed-truth coverage map and single-host selection-conditioned SBC control have now run; see [fixed-truth coverage-map results](docs/fixed_truth_coverage_map_results.md) and [single-host SBC results](docs/single_host_selection_sbc_results.md). The adaptive-quadrature normalization check is documented in [selection quadrature audit results](docs/selection_quadrature_audit_results.md). The multi-resolution comparison is implemented in [the grid-refinement script](scripts/run_selection_grid_refinement.py), with results emitted as a schema-checked CSV and tested in [the grid-refinement workflow](.github/workflows/selection-grid-refinement.yml). and automated by [the quadrature audit workflow](.github/workflows/selection-quadrature-audit.yml). These are controlled toy-model audits, not a complete observational dark-siren inference pipeline or a real detector model.


## Prior-support × detection-threshold factorial experiment

A separate paired experiment varies three prior intervals, four assumed detection thresholds, three fixed synthetic truths, and corrected versus uncorrected inference. Read the [prespecified protocol](docs/prior_threshold_factorial_protocol.md) before interpreting any output. The full protocol run uses 300 replicates per truth, a 4,001-point posterior grid, and seed 20261028; it is distinct from the older factorial-sensitivity workflow above.

To launch the full factorial run, open the [full prior-threshold factorial workflow](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/workflows/full-prior-threshold-factorial.yml), choose **Run workflow**, and use the protocol settings. The workflow is manually dispatched; editing or committing code does not launch the full experiment. A successful run should generate 21,600 raw rows, 72 summary rows, and 396 paired-contrast rows, plus `experiment_manifest.json`. These are expected counts, not claims that a full run has already completed. Inspect the run conclusion and download the artifacts before reporting results.

## Posterior-grid numerical sensitivity

A reproducible grid-resolution audit now compares the same simulated observations at 2,001-, 4,001-, and 8,001-point (H_0) grids for all four distance/host-mixture arms. It reports posterior median and interval-endpoint differences relative to the 8,001-point reference grid. The workflow runs 30 replicates per design cell and archives the raw posterior summaries and comparison table. This tests numerical discretization sensitivity only; it does not add observational realism. See [the audit script](scripts/run_joint_grid_resolution_audit.py), [its tests](tests/test_joint_grid_resolution_audit.py), and [the automated workflow](.github/workflows/joint-grid-resolution-audit.yml). Results should be interpreted only after that workflow succeeds and its artifact is inspected.

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

### Full prior-threshold factorial: audited results (2026-10-24)

The full 300-replicate, 4,001-point prior-support × assumed-detection-threshold run completed successfully and its raw artifact was inspected. See [the artifact audit and quantitative results](docs/full_prior_threshold_factorial_results_2026-10-24.md) and [the successful Actions run](https://github.com/vedantkushwaha033-blip/dark-siren-cosmology-robustness/actions/runs/40098731627). The artifact contains 21,600 per-run rows, 72 summary rows, 396 paired contrasts, and a runtime manifest. Initial results show strong prior-boundary dependence and no uniform advantage for the selection-corrected method; the document reports reference-cell paired effects and caveats. These remain synthetic-model diagnostics, not observational H0 measurements or evidence that the toy model is astrophysically realistic.

### Cross-experiment scientific synthesis

The [cross-experiment synthesis of selection normalization and prior-boundary effects](docs/cross_experiment_selection_and_prior_synthesis_2026-10-24.md) reconciles the prior-predictive SBC controls with fixed-truth coverage and prior-support experiments, distinguishes supported toy-model conclusions from untested claims, and prioritizes the next validation experiments. This is a synthesis of existing results, not a new simulation or observational result.

### Known-host calibration surface (artifact-backed results and independent audits complete)

The [prespecified protocol](docs/known_host_calibration_surface_protocol.md), [known-host calibration script](scripts/run_known_host_calibration_surface.py), [independent raw-record auditor](scripts/audit_known_host_calibration_surface.py), and [focused tests](tests/test_known_host_calibration_surface.py) plus [auditor integration tests](tests/test_audit_known_host_calibration_surface.py) isolate a single supplied host and vary nine fixed H0 truths across four fractional distance-noise levels. The analysis reports fixed-truth coverage with Wilson intervals, posterior-median and posterior-mean errors, RMSE, interval widths, posterior pulls, and prior-edge diagnostics. The manually dispatchable [workflow](.github/workflows/known-host-calibration-surface.yml) checks out the simulator at a pinned commit, independently recomputes summaries from raw records, verifies paired noise draws and provenance, hashes the inputs, and uploads the audit report with the raw data. NumPy 2.0 or newer is required because the pinned simulator uses `numpy.trapezoid`. Default design: 300 replicates per truth/noise cell (10,800 raw rows; 36 summary cells). The truths at H0 = 50 and 90 lie exactly on the prior boundaries [50, 90], so boundary behavior must be interpreted separately from interior-truth calibration. This is a synthetic one-host diagnostic, not a real gravitational-wave measurement. The archived workflow completed with independent raw-record and paired-comparison audits; full metrics and prior-boundary caveats are recorded in the [artifact-backed run report](docs/known_host_calibration_run_20261025.md) and [calibration synthesis](docs/known_host_calibration_synthesis_2026-10-25.md). This remains a synthetic one-host diagnostic, not a real gravitational-wave measurement.

### Single-host fixed-truth selection control

The [protocol](docs/single_host_fixed_truth_control_protocol.md), [implementation](scripts/run_single_host_fixed_truth_control.py), [focused tests](tests/test_single_host_fixed_truth_control.py), and [artifact-backed results](docs/single_host_fixed_truth_control_results_2026-10-24.md) are available. The workflow completed successfully; its 2,700 paired observations show severe edge undercoverage even in a one-host toy model and truth-dependent effects of selection correction. See the linked report for scope limits and uncertainty. The result does not isolate a unique cause or establish realistic detector behavior.

### Numerical grid convergence audit

The [numerical-grid diagnostic audit](docs/numerical_grid_diagnostic_audit_2026-10-24.md) records the successful cross-version numerical-diagnostics workflow and its artifact inspection. It distinguishes computational-grid stability from changing the prior itself, and notes that tail quantiles can converge more slowly than posterior medians in a targeted edge case.

The [single-host fixed-truth control](scripts/run_single_host_fixed_truth_control.py) and [focused tests](tests/test_single_host_fixed_truth_control.py) implement the previously prespecified protocol in [the protocol document](docs/single_host_fixed_truth_control_protocol.md). The manually dispatchable [GitHub Actions workflow](.github/workflows/single-host-fixed-truth-control.yml) runs 300 observations at each of nine fixed truths, compares paired corrected and uncorrected inference, and emits raw rows, coverage summaries with Wilson intervals, paired absolute-error bootstrap intervals, and a runtime manifest. Results are not claimed until that workflow completes and its artifact is checked.


The single-host fixed-truth control has completed successfully: [results and paired uncertainty analysis](docs/single_host_fixed_truth_control_results_2026-10-24.md), [workflow run](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40099867069), and artifact digest recorded in the report. Results reinforce that selection correction is truth-dependent and that severe edge undercoverage persists in the one-host control; this does not establish a unique cause or realistic detector behavior.


### Matched host-structure control

The [matched host-structure protocol](docs/matched_host_structure_control_protocol.md) specifies a paired experiment in which detected distances are generated from a two-host mixture and the same observations are analysed using the two-host model and a collapsed weighted-mean-redshift model, with and without the toy selection normalization. The implementation, focused tests, and [GitHub Actions workflow](.github/workflows/matched-host-structure-control.yml) produce per-run observations, coverage summaries, and paired bootstrap uncertainty. Results should be cited only after the workflow completes and its uploaded artifact is inspected. This is a synthetic model-compression diagnostic, not an observational H0 measurement.


### Host-mixture × detection-threshold sensitivity factorial (implementation added)

The prespecified [protocol](docs/host_mixture_threshold_sensitivity_protocol.md) now has an implementation in [`scripts/run_host_mixture_threshold_sensitivity.py`](scripts/run_host_mixture_threshold_sensitivity.py), focused tests in [`tests/test_host_mixture_threshold_sensitivity.py`](tests/test_host_mixture_threshold_sensitivity.py), and a [GitHub Actions workflow](.github/workflows/host-mixture-threshold-sensitivity.yml). The design crosses four fixed truths, three host-redshift pairs, three intrinsic-weight pairs, and three upper detection thresholds (108 cells; 300 detected observations per cell). The [completed results report](docs/host_mixture_threshold_sensitivity_results_2026-10-24.md) records the independent audit, equal-cell summaries, paired contrasts, grid-sensitivity check, and limitations. The [factorial workflow](https://github.com/vedantkushwaha033-blip/dark-siren-cosmology-robustness/actions/runs/40101614384) and [Python 3.10/3.11/3.12 test suite](https://github.com/vedantkushwaha033-blip/dark-siren-cosmology-robustness/actions/runs/40101614170) both succeeded.


## Cosmological distance-model robustness

The first paired distance-model experiment is complete and its artifact has been independently audited. It generates synthetic one-host distances under flat \u039bCDM with fixed \u03a9m = 0.3 and compares matched flat-\u039bCDM inference against the deliberately mismatched low-redshift approximation. The design uses three fixed H0 truths, six redshifts, 200 paired replicates per truth/redshift cell, a fixed [40, 100] prior, and a documented Gaussian distance-noise rule. The experiment produced 7,200 posterior rows, 36 summary rows, and 18 paired contrasts; all summary metrics and all 5,000-draw bootstrap intervals were independently recomputed from the raw artifact and matched to floating-point precision.

See the [prespecified protocol](docs/cosmological_model_robustness_protocol.md), [implementation](scripts/run_cosmological_model_robustness.py), [cosmology-aware inference primitives](src/dark_siren_robustness/cosmological_inference.py), and [artifact-backed results and independent audit](docs/cosmological_model_robustness_results_2026-10-24.md). The successful [Cosmological model robustness workflow](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40115709336) uploaded the artifact with digest `sha256:0cb8c5fbfe234ce75d6752460e7708b73f8ae8708d5a29baccc5aa66299b8fde`.

In this simulator, low-redshift inference increasingly underestimates H0 as redshift increases, and its interval coverage degrades severely at moderate/high redshift; matched flat-\u039bCDM inference remains substantially better calibrated across the tested range. These findings are specific to a simplified Gaussian-distance diagnostic, not a real gravitational-wave measurement or a complete hierarchical dark-siren analysis. The experiment does not include realistic detector selection, catalogue incompleteness, host-population uncertainty, or marginalization over \u03a9m. The full Python 3.10/3.11/3.12 matrix on the pre-cleanup PR #7 head `b25cea2d13c8aa75c20f309c7b7f5e3b505507cf` passed in [run #40822714311](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40822714311). That historical success does not certify subsequent commits; the new documentation commit must pass its own exact-head checks.


### Joint distance-model × host-mixture experiment

The paired 27-cell experiment is complete: 8,100 shared synthetic observations were analyzed by four inference arms (32,400 posterior rows). It reports paired contrasts for distance-model mismatch and host-mixture collapse, plus a difference-in-differences interaction contrast. The raw artifact was independently audited; see the [completed results and audit report](docs/joint_distance_model_host_mixture_results_2026-10-24.md), [protocol](docs/joint_distance_model_host_mixture_protocol.md), [runner](scripts/run_joint_distance_model_host_mixture.py), [tests](tests/test_joint_distance_model_host_mixture.py), and [workflow](https://github.com/vedantkushwaha033-blip/dark-siren-cosmology-robustness/actions/runs/40116696053). The findings remain specific to this simplified synthetic model.


## Completed factorial sensitivity experiment

The nine-scenario factorial study and all three oracle controls completed successfully on 2026-10-24. The final artifact contains 4,500 main factorial records, 3,000 missing-host control records, 1,000 full-population oracle records, and 500 known-host oracle records plus paired bootstrap contrasts. Independent post-download checks reproduced core summary statistics to floating-point precision and verified exact event/posterior equality across the paired full-population masks. The main result is substantial overcoverage (nominal 68% intervals covered 93.2–95.6%; nominal 95% intervals covered 99.2–100%). In a paired 500-replicate comparison, known-host minus full-population 68% coverage was −19.4 percentage points (paired bootstrap 95% interval −23.4 to −15.4 points); 95% coverage was −1.4 points (−2.6 to −0.4). Mean 68% and 95% widths narrowed by 5.38 and 6.04 km s⁻¹ Mpc⁻¹, while mean absolute posterior-median error increased by 0.80 (0.45 to 1.16). This indicates candidate-host mixture ambiguity is a major source of overcoverage, while residual overcoverage remains even with the true host supplied. Treat these outputs as a calibration warning, not evidence of calibrated precision. See the [audited results report](docs/factorial_sensitivity_results_2026-10-24.md) and [successful workflow](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40142824335).

The oracle runner's `median_bias` definition has since been standardized to match the main factorial's median-of-errors definition, with regression tests added. The historical artifact remains unchanged; rerun the workflow before using corrected oracle median-bias values. This is a synthetic model study, not an observational H0 measurement.

## Latest validation milestone

The full test suite passed on Python 3.10, 3.11, and 3.12. The joint distance-model × host-mixture experiment passed its independent output audit, and the posterior-grid resolution audit completed successfully. Read the [validation milestone report](docs/validation_milestone_2026-10-24.md) for run links, counts, numerical findings, and limitations. These results are synthetic diagnostics, not an observational H0 measurement.


## Selected-likelihood calibration grid sensitivity

The selected-likelihood PIT calibration runner accepts `--grid-points` so the same seeded synthetic observations can be evaluated on different posterior-grid resolutions. The new [grid-sensitivity protocol](docs/calibration_grid_sensitivity_protocol.md) compares 1,201, 2,401, and 4,801 grid points using 300 replicates per grid in CI. Each resolution's raw records and summary are independently audited before a comparison CSV is written.

To reproduce the screening experiment locally:

```bash
python scripts/run_calibration_grid_sensitivity.py \
  --replicates 300 \
  --seed 20261025 \
  --grid-points 1201,2401,4801 \
  --output-dir artifacts/calibration_grid_sensitivity
```

The output includes `grid_sensitivity.csv` and a separate `pit_replicates.csv`, `pit_summary.json`, and `audit_report.json` for each grid. The screening run is intended to detect numerical sensitivity, not establish formal equivalence or exact calibration. Any publication-facing conclusion should use more replicates, additional seeds, prespecified tolerances, and realistic selection-model validation.

### Paired host-misspecification experiment

The [paired host-misspecification protocol](docs/paired_host_misspecification_protocol.md), [runner](scripts/run_paired_host_misspecification.py), [independent auditor](scripts/audit_paired_host_misspecification.py), and [adversarial regression tests](tests/test_paired_host_misspecification.py) add a stronger comparison to the earlier unpaired host stress test. Each replicate now generates one synthetic observation and evaluates all four inference assumptions on that identical observation. The experiment reports paired contrasts for reversed versus nominal two-host weights and omitted versus modelled hidden host, with stratified-by-replicate percentile bootstrap intervals. The nominal “correct weights” arm still omits a component of the generating population, so it is not an oracle model. This remains a single-truth, single-noise, single-prior toy diagnostic, not real-data validation.

The full 1,000-observation workflow completed successfully: [workflow and audit log](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40280238225), [raw artifact](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40280238225/artifacts/12631703296), and [audited results and interpretation](docs/paired_host_misspecification_results_2026-10-26.md). The independent audit reconstructed all 4,000 posterior fits and all 12 paired bootstrap contrast rows. The uploaded artifact SHA-256 is `abc6f451311b921a0370e94356848a65734ecfb554f8443b43b740bf2b9f0d47`.

### Paired host-misspecification sensitivity surface: audited results

The [prespecified 108-cell protocol](docs/host_misspecification_sensitivity_surface_protocol.md), [factorial runner](scripts/run_host_misspecification_sensitivity_surface.py), [independent auditor](scripts/audit_host_misspecification_sensitivity_surface.py), [focused tests](tests/test_host_misspecification_sensitivity_surface.py), and [full workflow](.github/workflows/host-misspecification-sensitivity-surface.yml) vary three injected H0 values, three noise scales, four hidden-host fractions, and three hidden-host redshifts. The 21,600-observation experiment (86,400 posterior fits) completed successfully and passed independent reconstruction: [workflow and audit](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40280537559), [artifact](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40280537559/artifacts/12632197386), and [quantitative results and limitations](docs/host_misspecification_sensitivity_surface_results_2026-10-26.md). The zero-hidden-fraction negative control passed exactly. Across the tested design, omitting the hidden component increased RMSE and reduced 68%/95% coverage as its true population fraction grew. These are conditional synthetic diagnostics, not real-data validation.

The [equal-cell marginal summary script](scripts/summarize_host_misspecification_sensitivity_surface.py) verifies the raw-to-summary reconstruction and produces equal-cell arm summaries plus stratified paired-bootstrap intervals across hidden-fraction strata.



## Research-readiness audits (2026-10-28)

Two new source-level / artifact-level analyses document limits that matter before any publication-facing claim:

- [Mixture, selection-normalization, and boundary source audit](docs/mixture_selection_boundary_source_audit_2026-10-28.md): reviews the current synthetic mixture weighting, log-space arithmetic, toy detection efficiency beta(H0), and prior-boundary handling. It found no source-level evidence of the historical missing-mass double-counting bug in the inspected functions, but does not validate a real detector-selection model.
- [Fixed-truth coverage reference-interval audit](docs/fixed_truth_coverage_reference_interval_audit_2026-10-28.md): compares Bayesian equal-tail credible-interval coverage with an inverted Gaussian distance confidence interval in the archived Arm A synthetic artifact. This distinguishes two coverage estimands; it does not establish universal posterior calibration or astrophysical validity.

These audits are included in PR #7, which is open and ready for independent review. The pre-cleanup PR #7 head passed the three-version matrix in [run #40822714311](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40822714311); the current documentation cleanup needs its own exact-head CI result. These are source/artifact audits, not proof of real detector validity. The toy low-redshift relation, finite distance-window selection rule, prior sensitivity, and grid-resolution stress cases remain explicit limitations.


### Dataset-level likelihood scope audit

The [dataset-level likelihood grid and scope audit](docs/dataset_likelihood_grid_and_scope_audit_2026-10-28.md) reviews the new independent-event composition prototype. It highlights that the code does not require the numerical grid to include the declared prior boundaries, and that shared completeness/nuisance parameters and selection conditioning remain outside its scope. These are source-review findings and recommendations, not evidence of a runtime failure. Resolve the prior-grid convention and add convergence tests before using this prototype for confirmatory calibration.


### Current literature and novelty boundary (30 October 2026)

A focused literature audit now reviews recent work on GP reconstruction of incomplete line-of-sight catalogues, spatial/magnitude-aware catalogue completion, empirical completeness testing, joint host-population inference, and large-scale-structure methods. See [the literature and novelty audit](docs/literature_and_novelty_audit_2026-10-30.md).

The audit concludes that catalogue incompleteness, host weighting, and flexible reconstruction are not novel topics. The defensible scope of this repository is a controlled statistical benchmark of external parent-frame validation, estimated inverse-inclusion weighting, finite validation-sample sensitivity, positivity diagnostics, and selection-model misspecification. It does not establish a new astrophysical estimator or superiority to GP/LIGHT, variance-completion, CHIMERA, or gwcosmo methods.

### Pilot manuscript note

- **Linear finite-population bootstrap pilot manuscript:** [reproducible methods, audited pilot results, interpretation, and limitations](manuscript/linear_finite_population_bootstrap_pilot_note.md). This is an exploratory simulation note, not a confirmatory or submission-ready paper. Its audit qualification and v2 requirements are stated explicitly.
