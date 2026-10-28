# Research-readiness audit — 2026-10-28

## Purpose and scope

This is a source-level audit of the repository branch `research/distance-uncertainty-likelihoods` as inspected on 2026-10-28. It is deliberately conservative. It is not a substitute for executing the full test suite, downloading and auditing current CI artifacts, or independent review by a statistician or gravitational-wave cosmologist.

## Executive assessment

The repository contains meaningful research-engineering work: a documented synthetic inference model, tests for mixture behavior and numerical edge cases, prespecified simulation protocols, and independent artifact-audit scripts. The strongest current contribution is a carefully bounded simulation/software-validation study.

The current evidence does **not** yet justify calling the project a validated method for real dark-siren cosmology, a new H0 measurement, or a demonstrated novel diagnostic. External CHIMERA/reference-pipeline validation remains a distinct open gate in the project documentation.

## What was inspected

- `README.md`
- `pyproject.toml`
- `src/dark_siren_robustness/hardened_inference.py`
- `tests/test_hardened_inference.py`
- `docs/prospective_diagnostic_comparison_protocol_2026-10-27.md`
- `docs/distance_uncertainty_model_sensitivity_interpretation.md`
- `docs/missing_host_mixture_likelihood_impact_audit_2026-10-27.md`
- Open pull request #3 and its description
- Pull-request review-thread and review-submission queries for PR #3 (no review threads or submitted reviews were returned by the connector query)

Repository links:
- [Repository](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness)
- [PR #3](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/pull/3)
- [Inference implementation](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/research/distance-uncertainty-likelihoods/src/dark_siren_robustness/hardened_inference.py)
- [Prospective diagnostic protocol](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/blob/research/distance-uncertainty-likelihoods/docs/prospective_diagnostic_comparison_protocol_2026-10-27.md)

## Evidence and status

| Area | Source-level finding | Status | Next requirement |
|---|---|---|---|
| Missing-host mixture | Implementation uses log-space mixture assembly; tests compare it with an explicitly concatenated population for several nonzero missing fractions and test zero/tiny-mass limits. | Source/test coverage inspected; latest test execution not performed in this audit. | Run current tests and trace all direct and indirect callers. |
| Posterior normalization | Posterior normalization is restricted to points inside the declared prior interval, avoiding artificial trapezoidal mass across a discontinuous prior boundary. | Implementation inspected. | Add/retain tests for irregular grids, extreme likelihood ratios, and posterior summaries. |
| Posterior summaries | An independent reference-review handoff reports that the suspected `density > 0` tail-underflow concern was tested and found not to cause material summary error in the tested regime. The more consequential demonstrated limitation was grid resolution: when the posterior width is only a few grid steps, interpolated quantiles and SD can be inaccurate. Reported examples on a 0.025 grid include 68% interval error around 0.006 at sigma=0.02 and SD 0.0071 instead of 0.01 at sigma=0.01. These are handoff-reported reference results; the package/script is not included in this repository and production behavior was not independently rerun in this audit. | Revise the earlier underflow concern downward; grid-resolution limitation is plausible and quantified in the separate reference report, not yet verified against the production code. | Preserve regression tests in the statistical implementation audit; decide whether to add an opt-in resolution warning only after defining a defensible resolution criterion and checking backwards compatibility. |
| Distance-noise models | Documentation correctly labels fixed-realized-sigma arm A as oracle-like and notes that arm C changes the generator away from reference H0=70. | Interpretation is appropriately qualified in inspected documentation. | Ensure every result table and figure preserves these distinctions. |
| Research question | The prospective protocol asks whether diagnostics distinguish host-population and distance-noise misspecification. | Proposed, not established as novel. | Complete literature comparison and ensure the chosen experiment actually varies both misspecification families jointly. |
| Current distance-uncertainty experiment | PR #3 describes a synthetic one-known-host low-redshift experiment. | Useful narrow diagnostic; not a test of full catalogue incompleteness and hidden-host interactions. | Do not use it alone to claim that diagnostics distinguish catalogue and noise misspecification. |
| Corrected hidden-host artifacts | The impact-audit document reports reruns and independent reconstruction for the named experiments on a corrected revision. | Historical evidence is documented, not independently re-downloaded or rerun during this audit. | Verify artifact identifiers/checksums and current workflows before manuscript freeze. |
| CI | GitHub Actions records were queried for the exact observed heads. PR #3 head `9fb48b4acb62d91ee9e8830635fc1785cbe75b5f` has successful tests and named experiment runs. PR #4 head `9300d11d74d758f6a5cf42987a787c7ac533a843` has a completed successful `tests` run (run 40615127431, run number 1591). | These recorded runs passed for those exact SHAs. They do not validate unpushed local commits or prove external validity, novelty, or publication readiness. | Re-query Actions after any new push/rebase; inspect relevant job logs/artifacts and verify their SHA matches the code being discussed. |
| External benchmark | Project documents state that exact CHIMERA v2.2 input fixtures/reference validation remain outstanding. | Not externally validated on the evidence inspected. | Locate exact inputs and execute a documented benchmark, or explicitly scope the paper as synthetic. |
| Literature novelty | Protocol lists relevant dark-siren completeness and calibration work and explicitly says the bibliography is not exhaustive. | Not established. | Read methods/results of closest papers and map overlap claim by claim. |
| Publication readiness | The repository has protocols, tests, and analysis artifacts, but external validity and novelty remain unresolved. | Not manuscript-ready as a claim of a new real-world method. A narrowly scoped synthetic methods/software note may become viable after validation. | Freeze primary estimand, verify results, complete literature audit, obtain domain-expert review. |

## Statistical and scientific cautions

### 1. Separate calibration from robustness

Simulation-based calibration (SBC) is most directly interpretable when parameters are drawn from the prior and the data are generated from the same model used for inference. It is not, by itself, proof that inference is robust under misspecification. Repeated-sampling coverage at fixed parameter values answers a different question and should be reported separately.

### 2. Avoid mixing multiple mechanisms in a contrast

A paired comparison can isolate an inference-model change only when both arms analyze the same simulated observations and the relevant generator is held fixed. If an arm also changes the data-generating noise scale, its contrast combines generator and inference changes. The inspected interpretation document already recognizes this for arm C; maintain that distinction throughout the project.

### 3. Treat oracle arms as diagnostic controls

An arm supplied with the true hidden-host fraction, redshift distribution, or true-H0-derived uncertainty is a useful upper/reference benchmark. It is not an operational estimator unless that information is actually available in a real analysis. Label oracle arms prominently in plots and tables.

### 4. Coverage requires context

Report coverage alongside interval width, bias, RMSE, and Monte Carlo uncertainty. Prior-boundary truths require separate treatment; average coverage can hide truth-dependent failures. Do not rank methods on tiny differences without uncertainty estimates and a prespecified primary contrast.

### 5. The model is not a full dark-siren pipeline

The inspected documentation states that the simplified inference uses the low-redshift approximation d_L approximately equal to c z/H0 and omits realistic gravitational-wave likelihoods, detector-selection normalization, and several population and catalogue uncertainties. Claims must remain conditional on this simplified setup until a defensible external benchmark is completed.

## Prioritized next actions

1. **Reconcile current evidence:** identify the exact current PR head, query current Actions runs, and verify the status of the latest full test matrix and each primary experiment.
2. **Close code correctness gates:** execute the test suite; inspect the full call graph of the corrected mixture API; add only targeted numerical tests that reveal a meaningful edge case.
3. **Freeze the primary study:** choose one primary estimand and a small number of prespecified contrasts. Ensure the design includes separate host-only, noise-only, and joint misspecification conditions if the central question is diagnostic distinguishability.
4. **Complete novelty assessment:** compare the exact proposed contribution against the closest dark-siren catalogue-completeness, Bayesian calibration, posterior predictive checking, and SBC literature.
5. **Resolve external-validation scope:** either reproduce the exact pinned CHIMERA benchmark with verified inputs or make the manuscript explicitly a synthetic-model study without external-validity claims.
6. **Obtain expert feedback:** seek a domain researcher to critique the model assumptions and the interpretation before describing the work as a publishable astrophysical methodology.

## Release gates

Do not claim a validated real-data method or new H0 measurement until the study has a realistic reference-pipeline comparison, verified data provenance, independent validation, and domain-expert scrutiny. Do not claim novelty until the literature review establishes a meaningful distinction from prior work. Do not claim the latest CI passes unless the current run has actually been inspected.

## Audit limitations

This audit inspected source and documentation through GitHub connector results and subsequently queried Actions records for the exact PR heads listed above. It did not execute production Python tests locally, download/reconstruct workflow artifacts, run simulations, or complete an exhaustive literature review. A separate reference-review handoff reports 15 successful tests against an independently written reference implementation (Python 3.12.8, NumPy 1.26.4, SciPy 1.13.1); those tests are not production-repository tests, and the referenced ZIP/script was not present among the attachments available for this update. Successful CI applies only to the exact recorded commit SHAs and workflows.


## Patch-bundle reconciliation update (2026-10-28)

The supplied `pr3_all_9_commits.mbox` was parsed directly: it contains **nine** format-patch messages touching **20 distinct repository paths**. The live PR #3 changed-file listing contains 39 paths. Only `README.md` and `docs/literature_and_novelty_audit_2026-10-26.md` overlap those path lists; the other 18 paths are not part of PR #3's current diff. Spot-checks of new report paths such as `docs/joint_selection_simulation_audit_2026-10-26.md` and `docs/diagnostic_power_results_2026-10-27.md` returned 404 at the observed PR #3 head. This is strong evidence that the nine-patch bundle is **not integrated as a whole**. It does not prove that none of its overlapping README/literature content is already present, and it does not reconcile the separate handoff's 13-commit series because that series is not among the available attachments. Do not apply the nine patches blindly or describe their results as validated by PR #3 CI.

The latest observed PR #4 test run for SHA `9300d11d74d758f6a5cf42987a787c7ac533a843` completed successfully (GitHub Actions run `40615127431`). This is evidence for that exact documentation-only revision; a later commit must be checked separately.


## Targeted grid-resolution stress test (2026-10-28)

A separate report, [posterior-grid resolution stress test](posterior_grid_resolution_stress_test_2026-10-28.md), documents a coarse-versus-fine comparison. In a deliberately narrow single-host likelihood, the coarse-grid posterior SD was 0.00091 versus 0.00654 on the fine grid at distance uncertainty 0.02 Mpc, and nearly zero versus 0.00327 at 0.01 Mpc. A direct-package audit (`scripts/audit_extreme_grid_resolution.py`), regression tests, and a CI step have now been added on this branch; the newest workflow must pass before calling them validated. This establishes a plausible failure mode under deliberately extreme precision, not that existing headline runs are wrong. No production inference behavior was changed.


## Exact-SHA artifact review update (2026-10-28)

The Actions run for commit `cbc1316b00c396e1795c2cc900a736e13e067d14` (run `40618560947`) subsequently completed successfully across Python 3.10, 3.11, and 3.12. I downloaded and inspected the Python 3.12 artifact. Its full distance-uncertainty experiment reports 32,400 raw records, zero failed runs, a successful independent formula reconstruction, and a PASS grid-stability check against a 3,201-point grid. The maximum absolute change in the 68% and 95% interval endpoints was below 0.00018 H0 units for this experiment; this supports numerical stability for that particular design and does not override the deliberately extreme stress test above.

A coverage stratification check of the same raw records shows why the pooled 68% coverage must not be interpreted as a calibration result: truths exactly at the prior boundaries H0=50 and H0=90 have zero central-68% interval coverage in this design. The pooled coverage is approximately 57.0–57.3% across the three arms, while restricting to the interior truth values 55–85 yields approximately 73.3–73.6%. These are descriptive re-aggregations of the artifact, not new simulations. The interior values are not nominal 68% either, so the result requires explanation and should not be presented as proof of calibration. Boundary and interior truth strata, interval width, bias/RMSE, and Monte Carlo uncertainty should be reported separately.

The artifact did **not** contain `artifacts/extreme_grid_resolution_audit.json`, even though the workflow ran the script: the upload-artifact path list omitted that JSON. This workflow omission was corrected in commit `f2ceab3fbc0c9c34dab1c4a54afbf1b2021534a8`. A fresh CI run for that commit is pending; its artifact must be inspected before claiming the JSON was archived successfully. The extreme-grid test and script were also independently exercised against a local copy of the package implementation: four of five stress cases were classified as under-resolved, and the deliberately broad `sigma_d=2` Mpc case was within the stated tolerance. This local check is supportive but does not replace the pending exact-SHA repository CI run.
