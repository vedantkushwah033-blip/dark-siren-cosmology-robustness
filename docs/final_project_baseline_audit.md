# Final project baseline audit — 2026-11-04

## Audit scope and evidence standard

This audit records the live GitHub state inspected on 2026-11-04 through the repository contents, Git refs, recursive tree, pull-request metadata, and GitHub Actions run/job endpoints. It distinguishes repository-visible evidence from progress reported in conversation but not yet present in the inspected commit. A file's presence is not evidence that its experiment ran; a protocol is not a result; a passing unit test is not scientific validation.

Repository: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness

## 1. Live source-control baseline

- Audited working branch: `research/bootstrap-failure-path-repairs-2026-11-02`.
- Exact branch HEAD at inspection: `80efba0614a9cdbb2a91e218aeff22ae53305008`.
- Latest commit message: `test: accept early status-ledger tamper rejection`.
- The commit is on the repair branch. No local working-tree state is visible through the GitHub API; therefore cleanliness of a local checkout is **not verified**.
- PR #15 is open and unmerged. Its head remains `research/estimated-inclusion-audit` at `4f700fcf120ce547de598e6199d3829bc88c8498` in the fetched PR metadata. PR #15 was not modified or merged.
- The baseline audit does not assume that the repair branch, PR #15, PR #16, PR #17, or default branch are interchangeable. Each is a separate ref and must be compared before integrating work.

## 2. Exact-head CI status

Workflow run `41749192480` is attached to exact HEAD `80efba0614a9cdbb2a91e218aeff22ae53305008`. At the time of inspection it was **in progress**, with three jobs:
- numerical-diagnostics (Python 3.10): in progress;
- numerical-diagnostics (Python 3.11): in progress;
- numerical-diagnostics (Python 3.12): in progress.

No completed green exact-head CI result was available at this audit timestamp. Earlier workflow failures and fixes are recorded in commit history, but cannot substitute for the current exact-head run. The current status must be refreshed before release.

## 3. Repository-visible numerical and uncertainty work

The inspected branch tree contains existing numerical and uncertainty-related components, including:
- `src/dark_siren_robustness/finite_parent_uncertainty.py`;
- `src/dark_siren_robustness/frame_precision_planning.py`;
- `src/dark_siren_robustness/uncertainty_likelihoods.py`;
- `src/dark_siren_robustness/uncertainty_sensitivity.py`;
- `scripts/run_targeted_grid_convergence.py`;
- `scripts/check_distance_uncertainty_grid_stability.py`;
- `scripts/run_selection_prior_sensitivity.py`;
- `scripts/audit_distance_uncertainty_sensitivity.py`;
- related tests and prior numerical/selection sensitivity reports.

However, the inspected branch tree did **not** contain the specifically requested files:
- `scripts/audit_h0_grid_convergence.py`;
- `scripts/audit_redshift_grid_convergence.py`;
- `scripts/audit_h0_prior_sensitivity.py`;
- `manuscript/dark_siren_robustness_research_paper.md`;
- `docs/final_project_baseline_audit.md` (this file is being added now).

This does not prove the files do not exist in an unpushed local checkout or another branch. It means they were not present at the audited branch HEAD. The claimed local run outcomes (H0-grid PASS, redshift-grid PASS, prior-domain FAIL; 12-replicate H0 grid; byte-identical deterministic replay; 1,105 passed/2 xfailed) are treated as **user-reported, not independently verified**, until their exact source commit, scripts, immutable artifacts, manifests, and audit outputs are available in a repository ref or supplied as files.

## 4. Existing scientific documentation and model scope

The repository already contains substantial source-level and design documentation, including:
- `docs/independent_likelihood_equation_audit_2026-10-29.md`;
- `docs/selection_mechanism_interpretation_audit.md`;
- `docs/model_consistency_audit.md`;
- `docs/estimated_inclusion_nonlinear_stress_protocol_2026-10-30.md`;
- `docs/linear_to_nonlinear_validation_bridge_2026-10-31.md`;
- `docs/nested_bootstrap_failure_accounting_repair_spec_2026-11-04.md`;
- `docs/bootstrap_failure_path_repair_implementation_record_2026-11-02.md`;
- literature positioning, novelty, numerical validation, calibration, prior sensitivity, and experiment protocol documents.

The likelihood audit describes a discrete-host mixture, optional positive-truncated distance noise, a selection efficiency term, and a fixed-detected-count likelihood of the form
`sum_i log L_i(d_i | H) - N log beta(H)`. It explicitly limits this to the stipulated synthetic fixed-count model; it is not a realistic detector-level likelihood or event-rate model.

The repository's selection-mechanism audit describes simple synthetic retention mechanisms, including exact-size random thinning and inverse-redshift-weighted thinning. These are not calibrated survey selection functions.

The estimated-inclusion protocol distinguishes catalogue-only, oracle diagnostic, and validation-estimated weighting. It explicitly states that an independent validation sample under the same latent population/inclusion law is not a transportability test. The protocol is not itself a confirmatory result.

The README contains older headline results and references historical CI runs. These claims must be re-linked to their own immutable artifacts and source SHAs before they are reused in the new manuscript. A historical README claim is not evidence that the current repair HEAD regenerated the result.

## 5. Bootstrap accounting repair status

The branch includes a dated specification and implementation record for separating outer procedure outcomes from inner resample outcomes. The intended statuses are success, non-evaluable, failed, and not-run; wrapper-level failure must not be assigned a fabricated number of inner attempts; and the FPC-scaled arm is a derived diagnostic rather than a separately sampled procedure.

At HEAD `80efba0...`, the most recent test changes address status-ledger reconciliation and replay-auditor exception typing. Exact-head CI was still running. The smoke artifact and independent replay gate are therefore **not passed on the evidence inspected here**.

## 6. Scientific risks that remain release-critical

1. **Catalogue inclusion vs GW detection:** these are distinct mechanisms and must have separate probabilities and normalizations. A correction for galaxy-catalogue retention does not automatically correct gravitational-wave detector selection.
2. **Toy selection direction:** one historical distance-threshold toy accepts events above a distance threshold. This is a lower-truncated observed-distance construction, not a realistic statement of GW detector sensitivity. Keep it explicitly labelled and do not silently reinterpret it.
3. **Positivity:** inverse-probability weighting requires support/positivity. Report probability minima/quantiles, near-zero fractions, maximum normalized weight, top-weight mass, ESS, and failure rates. Truncating weights changes the estimand and must be predeclared.
4. **Estimated propensity uncertainty:** plug-in logistic estimates condition on fitted coefficients. Validation-sample uncertainty must be propagated through weights, host-density estimation, and posterior summaries using a design justified for the actual validation sampling mechanism.
5. **Nonlinear bootstrap validity:** a finite-population correction or a successful linear-mean bootstrap benchmark does not validate a bootstrap for logistic fit → inverse weighting → KDE/host density → nonlinear posterior functionals.
6. **Posterior vs repeated-sampling uncertainty:** conditional posterior spread, Monte Carlo error, finite-parent variation, validation-sample variation, and repeated-sampling variability of posterior summaries are different quantities.
7. **Selection normalization:** the beta(H) term must match the exact event-generation/conditioning mechanism and population measure. Fixed-count conditioning discards rate/count information.
8. **Numerical convergence vs validity:** grid refinement and prior-domain sensitivity test numerical approximation/support dependence; neither validates the generative model.
9. **Monte Carlo precision:** 2,000 replicates are not automatically sufficient. Set targets for the paired estimand, coverage uncertainty, failure rate, and subgroup/cell comparisons; report achieved precision.
10. **Independent audit independence:** an auditor that calls the same production summary implementation can detect artifact tampering but is not a fully independent reconstruction of the statistic. Keep independent equations and code paths for likelihood, posterior summaries, paired contrasts, and provenance.

## 7. Manuscript and artifact baseline

At the audited HEAD, `manuscript/` contains existing drafts (`finite_validation_bootstrap_draft.md` and `linear_finite_population_bootstrap_pilot_note.md`), but not the requested integrated journal manuscript. There is no evidence in the inspected tree that all new numerical audit outputs, uncertainty propagation outputs, and q=1.2 pilot failures have been persisted and independently audited together.

The manuscript must therefore use explicit placeholders for unavailable results. Do not populate Results from conversational interim numbers alone. In particular, the reported two failed q=1.2 pilot rows out of 840 require the raw records, optimiser diagnostics, and independent replay before they can be labelled numerical optimiser failures or genuine separation.

## 8. Status matrix

| Workstream | Repository-visible implementation | Tests/CI | Executed evidence | Independent verification | Gate |
|---|---|---|---|---|---|
| Existing synthetic likelihood and selection model | Yes, with extensive prior audits | Historical tests/docs; exact current CI pending | Historical artifacts/docs exist | Some prior equation/artifact audits; scope-specific | Partial |
| Estimated inclusion pilot | Yes | Existing tests and historical CI records; current exact-head status pending | Historical pilot artifacts/docs; current new pilot claims not yet located | Existing pilot auditor; new q=1.2 claim unverified | Partial |
| Nested bootstrap failure accounting | Yes, repair branch | Exact-head matrix in progress | No new smoke artifact verified at this HEAD | Replay auditor implemented; current smoke gate not passed | Not passed |
| Requested three numerical audit scripts | Not present under those names at audited HEAD; related tools exist | Related tests exist | New user-reported runs not tied to visible artifacts | Byte-identical replay claim not verified here | Not passed |
| New uncertainty-propagation module | Related modules exist; user-reported new module/tests not identified at HEAD | Not verified for claimed 7 tests | New propagation run not visible | Not verified | Not passed |
| Final integrated manuscript | No | N/A | N/A | N/A | Not passed |
| Final referee report | No integrated final report identified | N/A | N/A | N/A | Not passed |

## 9. Next actions in dependency order

1. Refresh run `41749192480` and obtain all three job conclusions and logs.
2. Identify the exact branch/commit containing the reported local numerical scripts and uncertainty module. If they are only local, obtain the files and their source artifacts before claiming completion; do not recreate them blindly.
3. Persist the numerical audit configuration, matched-input fingerprints, posterior reconstruction outputs, tolerance definitions, and per-replicate status ledgers.
4. Resolve and independently audit the q=1.2 failed rows before deciding whether the cause is optimizer instability, separation, or another model/data issue.
5. Finish the uncertainty study on immutable pilot inputs and independently recompute its summaries.
6. Freeze a confirmatory design only after planning-pilot precision, failure rates, and all audit gates are satisfactory.
7. Write the integrated manuscript and referee report from verified artifacts only.

## 10. Overall assessment

This is a substantial synthetic-methodology codebase with extensive prior experiments and audits, but the complete research mission is **not complete**. The current branch is still undergoing software failure-accounting verification. Numerical/prior sensitivity and uncertainty-propagation claims supplied in the latest progress message are not yet tied to repository-visible artifacts at this HEAD. The defensible publication path is a carefully bounded synthetic simulation study; real LVK inference, observational H0 measurement, real-survey transportability, and realistic full detector-selection modelling remain outside scope.
