# Research-readiness report — 29 October 2026

## Decision

**Current disposition: GO for a carefully scoped synthetic-methodology manuscript draft and independent review; NO-GO for claims about real GW catalogues, real H0, or external-method superiority.** PR #7 is open and ready for independent review, not merged; PR #8 remains a separate draft.

## Completed and artifact-verified

1. **KDE bandwidth convention:** replaced the ambiguous SciPy covariance-factor usage with a direct reflected weighted-kernel sum; explicit factor semantics now use the weighted SD of the original sample, with a direct kernel-sum regression test.
2. **Catalogue incompleteness benchmark:** a locked 500-replicate Stage B experiment completed with six paired arms. Its independent posterior-grid and raw-table audit passed. Primary contrast: 5.9519 km s⁻¹ Mpc⁻¹ reduction in absolute posterior-mean error, bootstrap 95% interval [5.5410, 6.3573].
3. **H0-dependent-noise model:** a separate 1,200-row synthetic calibration artifact passed independent generator replay and posterior reconstruction to floating-point precision. A follow-up audit on commit `9c662b6` also independently reaggregated every headline summary from the raw CSV and passed. Prior-predictive mean rank 0.514; ten-bin descriptive chi-square p=0.450.
4. **Transparent documentation:** added result reports, a literature-positioning audit with an assumption-by-assumption comparison table, and an integrated manuscript draft. The manuscript distinguishes three separate synthetic studies and states non-claims.
5. **Reproducibility plumbing:** raw replicate records, seed schemes, posterior-grid outputs, SHA-256 manifests, auditors, adversarial tampering tests and CI workflows are in the repository.

6. **Configurable selection-function replication:** three predeclared conditions completed 250 paired replicates each (1,500 rows per condition across six arms). All three independent auditors reconstructed every saved posterior grid and aggregate arm; all manifest hashes/byte counts passed. I independently rederived bias, RMSE, coverage counts, paired primary contrasts and the 10,000-resample bootstrap from the downloaded raw CSVs; summaries matched to floating-point precision. Paired absolute-error reductions were 7.436 [6.798, 8.053], 2.246 [1.711, 2.764], and 2.354 [1.779, 2.923] km s⁻¹ Mpc⁻¹ for the low-floor, high-floor and broad-transition conditions. This strengthens the within-model synthetic result, not a real-data claim. See [confirmatory report](catalogue_selection_confirmatory_replication_2026-10-29.md) and [workflow #40819870231](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40819870231).

## Still open

- The expanded host-density screen contains 13 configurations × 40 paired replicates. All 13 configurations passed local per-configuration independent replay/posterior reconstruction (1,560 fit rows), and their summary statistics were reaggregated independently. The combined clean campaign and one-shot full audit passed at [#40802415000](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40802415000). Artifact `12866753211` (ZIP digest `sha256:0ec51d1c96e9f86d7f1ef2a6dbcb9e57a0c67fa5be2459617fcf51ccbc5565fe`) independently reconstructed all 1,560 rows, reaggregated 403 summary fields with zero discrepancy, and verified 14 manifest-listed files. One interval-coverage flag disagreement was boundary-sensitive (endpoint 9.37×10⁻⁵ from truth) and is disclosed in the audit; it did not fail the audit. The screen remains exploratory with imprecise coverage.
- **Exact-head test matrix (updated 29 October 2026):** PR #7 head `9bf7508d6b0e78f504139571cb9c4c6d7cb118ff` passed the PR-triggered numerical-diagnostics workflow [#40874302252](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40874302252). All three jobs completed successfully on Python 3.10, 3.11, and 3.12; the Python 3.12 job also completed the full-run diagnostics and independent artifact checks, while Python 3.10/3.11 skipped the explicitly Python-3.12-only full-run steps. This establishes test success for this exact commit, not scientific correctness or observational validity. The push-triggered matrix also passed at [#40874295080](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40874295080).
- The CHIMERA v2.2 external benchmark is not complete. The exact expected input dataset(s), compatibility with the tagged release, and execution output have not been verified. Do not substitute files from an older release.
- A focused comparison table is now included in the manuscript and literature audit for Gaussian-process reconstruction, spatial-correlation completion, clustering-preserving reconstruction (LIGHT/Ghosh–More), variance completion, realistic clustering simulations, the practical `gwcosmo` completeness test, host-population inference and CHIMERA. A systematic review and expert assessment of assumptions remain outstanding.
- No expert review, arXiv submission, peer review, or publication claim is made.

## Required gates before release

1. **Exact-head CI: PASS for current PR #7 head** `9bf7508d6b0e78f504139571cb9c4c6d7cb118ff`: PR-triggered run [#40874302252](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40874302252) and push-triggered run [#40874295080](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40874295080) both completed successfully across Python 3.10, 3.11, and 3.12. Python 3.12 ran the full-run diagnostics and independent artifact checks; steps marked Python-3.12-only were skipped on 3.10/3.11. CI pass is not evidence of observational validity.
2. The 13-configuration host-density screen and full independent reconstruction are complete: 1,560/1,560 row fits reconstructed, 403 aggregate fields rederived with zero discrepancy, 14 manifest files checked. The canonical report records one borderline coverage flag and remains explicitly exploratory.
3. The focused comparison table is now included in the manuscript and literature-positioning audit. A systematic literature review and external expert assessment remain outstanding.
4. **CHIMERA v2.2 external benchmark: NO-GO pending data.** The pinned release and official examples were inspected, but the exact compatible SNR-25 mock, completeness products and `pcat_incomplete.h5` were not verified as available. Do not substitute v2.0 inputs; reopen this gate only when exact, provenance-traceable inputs are obtained.
5. Re-read the manuscript tables against artifacts, have a knowledgeable researcher review the model assumptions, and retain all toy-model limitations in any abstract or outreach.

## Claim language

Acceptable: “A reproducible controlled synthetic benchmark quantifies the effect of a known catalogue-inclusion mechanism and several selection-law assumptions on a simplified dark-siren inference model.”

Unacceptable: “Solved the Hubble tension,” “measured a new H0,” “first method for incomplete galaxy catalogues,” “validated LVK performance,” or “inverse weighting universally fixes incompleteness.”

## Current follow-up documentation

- The host-density campaign's authoritative full results are in [host-density sensitivity screening results](host_density_sensitivity_screening_results_2026-10-29.md); the supplemental [follow-up validation note](followup_validation_host_density_and_sbc_2026-10-29.md) summarizes the result and cross-checks the baseline numbers.
- PR #8 remains a separate draft for paired grid/bandwidth sensitivity validation. PR #7 is open and unmerged; its exact head `9bf7508d6b0e78f504139571cb9c4c6d7cb118ff` passed both CI runs [#40874302252](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40874302252) and [#40874295080](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40874295080). The PR is a large review surface (173 commits ahead and 82 changed files), has no recorded submitted reviews or inline review threads, and targets `experiment/catalogue-density-estimators` rather than `main`; obtain maintainer review of scope/base-branch choice before merging. Scientific review and the external-pipeline/data gate remain outstanding.

## Final project scope

A useful and defensible outcome is a reproducible scientific-computing / statistical-methods study with explicit controls and independent audits. A real astrophysical result requires external data, a realistic detector/selection likelihood, appropriate population treatment and independent subject-matter review. The current simulations do not cross that boundary.


## Completed screening campaign — evidence boundary (added 2026-10-29)

- Run: [host-density validation campaign #40797592794](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40797592794), source commit `5bbe5a434878af396e8c58e3ddfe7af0356e29ed`, artifact ID `12865051001`.
- Configuration: 13 settings × 40 paired replicates, fixed H0=70, 6,000 parent galaxies, 40 detected events, threshold 250 Mpc, distance noise 140 Mpc. This is screening, not a confirmatory precision run.
- Integrity: every one of the 14 manifest-listed campaign files matched its SHA-256 and byte count after download. The unified independent audit in run #40802415000 subsequently reconstructed all 1,560 fit rows, rederived all 403 aggregate fields with zero discrepancy, and verified all 14 manifest-listed files. One interval-coverage flag differed only at a boundary-sensitive endpoint 9.37×10⁻⁵ from truth; this is disclosed in the canonical report. The campaign remains an exploratory screen, not a confirmatory precision run.
- Baseline: catalogue-only bias −8.634, RMSE 9.165, 68% coverage 10.0%, density L1 0.537; known inclusion correction bias +1.180, RMSE 3.896, coverage 55.0%, density L1 0.082. Paired known-inclusion versus oracle absolute-error difference −0.057, bootstrap interval [−0.199, +0.083] (inconclusive).
- Sensitivity: tested grid changes moved posterior means by at most 8.51×10⁻⁴ (coarse) and 1.03×10⁻⁴ (fine) km s⁻¹ Mpc⁻¹. Bandwidth factors 0.5–1.2 materially changed density and posterior summaries; the apparent improvement in catalogue-only H0 bias at high bandwidth is bias compensation, not reliable density recovery. Selection-curve assumptions also changed bias and coverage.
- Oracle screening: the same workflow ran 300 prior/fixed-truth replicates on a coarser grid and reported aggregate summaries in artifact `12865285154`. That workflow does not save raw rows/seeds, so it is not a row-level independently replayed artifact; it is weaker evidence than the 1,200-row H0-dependent-noise SBC artifact.

The dedicated [screening report](host_density_sensitivity_screening_results_2026-10-29.md) contains full tables and archive provenance. CHIMERA v2.2 remains unrun because exact target inputs have not been verified; the literature review remains focused rather than systematic; no domain-expert review is claimed.


## Exact-density oracle artifact (latest CI job)

At [run #40802415000](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40802415000), the exact-density continuous-target oracle job completed successfully. Its independent audit replayed all 1,200 truth/distance rows, reconstructed all posterior summaries and aggregates, verified the manifest, and passed with summary errors below 1.9×10⁻¹³. This is a fixed-noise low-redshift toy distinct from the H₀-dependent-noise experiment; neither validates a real detector pipeline.
