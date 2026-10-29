# Research-readiness report — 29 October 2026

## Decision

**Current disposition: GO for a carefully scoped synthetic-methodology manuscript draft; NO-GO for claims about real GW catalogues, real H0, or external-method superiority.** PR #7 remains a draft and is not merged.

## Completed and artifact-verified

1. **KDE bandwidth convention:** replaced the ambiguous SciPy covariance-factor usage with a direct reflected weighted-kernel sum; explicit factor semantics now use the weighted SD of the original sample, with a direct kernel-sum regression test.
2. **Catalogue incompleteness benchmark:** a locked 500-replicate Stage B experiment completed with six paired arms. Its independent posterior-grid and raw-table audit passed. Primary contrast: 5.9519 km s⁻¹ Mpc⁻¹ reduction in absolute posterior-mean error, bootstrap 95% interval [5.5410, 6.3573].
3. **H0-dependent-noise model:** a separate 1,200-row synthetic calibration artifact passed independent generator replay and posterior reconstruction to floating-point precision. A follow-up audit on commit `9c662b6` also independently reaggregated every headline summary from the raw CSV and passed. Prior-predictive mean rank 0.514; ten-bin descriptive chi-square p=0.450.
4. **Transparent documentation:** added result reports, a literature-positioning audit with an assumption-by-assumption comparison table, and an integrated manuscript draft. The manuscript clearly distinguishes the two toy models and states non-claims.
5. **Reproducibility plumbing:** raw replicate records, seed schemes, posterior-grid outputs, SHA-256 manifests, auditors, adversarial tampering tests and CI workflows are in the repository.

## Still open

- The expanded host-density screen contains 13 configurations × 40 paired replicates. All 13 configurations passed local per-configuration independent replay/posterior reconstruction (1,560 fit rows), and their summary statistics were reaggregated independently. The combined clean campaign and one-shot full audit are running at [#40802415000](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40802415000); check its final status and audit JSON before claiming exact-head CI PASS. The screen remains exploratory with imprecise coverage.
- **Latest full test matrix passed on code head `6d5f999e`.** The Python 3.10/3.11/3.12 matrix completed successfully in [run #40802418966](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40802418966). This establishes the test suite for that code revision, not physical validity of a real-data analysis. The new unified host-density campaign/audit remains in progress in [run #40802415000](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40802415000). The specialized H0-dependent-noise SBC workflow also passed on source commit `9c662b657c4472b9614bf013625af8647e1a2d5d`, with 1,200 rows and an independent auditor reconstructing saved per-replicate posteriors and summary aggregates.
- The CHIMERA v2.2 external benchmark is not complete. The exact expected input dataset(s), compatibility with the tagged release, and execution output have not been verified. Do not substitute files from an older release.
- A full reading and assumption-by-assumption comparison of the closest 2026 completeness/GP methods is outstanding. The repository includes a focused literature/positioning audit, not a systematic review.
- No expert review, arXiv submission, peer review, or publication claim is made.

## Required gates before release

1. **Completed for current branch head:** preserve and cite the passing exact-head Python 3.10/3.11/3.12 matrix [run 40800561553](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40800561553). If additional code or documentation commits are made, rerun and verify CI for that new head before release.
2. The 13-configuration host-density screen is now complete and its canonical artifact/report are checked. The full 40-replicate-per-setting campaign has not been reconstructed row-by-row by an independent posterior implementation; the independent posterior smoke audit covers two replicates only. Treat the screen as exploratory and do not claim full independent replay.
3. Add a concise comparison table distinguishing assumptions in this project, the 2026 GP reconstruction study, the 2026 gwcosmo completeness test, and the host-population inference preprint.
4. Inspect CHIMERA v2.2 data availability and run feasibility from the pinned release; document a no-go if exact inputs are unavailable.
5. Re-read the manuscript tables against artifacts, have a knowledgeable researcher review the model assumptions, and retain all toy-model limitations in any abstract or outreach.

## Claim language

Acceptable: “A reproducible controlled synthetic benchmark quantifies the effect of a known catalogue-inclusion mechanism and several selection-law assumptions on a simplified dark-siren inference model.”

Unacceptable: “Solved the Hubble tension,” “measured a new H0,” “first method for incomplete galaxy catalogues,” “validated LVK performance,” or “inverse weighting universally fixes incompleteness.”

## Current follow-up documentation

- The host-density campaign's authoritative full results are in [host-density sensitivity screening results](host_density_sensitivity_screening_results_2026-10-29.md); the supplemental [follow-up validation note](followup_validation_host_density_and_sbc_2026-10-29.md) summarizes the result and cross-checks the baseline numbers.
- Draft PR #8 contains paired grid/bandwidth sensitivity validation and remains separate from the main draft PR. Neither PR should be merged automatically; scientific review and the current-head test matrix remain gates.

## Final project scope

A useful and defensible outcome is a reproducible scientific-computing / statistical-methods study with explicit controls and independent audits. A real astrophysical result requires external data, a realistic detector/selection likelihood, appropriate population treatment and independent subject-matter review. The current simulations do not cross that boundary.


## Completed screening campaign — evidence boundary (added 2026-10-29)

- Run: [host-density validation campaign #40797592794](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40797592794), source commit `5bbe5a434878af396e8c58e3ddfe7af0356e29ed`, artifact ID `12865051001`.
- Configuration: 13 settings × 40 paired replicates, fixed H0=70, 6,000 parent galaxies, 40 detected events, threshold 250 Mpc, distance noise 140 Mpc. This is screening, not a confirmatory precision run.
- Integrity: every one of the 14 manifest-listed campaign files matched its SHA-256 and byte count after download. The small independent-audit data directory also passed its two-file hash/byte check. Its independent posterior reconstruction covered two smoke replicates only; the full campaign was not regenerated record-by-record by a separate posterior implementation.
- Baseline: catalogue-only bias −8.634, RMSE 9.165, 68% coverage 10.0%, density L1 0.537; known inclusion correction bias +1.180, RMSE 3.896, coverage 55.0%, density L1 0.082. Paired known-inclusion versus oracle absolute-error difference −0.057, bootstrap interval [−0.199, +0.083] (inconclusive).
- Sensitivity: tested grid changes moved posterior means by at most 8.51×10⁻⁴ (coarse) and 1.03×10⁻⁴ (fine) km s⁻¹ Mpc⁻¹. Bandwidth factors 0.5–1.2 materially changed density and posterior summaries; the apparent improvement in catalogue-only H0 bias at high bandwidth is bias compensation, not reliable density recovery. Selection-curve assumptions also changed bias and coverage.
- Oracle screening: the same workflow ran 300 prior/fixed-truth replicates on a coarser grid and reported aggregate summaries in artifact `12865285154`. That workflow does not save raw rows/seeds, so it is not a row-level independently replayed artifact; it is weaker evidence than the 1,200-row H0-dependent-noise SBC artifact.

The dedicated [screening report](host_density_sensitivity_screening_results_2026-10-29.md) contains full tables and archive provenance. CHIMERA v2.2 remains unrun because exact target inputs have not been verified; the literature review remains focused rather than systematic; no domain-expert review is claimed.


## Exact-density oracle artifact (latest CI job)

At [run #40802415000](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40802415000), the exact-density continuous-target oracle job completed successfully. Its independent audit replayed all 1,200 truth/distance rows, reconstructed all posterior summaries and aggregates, verified the manifest, and passed with summary errors below 1.9×10⁻¹³. This is a fixed-noise low-redshift toy distinct from the H₀-dependent-noise experiment; neither validates a real detector pipeline.
