# Research-readiness report — 29 October 2026

## Decision

**Current disposition: GO for a carefully scoped synthetic-methodology manuscript draft; NO-GO for claims about real GW catalogues, real H0, or external-method superiority.** PR #7 remains a draft and is not merged.

## Completed and artifact-verified

1. **KDE bandwidth convention:** replaced the ambiguous SciPy covariance-factor usage with a direct reflected weighted-kernel sum; explicit factor semantics now use the weighted SD of the original sample, with a direct kernel-sum regression test.
2. **Catalogue incompleteness benchmark:** a locked 500-replicate Stage B experiment completed with six paired arms. Its independent posterior-grid and raw-table audit passed. Primary contrast: 5.9519 km s⁻¹ Mpc⁻¹ reduction in absolute posterior-mean error, bootstrap 95% interval [5.5410, 6.3573].
3. **H0-dependent-noise model:** a separate 1,200-row synthetic calibration artifact passed independent generator replay and posterior reconstruction to floating-point precision. A follow-up audit on commit `9c662b6` also independently reaggregated every headline summary from the raw CSV and passed. Prior-predictive mean rank 0.514; ten-bin descriptive chi-square p=0.450.
4. **Transparent documentation:** added result reports, protocol, literature-positioning audit, and a first integrated manuscript draft. The manuscript clearly distinguishes the two toy models and states non-claims.
5. **Reproducibility plumbing:** raw replicate records, seed schemes, posterior-grid outputs, SHA-256 manifests, auditors, adversarial tampering tests and CI workflows are in the repository.

## Still open

- The expanded host-density sensitivity campaign must complete. Its current expedited run uses 40 replicates per sensitivity configuration and 300 exact-density replicates per truth regime. These are exploratory screening designs, not a high-precision confirmatory campaign.
- The general full-suite CI must pass on the latest branch head. The specialized H0-dependent-noise SBC workflow has passed on source commit `9c662b657c4472b9614bf013625af8647e1a2d5d`, with 1,200 rows and the enhanced auditor independently reconstructing all saved per-replicate posteriors and all summary aggregates. That specialized workflow is not a substitute for the full Python 3.10/3.11/3.12 suite; the latest exact-head full suite is still pending/queued.
- The CHIMERA v2.2 external benchmark is not complete. The exact expected input dataset(s), compatibility with the tagged release, and execution output have not been verified. Do not substitute files from an older release.
- A full reading and assumption-by-assumption comparison of the closest 2026 completeness/GP methods is outstanding. The repository includes a focused literature/positioning audit, not a systematic review.
- No expert review, arXiv submission, peer review, or publication claim is made.

## Required gates before release

1. Obtain a passing exact-head Python 3.10/3.11/3.12 test matrix and preserve its run URL.
2. Finish the host-density campaign, inspect raw tables and manifest, and require independent reconstruction checks to pass. Update the local-snapshot report using the regenerated artifact rather than mixing results across code versions.
3. Add a concise comparison table distinguishing assumptions in this project, the 2026 GP reconstruction study, the 2026 gwcosmo completeness test, and the host-population inference preprint.
4. Inspect CHIMERA v2.2 data availability and run feasibility from the pinned release; document a no-go if exact inputs are unavailable.
5. Re-read the manuscript tables against artifacts, have a knowledgeable researcher review the model assumptions, and retain all toy-model limitations in any abstract or outreach.

## Claim language

Acceptable: “A reproducible controlled synthetic benchmark quantifies the effect of a known catalogue-inclusion mechanism and several selection-law assumptions on a simplified dark-siren inference model.”

Unacceptable: “Solved the Hubble tension,” “measured a new H0,” “first method for incomplete galaxy catalogues,” “validated LVK performance,” or “inverse weighting universally fixes incompleteness.”

## Final project scope

A useful and defensible outcome is a reproducible scientific-computing / statistical-methods study with explicit controls and independent audits. A real astrophysical result requires external data, a realistic detector/selection likelihood, appropriate population treatment and independent subject-matter review. The current simulations do not cross that boundary.
