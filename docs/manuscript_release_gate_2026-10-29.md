# Manuscript release gate and final audit checklist — 29 October 2026

## Current disposition

**Suitable for internal circulation as a controlled synthetic-methodology draft; not yet ready to be described as externally validated, peer reviewed, or publication-ready.** The exact CHIMERA v2.2 external benchmark remains NO-GO pending exact, provenance-traceable input data. This checklist is a gate, not a claim that every item has passed.

## Claim-to-artifact audit

| Manuscript claim | Supporting artifact/evidence | Release treatment |
|---|---|---|
| Stage B selection benchmark: 500 paired replicates and six arms | Stage B raw artifact, independent posterior-grid/raw-table auditor, reported CI workflow | Keep as synthetic and conditional on the stated simulator. |
| Stage B primary paired absolute-error reduction 5.9519 km s⁻¹ Mpc⁻¹; bootstrap 95% interval [5.5410, 6.3573] | Raw-table reaggregation and independent posterior reconstruction documented in readiness report | Preserve the paired estimand and interval; do not translate into real-data performance. |
| Host-density sensitivity: 13 settings × 40 paired replicates × three arms | Campaign artifact ID `12866753211`; 1,560 fits reconstructed, 403 summaries reaggregated, 14 manifest files checked | Label exploratory; disclose the single boundary-sensitive coverage flag and bandwidth dependence. |
| H0-dependent-noise calibration: 1,200 synthetic datasets | Seed replay, posterior reconstruction and summary audit documented in the dedicated report | Separate prior-predictive SBC from conditional fixed-truth diagnostics; do not call it detector validation. |
| Selection-law replication: three 250-replicate conditions | Workflow `40819870231`; condition-specific posterior-grid/aggregate audits | Retain the shared-seed qualification; do not call the three conditions independent replications. |
| Comparison to other methods | Focused literature-positioning audit | Describe as focused, not systematic; no superiority claim without same-data, same-model comparison. |
| CHIMERA v2.2 reproduction | Exact required inputs have not been verified available or run | **NO-GO.** Do not claim reproduction or substitute a different-version mock. |

## Final checks required before any public submission

### Scientific
- [ ] Independent researcher reviews the generative model, likelihood normalization, selection conditioning, KDE reflection/bandwidth convention, and posterior-grid numerical accuracy.
- [ ] Verify that each manuscript table and abstract number maps to the intended artifact, configuration, seed range, and audit report.
- [ ] Retain all negative controls, inconclusive contrasts, boundary-sensitive discrepancies, shared-seed qualifications, and model-form limitations.
- [ ] Avoid language implying a real (H_0) measurement, LVK/GWTC validation, a solution to the Hubble tension, universal weighting superiority, or novelty priority.

### Reproducibility
- [ ] Run the full test suite and numerical diagnostics on the exact final code SHA and record run URLs/results.
- [ ] Ensure all raw outputs used for headline results have a manifest with byte counts and SHA-256 checksums.
- [ ] Verify clean-environment install/run instructions and dependency/version pins.
- [ ] Confirm that each figure/table can be regenerated from committed scripts and documented artifacts.
- [ ] Preserve a record of all workflow artifacts and the exact source commit used for each result.

### Manuscript and release package
- [ ] Proofread equations, notation, units, table captions, references, and internal relative links.
- [ ] Ensure the data/code availability statement matches actual access and licensing; do not promise public redistribution of restricted external data.
- [ ] Add author contributions, funding, conflicts-of-interest and acknowledgements only when accurate; do not invent declarations.
- [ ] Select a venue only after reading its official scope, data/code policy, and submission requirements.
- [ ] Obtain domain-expert feedback before submission.
- [ ] Do not label the draft submitted, accepted, peer reviewed, or published unless that event actually occurs.

## External benchmark recovery gate

Reopen the CHIMERA v2.2 gate only when every required input is available from an authoritative source with access/licence terms and SHA-256 checksums, schemas match the pinned notebook, the exact pinned environment and commands are recorded, and an independent audit can reconstruct the run. Until then, an honest synthetic-only manuscript is preferable to a misleading external-validation claim.

## Status log

- 2026-10-29: Existing synthetic campaigns and selected independent artifact audits are documented in the repository.
- 2026-10-29: Exact CHIMERA v2.2 data products remain unverified; an unsent maintainer request draft is available at `docs/chimera_v2_2_data_request_draft.md`.
- 2026-10-29: Current manuscript remains internal draft v0.3. No expert review, submission, peer review, or publication is claimed.
