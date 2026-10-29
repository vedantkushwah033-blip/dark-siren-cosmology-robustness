# Release-gate update — 29 October 2026

## Exact-head automated validation

Working branch: `research/audit-posterior-grid-integrity`  
Current head at this update: `20b4a17e216e9de8d76e7b4bfbe5a6f6b03e9dca`  
Exact-head Actions run: [#40854750158](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40854750158); independent duplicate exact-head run: [#40854458789](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40854458789). Earlier code-head run: [#40832784066](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40832784066); earlier pull-request merge-ref run: [#40832806297](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40832806297), testing merge ref `29308d22e59d4e1c8f93538953e92de67fe27881`.

Both exact-head runs completed successfully; each passed the numerical-diagnostics matrix on Python 3.10, 3.11, and 3.12. Jobs exercise likelihood composition, H0-dependent noise and selection-conditioned inference, independent artifact auditors, catalogue thinning, host-density estimation, selection-function misspecification, distance-uncertainty sensitivity, posterior-grid stability, and associated diagnostics. These results support software reproducibility on the tested environments, not realism of the synthetic assumptions. The cancelled duplicate run #40854450158 is not used as validation evidence.

## External CHIMERA v2.2 benchmark gate

**Status remains NO-GO until the exact pinned v2.2 inputs are obtained and independently verified.**

The official v2.2 release page exposes source archives but no verified compatible mock-data bundle as a release asset. The pinned incomplete-catalogue notebook at commit `fd35f897e36850d865c760c13325c5a46ebe94db` expects mock-specific products including a 100-event pixelated SNR-25 mock and `pcat_incomplete.h5`, plus completeness and selection products. Those exact compatible inputs have not been located and verified. The documented v2.0 fixtures are not an acceptable substitute for a v2.2 benchmark. Notebook saved outputs alone do not count as a reproduced benchmark.

Official references:
- [CHIMERA v2.2 release](https://github.com/CosmoStatGW/CHIMERA/releases/tag/v2.2)
- [Pinned incomplete-catalogue notebook](https://github.com/CosmoStatGW/CHIMERA/blob/fd35f897e36850d865c760c13325c5a46ebe94db/examples/incomplete_galaxy_catalog.ipynb)
- [Existing execution contract](chimera_v2_2_notebook_execution_contract_2026-10-26.md)

## Scientific disposition

The controlled synthetic results can support a carefully scoped methods manuscript draft and external review. They do not establish a real-data H0 estimate, a result from LVK/GWTC observations, superiority to CHIMERA or other external methods, or completed peer review/publication. Preserve these release gates in the README and manuscript until the missing benchmark inputs are recovered and the pinned pipeline is executed with auditable outputs.


## Confirmatory seed dependence and bias-compensation follow-up

The 750-replicate confirmatory study has overlapping seed ranges across its three condition-specific jobs: 150 seeds are shared by low-floor/high-floor, 150 by high-floor/broad-transition, and 50 by low-floor/broad-transition. Within-condition paired contrasts remain valid, but cross-condition estimates must not be described as independent replications. The four-cell shallow-transition follow-up shows an RMSE advantage at true H0=60, but that advantage disappears or reverses at H0=80; an interval spanning zero is inconclusive, not affirmative evidence for either arm.


## Status correction — 29 October 2026, latest exact-head check

The latest PR-triggered numerical-diagnostics run checked during this update is [run #40874302252](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40874302252). All three matrix jobs completed successfully: Python 3.10, 3.11, and 3.12; no failed steps were reported. The workflow exercised the test suite's synthetic likelihood, selection normalization, artifact-auditor and numerical-diagnostic checks. This is software-level evidence only.

The prior documentation/seed-overlap commit's workflow [#40858566961](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40858566961) is **cancelled** and is not counted as a passing run. The newer successful matrix above is the current passing evidence identified in this check.

The independent equation audit at [docs/independent_likelihood_equation_audit_2026-10-29.md](independent_likelihood_equation_audit_2026-10-29.md) documents the implemented discrete-host mixture likelihood, H0-dependent noise, fixed-detected-count selection correction, positive-truncation normalization, and endpoint-inclusive prior normalization. It is an equation-level review, not an independent numerical rerun or detector-level validation.

**Pull-request integration remains a separate gate:** PR #7 is open and unmerged, has 82 changed files and 173 commits ahead according to the PR description, and targets `experiment/catalogue-density-estimators` rather than `main`. The base branch should be confirmed as intentional; the large accumulated diff deserves maintainer review and may warrant splitting before merge. No automatic merge or retargeting was performed.

**Still NO-GO:** CHIMERA v2.2 exact-input benchmark, real LVK/GWTC observational validation, external domain-expert sign-off, and publication/peer-review claims. The internal manuscript remains a synthetic-methodology draft.
