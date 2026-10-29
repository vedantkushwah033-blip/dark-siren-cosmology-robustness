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
