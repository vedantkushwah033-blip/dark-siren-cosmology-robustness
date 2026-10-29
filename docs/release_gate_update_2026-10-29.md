# Release-gate update — 29 October 2026

## Exact-head automated validation

Working branch: `research/audit-posterior-grid-integrity`  
Verified head recorded for this check: `228f1729fbd932d0d4dcdaa072f55c8ae83cb8bc`  
GitHub Actions run: [#40824083391](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40824083391)

The numerical-diagnostics matrix completed successfully on Python 3.10, 3.11, and 3.12. The jobs exercised likelihood composition, H0-dependent noise and selection-conditioned inference, independent artifact auditors, catalogue thinning, host-density estimation, selection-function misspecification, distance-uncertainty sensitivity, posterior-grid stability, and the associated full-run diagnostic/audit steps. This is evidence for reproducibility and software validation on the tested environments; it is not evidence that the synthetic assumptions describe real gravitational-wave observations.

## External CHIMERA v2.2 benchmark gate

**Status remains NO-GO until the exact pinned v2.2 inputs are obtained and independently verified.**

The official GitHub release API response for tag `v2.2` reports an empty release-assets list (`assets: []`). The pinned incomplete-catalogue notebook at commit `fd35f897e36850d865c760c13325c5a46ebe94db` expects mock-specific products including a 100-event pixelated SNR-25 mock and `pcat_incomplete.h5`, plus completeness and selection products. The documented v2.0 fixture is not an acceptable substitute. Notebook saved outputs alone do not count as a reproduced benchmark.

Official references:
- [CHIMERA v2.2 release](https://github.com/CosmoStatGW/CHIMERA/releases/tag/v2.2)
- [Pinned incomplete-catalogue notebook](https://github.com/CosmoStatGW/CHIMERA/blob/fd35f897e36850d865c760c13325c5a46ebe94db/examples/incomplete_galaxy_catalog.ipynb)
- [Existing execution contract](chimera_v2_2_notebook_execution_contract_2026-10-26.md)

## Scientific disposition

The controlled synthetic results can support a carefully scoped methods manuscript draft and external review. They do not establish a real-data H0 estimate, a result from LVK/GWTC observations, superiority to CHIMERA or other external methods, or completed peer review/publication. Preserve these release gates in the README and manuscript until the missing benchmark inputs are recovered and the pinned pipeline is executed with auditable outputs.
