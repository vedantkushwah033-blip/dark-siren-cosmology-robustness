# CHIMERA benchmark target/version audit (2026-10-26)

## Decision

**Do not treat the current CHIMERA v2.0 Zenodo fixture as a complete, ready-to-run reproduction of the incomplete-catalogue tutorial.** Inspection of the pinned source tree found a material version/data mismatch that must be resolved before an external benchmark can be claimed.

## What was checked

- Reference code commit: `080a1ad716b02867a8d1dba27b3724a495ee3990`.
- The pinned commit's recursive Git tree contains `docs/user_guide/tutorials/empty_galaxy_catalog.ipynb` and `docs/user_guide/tutorials/complete_galaxy_catalog.ipynb`, but **does not contain** `examples/incomplete_galaxy_catalog.ipynb`.
- The pinned `empty_galaxy_catalog.ipynb` is a GW-only/empty-catalog example. It loads non-pixelated posterior samples and injections, defines a population and selection function, and computes an H0 posterior. It does not establish a test of the incomplete-catalogue correction.
- The incomplete-catalogue example exists on the moving `main` branch, but it is not part of the pinned v2.0 commit. It must not be silently combined with the v2.0 code and described as a single pinned reference.
- The current-main notebook references these inputs/products:
  - `O5snr25_100evs_mock1_pixelated_nsidelist512-256-128-64-32-16-8_meanpixels15_skyconf0.9_nestFalse.h5`
  - `injections_Ninj_2e7_O5Like_snr20.h5`
  - `completeness_masks_5x_pix32False.h5`
  - `completeness_masks_5x_pix32False_mock1_mcut10.h5`
  - `temp_gw_mask_idxs_mock1_mcut10.h5`
  - `pcat_incomplete.h5`
- Of these, the v2.0 Zenodo fixture currently prepared by this repository includes the injection file and the parent catalogue and pixelated PE file listed in that record. That is **not enough by itself** to reproduce the current-main notebook: the notebook uses a different mock-event file and several additional derived completeness/catalogue products. Their provenance and exact generation steps have not yet been established.
- The upstream README identifies CHIMERA as GPL-3.0 licensed. Keep the benchmark as a separately installed dependency and do not copy upstream source into this repository without checking licence obligations.

## Sources inspected

- [Pinned CHIMERA v2.0 source tree](https://github.com/CosmoStatGW/CHIMERA/tree/080a1ad716b02867a8d1dba27b3724a495ee3990)
- [Pinned empty-catalogue tutorial](https://github.com/CosmoStatGW/CHIMERA/blob/080a1ad716b02867a8d1dba27b3724a495ee3990/docs/user_guide/tutorials/empty_galaxy_catalog.ipynb)
- [Current-main incomplete-catalogue example](https://github.com/CosmoStatGW/CHIMERA/blob/main/examples/incomplete_galaxy_catalog.ipynb) — useful for discovery only, not a pinned reproducibility target
- [CHIMERA README and licence statement](https://github.com/CosmoStatGW/CHIMERA/blob/080a1ad716b02867a8d1dba27b3724a495ee3990/README.md)
- [CHIMERA v2.0 mock-data record](https://zenodo.org/records/17143346)

## Required go/no-go steps

1. Choose one consistent reference: either a release/commit that contains the desired incomplete-catalogue example, or a different v2.0-supported scientific question.
2. Pin the exact code commit and data-record versions together.
3. Trace the origin and generation commands for the mock event file, masks, completeness interpolants, temporary mask-index file, and `pcat_incomplete.h5`.
4. Establish which inputs are downloadable, licensed for the intended use, and reproducible from upstream source.
5. Run a small upstream-supported end-to-end smoke case before comparing posterior outputs.
6. Only then compare the reference and this repository, after matching estimands, priors, selection treatment, and population assumptions.

Until these gates pass, the external pipeline benchmark is **prepared but not executed**, and real-data/real-pipeline validation remains incomplete.
