# CHIMERA benchmark target/version audit (2026-10-26)

## Decision

**Do not treat the current CHIMERA v2.0 Zenodo fixture as a complete, ready-to-run reproduction of the incomplete-catalogue tutorial.** Inspection of the pinned source tree found a material version/data mismatch that must be resolved before an external benchmark can be claimed.

## Updated candidate choice after checking official releases

The official GitHub release list contains **CHIMERA v2.2**, published 2026-06-09, at commit `fd35f897e36850d865c760c13325c5a46ebe94db`. Its release notes explicitly list `examples/incomplete_galaxy_catalog.ipynb` and describe API-breaking changes relative to earlier releases. Therefore, if the scientific target is specifically incomplete-catalogue inference, **v2.2 is a more coherent code target than v2.0**.

This does not yet make the benchmark reproducible: the v2.2 notebook's 100-event SNR-25 mock posterior file and saved derived products are not listed in the v2.0 Zenodo data record. The official v2.1 Zenodo record is a software archive, not evidence that these exact notebook inputs are available there. The v2.2 source pin and every data product must be tracked separately, with provenance and generation settings established before execution.

Candidate reference:
- Code release: [CHIMERA v2.2](https://github.com/CosmoStatGW/CHIMERA/releases/tag/v2.2)
- Exact commit: `fd35f897e36850d865c760c13325c5a46ebe94db`
- Notebook: [incomplete-galaxy-catalog example at the pinned v2.2 commit](https://github.com/CosmoStatGW/CHIMERA/blob/fd35f897e36850d865c760c13325c5a46ebe94db/examples/incomplete_galaxy_catalog.ipynb)

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

## Additional preprocessing path confirmed from current-main source

Source inspection confirms that some derived products can be regenerated in principle, provided the correct underlying catalogue and GW-event inputs are obtained:

- `mask_completeness.get_masks(..., dir_masks=None)` computes sky masks from galaxy right ascension and declination and saves the mask products.
- `mask_completeness.get_completeness_interpolants(dir_compl=None)` computes and saves completeness-vs-redshift interpolants, requiring the relevant galaxy redshifts, theoretical galaxy-density function, and matching completeness settings.
- `set_gw_mask_idxs(pix_gw, nside_gw, nest_gw)` maps event sky pixels onto the completeness masks. The notebook currently loads these indices from `temp_gw_mask_idxs_mock1_mcut10.h5`.
- `precompute_p_cat(completeness, cosmo, z_grids, fname_data_gal, data_gw_pixelated, z_err=...)` can create the `p_cat`, `N_gal`, and `P_compl` arrays; `save_p_cat`/ `load_p_cat` persist and reload them. This is described upstream as an expensive one-time computation.

This narrows the unknowns but does **not** resolve them: the exact source catalogue used by the current-main notebook, its relationship to the v2.0 parent catalogue, the 100-event SNR-25 mock data, the (m_{cut}=10) galaxy selection, the redshift-error convention, and the provenance of the saved `pcat_incomplete.h5` still need to be pinned and checked. Recomputing derived files with different inputs or settings would create a new benchmark, not reproduce the notebook's original result.

## Required go/no-go steps

1. Choose one consistent reference: either a release/commit that contains the desired incomplete-catalogue example, or a different v2.0-supported scientific question.
2. Pin the exact code commit and data-record versions together.
3. Trace the origin and generation commands for the mock event file, masks, completeness interpolants, temporary mask-index file, and `pcat_incomplete.h5`.
4. Establish which inputs are downloadable, licensed for the intended use, and reproducible from upstream source.
5. Run a small upstream-supported end-to-end smoke case before comparing posterior outputs.
6. Only then compare the reference and this repository, after matching estimands, priors, selection treatment, and population assumptions.

Until these gates pass, the external pipeline benchmark is **prepared but not executed**, and real-data/real-pipeline validation remains incomplete.

## Follow-up source and reproducibility check (2026-10-26)

### Public repository does not itself supply the tutorial's data directory

The v2.2 notebook contains executed output showing that it loaded the five masks, five completeness interpolants, and a precomputed catalogue likelihood file, and then evaluated 100 events over 50 H0 grid points. But its code refers to local relative paths under `./data/`; the notebook does not provide download URLs or a script for obtaining the exact six files listed above. The line `# precomputed elesewhere due to very big gaalxy catalog` explicitly identifies `pcat_incomplete.h5` as an externally precomputed product rather than something generated by the notebook itself.

The CHIMERA issue tracker also contains a 2024 report that a tutorial's required `data` directory was absent from the public checkout. A maintainer replied that a download script was being worked on. This is historical context only; it does not prove the status of the current v2.2 files, but it reinforces that notebook presence and data availability must be verified separately:
- [Issue #2: 'data' file not present](https://github.com/CosmoStatGW/CHIMERA/issues/2)

### Packaging metadata has a licence inconsistency

At the pinned v2.2 commit, the repository README badge says GPLv3 and the `LICENSE` file contains the GNU General Public License v3 text, while `pyproject.toml` declares `license = "MIT"`. This is an upstream metadata inconsistency; this project should not infer that the code is MIT-licensed from the packaging field alone. The benchmark plan remains to install the pinned upstream dependency separately and avoid copying CHIMERA source into this repository. If redistribution or source incorporation becomes necessary, ask upstream maintainers to clarify the intended licence and follow the authoritative licence file/legal guidance.

Sources inspected:
- [v2.2 README at the pinned commit](https://github.com/CosmoStatGW/CHIMERA/blob/fd35f897e36850d865c760c13325c5a46ebe94db/README.md)
- [v2.2 LICENSE at the pinned commit](https://github.com/CosmoStatGW/CHIMERA/blob/fd35f897e36850d865c760c13325c5a46ebe94db/LICENSE)
- [v2.2 pyproject.toml at the pinned commit](https://github.com/CosmoStatGW/CHIMERA/blob/fd35f897e36850d865c760c13325c5a46ebe94db/pyproject.toml)

### Current decision

The source audit is more precise, but the benchmark remains **NO-GO for exact end-to-end reproduction** until the exact 100-event SNR-25 event file and the underlying mock galaxy catalogue / precomputed `pcat_incomplete.h5` are located with documented provenance, or regenerated from verified source data with a fully specified recipe. Do not treat the saved notebook output as an independently rerun result. Do not substitute the v2.0 300-event SNR-20 sample and label it a v2.2 replication.

### Additional upstream release/archive inventory (checked 2026-10-26)

A recursive tree inspection of the pinned v2.2 commit found the incomplete-catalogue notebook and CHIMERA's Python module `CHIMERA/data.py`, but no tracked `data/` directory or dedicated download helper for the notebook's six required products. The GitHub Releases API lists v2.2, v2.1, v2.0, and v1.0.0 releases, each with **no attached release assets**. Thus there is no release-attached archive that resolves the missing v2.2 products.

The public Zenodo record linked to CHIMERA v2.0 lists four data files: the 20-million-injection sample, a 1.6-million-galaxy parent catalogue, and pixelated/non-pixelated SNR-20 PE samples. It does not list the v2.2 notebook's 100-event SNR-25 mock, completeness masks/interpolants, GW-mask cache, or `pcat_incomplete.h5`. This confirms that the v2.0 record is a useful source for a separate v2.0 fixture, but not a complete source for the matched v2.2 benchmark.

Evidence:
- [Pinned v2.2 repository tree](https://github.com/CosmoStatGW/CHIMERA/tree/fd35f897e36850d865c760c13325c5a46ebe94db)
- [CHIMERA GitHub releases](https://github.com/CosmoStatGW/CHIMERA/releases)
- [CHIMERA v2.0 Zenodo archive and file inventory](https://zenodo.org/records/17143346)

**Decision unchanged:** do not download the v2.0 fixture and present it as v2.2. The next defensible step is to request/locate the missing data provenance or explicitly design a separately labelled v2.0 compatibility benchmark.


## Licence metadata discrepancy — do not infer redistribution rights

A separate inspection of the pinned v2.2 source metadata found an unresolved inconsistency: the repository README and top-level `LICENSE` identify GPL-3.0, while `pyproject.toml` declares MIT. These declarations are not interchangeable. Until the upstream maintainers clarify which licence applies to the distributed code, this project should treat CHIMERA as an external, separately installed dependency and should not copy its source files or assume that the MIT declaration permits redistribution. The licence of the scientific data files is a separate question and must be checked independently; a code licence does not automatically establish the data licence.

This is a provenance/compliance blocker, not evidence that the numerical pipeline is incorrect. It does not prevent a local, properly attributed benchmark in principle, but the exact code version, data terms, and redistribution permissions must be recorded before sharing a bundled benchmark artifact.
