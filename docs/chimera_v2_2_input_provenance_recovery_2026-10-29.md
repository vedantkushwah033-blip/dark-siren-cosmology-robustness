# CHIMERA v2.2 input provenance and recovery checklist — 29 October 2026

## Decision

**The pinned incomplete-catalogue benchmark is blocked by an input-provenance gap, not merely by installation or code execution.** Do not report it as reproduced until the exact notebook inputs are obtained, checksummed, and used in a successful run.

Pinned source: CHIMERA commit `fd35f897e36850d865c760c13325c5a46ebe94db`, notebook `examples/incomplete_galaxy_catalog.ipynb`.

## Direct source inspection

The pinned notebook's first data-setup cell sets:

```python
dir_data = "/leonardo_work/IscrC_LIGEA3/O5_B25/data/"
mock_id = 1
mcut = 10
file_gw_ev_pixelated = "./data/O5snr25_100evs_mock1_pixelated_nsidelist512-256-128-64-32-16-8_meanpixels15_skyconf0.9_nestFalse.h5"
file_inj = "./data/injections_Ninj_2e7_O5Like_snr20.h5"
file_compl_mask = "./data/completeness_masks_5x_pix32False.h5"
file_compl_intp = "./data/completeness_masks_5x_pix32False_mock1_mcut10.h5"
dir_temp_mask_idx = "./data/temp_gw_mask_idxs_mock1_mcut10.h5"
```

The notebook then loads the precomputed incomplete-catalogue products from `./data/pcat_incomplete.h5`. It applies an SNR threshold of 25 to the injection data. The hard-coded `dir_data` points to a project-specific HPC filesystem path; it is not a portable public download instruction.

## Public-source checks performed

1. The official GitHub v2.2 release API response has an empty `assets` list. The release notes describe the software release and example notebooks, but do not provide a bundled dataset archive.
2. The repository tree at the pinned v2.2 commit contains the notebook and software, but no matching HDF5 data products under the repository tree.
3. The public v2.0 Zenodo fixture previously inspected by this project does not match the notebook's required 100-event pixelated SNR-25 mock. It is not a valid replacement for this benchmark.
4. Exact-name searches for `pcat_incomplete.h5` and `O5snr25_100evs_mock1_pixelated` did not identify a publicly indexed copy. Search absence is not proof that no copy exists elsewhere.

Official references:
- Release: https://github.com/CosmoStatGW/CHIMERA/releases/tag/v2.2
- Pinned notebook: https://github.com/CosmoStatGW/CHIMERA/blob/fd35f897e36850d865c760c13325c5a46ebe94db/examples/incomplete_galaxy_catalog.ipynb
- Pinned repository tree: https://github.com/CosmoStatGW/CHIMERA/tree/fd35f897e36850d865c760c13325c5a46ebe94db

## Exact input manifest required before running

Recover and record all of the following, with source, license/access conditions, file size, and SHA-256:

- [ ] `O5snr25_100evs_mock1_pixelated_nsidelist512-256-128-64-32-16-8_meanpixels15_skyconf0.9_nestFalse.h5`
- [ ] `injections_Ninj_2e7_O5Like_snr20.h5` (verify that this is the intended file and document the notebook's SNR-25 filtering)
- [ ] `completeness_masks_5x_pix32False.h5`
- [ ] `completeness_masks_5x_pix32False_mock1_mcut10.h5`
- [ ] `temp_gw_mask_idxs_mock1_mcut10.h5`
- [ ] `pcat_incomplete.h5`
- [ ] Any additional files generated or consumed by later notebook cells, including cached outputs that are inputs to the benchmark rather than mere display artifacts.

For each file, confirm schema/dataset names, mock ID, magnitude cut, HEALPix nside/order, sky-confidence convention, SNR selection, and consistency with the pinned notebook. Preserve the original files read-only and log the exact code commit, package versions, command, exit status, and output artifact hashes.

## Recovery path

The strongest next action is to request the exact data bundle or an authoritative archival DOI from the CHIMERA maintainers/authors of the v2.2 incomplete-catalogue example. Ask specifically for the six named products above, their checksums, redistribution permission, and any setup script that generates them from parent products. Do not imply that such contact has already occurred.

If the original products cannot be obtained, use one of two clearly separated options:

- **Option A — no external benchmark claim:** finalize the synthetic-methodology paper and state that the external comparison remains blocked.
- **Option B — new, version-matched reproduction:** construct an independently documented dataset and run a separately labelled CHIMERA-compatible experiment. This is a new benchmark, not reproduction of the pinned notebook; all differences and deviations must be prespecified.

## Release gate

The benchmark may be marked PASS only when the exact input manifest is complete, hashes and provenance are independently checked, the pinned notebook or a fully logged equivalent completes from a clean environment, and its numerical outputs are independently reconstructed. Until then the project is **NO-GO for claims of CHIMERA v2.2 validation or external-method superiority**.
