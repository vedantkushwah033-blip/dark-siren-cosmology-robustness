# CHIMERA v2.2 benchmark input provenance audit (2026-10-29)

## Decision

**External benchmark status: NO-GO for the planned incomplete-catalogue replication.** The pinned CHIMERA v2.2 source and notebook can be inspected, but the exact version-compatible input bundle required by that notebook was not found in the tagged release assets or the repository tree during this audit. Do not substitute a different-version data bundle and call it a v2.2 reproduction.

## Source inspection

- Pinned source release: [CosmoStatGW/CHIMERA v2.2 release](https://github.com/CosmoStatGW/CHIMERA/releases/tag/v2.2), published 9 June 2026.
- Release API metadata reports no downloadable release assets.
- The v2.2 source tree contains [examples/incomplete_galaxy_catalog.ipynb](https://github.com/CosmoStatGW/CHIMERA/blob/v2.2/examples/incomplete_galaxy_catalog.ipynb). Its first data-loading cell expects local files including:
  - `./data/O5snr25_100evs_mock1_pixelated_nsidelist512-256-128-64-32-16-8_meanpixels15_skyconf0.9_nestFalse.h5`
  - `./data/injections_Ninj_2e7_O5Like_snr20.h5`
  - `./data/completeness_masks_5x_pix32False.h5`
  - `./data/completeness_masks_5x_pix32False_mock1_mcut10.h5`
  - `./data/temp_gw_mask_idxs_mock1_mcut10.h5`
- Those expected data files are not present in the tagged repository tree; the notebook refers to a local `./data/` directory.
- The v2.2 release notes explicitly state that legacy example notebooks and `galcat_test.h5` were removed, and describe API-breaking changes.

Relevant pinned source paths:
- [v2.2 incomplete-catalogue notebook](https://github.com/CosmoStatGW/CHIMERA/blob/v2.2/examples/incomplete_galaxy_catalog.ipynb)
- [v2.2 completeness implementation](https://github.com/CosmoStatGW/CHIMERA/blob/v2.2/CHIMERA/catalog/completeness.py)
- [v2.2 release notes](https://github.com/CosmoStatGW/CHIMERA/releases/tag/v2.2)

## Why the older public bundle is not an acceptable substitute

The public [CHIMERA v2.0 Zenodo record](https://zenodo.org/records/17143346) lists an older software/data bundle with `PE_O5Like_snr20.h5`, `injections_Ninj_2e7_O5Like_snr20.h5`, `parent_galaxy_catalog_1.6M_dVdz.h5`, and an SNR-20 pixelated PE file. It does not list the SNR-25, 100-event mock or the v2.2 incomplete-catalogue completeness products required above. It is a different version and a different data design, so it cannot silently stand in for the pinned v2.2 benchmark.

This audit used release metadata, the versioned repository tree, the notebook's explicit file paths, and the public Zenodo record listing. It did **not** download or execute the CHIMERA notebook, and it does not claim a benchmark run.

## Re-entry criteria

The benchmark can move from NO-GO only after the exact required files (or an authoritative version-matched source for them) are located and their provenance, checksums, sizes, licensing/access conditions, and compatibility with the v2.2 notebook are recorded. Then:
1. pin the v2.2 source commit and environment;
2. verify every input file by checksum and expected HDF5 keys/shapes;
3. run the documented notebook or a minimal, explicitly equivalent script end to end;
4. preserve stdout/stderr, runtime, software versions, outputs, and checksums;
5. independently inspect the resulting likelihood/posterior and compare against a declared reference.

Until these gates pass, the project's defensible endpoint remains a synthetic methodological study with an unresolved external-validation gate—not a real-data replication.
