# CHIMERA v2.2 incomplete-catalogue notebook: execution contract (2026-10-26)

## Purpose and status

This document records a source-level inspection of the incomplete-catalogue example at the immutable CHIMERA v2.2 commit `fd35f897e36850d865c760c13325c5a46ebe94db`.

**Status: source inspection complete; benchmark not executed by this repository.** The notebook contains saved outputs from an upstream run, but those outputs are not a substitute for reproducing that run from verified inputs in our own environment.

Reference:
- Repository: https://github.com/CosmoStatGW/CHIMERA
- Pinned notebook: https://github.com/CosmoStatGW/CHIMERA/blob/fd35f897e36850d865c760c13325c5a46ebe94db/examples/incomplete_galaxy_catalog.ipynb
- Release: https://github.com/CosmoStatGW/CHIMERA/releases/tag/v2.2

## What the pinned notebook actually does

The notebook's key configuration is:

| Item | Value read from notebook | Implication |
|---|---|---|
| Mock identifier | `mock_id = 1` | All mock-specific products must refer to mock 1. |
| Magnitude cut | `mcut = 10` | Completeness/interpolation products must use this cut. |
| GW events | 100-event pixelated mock | This is not the 300-event SNR-20 pixelated fixture listed in the v2.0 Zenodo record. |
| Detection threshold | `snr_cut=25` | Selection injections are filtered at SNR 25. |
| Fiducial cosmology | `H0=70`, `Om0=0.25` | The displayed posterior includes a fiducial H0 of 70. |
| H0 evaluation grid | 50 points from 50 to 90 | This is the notebook's reported grid, not a general recommendation for inference. |
| Selection injections | `N_inj=20*1e6` | The pipeline expects a 20-million-injection selection model. |
| Completeness masks | 5 masks, nside 32, `nest=False` | Pixel ordering and mask configuration must be consistent. |
| Population | conditioned PLP mass model and Madau–Dickinson rate | A comparison must match these assumptions or explicitly label the difference. |

The notebook builds the masked completeness object, loads saved GW-mask indices, loads precomputed `p_cat`, `N_gal`, `P_compl`, and redshift grids from `pcat_incomplete.h5`, then evaluates the CHIMERA hyperlikelihood and its selection function.

## Required inputs and derived products

| Input/product | Notebook filename | State in the currently documented v2.0 fixture |
|---|---|---|
| Pixelated GW posterior samples | `O5snr25_100evs_mock1_pixelated_nsidelist512-256-128-64-32-16-8_meanpixels15_skyconf0.9_nestFalse.h5` | **Missing / not matched**: the v2.0 fixture lists a different, 300-event SNR-20 file. |
| Selection injections | `injections_Ninj_2e7_O5Like_snr20.h5` | Listed in v2.0 Zenodo fixture; must still be checksum-verified and filtered with the notebook's SNR threshold. |
| Completeness masks | `completeness_masks_5x_pix32False.h5` | Not listed in v2.0 fixture. |
| Mock-specific completeness interpolants | `completeness_masks_5x_pix32False_mock1_mcut10.h5` | Not listed in v2.0 fixture; requires matching mask, galaxy-density model, redshift range, and magnitude cut. |
| GW-to-mask index cache | `temp_gw_mask_idxs_mock1_mcut10.h5` | Not listed in v2.0 fixture; can potentially be regenerated from the pinned code and matching event pixels/masks. |
| Precomputed catalogue likelihood terms | `pcat_incomplete.h5` | Not listed in v2.0 fixture; described by the notebook as precomputed elsewhere because the galaxy catalogue is very large. |
| Underlying mock galaxy data | Not directly loaded in the shown notebook cells | Provenance, exact source path, generation recipe, and redistribution availability remain unresolved. |

The notebook contains a local/HPC directory string (`/leonardo_work/IscrC_LIGEA3/O5_B25/data/`), but the shown cells use relative `./data/` paths. The directory string alone is not a public data source or a reproducible acquisition instruction.

The saved notebook output reports very large mask galaxy counts (including tens to hundreds of millions in some masks). This is an additional reason not to assume that a small parent-catalogue fixture is the same underlying dataset used to create these products. Counts in saved output are clues for provenance investigation, not independently verified properties of an acquired file.

## Scientific comparison contract

Before comparing our model with CHIMERA, record and align at minimum:

1. The exact CHIMERA commit and dependency/environment versions.
2. The exact event list, event count, PE prior correction, pixelization parameters, and SNR threshold.
3. The exact injection set and selection-function normalization.
4. The underlying galaxy catalogue and its generation/provenance, magnitude cut, masks, and completeness interpolants.
5. The cosmology, population/rate/mass assumptions, H0 prior and evaluation grid.
6. Whether the comparison targets per-event likelihoods, selection normalization, or the final H0 posterior.
7. Numerical tolerances and an acceptance rule chosen before inspecting comparison results.

A posterior-curve resemblance alone is not sufficient. If the implementations use different inputs or assumptions, report those differences and do not call the result a matched benchmark.

## Source-level likelihood compatibility audit

A source inspection of the pinned CHIMERA v2.2 commit `fd35f897e36850d865c760c13325c5a46ebe94db` shows why the current finite-host Gaussian toy likelihood is a control model, not an end-to-end substitute for CHIMERA.

1. **Catalogue plus missing/background population.** In `CHIMERA/catalog/catalog.py`, `pixelated_catalog.p_gal` constructs
   `p_gal = fR * p_cat + (1 - P_compl) * p_bkg`.
   Here `p_cat` is a precomputed event × pixel × redshift catalogue term; `P_compl` is the completeness evaluated on the analysis grid; and `p_bkg` is a normalized background redshift density. The completeness normalization `fR` is computed from the completeness function and background density. This is not equivalent, in general, to a single fixed hidden-host redshift with a scalar missing fraction.
2. **Event numerator.** In `CHIMERA/likelihood.py`, the pixelated numerator combines the 3D GW posterior density, the CBC redshift/population density, and the cosmological Jacobian `dDL/dz * (1+z)^2`; it integrates over redshift and sums contributions over sky pixels.
3. **Detection selection.** The hyperlikelihood calls the selection model's `N_exp`. In scale-free mode it subtracts `log(N_exp)` from the event log-likelihood; in the non-scale-free mode it includes the rate/exposure factor and subtracts `N_exp / N_events`. The notebook configures an injection file and `N_inj = 20e6`, with the injection data filtered at `snr_cut=25`.
4. **Current repository model boundary.** The finite-host robustness experiments use a one-dimensional synthetic distance observation and discrete host redshifts. They do not currently ingest the notebook's pixelated GW posterior samples, reconstruct `p_cat` from the underlying galaxy catalogue, or reproduce CHIMERA's event numerator and selection normalization as a matched likelihood. A visually similar H0 curve would not close these gaps.

**Consequence:** external validation must be staged. First establish input provenance and API/environment compatibility; then reproduce one documented CHIMERA likelihood calculation on a small verified fixture; only after that should a matched end-to-end comparison be attempted. If the upstream event file and `pcat_incomplete.h5`/underlying mock-galaxy source remain unavailable, report the benchmark as blocked rather than replacing them with a different fixture.

## Next go/no-go gate

**Go** only when the event file and all essential inputs have a documented source, matching configuration, and verified checksums—or when every derived product can be deterministically regenerated from verified source data using pinned code.

**No-go for end-to-end comparison** while the 100-event SNR-25 event file and the underlying mock-galaxy data / `pcat_incomplete.h5` provenance remain unresolved. Continue with provenance tracing and lightweight source/API compatibility checks, but do not substitute the v2.0 fixture silently or describe the benchmark as run.

## Provenance boundary

This is a source-level execution contract, not an executed reproduction, independent validation of CHIMERA, or evidence that our current model is astrophysically valid. Any future benchmark report must be generated from actual run logs and verified artifacts.


## Public-release provenance recheck (2026-10-28)

A follow-up check of the official [CHIMERA GitHub releases](https://github.com/CosmoStatGW/CHIMERA/releases) and the public [CHIMERA v2.0 Zenodo record](https://zenodo.org/records/17143346) does not close the v2.2 input gap.

- The official GitHub release page identifies v2.2 as a software release at commit `fd35f89`; the release notes describe software/API changes and the new example notebooks, but do not document a downloadable archive containing the exact 100-event SNR-25 mock and the v2.2-derived incomplete-catalogue products required above.
- The public v2.0 Zenodo record lists an SNR-20 pixelated event file named `PE_O5Like_snr20_pixelated_nsidelist512-256-128-64-32-16-8_meanpixels15_skyconf0.9_nestFalse.h5`, plus an injection file and parent galaxy catalogue. This is not the v2.2 notebook's required `O5snr25_100evs_mock1_pixelated_...` file. Similar filenames or shared injection products do not establish matching event samples, selection cuts, masks, or precomputed catalogue likelihood terms.
- The search did not identify a public source for `pcat_incomplete.h5` or the exact mock-1, magnitude-cut-10 completeness products. Absence from these checked public records is not proof that the files are unavailable elsewhere; it means their source and checksums remain unverified for this project.

**Decision unchanged: NO-GO for a matched end-to-end CHIMERA v2.2 benchmark.** The v2.0 fixture may be used only for a separately named, version-matched compatibility exercise if its full assumptions are documented; it must not be presented as reproducing the v2.2 incomplete-catalogue notebook.
