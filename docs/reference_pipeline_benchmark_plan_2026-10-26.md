# Published-pipeline benchmark: feasibility and go/no-go plan (2026-10-26)

## Decision

The next scientifically useful milestone is a benchmark against an established completeness-aware pipeline. Do not present the current toy-model experiments as a substitute for this step.

Two realistic reference options were checked:

1. **CHIMERA** — public source code and tutorials cover complete and incomplete galaxy catalogues. The required O5-like mock inputs are available from the public [CHIMERA v2.0 Zenodo record](https://zenodo.org/records/17143346), which lists 662.3 MB of files. Relevant listed files include the 156.1 MB 20-million-injection HDF5 file, a 133.8 MB parent-galaxy catalogue, and 180.1/192.2 MB GW posterior-sample files (the latter is the pixelated version used by the incomplete-catalogue notebook). This makes CHIMERA a plausible benchmark candidate, not an impossible one. However, the incomplete-catalogue notebook also references completeness-mask/interpolation products and temporary mask indices; their generation and provenance must be traced before the example is reproducible end-to-end.
2. **gwcosmo** — mature public dark-siren pipeline with a documented command-line interface. Its dark-siren method requires GW posterior samples, sky maps, a selection-effect injection set, and a precomputed line-of-sight redshift-prior HDF5 file. The official documentation notes that injections and LOS priors are not bundled with the package by default. A versioning issue must also be resolved: the official documentation inspected here identifies itself as 3.0.1.dev2, while PyPI lists a 3.1.0 release uploaded on 22 July 2026. The benchmark must pin one exact release/commit and confirm that the data/API instructions match that version. This creates a data-provenance and compute prerequisite, not merely a Python dependency.

Sources:
- [CHIMERA public repository](https://github.com/CosmoStatGW/CHIMERA)
- [CHIMERA incomplete-galaxy-catalogue example](https://github.com/CosmoStatGW/CHIMERA/blob/main/examples/incomplete_galaxy_catalog.ipynb)
- [CHIMERA v2.0 mock-data record (Zenodo 17143346)](https://zenodo.org/records/17143346)
- [CHIMERA documentation](https://chimera-gw.readthedocs.io/latest/)
- [CHIMERA completeness implementation (mask/interpolant generation and loading)](https://github.com/CosmoStatGW/CHIMERA/blob/main/CHIMERA/catalog/completeness.py)
- [gwcosmo data requirements](https://lscsoft.docs.ligo.org/gwcosmo/data.html)
- [gwcosmo inference command and input requirements](https://lscsoft.docs.ligo.org/gwcosmo/cosmological-inference.html)
- [gwcosmo installation and source](https://lscsoft.docs.ligo.org/gwcosmo/installation.html)
- [gwcosmo PyPI release history](https://pypi.org/project/gwcosmo/)
- [gwcosmo pipeline validation and systematics paper (2026)](https://doi.org/10.1103/53c5-mtk5)

## Candidate CHIMERA v2.0 input manifest

The Zenodo record lists the following upstream checksums (MD5) and sizes. MD5 is recorded here for source verification; after downloading, compute and store SHA-256 as the repository's stronger local provenance identifier.

| File | Size | Zenodo-listed MD5 | Role in candidate benchmark |
|---|---:|---|---|
| [CHIMERA v2.0 software archive](https://zenodo.org/records/17143872) | 9.7 MB | 6b11c8659c0627dcfcf85f35a2473fb2 | Versioned upstream code (separate Zenodo software record 17143872) |
| [injections_Ninj_2e7_O5Like_snr20.h5](https://zenodo.org/records/17143346/files/injections_Ninj_2e7_O5Like_snr20.h5) | 156.1 MB | 213ab2655b9fac5fb8cdc022fe8bd555 | Selection-effect injections |
| [parent_galaxy_catalog_1.6M_dVdz.h5](https://zenodo.org/records/17143346/files/parent_galaxy_catalog_1.6M_dVdz.h5) | 133.8 MB | c59d9feece4e0deab2b2b544cdb86303 | Parent galaxy sample |
| [PE_O5Like_snr20.h5](https://zenodo.org/records/17143346/files/PE_O5Like_snr20.h5) | 180.1 MB | f887599e7e11b19ce21039ecb51616c7 | Non-pixelated event posterior samples |
| [PE_O5Like_snr20_pixelated_nsidelist512-256-128-64-32-16-8_meanpixels15_skyconf0.9_nestFalse.h5](https://zenodo.org/records/17143346/files/PE_O5Like_snr20_pixelated_nsidelist512-256-128-64-32-16-8_meanpixels15_skyconf0.9_nestFalse.h5) | 192.2 MB | 165914665a87ebb971fad2ed3bc66f15 | Pixelated event samples used by the incomplete-catalogue example |

The three files directly referenced by the incomplete-catalogue example total approximately 482.1 MB (injections, parent galaxy catalogue, pixelated posterior samples), excluding generated completeness masks, interpolation products, temporary mask indices, and software dependencies. These data are publicly listed in [Zenodo record 17143346](https://zenodo.org/records/17143346); the versioned v2.0 code is in the separate [CHIMERA software release record 17143872](https://zenodo.org/records/17143872). Do not mistake the small 63.5 kB archive listed alongside the data for the full software source. An opt-in helper, [prepare_chimera_v2_benchmark_data.py](../scripts/prepare_chimera_v2_benchmark_data.py), streams the three-file fixture, checks the Zenodo-listed MD5 before publishing each file, computes a local SHA-256, and writes a provenance manifest. Its test suite uses mocked responses and never downloads the large files.

For a local feasibility check, first run:
```bash
python scripts/prepare_chimera_v2_benchmark_data.py --dry-run
```
After reviewing the upstream licence and confirming sufficient storage, run:
```bash
python scripts/prepare_chimera_v2_benchmark_data.py --output-dir external_data/chimera_v2
```
Do not commit these large files or add them to ordinary CI. Verify the CHIMERA v2.0 source archive, generate/recover the mask/interpolant products, and first run a small upstream-supported configuration before attempting the full example.

## Feasibility comparison

| Requirement | CHIMERA | gwcosmo |
|---|---|---|
| Published incompleteness-aware method | Yes | Yes |
| Public source and documentation | Yes | Yes |
| Small self-contained incomplete-catalogue fixture confirmed in source reviewed here | No; notebook references external O5 mock/injection/mask files | No; LOS-prior HDF5 and injections are external to the package by default |
| Selection-effect modelling | Explicit population and selection-function components | Explicit injections and selection-effect requirements |
| Direct drop-in replacement for this repository's toy likelihood | No | No |
| Immediate low-cost CI benchmark justified from current evidence | No | No |

The current evidence supports CHIMERA as the first candidate for a bounded benchmark because its mock events and injections have a public Zenodo record. Its completeness implementation supports computing masks and completeness interpolants from scratch when the corresponding file paths are set to None; however, the example notebook instead loads precomputed mask/interpolant files and a temporary GW-mask-index file. A reproduction must either recover those exact derived files or explicitly run and validate the upstream preprocessing path from the parent galaxy catalogue and pixelated GW samples. gwcosmo remains a second candidate; its current package release and the older documentation version inspected must be reconciled before a pinned run. Neither pipeline has yet been executed in this repository.

## Staged benchmark plan

### Gate A — freeze a reference target

Choose one specific upstream release or commit, cite the corresponding method paper, and record its licence, environment, and the exact inference mode. Do not track a moving main branch as the scientific reference.

### Gate B — secure a small, redistributable benchmark fixture

The fixture must include or point to stable public files for:
- GW posterior samples and sky map (or the exact upstream mock-event equivalents);
- the line-of-sight catalogue/redshift prior;
- selection-effect injections compatible with the selected detection threshold;
- galaxy catalogue/completeness metadata and any mask files;
- the data-generation truth and configuration.

Record upstream URLs, licences, file sizes, SHA-256 hashes, and any data that cannot be mirrored in the repository. A benchmark is not reproducible if its crucial inputs are implicit paths on a lab cluster.

### Gate C — run the unmodified upstream baseline

Run the reference pipeline without local changes. Save the exact command/configuration, environment versions, logs, output posterior, runtime, and all input hashes. First reproduce the upstream example or its published validation before changing a model assumption.

### Gate D — define a fair comparison

The current simulator uses a low-redshift, discrete-host Gaussian distance likelihood with fixed host redshifts and no realistic sky localization or detector-selection model. It is not mathematically equivalent to either full pipeline. Comparisons must therefore use a clearly shared estimand and a matched generative process. Do not compare absolute H0 posterior widths from incompatible models and call the difference a method improvement.

A defensible bridge experiment would be:
1. use an upstream-supported synthetic/mock configuration with known truth;
2. verify the unmodified pipeline recovers the truth under its correctly specified assumptions;
3. introduce one controlled catalogue-completeness or host-weight perturbation while holding the event realization fixed;
4. compare bias, RMSE, interval coverage, and posterior width with paired Monte Carlo uncertainty;
5. independently validate the data-generation protocol and output summaries.

### Gate E — decide whether integration belongs in this repository

Only add a CI job if the fixture is legally redistributable, small enough for routine CI, and reproducible on supported Python/platform versions. Otherwise, maintain a documented external benchmark script with pinned versions and immutable inputs, and keep the existing lightweight toy-model tests as unit/regression controls.

## Go/no-go criteria

**Go** when the upstream reference can be installed or run in a pinned environment, the input data have traceable provenance and checksums, and the baseline result can be independently reproduced.

**No-go for scientific claims** if the required LOS prior, injections, catalogue masks, or event files are unavailable or unverified; if the baseline cannot be reproduced; or if the toy likelihood and reference likelihood do not share a valid comparison estimand.

## Current status

This is a feasibility review, not a completed upstream benchmark. No CHIMERA or gwcosmo analysis has been run as part of this repository, and no real gravitational-wave/catalogue result is claimed. The next action should be selecting and verifying the exact small upstream fixture, not adding another toy-model experiment merely to increase the number of simulations.
