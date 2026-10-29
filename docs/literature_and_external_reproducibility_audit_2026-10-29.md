# Literature and External-Reproducibility Audit

**Audit date:** 2026-10-29  
**Branch:** `research/literature-audit-2026-10-29`  
**Scope:** Publicly accessible papers, publisher pages, arXiv records, official CHIMERA repository/release/notebook, official CHIMERA documentation, and official release-policy documentation. This is a source-backed desk audit, not a reproduction of any paper or CHIMERA run.

## Executive decision

The defensible framing is a **controlled synthetic statistical-methods and calibration study** of how host-density estimation, catalogue incompleteness, weighting and selection-model misspecification interact. It is not a new real-data (H_0) measurement, not a validated replacement for CHIMERA, and not evidence that inverse-inclusion weighting is universally corrective.

The literature is now particularly close to this project's subject. The relevant work includes:
- direction-dependent completeness in the galaxy-catalogue dark-siren method (Gray, Messenger & Veitch 2022);
- a methodological guide that stresses correct hierarchical likelihoods and selection treatment (Gair et al. 2023);
- mock-catalogue studies of incorrect host weights and their biases (Hanselman et al. 2025; Alfradique, Bom & Castro, arXiv:2503.18887);
- a CHIMERA implementation of incompleteness and stellar-mass weighting (Borghi et al. 2026);
- Gaussian-process reconstruction of missing-galaxy clustering, published in *Physical Review D* in October 2026 (Tagliazucchi et al. 2026);
- clustering-aware simulation studies showing that clustering and localization can change how much catalogue completeness is needed (Kalomenopoulos et al. 2026; VanWyngarden et al., arXiv:2511.04786).

Therefore, broad novelty claims about catalogue incompleteness, host weighting, missing-host mixtures, clustering reconstruction, or selection effects are not defensible. A narrower possible contribution is a transparent, paired and independently reconstructed map of **failure regimes** under a clearly defined synthetic data-generating process, including cases where common diagnostics do or do not reveal miscalibration. Novelty remains provisional until the literature search and method comparison are broadened and the exact estimand and implementation differences are specified.

### Status summary

| Claim or deliverable | Assessment | Reason |
|---|---|---|
| Controlled synthetic-methods study | **GO** | Matches the current scope and can be made fully reproducible without external data. |
| New real-data (H_0) result | **NO-GO** | No real LVK/GWTC analysis is established by this project. |
| Broadly novel catalogue-incompleteness method | **NO-GO** | Directly related methods and empirical/simulation studies already exist. |
| Narrow failure-regime/calibration study | **CONDITIONAL GO** | Potentially useful if the factorial question, controls, uncertainty, and comparison against prior methods are precise. |
| Exact reproduction of CHIMERA v2.2 incomplete-catalogue notebook | **NO-GO for now** | The pinned notebook names six external data products; the release API shows no attached release assets, and the pinned repository tree contains the notebook but not those exact named data products. Their possible existence elsewhere is not ruled out. The notebook was not executed in this audit. |
| Small component benchmark with public data | **CONDITIONAL GO, as a different benchmark** | A separate, explicitly scoped component check may be possible; it must not be described as CHIMERA v2.2 reproduction. |
| Independent reconstruction as sole novelty claim | **NO-GO** | Strong reproducibility practice, but not by itself a new cosmological inference method. |

## 1. Literature map and annotated bibliography

### 1.1 Core statistical framework

**Gair, J. R. et al. (2023), “The Hitchhiker’s Guide to the Galaxy Catalog Approach for Dark Siren Gravitational-wave Cosmology,” *The Astronomical Journal*, 166(1), 22. DOI: [10.3847/1538-3881/acca78](https://doi.org/10.3847/1538-3881/acca78); arXiv: [2212.08694](https://arxiv.org/abs/2212.08694). Peer-reviewed journal article.**

This is a foundational methods reference for the galaxy-catalogue approach. It discusses common implementation errors, including selection effects and inconsistent Bayesian implementations, and explains why apparently biased results can arise from an incorrectly implemented likelihood rather than from the method itself. It should ground the manuscript's likelihood notation, conditioning assumptions, selection normalization, and terminology. The current project is a much smaller controlled simulator; it should not imply that its Gaussian toy likelihood represents the full likelihood in this paper.

**Gray, R., Messenger, C. & Veitch, J. (2022), “A Pixelated Approach to Galaxy Catalogue Incompleteness: Improving the Dark Siren Measurement of the Hubble Constant,” *Monthly Notices of the Royal Astronomical Society*, 512(1), 1127–1140. DOI: [10.1093/mnras/stac366](https://doi.org/10.1093/mnras/stac366); arXiv: [2111.04629](https://arxiv.org/abs/2111.04629). Peer-reviewed journal article.**

This work models direction-dependent catalogue completeness by pixelating the sky and combining pixel-level inference. It directly rules out framing catalogue incompleteness itself as a new topic. The distinction from the present repository must be methodological and explicit: this repository's current controlled synthetic estimators and misspecification contrasts are not a replacement for a sky-pixelated observational likelihood.

**Finke, A. et al. (2021), “Cosmology with LIGO/Virgo dark sirens: Hubble parameter and modified gravitational wave propagation,” *JCAP* 2021(08), 026. DOI: [10.1088/1475-7516/2021/08/026](https://doi.org/10.1088/1475-7516/2021/08/026); arXiv: [2101.12660](https://arxiv.org/abs/2101.12660).** This work applies a galaxy-catalogue framework to GWTC-2 and GLADE and studies systematic effects. It is useful for positioning catalogue-based inference in the broader dark-siren pipeline landscape, rather than presenting the repository as the first systematic audit.

### 1.2 Host weighting and misspecification

**Hanselman, A. G., Vijaykumar, A., Fishbach, M. & Holz, D. E. (2025), “Gravitational-wave dark siren cosmology systematics from galaxy weighting,” *The Astrophysical Journal*, 979, 9. DOI: [10.3847/1538-4357/ad9393](https://doi.org/10.3847/1538-4357/ad9393); arXiv: [2405.14818](https://arxiv.org/abs/2405.14818). Peer-reviewed journal article.**

Using mock galaxy catalogues, the authors find that incorrect host weighting can bias (H_0) through both a mismodelled galaxy redshift distribution and preferential weighting of incorrect candidate hosts. The impact depends on line-of-sight galaxy counts, distance uncertainty, and correlations among galaxy properties. They also discuss hierarchical inference as a diagnostic for incorrect weights. This is a particularly close conceptual overlap. A defensible distinction would require a precisely identified additional factor or failure mode, matched controls, and evidence that the repository answers a question not already answered by these experiments.

**Alfradique, V., Bom, C. R. & Castro, T. (2025), “Systematic bias in dark siren statistical methods and its impact on Hubble constant measurement,” arXiv:2503.18887.** [arXiv record](https://arxiv.org/abs/2503.18887). The source inspected is an arXiv manuscript; a peer-reviewed journal status was **not verified in this audit**. It uses end-to-end simulations and examines catalogue incompleteness, host-weight choices, event localization and distance uncertainty. It reports that averaging can obscure event-level posterior offsets. This is a direct warning against relying only on an overall mean error or a pooled coverage number. The repository should report per-regime results and finite-simulation uncertainty, and should not describe this source as peer-reviewed without rechecking its current record.

### 1.3 Incompleteness, clustering and reconstruction

**Borghi, N., Moresco, M., Tagliazucchi, M. & Cuomo, G. (2026), “Echoes from the Dark: Galaxy Catalog Incompleteness in Standard Siren Cosmology,” *Astronomy & Astrophysics*, 706, A199. arXiv: [2509.18243](https://arxiv.org/abs/2509.18243).** The arXiv full text and official CHIMERA release notes were inspected. The work implements incompleteness and host weighting in CHIMERA and studies 100 O5-like binary-black-hole events across catalogue-completeness and stellar-mass-weighting choices. It reports percent-level forecast constraints in the specified mock setup and sensitivity of the (H_0)-population-mass correlation to completeness and weighting. These are forecasts under stated simulation assumptions, not current observed (H_0) precision. This paper is a central related-work citation and strongly narrows any broad novelty claim.

**Tagliazucchi, M., Gair, J., Barbieri, R. & Moresco, M. (2026), “Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs,” *Physical Review D* 114, 083523. DOI: [10.1103/vmzf-1wxq](https://doi.org/10.1103/vmzf-1wxq); arXiv: [2608.24992](https://arxiv.org/abs/2608.24992). Published 13 October 2026.**

This recent paper proposes Gaussian-process reconstruction of the line-of-sight galaxy redshift distribution while accounting for survey selection, preserving clustering information that homogeneous completion discards. It reports improved forecast precision over homogeneous completion in its mock scenarios. The current repository does not implement this GP method, so it must not imply that its KDE estimator reconstructs missing-galaxy clustering in the same way. The paper also means that a generic proposal to “try GP reconstruction for incomplete catalogues” is no longer a distinctive next step; a more defensible project would compare failure regimes or calibration under matched assumptions.

**Kalomenopoulos, M. et al. (2026), “Clustering effects on the dark siren determination of (H_0): a simulation study,” *Monthly Notices of the Royal Astronomical Society*, 552(1), stag1662. DOI: [10.1093/mnras/stag1662](https://doi.org/10.1093/mnras/stag1662). Published 3 September 2026.**

This work uses high-resolution cosmological simulations with clustering and studies completeness, localization and event counts. It reports that clustered catalogues can remain informative at lower completeness than a uniform catalogue in its setup. This is relevant because catalogue completeness cannot be treated as a single scalar independently of the spatial structure and localization regime. The repository's synthetic setup should clearly state whether clustering is absent, simplified, or conditioned on, and should not generalize scalar-completeness findings to realistic surveys.

**VanWyngarden, M. et al. (2025), “How Low Can You Go: Constraining the Effects of Catalog Incompleteness on Dark Siren Cosmology,” arXiv:2511.04786.** [arXiv record](https://arxiv.org/abs/2511.04786). A preprint record was inspected; journal status was not verified here. It studies bright-galaxy-only catalogues and shows in its mocks that clustering around bright galaxies can preserve useful information for well-localized events, while unclustered catalogues or broad localization can produce bias. It supports stratifying conclusions by clustering and localization, rather than treating incompleteness fraction alone as the causal factor.

### 1.4 Calibration and computation

**Talts, S., Betancourt, M., Simpson, D., Vehtari, A. & Gelman, A. (2020), “Validating Bayesian Inference Algorithms with Simulation-Based Calibration.” arXiv: [1804.06788](https://arxiv.org/abs/1804.06788).**

This is a foundational reference for simulation-based calibration (SBC). SBC assesses whether an implementation's posterior computations are consistent with the generative model under prior-predictive simulations. It does not establish that the generative model matches reality. The manuscript must distinguish algorithmic SBC from fixed-truth frequentist coverage, posterior interval coverage at selected truths, and robustness under misspecification. Passing SBC under a self-consistent simulator is not proof of scientific adequacy under a wrong model.

## 2. Novelty assessment

### Already established in the literature

1. Galaxy-catalogue dark-siren inference and its hierarchical likelihood are established.
2. Catalogue incompleteness and direction-dependent completeness modelling are established.
3. Host-galaxy weighting and bias from incorrect weighting are established.
4. Selection effects and the need for consistent normalization are established.
5. Homogeneous missing-galaxy completion, clustering-aware approaches, and GP reconstruction are established.
6. Simulation-based calibration is established as a computational validation technique.

### Potentially useful but not yet demonstrated as novel

A reproducible comparison of **joint** host-density and selection-function misspecification, with paired experiments, oracle controls, calibration diagnostics, and independent reconstruction could be useful. However, “joint misspecification” is not automatically novel simply because a repository combines several existing ingredients. To make a contribution, the paper must state:
- the exact estimand and data-generating model;
- which factor(s) prior studies did not vary jointly or did not diagnose;
- the conditions under which a method succeeds or fails;
- why the result is not a trivial consequence of the toy model;
- how Monte Carlo uncertainty and multiple comparisons are handled;
- which results are exploratory versus confirmatory.

Independent artifact reconstruction is a meaningful reproducibility safeguard. It becomes a stronger research-software contribution if the study systematically injects known faults (e.g., wrong selection normalization, altered weights, corrupted summaries, or a grid too coarse to resolve a narrow posterior) and reports which audit layers detect each fault. Reconstruction alone does not establish a new astrophysical inference method.

### Claims to avoid

- “First,” “novel method,” or “state of the art” without a systematic comparison and explicit scope.
- “Inverse-inclusion weighting corrects incompleteness” as a universal statement.
- “The method is calibrated” without specifying the calibration experiment and its assumptions.
- “The model is validated” merely because unit tests or CI pass.
- “CHIMERA validated/reproduced” without running the exact pinned software, inputs, and configuration.
- Any suggestion that synthetic (H_0) results are measurements of the real Universe.

### Recommended framing

> A reproducible synthetic study mapping the calibration and bias failure regimes of KDE-based host-density estimation and inverse-inclusion weighting under explicitly controlled catalogue-selection and host-population misspecification.

Treat this as a working framing, not a priority claim.

## 3. Three ranked next research questions

### Rank 1 — Diagnostic failure map under joint misspecification (best near-term fit)

**Question:** Under which combinations of selection-function misspecification, host-density error, KDE bandwidth and localization/noise regime do posterior summaries become materially biased or miscalibrated, and which diagnostics detect those failures?

**Method:** Prespecify a small factorial design; separate calibration from evaluation simulations; compare a known-density oracle, the current KDE method, and a clearly defined weighting baseline; use paired random seeds where scientifically valid; record bias, RMSE, interval coverage, interval width, SBC ranks where appropriate, and Monte Carlo uncertainty.

**Controls:** matched data-generating and inference models; correct-selection negative control; deliberate selection-normalization defect; truth strata including prior boundaries and interior; independent raw-artifact reconstruction; predeclared thresholds.

**Useful result even if the proposed correction fails:** a map of conditions where common diagnostics miss a known failure, or where a seemingly improved RMSE comes at the cost of undercoverage or dependence on boundary placement.

**Feasibility:** High; synthetic-only and independently testable. **Novelty status:** plausible but unproven; compare against Hanselman et al. and Alfradique et al. before claiming a gap.

### Rank 2 — Does clustering change the ranking of completion/weighting methods?

**Question:** Does the ranking of homogeneous completion, KDE-based estimation and weighting change when the same completeness fraction is applied to clustered versus unclustered host populations?

**Method:** Matched clustered and unclustered synthetic catalogues with identical marginal redshift distributions and matched event/noise realizations. Vary completeness and localization separately.

**Controls:** complete-catalogue oracle; matched marginal distributions; identical selection normalization; repeated catalogues rather than a single fixed realization; external method comparisons.

**Useful negative result:** show that apparent improvements from weighting disappear or reverse once clustering and localization are matched.

**Feasibility:** Moderate; current simulator may need substantial extension and the question overlaps strongly with recent clustering literature. **Novelty status:** only if a specific unaddressed interaction or diagnostic is identified.

### Rank 3 — Fault-injection benchmark for independent reconstruction

**Question:** Which classes of scientific-pipeline error are detected by unit tests, SBC, raw-data reconstruction, convergence tests, and checksum/manifest audits—and which escape all of them?

**Method:** Inject controlled defects into a toy pipeline and build a detection matrix across audit layers.

**Controls:** known-good implementation, seeded defects, blinded audit labels if feasible, fixed artifacts, deterministic reruns, and explicit false-positive/false-negative accounting.

**Useful result:** a reproducible map of audit coverage and blind spots, even if no single audit detects every defect.

**Feasibility:** High. It is more naturally a computational-reproducibility/software-methods contribution than a cosmology-method paper. **Novelty status:** requires comparison with existing scientific-software testing literature.

## 4. CHIMERA v2.2 external-benchmark provenance

### Exact version inspected

- Repository: [CosmoStatGW/CHIMERA](https://github.com/CosmoStatGW/CHIMERA)
- Release: [v2.2](https://github.com/CosmoStatGW/CHIMERA/releases/tag/v2.2)
- Pinned notebook: [`examples/incomplete_galaxy_catalog.ipynb` at `fd35f897e36850d865c760c13325c5a46ebe94db`](https://github.com/CosmoStatGW/CHIMERA/blob/fd35f897e36850d865c760c13325c5a46ebe94db/examples/incomplete_galaxy_catalog.ipynb)
- Notebook Git blob SHA reported by GitHub: `20934650723fdd35feabca50a5fc9d927edde62a`
- Release API: tag `v2.2`; the GitHub release's `assets` array was empty when inspected on 2026-10-29. GitHub release page notes the v2.2 API restructuring and lists the incomplete-catalogue notebook as a new example.
- CHIMERA repository licence shown on the official repository page: GPL-3.0. Check the exact pinned licence before redistribution or combining code.

### Data products explicitly referenced by the pinned notebook

The notebook code sets a local `./data/` directory and names these inputs:

1. `O5snr25_100evs_mock1_pixelated_nsidelist512-256-128-64-32-16-8_meanpixels15_skyconf0.9_nestFalse.h5`
2. `injections_Ninj_2e7_O5Like_snr20.h5`
3. `completeness_masks_5x_pix32False.h5`
4. `completeness_masks_5x_pix32False_mock1_mcut10.h5`
5. `temp_gw_mask_idxs_mock1_mcut10.h5`
6. `pcat_incomplete.h5`

It also sets an HPC-specific path, `/leonardo_work/IscrC_LIGEA3/O5_B25/data/`, and comments that the precomputed catalogue term is generated elsewhere because of the very large galaxy catalogue. The notebook calls the selection function with `N_inj=20*1e6`, so the injection file's provenance and exact compatibility matter. The pinned repository tree contains the notebook but did not contain any of the six exact named data products at the paths searched. The release's attached-asset list was empty. This is evidence that the required inputs are **not included in the inspected tracked tree or attached release assets**, not proof that they are unavailable elsewhere.

### Data-status classification

| Item | Status as of audit | What is established |
|---|---|---|
| Pinned notebook | **Verified publicly available** | Raw notebook was retrieved from the exact commit. |
| Six named HDF5 inputs | **Not located in inspected pinned tree/release assets** | Exact names occur in the notebook; those files were not found in the inspected tree and the release has no attached assets. |
| Official external archive/DOI for all six exact inputs | **Not verified** | No authoritative archive linking all exact files and checksums was verified in this audit. |
| Exact generation scripts for all six inputs | **Not verified** | Notebook comments identify some precomputation but do not by themselves establish a public end-to-end generation path. |
| Software environment and compute requirements | **Partially documented** | The notebook sets `CHIMERA_ENABLE_GPU=False`, imports the v2.2 API, and requests 20 million injections; a full execution environment and data acquisition route were not established. |
| Notebook execution | **Not performed** | No claim of successful execution is made. |

**Decision:** exact CHIMERA v2.2 reproduction remains **NO-GO until every required file, checksum/version, licence/access condition and generation path is verified**. Do not substitute v2.0 or a different mock and label it v2.2 reproduction. If a smaller official tutorial is run successfully with its own available inputs, report it as a separate software smoke test, not as reproduction of the incomplete-catalogue benchmark.

The official [CHIMERA Catalog API documentation](https://chimera-gw.readthedocs.io/latest/api/catalog/) confirms that catalogue calculations use galaxy-property files, pixelated GW data, a completeness object, cosmology, interpolation products and host weights. This supports the conclusion that input provenance is a scientific requirement, not merely a download inconvenience.

## 5. Alternative external benchmarks

| Candidate | What it can test | Prerequisites / cost | Permissible claim |
|---|---|---|---|
| CHIMERA v2.2 incomplete-catalogue notebook | End-to-end compatibility with that exact workflow | All six named files, version-matched software, data rights, substantial storage/compute and exact configuration | Exact reproduction only after matching inputs, versions and outputs; currently **NO-GO**. |
| CHIMERA complete-catalogue example with its own documented inputs | API integration and one complete-catalogue path | Verify exact input availability and version; likely easier than incomplete-catalogue example | “Complete-catalogue software smoke test,” not incomplete-catalogue validation. |
| Public paper/mock outputs from a related study | A narrow numerical cross-check if the published artifact includes sufficient inputs and metadata | Must verify source, licence, definitions and matching estimand | “Comparison against published mock results” only if independently reconstructed; otherwise contextual comparison. |
| Public GW event products plus a galaxy catalogue | Data-ingestion and limited component demonstrations | Need selection-function treatment, catalogue completeness characterization, event-product provenance, licences and appropriate likelihood | Component demonstration only; not a valid (H_0) result merely from cross-matching. |
| Synthetic clustered/unclustered catalogue pair | Controlled causal comparison of clustering, completeness and localization | Build and validate the generator; modest-to-high compute depending on repeated event simulations | A synthetic-methods result under the stated generator, not external validation. |

## 6. Release and publication requirements

- **MNRAS official instructions:** require a Data Availability Statement and encourage access to data and software underlying conclusions; support full data/software citations. The journal also states that AI use in drafting, code, analysis or translation should be disclosed in the cover letter and manuscript Methods or Acknowledgements, and that authors remain responsible for the content. See [MNRAS Instructions to Authors](https://academic.oup.com/mnras/pages/General_Instructions).
- **arXiv official submission guidance:** submissions should be topical, refereeable scientific contributions and are subject to moderation. First submissions or new categories may require endorsement; endorsement is not peer review. See [arXiv submission guidelines](https://info.arxiv.org/help/submit/index.html) and [endorsement policy](https://info.arxiv.org/help/endorsement.html).
- **Zenodo official GitHub integration:** when a repository is connected and configured, new GitHub releases can be ingested and archived. This creates a versioned archival record; it does not independently certify correctness. See [Zenodo: enable a repository](https://help.zenodo.org/docs/github/enable-repository/).
- **CHIMERA licence:** the official repository page identifies GPL-3.0. If code is copied, adapted, bundled or redistributed, review GPL obligations and dependencies; do not assume the repository's own licence can override upstream terms. Cite CHIMERA and its version.
- **Third-party data:** do not redistribute HDF5 inputs, mock catalogues or GW/galaxy data unless their licence and access terms permit it. Prefer links and checksums to authorized sources when redistribution is not allowed.
- **Authorship and AI assistance:** record actual human contributions, code review, analysis, and writing. AI systems should not be listed as authors. Disclose AI assistance in line with the chosen venue's policy and ensure a human verifies every scientific claim and citation.

No journal is recommended as guaranteed suitable, and no acceptance prediction is made.

## 7. Prioritized release roadmap

### Priority 1 — Prevent incorrect claims

1. Replace any broad “novel method” or “first” wording with a precise research question and a comparison against the closest papers above.
2. Label every headline result as synthetic and state the data-generating assumptions beside it.
3. Distinguish algorithmic SBC, fixed-truth coverage, posterior interval width, and robustness under misspecification.
4. Separate oracle information from information available to a real analysis.
5. Report truth-boundary and interior results separately where the prior support makes boundary coverage special.
6. Do not claim CHIMERA v2.2 reproduction until all six inputs and the execution path are verified.

**Evidence/output:** revised claims table and manuscript language; checked citations; explicit GO/NO-GO list.

### Priority 2 — Make the synthetic result independently reproducible

1. Pin a commit, dependency lockfile and Python versions.
2. Store configurations, seeds, raw per-replicate output, and machine-readable summary tables.
3. Include SHA-256 checksums for all generated artifacts and a script that rebuilds each headline figure/table from raw data.
4. Use an independent reconstruction path that does not import the experiment's own summary functions.
5. Add known-answer and fault-injection tests for likelihood normalization, weight handling, prior boundaries, duplicated/missing rows, and coarse-grid posteriors.
6. Record Monte Carlo uncertainty for coverage, paired contrasts and bootstrap intervals; preserve failed or exploratory runs with clear labels.
7. Run multi-resolution convergence checks and retain diagnostics even when a run fails a convergence threshold.

**Evidence/output:** clean-environment reproduction log, raw-to-summary audit, checksum manifest, and release checklist.

### Priority 3 — External benchmark

1. Retrieve or obtain an authoritative data manifest for all six CHIMERA inputs.
2. Verify exact SHA-256 checksums, size, source URL/DOI, licence and access conditions.
3. Confirm the notebook's commit, package version, configuration, Python/JAX environment and injection count.
4. Determine whether each input can be obtained legally and whether any generation scripts are public.
5. Run only after the data and compute path is clear; preserve the complete log and output manifest.
6. If any input is not available, record the exact blocker and either postpone the benchmark or choose a different benchmark with a different claim.

**Evidence/output:** per-file provenance table, environment lock, run log, independent result comparison. No contacting maintainers or accepting terms is done by this audit.

### Priority 4 — Manuscript and release polish

1. Add a concise methods diagram and equations mapped to implementation functions.
2. Include a limitation table: no real LVK analysis, simplified catalogue structure, no validated clustering reconstruction unless actually implemented, and sensitivity to prior/noise assumptions.
3. Provide a Data Availability Statement and software citation.
4. Tag a release only after the final tests pass; archive that release through Zenodo if enabled.
5. Include a contribution statement and AI-assistance disclosure consistent with the chosen venue.
6. Provide all figures as scripts plus data, rather than relying on manually edited outputs.

## 8. Source register and search limitations

This audit used primary/publisher sources, arXiv records, official CHIMERA GitHub content and documentation, and official policy pages. It did not perform an exhaustive systematic review of every dark-siren, survey-completeness, inverse-probability-weighting or scientific-software paper. The literature matrix is therefore a **verified core set, not an exhaustive bibliography**.

Specific limits:
- The exact CHIMERA notebook was inspected, but it was not executed.
- An empty release-asset list and absence from the inspected Git tree do not establish that files are unavailable from all external archives or by an official access process.
- The journal status of arXiv:2503.18887 and arXiv:2511.04786 was not verified in the sources inspected here.
- No independent reproduction of the papers' headline numerical results was performed.
- The recent GP and clustering papers materially update the novelty assessment; citations should be rechecked immediately before manuscript submission.

## Official and primary-source links

1. [Gair et al. 2023, arXiv:2212.08694](https://arxiv.org/abs/2212.08694) — [DOI](https://doi.org/10.3847/1538-3881/acca78)
2. [Gray, Messenger & Veitch 2022, arXiv:2111.04629](https://arxiv.org/abs/2111.04629) — [DOI](https://doi.org/10.1093/mnras/stac366)
3. [Finke et al. 2021, arXiv:2101.12660](https://arxiv.org/abs/2101.12660) — [DOI](https://doi.org/10.1088/1475-7516/2021/08/026)
4. [Hanselman et al. 2025, arXiv:2405.14818](https://arxiv.org/abs/2405.14818) — [DOI](https://doi.org/10.3847/1538-4357/ad9393)
5. [Alfradique, Bom & Castro, arXiv:2503.18887](https://arxiv.org/abs/2503.18887)
6. [Borghi et al., arXiv:2509.18243](https://arxiv.org/abs/2509.18243)
7. [Tagliazucchi et al. 2026, arXiv:2608.24992](https://arxiv.org/abs/2608.24992) — [DOI](https://doi.org/10.1103/vmzf-1wxq)
8. [Kalomenopoulos et al. 2026, MNRAS DOI](https://doi.org/10.1093/mnras/stag1662)
9. [VanWyngarden et al., arXiv:2511.04786](https://arxiv.org/abs/2511.04786)
10. [Talts et al., arXiv:1804.06788](https://arxiv.org/abs/1804.06788)
11. [CHIMERA v2.2 release](https://github.com/CosmoStatGW/CHIMERA/releases/tag/v2.2)
12. [Pinned CHIMERA incomplete-catalogue notebook](https://github.com/CosmoStatGW/CHIMERA/blob/fd35f897e36850d865c760c13325c5a46ebe94db/examples/incomplete_galaxy_catalog.ipynb)
13. [CHIMERA Catalog API documentation](https://chimera-gw.readthedocs.io/latest/api/catalog/)
14. [MNRAS author instructions](https://academic.oup.com/mnras/pages/General_Instructions)
15. [arXiv submission guidance](https://info.arxiv.org/help/submit/index.html) and [endorsement policy](https://info.arxiv.org/help/endorsement.html)
16. [Zenodo GitHub integration](https://help.zenodo.org/docs/github/enable-repository/)
