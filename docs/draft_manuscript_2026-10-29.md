# Selection-Function and Catalogue-Incompleteness Audits for Dark-Siren H0 Inference: Controlled Synthetic Benchmarks and Reproducible Validation

**Draft manuscript v0.1 — 29 October 2026**  
**Status:** internal research draft; not submitted or peer reviewed.

## Abstract

Dark-siren cosmology infers cosmological parameters from gravitational-wave distance information without requiring a uniquely identified electromagnetic counterpart, and its performance can depend on how the potential host-galaxy population and observational selection are represented. We present two controlled synthetic software-and-methodology audits designed to make distinct failure modes observable and reproducible. The first compares six host-density treatments in a simplified low-redshift model with a deliberately incomplete galaxy catalogue. A 500-replicate paired experiment finds that ignoring the declared catalogue selection produces a mean posterior-mean bias of −8.85 km s⁻¹ Mpc⁻¹ and nominal 68% interval coverage of 9.6%, whereas inverse-inclusion weighting under the exactly known toy selection produces bias +1.02 km s⁻¹ Mpc⁻¹ and coverage 69.6%. The prespecified paired reduction in absolute posterior-mean error is 5.95 km s⁻¹ Mpc⁻¹ (percentile-bootstrap 95% interval 5.54–6.36). The second experiment introduces an H0-dependent distance-noise scale and a selection normalization computed from that same noise model. Across 1,200 synthetic datasets, an independent audit replayed each random seed and reconstructed all saved distances and posterior summaries to floating-point precision; the prior-predictive posterior-CDF rank mean was 0.514 and its descriptive decile-bin chi-square p-value was 0.450. Fixed-truth runs display conditional shrinkage toward the prior center. Both experiments are synthetic, differ in their likelihood assumptions, and test internal behavior rather than detector realism. They do not constitute a gravitational-wave catalogue analysis or a measurement of the real Hubble constant.

**Keywords:** dark sirens; Hubble constant; catalogue incompleteness; selection function; inverse-inclusion weighting; posterior calibration; reproducible simulation.

## 1. Introduction

Gravitational-wave standard sirens provide luminosity-distance information that can be combined with redshift information to constrain the cosmic expansion rate. For dark sirens, where no unique electromagnetic host is identified, an inference procedure can marginalize over a population of candidate host galaxies. The resulting likelihood is affected by the representation of that population and by how the detection process and catalogue selection are modeled.

Catalogue incompleteness and selection-aware host-population reconstruction are established topics, not new problems introduced by this work. Recent research has considered Gaussian-process reconstruction of incomplete galaxy catalogues, robust catalogue-completeness tests in practical dark-siren pipelines, and inference for populations of gravitational-wave-hosting galaxies. Our goal is narrower: to construct controlled simulations where individual statistical and numerical assumptions can be changed independently, and where the resulting artifacts can be re-derived from their seeds and raw records.

This study reports two related but technically distinct audits. Study A tests host-density estimation when the catalogue inclusion law is known or deliberately misspecified. Study B tests a separate discrete-host likelihood whose distance-noise scale depends on H0 and whose detection efficiency uses that same heteroscedastic noise model. Study A is a repeated-sampling comparison of host-density arms under a fixed-noise toy model. Study B is a simulation-based calibration check of a new likelihood component. Results from these studies are not pooled because their generative and likelihood assumptions differ.

The contribution is primarily methodological and reproducibility-focused. We report a corrected physical bandwidth convention for a reflected weighted KDE, paired estimands, raw replicate tables, checksum manifests, and independent auditors that reconstruct derived quantities rather than simply comparing an aggregate with its own saved summary.

## 2. Scope and scientific claims

The work is an exploratory synthetic-methodology study. All quantities labelled H0 are parameters within the defined toy models. No LVK event data, GWTC host catalogue, or real detector-response likelihood is analyzed here. The study does not resolve the Hubble tension, estimate a new real-world H0, or demonstrate that inverse-probability weighting is universally reliable.

The two experiments have different purposes:

- **Study A:** quantify the effect of a known toy catalogue-inclusion mechanism and a small set of assumed-selection misspecifications on host-density and posterior summaries.
- **Study B:** verify that event likelihoods, H0-dependent measurement uncertainty, the detection-efficiency normalization, synthetic data generation, and posterior summaries are numerically consistent within a discrete-host toy model.

The scope boundary is important. Software agreement and calibration in a deliberately specified synthetic model are necessary checks, but they do not establish physical realism or validity for real gravitational-wave data.

## 3. Study A: host-density and catalogue-selection benchmark

### 3.1 Generative setup

The primary Stage B experiment uses 500 paired replicate datasets. Each replicate begins with 4,000 parent galaxies and yields 30 detected events. The true parameter is H0 = 70 km s⁻¹ Mpc⁻¹. The inference grid contains 81 equally spaced H0 values over [50, 90] km s⁻¹ Mpc⁻¹, and the redshift grid contains 401 points over [0.01, 0.40]. The toy distance model uses a fixed Gaussian noise scale of 140 Mpc and a detection threshold of 250 Mpc.

The catalogue inclusion probability follows a declared logistic transition with floor 0.1, true midpoint 0.22, and true width 0.05. Catalogue selection is applied to the simulated parent population. The event-host density includes the declared merger-weight component. The study therefore treats the toy parent/host law and the inclusion law as known generative quantities, even though some inference arms deliberately receive only the incomplete catalogue or an incorrect inclusion curve.

### 3.2 Compared arms

Six arms are paired within every replicate:

1. **Analytic truth:** direct use of the declared normalized target density on the redshift grid; an inaccessible diagnostic control, not an estimator available in a real analysis.
2. **Oracle KDE:** reflected KDE estimated from the complete simulated parent population.
3. **Catalogue only:** reflected KDE estimated from the incomplete catalogue without inclusion correction.
4. **Correct weighting:** inverse-inclusion weighting with the exactly known toy inclusion probabilities.
5. **Shifted transition:** the inclusion midpoint supplied to inference is 0.27 rather than the true 0.22.
6. **Shallow transition:** the supplied transition width is 0.075 rather than the true 0.05.

The corrected estimator uses weighted Scott bandwidth by default, h = weighted SD × neff^(−1/5). An explicit factor f means physical bandwidth h = f × weighted SD of the original sample. This differs from passing a value directly as SciPy’s gaussian_kde covariance-factor setting, especially after reflected-sample augmentation. A direct kernel-sum regression test was added to guard this distinction.

### 3.3 Primary estimand and uncertainty

For replicate r, define

Dᵣ = |posterior mean H0 from catalogue-onlyᵣ − H0,true| − |posterior mean H0 from correct-weightingᵣ − H0,true|.

Positive D means correct weighting reduced absolute posterior-mean error. This is the preregistered primary paired contrast. Its uncertainty is reported using a 10,000-resample percentile bootstrap; the paired standard error and normal-approximation interval are diagnostics. Per-arm coverage uncertainty is summarized with Wilson 95% intervals. The raw records and saved posterior grids are checked by a separate auditor.

Stage A used 200 paired replicates to estimate the primary contrast’s standard deviation. Its observed mean was 6.536 km s⁻¹ Mpc⁻¹ and sample standard deviation 4.322. The prespecified target half-width of 0.5 implies N = ceil[(1.96 × 4.322 / 0.5)²] = 288; the protocol’s minimum set the confirmation at 500 replicates. Stage B used the disjoint seeds 20262030–20262529, so the confirmatory run did not reuse Stage A replicate seeds.

### 3.4 Stage B results

**Table 1. Fixed-H0 Stage B posterior summaries.** Bias is mean posterior mean minus truth. Coverage is for nominal 68% equal-tail intervals; coverage intervals are Wilson 95% intervals.

| Arm | Bias ± Monte Carlo SE (km s⁻¹ Mpc⁻¹) | RMSE | 68% coverage | Wilson 95% interval | Mean interval width | Mean density ISE |
|---|---:|---:|---:|---:|---:|---:|
| Analytic truth | +0.612 ± 0.161 | 3.647 | 72.0% | 67.9–75.8% | 7.351 | 0.0000 |
| Oracle KDE | +0.918 ± 0.160 | 3.695 | 72.8% | 68.7–76.5% | 7.376 | 0.0200 |
| Catalogue only | −8.847 ± 0.138 | 9.369 | 9.6% | 7.3–12.5% | 8.745 | 0.8624 |
| Correct weighting | +1.017 ± 0.163 | 3.787 | 69.6% | 65.4–73.5% | 7.382 | 0.0386 |
| Shifted transition | +0.253 ± 0.161 | 3.615 | 72.4% | 68.3–76.1% | 7.761 | 0.1280 |
| Shallow slope | −1.177 ± 0.157 | 3.709 | 70.4% | 66.3–74.2% | 7.843 | 0.0992 |

The primary paired error reduction was 5.9519 km s⁻¹ Mpc⁻¹ (paired SD 4.6784; standard error 0.2092; normal 95% interval 5.5418–6.3620; percentile-bootstrap 95% interval 5.5410–6.3573). Correct weighting produced lower absolute error than catalogue-only in 442 of 500 replicates (88.4%; Wilson 95% interval 85.3–90.9%).

The dominant finding is specific and straightforward: in this toy setup, ignoring the declared catalogue incompleteness strongly distorts the inferred host density and substantially degrades H0 error and interval coverage. Inverse-inclusion weighting under the exactly known inclusion law substantially reduces that particular error. This does not demonstrate that the selection law would be known in an actual survey, nor that weighting is generally robust under low inclusion probabilities or misspecified host populations.

The misspecified arms deserve cautious interpretation. The shifted midpoint arm had slightly lower absolute posterior-mean error than correct weighting in this particular experiment, which can arise from bias compensation and is not evidence that a misspecified selection law is safer. The shallow-width contrast versus correct weighting had a bootstrap interval spanning zero in Stage B, so the apparent larger degradation in Stage A did not become a robust finding. The paper therefore does not claim a general penalty for all selection-curve misspecification.

### 3.5 Stage B integrity audit

The simulation, independent posterior-grid auditor and supplementary raw-table analysis completed successfully for the locked Stage B artifact. The saved artifact contains 3,000 rows (500 replicates × six arms), with paired seeds and one true H0 per replicate. Checksums were verified, and arm metrics and the primary contrast were independently recomputed from the raw CSV. The workflow and artifact provenance are documented in the repository’s Stage B report.

## 4. Study B: H0-dependent noise and selection-conditioned likelihood

### 4.1 Model

Study B is a separate distance-only discrete-host model. For a host at redshift z and candidate H0,

μ(z,H0) = c z/H0,  
σ(z,H0) = sqrt(σfloor² + [fσ μ(z,H0)]²),

with fractional noise fσ = 0.12 and distance-noise floor σfloor = 8 Mpc. Event likelihoods marginalize over 12 discrete candidate redshifts with positive normalized weights. The population detection probability is computed from the same H0-dependent noise law for threshold dmin = 300 Mpc. For a fixed number Ndet of detected events,

log Lselected(H0) = Σᵢ log Lᵢ(H0) − Ndet log β(H0).

The likelihood conditions on the observed count; it does not infer an event rate or the distribution of the number of detections. Each dataset contains three detected events. Inference uses a uniform prior over [50, 90] km s⁻¹ Mpc⁻¹ and 161 H0-grid points.

### 4.2 Simulation and audit

The experiment used 300 replicate datasets for each of four truth regimes: H0 drawn from the uniform prior, and fixed H0 = 60, 70, 80. It preserves each seed, true value and detected-distance vector. The independent auditor replays the generator and independently recomputes the host-mixture likelihood, H0-dependent detection efficiency, selected likelihood and posterior summaries.

The latest successful artifact is from source commit `e9a512d46aeeb7070649ac8f851c1d86103f4379`, [GitHub Actions run #40797963979](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40797963979), artifact ID `12864343125`. It contains 1,200 records (300 per truth regime). The raw CSV is 278,551 bytes, SHA-256 `e239b18f452084790b9fa8aee2c8305f0d5dfef79c4555358f197f096626576d`; the summary is 5,464 bytes, SHA-256 `ba32c697119365a37b8dc55da648fddab235e8a85280cdc022a39ab02dfff38d`. The strengthened independent audit passed with `summary_aggregates_reconstructed: true`, replaying all 1,200 seeds and detected distances and independently recomputing posterior summaries and headline aggregate statistics. Maximum absolute posterior-mean reconstruction error was 4.3×10⁻¹⁴; medians and interval endpoints differed by at most 2.9×10⁻¹⁴; CDF ranks differed by 3.4×10⁻¹⁶. The raw CSV is identical to the prior artifact because the seeded simulation is unchanged; the summary hash changed because the recorded source commit changed.

### 4.3 Calibration results

**Table 2. H0-dependent-noise toy calibration.** Bias is posterior mean minus truth. Monte Carlo SE refers to the estimated mean bias. Coverage is compared with nominal equal-tail posterior intervals.

| Truth regime | Bias ± SE (km s⁻¹ Mpc⁻¹) | RMSE | 68% coverage (Wilson 95% interval) | 90% coverage (Wilson 95% interval) | Mean posterior-CDF rank |
|---|---:|---:|---:|---:|---:|
| Prior-drawn H0~Uniform(50,90) | −0.664 ± 0.442 | 7.678 | 72.7% (67.4–77.4%) | 93.3% (89.9–95.6%) | 0.514 ± 0.016 |
| Fixed H0=60 | +4.901 ± 0.420 | 8.758 | 67.3% (61.8–72.4%) | 92.3% (88.8–94.8%) | 0.391 ± 0.016 |
| Fixed H0=70 | +0.949 ± 0.364 | 6.366 | 74.3% (69.1–78.9%) | 96.3% (93.6–97.9%) | 0.470 ± 0.015 |
| Fixed H0=80 | −4.368 ± 0.223 | 5.829 | 88.3% (84.2–91.5%) | 97.0% (94.4–98.4%) | 0.604 ± 0.010 |

For the prior-drawn mode, the ten rank-decile counts were [18,31,36,30,27,35,33,34,31,25]; a descriptive chi-square statistic was 8.87 on nine degrees of freedom (p=0.450). The mean rank was 0.514. These results are broadly compatible with prior-predictive rank calibration at this Monte Carlo resolution, although empirical interval coverage is somewhat above nominal.

At fixed H0, the posterior mean shows shrinkage toward the center of the uniform prior: it tends to be high at H0=60 and low at H0=80. That conditional pattern is not inconsistent with prior-predictive calibration. With only three events, posterior estimates remain prior-sensitive. Fixed-truth rank distributions are not expected to satisfy the same uniformity diagnostic as ranks pooled over a parameter drawn from the prior; these fixed-truth ranks are reported descriptively rather than used as a prior-predictive SBC pass/fail test.

## 5. Numerical and software validation

The KDE correction was motivated by a concrete interpretation issue: SciPy’s gaussian_kde bw_method is a covariance-factor convention, not a direct physical bandwidth multiplier. Reflection augments the sample, so directly using that factor can cause bandwidth behavior different from the intended rule based on the original catalogue. The revised estimator directly sums reflected weighted Gaussian kernels, defines the default Scott rule using weighted effective sample size, and normalizes on the declared finite redshift grid.

For the two completed studies, the validation record includes:
- direct-sum regression tests for the corrected KDE bandwidth semantics;
- finite-domain density normalization;
- endpoint-aware posterior normalization tests;
- input-domain checks for H0 grids, distances, thresholds and noise parameters;
- posterior-grid reconstruction and rehashed-tamper negative controls;
- SHA-256 manifests for raw replicate records and summaries;
- seed-level data-generation replay and independently coded likelihood/posterior reconstruction.

These checks establish reproducible software behavior for the declared models. They are not a substitute for external code validation, expert review, or physical validation against a realistic data/selection pipeline.

## 6. Relationship to existing research

Recent work has used Gaussian-process reconstruction of incomplete galaxy catalogues to improve dark-siren inference, proposed robust catalogue-completeness tests in practical gwcosmo workflows, and extended inference to the population properties of gravitational-wave-hosting galaxies. These results show that catalogue completeness and host-population modeling are active research areas.

The project should therefore be positioned as a reproducible, controlled robustness benchmark, not as the first treatment of catalogue incompleteness or as a replacement for full hierarchical inference. Its value is the ability to compare isolated assumptions on matched synthetic data with explicit estimands, independently reconstructed artifacts and transparent failure modes.

The closest work differs from this project in both method and evidential scope:

| Work | Approach and evidence scope | Difference from this project | What this project does **not** claim |
|---|---|---|---|
| Tagliazucchi et al. (2026), *Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs*, [Phys. Rev. D 114, 083523](https://doi.org/10.1103/vmzf-1wxq), [arXiv:2608.24992](https://arxiv.org/abs/2608.24992) | Gaussian-process reconstruction of a line-of-sight galaxy-redshift distribution in a selection-aware hierarchical dark-siren likelihood; mock studies compare reconstructed distributions against homogeneous completion. | Here the controlled contrast is KDE bandwidth, inverse-inclusion weighting under a declared selection law, deliberate inclusion-curve changes, and artifact auditability. The estimands and models are not matched for a direct performance ranking. | No claim to be the first to reconstruct incomplete galaxy populations or to outperform GP reconstruction. |
| Datrier & Hendry (2026), *Implementing a robust test of galaxy catalogue completeness for dark siren measurements of the Hubble constant*, [MNRAS 549, stag828](https://academic.oup.com/mnras/article/549/1/stag828/8666390) | An apparent-magnitude completeness-limit diagnostic integrated into the practical `gwcosmo` workflow and applied to GWTC catalogues. | This project stipulates a known logistic inclusion probability and tests applying/misspecifying that law; it does not implement their magnitude-limit method or analyze GWTC data. | No claim to validate their completeness test or establish real-catalogue performance. |
| Gray, Williams & Papadopoulos (2026), *Expanding the scope of dark siren cosmology: Inferring the population properties of gravitational wave-hosting galaxies*, [arXiv:2609.04991](https://arxiv.org/abs/2609.04991) | A broader host-population inference approach intended to model host-population properties jointly with cosmological parameters in a practical dark-siren framework. | This project uses a small declared synthetic host law, fixed toy H0 designs, and a known-inclusion control; it does not fit a flexible host population jointly with cosmology. | No claim that this simplified benchmark competes with or replaces joint host-population inference. |
| CHIMERA v2.2, [tagged release](https://github.com/CosmoStatGW/CHIMERA/releases/tag/v2.2) | Public reference software relevant to hierarchical dark-siren and incomplete-catalogue workflows. Reproduction requires a compatible pinned code and data version. | The exact v2.2 mock-event and derived completeness products required by the relevant notebook have not been verified here. | No external CHIMERA replication or validation claim. |

This is a focused positioning exercise, not a systematic review. The defensible contribution is the controlled benchmark, corrected physical KDE-bandwidth convention, paired estimands and reproducible auditing—not superiority over these existing analyses.

CHIMERA v2.2 is a relevant potential external benchmark, but the exact previously requested large mock input and completeness products have not been verified in this work. A real external benchmark remains pending. It would require the exact version, data files, compatible interface, command, environment and output artifacts to be pinned and audited. A different-version data file must not be substituted silently.

## 7. Limitations

The experiments have intentionally simplified assumptions. They use toy redshift populations, fixed event counts, simple distance likelihoods and declared selection functions. They omit realistic detector response and network selection, waveform and calibration uncertainty, sky position and inclination marginalization, event-rate marginalization, peculiar velocities, redshift uncertainty, realistic cosmological distance-redshift relations and survey-level hierarchical galaxy-population inference. Study A’s correct-weighting arm assumes exactly known inclusion probabilities; this is a strong control rarely available in a real survey. Study B’s independent audit proves agreement between independently coded calculations of its own toy specification, not realism of that specification.

The expanded host-density screen subsequently completed in GitHub Actions run [#40797592794](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40797592794): 13 configurations × 40 paired replicates. All 14 manifest-listed campaign files passed SHA-256 and byte-count checks, and the saved per-replicate CSV metrics were independently re-aggregated. The workflow's independent posterior reconstruction covered a separate two-replicate benchmark only, not every row of all 13 campaign configurations. The screening indicates negligible H0-grid changes over tested resolutions but material bandwidth sensitivity. At large bandwidth, catalogue-only H0 bias sometimes appears smaller through bias compensation despite poor host-density estimation; this is not evidence that oversmoothing solves catalogue selection. See the dedicated [screening report](host_density_sensitivity_screening_results_2026-10-29.md).

## 8. Conclusions

In a declared fixed-H0 synthetic experiment, ignoring catalogue incompleteness caused a large downward bias and poor nominal 68% interval coverage; inverse-inclusion weighting with the exactly known toy inclusion law substantially reduced the error. The preregistered paired contrast was positive with a bootstrap interval excluding zero. A separate H0-dependent-noise selection-conditioned likelihood also passed seed-replay and posterior-reconstruction audits, while exposing expected conditional shrinkage under a weakly informative three-event experiment.

These results demonstrate the value of explicit controls and independently auditable simulation artifacts. They do not establish a real-data method, a new H0 measurement, or universal robustness to catalogue incompleteness. Any broader scientific claim should wait for the completed sensitivity campaign, a careful literature comparison and an exact, validated external benchmark.

## References

1. Tagliazucchi, Gair, Barbieri & Moresco, “Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs,” Physical Review D 114, 083523 (2026). [DOI: 10.1103/vmzf-1wxq](https://doi.org/10.1103/vmzf-1wxq); [arXiv:2608.24992](https://arxiv.org/abs/2608.24992).
2. “Implementing a robust test of galaxy catalogue completeness for dark siren measurements of the Hubble constant,” Monthly Notices of the Royal Astronomical Society 549(1), stag828 (2026). [Publisher article](https://academic.oup.com/mnras/article/549/1/stag828/8666390).
3. Gray, Williams & Papadopoulos, “Expanding the scope of dark siren cosmology: Inferring the population properties of gravitational wave-hosting galaxies” (2026). [arXiv:2609.04991](https://arxiv.org/abs/2609.04991).
4. CHIMERA v2.2 software release (9 June 2026). [GitHub release and versioned source](https://github.com/CosmoStatGW/CHIMERA/releases/tag/v2.2).
