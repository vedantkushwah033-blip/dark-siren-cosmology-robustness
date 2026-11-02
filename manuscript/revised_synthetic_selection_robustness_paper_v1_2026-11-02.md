# Catalogue-selection misspecification and host-density reconstruction in simplified dark-siren inference: a reproducible synthetic robustness study

**Manuscript status:** research draft v1.0, 2 November 2026. Internal working document; not submitted, not peer reviewed, and not yet cleared for external submission.

**Evidence basis:** the repository's locked Stage B report, Study B H0-dependent-noise report, Study C selection-condition reports, numerical validation records, and claim-to-evidence ledger. The linear finite-population bootstrap benchmark and the newest estimated-inclusion curvature experiments are deliberately not treated as evidence for the cosmological claims in this manuscript.

## Abstract

Dark-siren inference combines gravitational-wave luminosity-distance information with a statistical description of the potential host-galaxy population. When galaxy catalogues are incomplete, the inferred host-redshift distribution can be distorted by catalogue selection, and this distortion can propagate into posterior summaries for the Hubble constant, H0. We present a reproducible set of controlled synthetic experiments that isolate selected aspects of catalogue incompleteness, host-density estimation, selection-curve misspecification, and likelihood conditioning. The central fixed-H0 benchmark uses 500 paired synthetic datasets, each with 4,000 parent galaxies and 30 detected events, at injected H0 = 70 km s−1 Mpc−1. Under the declared logistic catalogue-retention law, catalogue-only inference had posterior-mean bias −8.847 km s−1 Mpc−1 and nominal 68% interval coverage 9.6%; inverse-inclusion weighting using the exactly known synthetic retention probabilities had bias +1.017 km s−1 Mpc−1 and coverage 69.6%. The prespecified paired reduction in absolute posterior-mean error was 5.9519 km s−1 Mpc−1 (10,000-resample percentile-bootstrap 95% interval 5.5410–6.3573). A separate selection-function experiment varied the true inclusion floor and transition width. In its 250-replicate-per-condition follow-up, the paired error reduction from catalogue-only to known-law weighting was 7.436, 2.246, and 2.354 km s−1 Mpc−1 across the low-floor, high-floor, and broad-transition conditions, respectively. Those conditions share some seed labels and are not independent replications. A distinct 1,200-dataset distance-only experiment examined an H0-dependent noise model and its associated detection-efficiency normalization; its independent auditor replayed the data-generation seeds and reconstructed the likelihood and posterior summaries, while fixed-truth results exposed prior-sensitive conditional behavior. These results support a narrow conclusion: within the declared one-dimensional synthetic designs, ignoring catalogue selection can substantially degrade inference, while the magnitude of the weighting benefit depends on the selection regime and estimator assumptions. They do not establish performance on real LVK data, a generally reliable estimator for real surveys, or a new measurement of H0. All numerical claims are tied to their specific synthetic design and audit scope.

**Keywords:** dark sirens; standard sirens; Hubble constant; galaxy catalogue incompleteness; inverse-inclusion weighting; selection functions; host-galaxy population; simulation-based validation; reproducibility.

---

## 1. Introduction

Gravitational-wave standard sirens provide luminosity-distance information that can be combined with redshift information to constrain cosmological parameters. For a dark siren, no unique electromagnetic counterpart identifies the source redshift. Instead, the analysis marginalizes over a population of potential host galaxies, using a galaxy catalogue or a model for the underlying host distribution. The inferred cosmological posterior can therefore depend on how the host population, catalogue incompleteness, and gravitational-wave detection process are represented.

Catalogue incompleteness is an established challenge in dark-siren cosmology. Existing approaches include direction-dependent completeness corrections, large-scale-structure-informed variance completion, spatially correlated catalogue reconstruction, Gaussian-process reconstruction of line-of-sight galaxy distributions, and statistical tests for magnitude-completeness limits. Recent work has also considered joint inference of host-galaxy population properties and cosmological parameters. The present study does not introduce catalogue incompleteness as a new problem and does not claim to replace these methods.

Our narrower objective is to build controlled simulations in which selected statistical assumptions can be changed independently and their consequences measured using paired estimands, explicit uncertainty summaries, and inspectable computational artifacts. We focus on a one-dimensional redshift-only setting, where a declared parent galaxy population is thinned by a known retention mechanism. This permits a clean comparison among catalogue-only inference, a diagnostic oracle using known inclusion probabilities, and deliberately misspecified selection inputs. It also allows numerical and software errors to be investigated separately from the assumptions of the data-generating model.

Three experiments are reported, and their results are not pooled:

1. **Study A — fixed-H0 host-density benchmark.** A 500-replicate paired experiment compares six host-density inputs under a fixed synthetic cosmology and a known logistic catalogue-selection law.
2. **Study B — H0-dependent-noise likelihood diagnostic.** A separate 1,200-dataset experiment tests a discrete-host likelihood with H0-dependent distance uncertainty and a detection-efficiency normalization derived from the same noise model.
3. **Study C — selection-regime sensitivity.** A separate experiment changes the true catalogue-inclusion floor or transition width and examines how the prespecified catalogue-only versus known-law-weighting contrast changes.

A follow-up host-density screen examines numerical bandwidth and integration-grid sensitivity. It is exploratory and is not treated as a powered confirmation.

The intended contribution is a reproducible robustness benchmark and an evidence trail: paired comparisons, a corrected physical interpretation of weighted KDE bandwidth, explicit separation of oracle and implementable arms, seed replay where available, independent reconstruction of derived outputs, and a clear accounting of what each audit does and does not establish.

## 2. Scope and claim boundaries

All results in this manuscript are synthetic. The H0 values are injected parameters in controlled simulations. No LVK/GWTC event data, real galaxy catalogue, detector-calibrated selection function, or observational H0 estimate is analyzed.

The study does not:
- resolve or adjudicate the Hubble tension;
- estimate a new real-world H0;
- validate a real survey's completeness function;
- demonstrate that inverse-inclusion weighting is universally reliable;
- compare performance directly with Gaussian-process, LIGHT, variance-completion, CHIMERA, or other modern catalogue-reconstruction pipelines;
- establish a realistic three-dimensional galaxy or detector-response model.

The oracle/known-inclusion arm is a diagnostic control: it is given the true synthetic inclusion probabilities and is not available in that form in a real survey. It must not be described as a practical estimator unless the inclusion probabilities are independently known or estimated with an explicitly validated procedure.

## 3. Statistical setup

### 3.1 Parent population and catalogue retention

Let a simulated parent frame contain galaxies indexed by j, with redshift z_j. A catalogue-retention indicator R_j is generated according to a declared probability law

[
R_j mid z_j sim operatorname{Bernoulli}{pi(z_j)},
]

where (pi(z)) is the synthetic inclusion probability. The observed catalogue consists of parent-frame galaxies with R_j = 1. The principal benchmark uses a logistic transition with an inclusion floor of 0.10, midpoint 0.22, and width 0.05. The simulation retains the complete parent frame so that the true host-density target and the oracle inclusion probabilities are known for diagnostic purposes.

For a retained galaxy, inverse-inclusion weighting uses a weight proportional to

[
w_j = rac{1}{pi(z_j)}.
]

The weighted empirical distribution is then smoothed to form a one-dimensional host-redshift density. This is a controlled implementation of inverse-probability weighting; it does not account for all sources of uncertainty in an estimated selection model and should not be interpreted as a full hierarchical treatment of catalogue completeness.

### 3.2 Host-density and posterior comparisons

Study A compares six paired arms within each synthetic replicate:

1. **Analytic truth:** the declared normalized synthetic target density; an inaccessible diagnostic control.
2. **Oracle KDE:** a reflected KDE estimated from the complete parent population.
3. **Catalogue only:** a reflected KDE estimated from the incomplete catalogue without correction.
4. **Known-law weighting:** inverse-inclusion weighting using the exact synthetic retention probabilities.
5. **Shifted transition:** weighting based on a supplied logistic midpoint of 0.27 rather than the true 0.22.
6. **Shallow transition:** weighting based on a supplied transition width of 0.075 rather than the true 0.05.

The primary paired estimand for replicate r is

[
D_r =
left|widehat H_{0,mathrm{catalogue},r}-H_{0,mathrm{true}}ight|
-
left|widehat H_{0,mathrm{weighted},r}-H_{0,mathrm{true}}ight|.
]

Positive D means that the known-law weighting arm has smaller absolute posterior-mean error than the catalogue-only arm. The estimand concerns posterior-mean absolute error only; it does not by itself establish better interval coverage, calibration, or general estimator superiority.

### 3.3 KDE bandwidth convention

The weighted KDE implementation uses reflected kernels to address the finite redshift domain. A software review identified that SciPy's `gaussian_kde` bandwidth argument is a covariance-factor convention rather than a direct physical bandwidth multiplier. The revised implementation defines the physical bandwidth explicitly, including a weighted Scott-style default based on the weighted effective sample size. The explicit bandwidth factor is interpreted relative to the weighted standard deviation of the original sample, rather than being passed uncritically as a covariance factor after reflection augmentation. Direct kernel-sum regression tests were added to protect this interpretation.

This correction is important because a numerically valid KDE can still implement a different smoothing rule than the one intended by the statistical specification.

### 3.4 H0-dependent distance noise and selection conditioning

Study B is a separate distance-only discrete-host model. For a host at redshift z and candidate H0,

[
mu(z,H_0)=rac{cz}{H_0}, qquad
sigma(z,H_0)=sqrt{sigma_{m floor}^2+[f_sigmamu(z,H_0)]^2},
]

with fractional noise (f_sigma=0.12) and noise floor (sigma_{m floor}=8) Mpc. Each event likelihood marginalizes over 12 discrete candidate redshifts with positive normalized host weights. The detection probability is computed using the same H0-dependent noise model and a distance threshold of 300 Mpc. For a fixed number of detected events, the implemented selected log-likelihood is

[
log mathcal L_{m selected}(H_0)
=
sum_{i=1}^{N_{m det}}log mathcal L_i(H_0)
-
N_{m det}log eta(H_0),
]

where (eta(H_0)) is the model's detection efficiency. This experiment conditions on the observed number of detections; it does not infer an event rate or model the detected-count distribution. Each dataset contains three detected events, and inference uses a uniform prior on [50, 90] km s−1 Mpc−1 with 161 H0 grid points.

The fixed-count selection-conditioned formulation is reported as implemented in this toy model. It should not be conflated with a full population-level gravitational-wave likelihood that marginalizes over event-rate and other population hyperparameters.

## 4. Simulation designs

### 4.1 Study A: fixed-H0 Stage B benchmark

Each replicate contains 4,000 parent galaxies and 30 detected events. The injected value is H0 = 70 km s−1 Mpc−1. The H0 grid has 81 equally spaced points over [50, 90] km s−1 Mpc−1, and the redshift grid has 401 points over [0.01, 0.40]. The toy distance model uses a fixed Gaussian noise scale of 140 Mpc and an observed-distance detection threshold of 250 Mpc. The true catalogue-retention curve has floor 0.10, midpoint 0.22, and width 0.05.

Stage A used 200 paired replicates to estimate variability of the primary contrast. The prespecified Stage B confirmation used 500 replicates with a disjoint seed range. The Stage B artifact contains 3,000 rows (500 replicates × six arms). The six arms are paired within each replicate.

### 4.2 Study C: selection-regime stress conditions

Study C retains the parent-population size, detected-event count, injected H0, and inference-grid sizes of Study A, but varies the true selection curve. The baseline has floor 0.10, midpoint 0.22, and width 0.05. The stress conditions use:
- low inclusion floor: 0.05;
- high inclusion floor: 0.25;
- broader transition: width 0.08.

The baseline pilot used 200 paired replicates. An initial stress screen used 50 replicates per condition. A larger follow-up used 250 paired replicates per condition. Each condition compares the same six conceptual arms as Study A.

**Dependence qualification:** the three 250-replicate condition seed ranges overlap: the low-floor and high-floor conditions share 150 integer seeds, high-floor and broad-transition share 150, and low-floor and broad-transition share 50. The within-condition paired contrast remains paired by construction. The conditions must not be described as independent replications, and cross-condition uncertainty must account for shared seed labels and possible dependence.

### 4.3 Study B: H0-dependent-noise calibration

Study B uses 300 datasets in each of four truth regimes: H0 drawn from Uniform(50, 90), and fixed H0 values of 60, 70, and 80 km s−1 Mpc−1. The total is 1,200 datasets. The auditor replays all data-generation seeds and detected distances, then independently recomputes the host-mixture likelihood, detection efficiency, selected likelihood, and posterior summaries.

### 4.4 Host-density numerical sensitivity screen

A separate exploratory campaign varied 13 configurations across bandwidth factors, integration-grid resolutions, and supplied selection-curve parameters, with 40 paired replicates per configuration. It used 6,000 parent galaxies and 40 detected events per replicate. Its purpose was to identify numerical sensitivities, not to establish a definitive estimator ranking.

## 5. Results

### 5.1 Study A: catalogue selection can dominate the fixed-H0 toy benchmark

The Stage B results are summarized in Table 1. Bias is the mean posterior mean minus the injected truth. Coverage refers to nominal 68% equal-tail posterior intervals; the reported uncertainty interval is a Wilson 95% interval.

| Arm | Bias ± Monte Carlo SE (km s−1 Mpc−1) | RMSE | 68% coverage | Wilson 95% interval | Mean interval width | Mean density ISE |
|---|---:|---:|---:|---:|---:|---:|
| Analytic truth | +0.612 ± 0.161 | 3.647 | 72.0% | 67.9–75.8% | 7.351 | 0.0000 |
| Oracle KDE | +0.918 ± 0.160 | 3.695 | 72.8% | 68.7–76.5% | 7.376 | 0.0200 |
| Catalogue only | −8.847 ± 0.138 | 9.369 | 9.6% | 7.3–12.5% | 8.745 | 0.8624 |
| Known-law weighting | +1.017 ± 0.163 | 3.787 | 69.6% | 65.4–73.5% | 7.382 | 0.0386 |
| Shifted transition | +0.253 ± 0.161 | 3.615 | 72.4% | 68.3–76.1% | 7.761 | 0.1280 |
| Shallow transition | −1.177 ± 0.157 | 3.709 | 70.4% | 66.3–74.2% | 7.843 | 0.0992 |

The primary paired reduction in absolute posterior-mean error was 5.9519 km s−1 Mpc−1. The paired standard deviation was 4.6784 and the Monte Carlo standard error was 0.2092. The normal-approximation 95% interval was [5.5418, 6.3620], and the 10,000-resample percentile-bootstrap 95% interval was [5.5410, 6.3573]. Known-law weighting had lower absolute error than catalogue-only inference in 442 of 500 replicates (88.4%; Wilson 95% interval 85.3–90.9%).

The dominant result is that ignoring the stipulated selection mechanism severely distorts the inferred host density in this particular design. The catalogue-only posterior mean has a large negative bias, and its nominal 68% interval covers the injected value in only 9.6% of replicates. Weighting with the exact synthetic inclusion probabilities greatly reduces this error. The result is conditional on the specified simulator and on access to the true inclusion law.

The misspecified arms do not show a simple monotonic degradation. The shifted-midpoint arm has slightly lower absolute posterior-mean error than the known-law arm in this run, plausibly through finite-sample bias compensation. The shallow-width contrast against known-law weighting is inconclusive in Stage B. These observations argue against a broad statement that every selection-curve misspecification necessarily worsens every chosen posterior summary.

### 5.2 Study C: the size of the weighting contrast varies by selection regime

In the 200-replicate baseline pilot, catalogue-only inference had bias −9.468 km s−1 Mpc−1, RMSE 9.984, and 68% coverage 7.5%. Known-law weighting had bias +0.351, RMSE 3.785, and coverage 71.0%. The paired absolute-error reduction was 6.536 km s−1 Mpc−1 (bootstrap 95% interval [5.920, 7.119]).

In the larger 250-replicate-per-condition follow-up, the prespecified paired reduction in absolute error from catalogue-only to known-law weighting was:

| True selection condition | Mean paired absolute-error reduction (km s−1 Mpc−1) | Bootstrap 95% interval | Catalogue-only bias | Known-law weighting bias | Known-law 68% coverage |
|---|---:|---:|---:|---:|---:|
| Inclusion floor 0.05 | 7.436 | [6.798, 8.053] | −10.617 | +1.185 | 65.6% |
| Inclusion floor 0.25 | 2.246 | [1.711, 2.764] | −4.717 | +0.938 | 69.6% |
| Transition width 0.08 | 2.354 | [1.779, 2.923] | −4.994 | +1.304 | 70.8% |

The within-condition paired estimates indicate that the benefit of the known selection correction is larger in the low-floor condition than in the high-floor and broad-transition conditions. Because seed labels overlap across conditions, these three rows should not be interpreted as independent replications or combined using an independence assumption. In addition, the misspecified shallow-transition arm had lower RMSE than the known-law arm in the high-floor and broad-transition conditions. This is consistent with finite-sample bias compensation and reinforces the need to report the declared primary estimand alongside RMSE, coverage, and model specification.

The audit reconstructed posterior grids and aggregate metrics from exported replicate-level records and checked artifact integrity. It did not independently regenerate the complete Study C simulations from their random seeds. This distinction is material: independent reconstruction from saved records is not the same as independent seed-level reproduction.

### 5.3 Study B: prior-predictive rank behavior and conditional fixed-truth behavior differ

Study B's 1,200-dataset experiment produced the following results:

| Truth regime | Bias ± SE (km s−1 Mpc−1) | RMSE | 68% coverage | 90% coverage | Mean posterior-CDF rank |
|---|---:|---:|---:|---:|---:|
| Prior-drawn H0 ~ Uniform(50,90) | −0.664 ± 0.442 | 7.678 | 72.7% | 93.3% | 0.514 ± 0.016 |
| Fixed H0 = 60 | +4.901 ± 0.420 | 8.758 | 67.3% | 92.3% | 0.391 ± 0.016 |
| Fixed H0 = 70 | +0.949 ± 0.364 | 6.366 | 74.3% | 96.3% | 0.470 ± 0.015 |
| Fixed H0 = 80 | −4.368 ± 0.223 | 5.829 | 88.3% | 97.0% | 0.604 ± 0.010 |

For prior-drawn H0, the ten rank-decile counts were [18, 31, 36, 30, 27, 35, 33, 34, 31, 25], with a descriptive chi-square statistic of 8.87 on nine degrees of freedom (p = 0.450). The mean rank was 0.514. These results are broadly compatible with prior-predictive rank calibration at the available Monte Carlo resolution, although empirical interval coverage is somewhat above nominal.

At fixed H0, posterior means show shrinkage toward the center of the uniform prior: estimates tend to be high at H0 = 60 and low at H0 = 80. This conditional behavior is not assessed by the same rank-uniformity criterion as prior-predictive simulation-based calibration. With only three detected events per dataset, posterior estimates remain sensitive to the prior. The fixed-truth results are therefore reported as conditional diagnostics, not as an SBC pass/fail result.

### 5.4 Numerical sensitivity and grid adequacy

The 13-configuration host-density screen found negligible posterior-mean differences between the tested coarse and fine integration grids at the tested settings, but material sensitivity to KDE bandwidth. At the weighted-Scott baseline, RMSE values were 3.936 for the complete-parent oracle KDE, 3.896 for known-completeness weighting, and 9.165 for catalogue-only inference. Increasing the physical bandwidth factor to 1.2 raised known-completeness RMSE to 5.330, while catalogue-only RMSE decreased to 4.238. The latter should not be interpreted as successful recovery of the true host density: oversmoothing can accidentally compensate bias for a particular posterior summary.

The campaign auditor replayed the 1,560 saved fits, recomputed 403 aggregate fields with zero discrepancy, and verified 14 manifest-listed files. One binary 68% coverage flag differed for an oracle-arm replicate whose independently reconstructed endpoint was within 9.37 × 10−5 of the injected truth, below the prespecified endpoint tolerance of 0.10. This was recorded as a boundary-sensitive classification difference, not a posterior reconstruction failure. The screen remains exploratory, with only 40 replicates per setting.

A separate extreme numerical stress test compared a 1,601-point H0 grid with a 400,001-point reference for progressively narrower distance-noise scales. The coarse grid met the prespecified 0.001 quantile tolerance only at σd = 2.0 Mpc; narrower cases were flagged as underresolved. This deliberately extreme test does not imply that the main experiment's settings were underresolved. It demonstrates that grid adequacy must be evaluated relative to the likelihood width.

## 6. Computational reproducibility and audit scope

The validation record includes direct-sum regression tests for the corrected KDE bandwidth semantics, finite-domain density normalization, endpoint-aware posterior normalization, input-domain checks, posterior-grid reconstruction, rehashed-tamper negative controls, SHA-256 manifests, and independent likelihood/posterior calculations.

The principal audits have different scopes:

- **Study A Stage B:** checksums and raw replicate aggregates were independently recomputed; posterior grids and the primary paired estimand were reconstructed by an independent auditor.
- **Study B:** the independent auditor replayed all 1,200 data-generation seeds and detected distances, then independently recomputed the host-mixture likelihood, detection efficiency, selected likelihood, posterior summaries, and headline aggregates.
- **Study C:** artifact checksums, saved posterior grids, replicate-level outputs, and aggregate metrics were independently checked/reconstructed. The full simulations were not independently regenerated from their seeds.
- **Host-density sensitivity screen:** the auditor replayed seeds and independently reconstructed the 1,560 saved fits and campaign aggregates, subject to the disclosed boundary-sensitive coverage classification difference.

A passing software test or artifact auditor establishes agreement with the declared computational specification within the audit's scope. It does not prove that the specification is a realistic physical model, that a selection law is correct for an observed survey, or that the method is calibrated outside the simulated regimes.

## 7. Relation to existing work

Dark-siren catalogue incompleteness and host weighting have been studied using a range of approaches. Pixelated catalogue methods model direction-dependent completeness; variance-completion methods incorporate large-scale-structure information; LIGHT reconstructs spatially correlated galaxy fields and magnitude distributions; Gaussian-process approaches reconstruct line-of-sight redshift distributions while accounting for survey selection; and robust completeness tests estimate apparent-magnitude limits from galaxy samples. Recent work has also considered joint inference of host-galaxy population properties.

The present one-dimensional benchmark differs in purpose and scope. It stipulates a known synthetic retention law and tests the behavior of catalogue-only and inverse-inclusion-weighted host-density estimators, with deliberate changes to the selection curve. It does not reconstruct angular clustering, infer galaxy magnitudes, use a cosmological simulation of the full galaxy field, or analyze real GW events. No direct performance ranking against these methods is warranted because the models, inputs, and estimands are not matched.

The closest literature includes:
- Tagliazucchi et al. (2026), Gaussian-process reconstruction of incomplete line-of-sight catalogues, *Physical Review D* 114, 083523, DOI: 10.1103/vmzf-1wxq.
- Leyde, Baker & Enzi (2026), spatially correlated galaxy-field and magnitude reconstruction, *JCAP* 01 (2026) 013, DOI: 10.1088/1475-7516/2026/01/013.
- Dalang, Fiorini & Baker (2026), large-scale-structure prior information and variance completion, *JCAP* 01 (2026) 034, DOI: 10.1088/1475-7516/2026/01/034.
- Datrier & Hendry (2026), robust magnitude-completeness testing in the `gwcosmo` pipeline, *MNRAS* 549, stag828, DOI: 10.1093/mnras/stag828.
- Alfradique, Bom & Castro (2025), systematic biases from catalogue incompleteness and host weighting, *Physical Review D* 112, 063561, DOI: 10.1103/vd36-3mys.
- Gray, Williams & Papadopoulos (2026), joint inference of host-galaxy population properties, arXiv:2609.04991.

This is a focused positioning summary, not a systematic review. A full literature review and a direct, matched external benchmark remain future work.

## 8. Limitations

The experiments use simplified redshift-only host populations, fixed event counts, simplified distance likelihoods, and declared selection functions. The principal known-law weighting arm receives the exact synthetic inclusion probabilities. This is an oracle condition, not a realistic assumption about a survey. The models omit realistic sky localization, angular clustering, a three-dimensional galaxy field, galaxy magnitudes and luminosities, redshift uncertainties, peculiar velocities, waveform and calibration uncertainty, inclination marginalization, detector-network response, realistic detection pipelines, and full event-rate marginalization.

The fixed-H0 benchmark uses one injected value, one parent-population design, one event count, and one main noise setting. Study C varies selected aspects of the retention law but does not cover the full space of selection-model misspecification. The Study B experiment uses only three detected events and a discrete host set; it demonstrates numerical consistency and exposes conditional prior sensitivity, not realistic detector calibration. Study C's seed overlap prevents treating its conditions as independent replications. The linear finite-population bootstrap benchmark is a separate statistical-software study and provides no direct evidence about the nonlinear dark-siren posterior.

The estimated-inclusion pilot and nonlinear-curvature studies are not included in the quantitative claims of this manuscript. They should be incorporated only after their exact source commits, manifests, completed runs, raw artifacts, independent audit reports, and analysis scripts are verified. A planned protocol or a user-reported result is not sufficient evidence for a manuscript result.

## 9. Conclusions

In a controlled one-dimensional synthetic benchmark, ignoring a stipulated galaxy-catalogue selection mechanism produced substantial bias and severe undercoverage, while inverse-inclusion weighting with the exact synthetic retention probabilities reduced the prespecified paired absolute posterior-mean error. The size of this contrast varied across the tested selection regimes. Deliberately misspecified selection curves sometimes produced apparent improvements in selected finite-sample metrics, consistent with bias compensation and demonstrating that method rankings depend on the estimand and simulation regime.

A separate H0-dependent-noise experiment showed broadly compatible prior-predictive rank behavior at the available Monte Carlo resolution while revealing conditional fixed-truth shrinkage toward the prior center. Numerical audits established reproducibility of the declared synthetic calculations at their stated audit levels, but they do not establish physical realism or external validity.

The defensible contribution is therefore a controlled and auditable synthetic robustness study, not a new cosmological measurement or a universally reliable correction method. The next scientific step is not to add more loosely connected exploratory simulations. It is to freeze a new, independently seeded protocol for estimated-inclusion uncertainty, define its primary estimands and failure accounting, independently regenerate the full artifacts, and validate the nonlinear posterior under a more realistic selection-conditioned generative model. An exact compatible external benchmark and domain-expert review are required before any claim of observational relevance.

## 10. Reproducibility and release checklist

Before submission or external circulation, complete all of the following:

- [ ] Freeze the exact manuscript source commit and hash.
- [ ] Download each artifact cited by a table or figure and verify byte count and SHA-256 against its manifest.
- [ ] Recompute every reported table from raw replicate-level records using a separately maintained analysis script.
- [ ] Record source commit, workflow run, artifact ID, configuration, seed range, denominator, and audit scope for every quantitative claim.
- [ ] Reconcile all Study C statements with the shared-seed qualification and the fact that its audit did not replay the full simulation from seeds.
- [ ] Verify the current manuscript's literature review and bibliographic metadata against publisher or arXiv records.
- [ ] Resolve the confirmatory design register for any new study before generating outcomes.
- [ ] Obtain independent domain-expert review of the likelihood conditioning, selection normalization, host-population assumptions, and physical interpretation.
- [ ] Keep the CHIMERA v2.2 external benchmark marked NO-GO until exact compatible inputs and a complete independently auditable run are available.
- [ ] Run exact-head CI and independently reproduce the final artifact set after all code and manuscript changes.

## References

1. Schutz, B. F. (1986). Determining the Hubble constant from gravitational wave observations. *Nature*, 323, 310–311.
2. Holz, D. E., & Hughes, S. A. (2005). Using gravitational-wave standard sirens. *The Astrophysical Journal*, 629, 15–22.
3. Chen, H.-Y., Fishbach, M., & Holz, D. E. (2018). A two per cent Hubble constant measurement from standard sirens within five years. *Nature*, 562, 545–547.
4. Fishbach, M., Holz, D. E., & Farr, W. M. (2019). Does the black hole merger GW170814 have an electromagnetic counterpart? *The Astrophysical Journal Letters*, 863, L41.
5. Gray, R., et al. (2020). Cosmological inference using gravitational wave standard sirens: A mock data analysis. *Physical Review D*, 101, 122001.
6. Palmese, A., et al. (2023). (Galaxy-catalogue completeness and dark-siren inference; verify exact bibliographic metadata before submission.)
7. Alfradique, V., Bom, C. R., & Castro, T. (2025). Systematic bias in dark siren statistical methods and its impact on Hubble constant measurement. *Physical Review D*, 112, 063561. DOI: 10.1103/vd36-3mys.
8. Leyde, K., Baker, T., & Enzi, W. (2026). Cosmic cartography. Part II. Completing galaxy catalogs for gravitational-wave cosmology. *Journal of Cosmology and Astroparticle Physics*, 01 (2026) 013. DOI: 10.1088/1475-7516/2026/01/013.
9. Dalang, C., Fiorini, B., & Baker, T. (2026). Large scale structure prior knowledge in the dark siren method. *Journal of Cosmology and Astroparticle Physics*, 01 (2026) 034. DOI: 10.1088/1475-7516/2026/01/034.
10. Datrier, L. E. H., & Hendry, M. A. (2026). Implementing a robust test of galaxy catalogue completeness for dark siren measurements of the Hubble constant. *Monthly Notices of the Royal Astronomical Society*, 549, stag828. DOI: 10.1093/mnras/stag828.
11. Tagliazucchi, M., Gair, J., Barbieri, R., & Moresco, M. (2026). Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs. *Physical Review D*, 114, 083523. DOI: 10.1103/vmzf-1wxq.
12. Gray, R., Williams, D., & Papadopoulos, A. (2026). Expanding the scope of dark siren cosmology: Inferring the population properties of gravitational wave-hosting galaxies. arXiv:2609.04991.

**Bibliographic note:** Reference 6 and the precise metadata for references 1–5 must be checked against the publisher/arXiv records before submission. The references above are a working bibliography, not a completed systematic literature review.
