# Repeated-sampling H0 calibration protocol for finite-validation inclusion uncertainty

**Protocol date:** 2026-11-08  
**Status:** preregistration/design proposal; the repeated-sampling runner and independent artifact auditor are not implemented.  
**Target branch:** `research/independent-posterior-audit-claim-fix`  
**Scope:** synthetic low-redshift Euclidean toy inference only.

## 1. Question and estimands

The next experiment asks whether uncertainty in a finite-validation estimate of galaxy-catalogue inclusion probabilities changes the repeated-sampling behaviour of selection-conditioned dark-siren (H_0) posteriors.

The primary scientific question is:

> Under a fixed synthetic parent population and a separate GW detection process, how do finite validation size and inclusion-model misspecification affect the bias, interval width, and repeated-sampling coverage of (H_0) posteriors obtained from plug-in and uncertainty-integrated IPW host densities?

Primary outcomes, reported separately by scenario and inference arm:

1. Mean posterior-mean error and RMSE for fixed (H_0^star).
2. Empirical 68% and 95% interval coverage, with binomial Monte Carlo standard errors and confidence intervals.
3. Mean 68% and 95% interval widths.
4. Posterior boundary mass and fraction of non-evaluable/failed attempts.
5. Per-replicate proposal ESS, H0-specific likelihood-contribution ESS, maximum contribution weight, and log-marginal-likelihood batch MCSE.
6. Paired replicate-level contrasts against catalogue-only, with complete-pair counts and explicit failure accounting.

Density IAE remains a separate upstream outcome. It must not be substituted for posterior calibration.

## 2. Synthetic data-generating process

Use the existing declared two-overdensity parent intensity on (z\in[0.01,0.40]):
[
f_{\rm parent}(z)\propto
1+2.2\exp\left[-\tfrac12\left(\frac{z-0.105}{0.025}\right)^2\right]
+1.6\exp\left[-\tfrac12\left(\frac{z-0.275}{0.040}\right)^2\right].
]

For each replicate, draw a finite parent frame of (N_{\rm parent}=4,000) units. Generate catalogue inclusion independently conditional on redshift:
[
R_i\mid z_i\sim\operatorname{Bernoulli}(\pi(z_i)),\qquad
\pi(z)=\operatorname{expit}\{\alpha-x+\gamma(x^2-c)\},
]
where (x=(z-0.01)/0.39), (c=E_{\rm parent}[x^2]), and (\alpha) is quadrature-calibrated to mean retention 0.25 under the declared parent intensity.

Primary curvature scenarios are (\gamma\in\{-1.2,0,+1.2\}). The γ=0, degree-1 fitted model is the matched inclusion-link scenario. The nonzero-curvature, degree-1 fitted model is deliberately misspecified. Degree-2 fitting is a correctly specified model-family sensitivity for all three scenarios; it is not a claim that the real world follows a quadratic logit.

Select validation units by simple random sampling without replacement from the full parent frame, independently of (R_i). Exclude those units from the analysis catalogue. The validation labels are observed for every selected validation unit. The planned validation sizes are (n_V\in\{100,300,1000\}), nested within each replicate to preserve pairing. Catalogue inclusion and GW detection must remain separate mechanisms.

## 3. GW event generation

For every replicate, generate one shared set of (N_{\rm det}=30) detected-event distances at fixed (H_0^\star=70\;\mathrm{km\,s^{-1}\,Mpc^{-1}}), using the declared merger-weighted parent redshift law and Gaussian distance noise with σ=140 Mpc. Apply the fixed observed-distance threshold (d\ge250\) Mpc by rejection sampling. Use a dedicated event RNG stream independent of parent, validation-split, and catalogue-retention streams.

The same parent, validation split, catalogue, and detected-event realization must be reused across all inference arms for a replicate. The event simulation is not conditioned on galaxy-catalogue inclusion: GW detection and catalogue inclusion are separate.

The event count is conditioned upon in the toy likelihood. Do not add rate inference or a real detector response.

## 4. Inference arms

All arms use the same predeclared redshift grid, (H_0) grid, merger-weight function, distance-noise model, threshold, and KDE bandwidth rule.

1. **Oracle parent reference:** use the known analytic parent intensity times merger weight. This is a scoring/reference arm, not an operational estimator.
2. **Catalogue-only:** estimate the host density from retained analysis-catalogue redshifts without inclusion correction.
3. **Known-inclusion IPW:** use the simulator's true π(z) for analysis-catalogue objects. This is an oracle diagnostic, never supplied to a fitted arm.
4. **Estimated plug-in IPW:** fit the declared unpenalized logistic model on validation labels only and use its MLE predictions. Retain fit failures as FAILED or NON_EVALUABLE; never replace them silently with the Bayesian fit.
5. **Estimated uncertainty-integrated IPW:** approximate the fixed-basis proper-prior validation-only posterior with the weighted Student-t importance sampler. For every coefficient proposal, reconstruct a host density, evaluate every event numerator and the GW detection efficiency from that same density, apply the selection correction within the proposal, and then marginalize the selected likelihood using the proposal's log importance weights.

The Bayesian arm is a modular/cut analysis: GW data do not update inclusion coefficients. State this explicitly in all reports. Do not average selection efficiencies across coefficient proposals before applying the selection correction. Do not average separately normalized H0 posteriors.

## 5. Prior and bandwidth decisions

The fixed-domain coefficient prior currently used by the importance sampler is independent Gaussian with mean zero and coefficient SD 1.0. Treat this as a candidate reference, not a proven robust choice. Run prior-scale sensitivity at 0.75, 1.0, and 1.5 on the same validation records and report posterior predictive calibration, importance ESS, H0-specific likelihood-contribution ESS, and H0 outcomes.

The Scott multiplier 0.15 is the current primary candidate from the density-recovery experiments, but bandwidth effects are substantial. Before the H0 confirmatory run, freeze a bandwidth rule using density-only training/holdout evidence and do not select it using H0 coverage or posterior width. Include 0.40 and 1.00 as declared sensitivity settings. Use the same bandwidth rule across coefficient proposals and inference arms.

## 6. Staged execution and compute gates

Do not launch the full confirmatory design before the following stages pass.

### Stage A — engineering pilot

- 3 curvature scenarios.
- Validation sizes 100 and 300.
- 5 independent replicate seeds per cell.
- 1,000 importance proposals per Bayesian fit.
- Coarse but predeclared (z) and (H_0) grids, followed by at least one finer-grid comparison.
- Full raw records and independent artifact reconstruction.

Purpose: verify runtime, failure statuses, artifact size, independent reconstruction, event-generation correctness, and whether proposal and likelihood-contribution ESS thresholds are achievable. These five-replicate results are engineering diagnostics only.

### Stage B — Monte Carlo and grid pilot

Using new, non-overlapping seed blocks, increase proposal draws (for example 2,000 and 4,000) and compare at least two redshift and (H_0) grid resolutions. Compare posterior means, quantiles, boundary mass, and log marginal likelihoods. The grid must resolve the posterior; an interval that is clipped by the grid is not a valid coverage success.

If the H0-specific contribution ESS is low at any material region of the grid, increase proposal quality/draw count or mark the replicate NON_EVALUABLE. Do not discard low-ESS replicates from the attempt ledger.

### Stage C — confirmatory repeated-sampling coverage

Only after Stages A and B establish stable computation, freeze:
- the DGP and scenario cells;
- the prior scale and sensitivity family;
- KDE bandwidth rule;
- proposal count and ESS gates;
- redshift/H0 grids and convergence tolerance;
- replicate count and seed list;
- primary paired contrasts and interval estimands;
- all status and exclusion rules.

Use at least 100 independent replicates per confirmatory cell as a starting point, and report the binomial uncertainty of coverage. If precision requirements imply more replicates, increase the count before inspecting confirmatory outcomes. Keep discovery/tuning seeds disjoint from confirmatory seeds.

## 7. Required raw artifacts

The runner must write immutable, machine-readable artifacts containing:

- exact config, software versions, source commit, seed, RNG stream identifiers, and a SHA-256 manifest;
- one replicate/scenario attempt row for every planned attempt and every inference arm;
- parent redshifts, retention labels, validation IDs, analysis IDs, true inclusion probabilities (scoring/auditor only), and observed GW distances;
- validation model degree, prior settings, optimizer/sampler status, proposal coefficients, log importance weights, normalized weights, and all ESS/MCSE diagnostics;
- host-density grids for all proposals or a lossless deterministic representation sufficient to reconstruct them;
- per-proposal event log likelihoods, per-proposal β_GW(H0), selected log likelihoods, and marginalized log likelihood;
- normalized H0 posterior density on the exact grid, posterior summaries, boundary mass, and coverage indicators;
- attempt lifecycle status: SUCCESS, NON_EVALUABLE, FAILED, or NOT_RUN, with a reason.

Never serialize only aggregate means. Preserve all planned attempts, including failed and non-evaluable cases. Hash exact raw file bytes and record array shape/dtype/endianness where binary arrays are used.

## 8. Independent audit requirements

The auditor must not merely read the runner's summary JSON and repeat its claims. It must independently:

1. Verify config, seed streams, parent/validation/analysis partition disjointness, and the catalogue-retention realization.
2. Reconstruct the true inclusion probabilities from the declared DGP and verify that they were not passed to fitted arms.
3. Recompute fitted-arm inclusion probabilities from stored coefficients and validation redshifts; verify coefficient/proposal dimensions and importance-weight normalization.
4. Reconstruct selected likelihoods from stored host densities and event distances using the slow reference likelihood path; verify per-draw β_GW and the order of selection correction.
5. Recompute the weighted log marginal likelihood and H0 posterior normalization/quantiles/coverage from raw arrays.
6. Recompute all aggregate bias, RMSE, coverage, interval-width, paired-contrast, Monte Carlo error, and failure-ledger totals.
7. Verify hashes and report any missing, extra, duplicated, malformed, or non-finite records.

The auditor may share SciPy's underlying KDE primitive with the production path, but it must use a distinct orchestration path and explicitly disclose that shared numerical primitive.

## 9. Calibration analyses are distinct

### Matched-model SBC

Run a separate simulation-based calibration experiment in which inclusion coefficients are drawn from the declared prior, validation labels are generated under the same model, and posterior ranks are computed for the coefficients and predeclared scalar functionals. Assess rank uniformity with uncertainty bands and repeat across independent seeds. This tests the validation-posterior computation under its matched model.

Do not label the deliberately misspecified degree-1 cases as SBC. They are fixed-truth misspecification/coverage experiments.

### Fixed-truth H0 coverage

For every DGP scenario and inference arm, evaluate empirical coverage of the declared 68% and 95% H0 credible intervals at the fixed (H_0^\star). Report binomial standard errors/intervals, interval widths, bias, RMSE, boundary mass, and all failed/non-evaluable attempts. Coverage is not inferred from a single end-to-end smoke test.

### Predictive diagnostics

Only after the repeated-sampling data exist, assess whether pre-posterior diagnostics predict H0 calibration failures. Train and evaluate any diagnostic prediction rule on disjoint seed blocks or held-out scenario cells. Do not use the same replicates for rule discovery and claimed predictive validation.

## 10. Claim boundary and stop rules

This protocol does not establish a result. It defines the next experiment.

Stop or mark a replicate NON_EVALUABLE if:
- the validation posterior sampler fails its declared ESS/MCSE gates;
- any H0 grid point has inadequate likelihood-contribution ESS;
- any required coefficient draw, host density, likelihood, or efficiency is invalid;
- grid refinement changes the declared posterior summaries beyond the frozen tolerance;
- the posterior has material boundary mass under the predeclared threshold;
- a required raw artifact or independent audit check is missing.

Do not silently drop or replace these replicates. Report unconditional planned-attempt denominators and conditional complete-pair analyses separately.

Even a successful study would support only a controlled synthetic-methodology conclusion under the declared parent population, inclusion law, event model, priors, bandwidth, and selection rule. It would not validate a real galaxy survey, LVK/GWTC event analysis, or observational (H_0) measurement, and would not establish universal estimator superiority.

## 11. Current implementation checkpoint

At PR head `dae374c17a462b65aec94c06ca05f730ac094e37`:

- The finite-validation density-recovery experiments are complete for the two declared designs and their raw artifacts passed independent audits.
- The fixed-basis weighted Bayesian inclusion importance sampler is implemented.
- The vectorized H0 propagation path and a slower reference path are implemented.
- The end-to-end parent-frame → disjoint validation → Bayesian inclusion posterior → detected-event distances → draw-specific selected likelihood → normalized H0 posterior smoke test is in the test suite.
- H0-specific likelihood-contribution ESS and batch MCSE diagnostics are implemented.
- The full numerical-diagnostics matrix and focused inclusion-posterior sampler matrix passed on Python 3.10, 3.11, and 3.12 (runs `42482040471` and `42482040468`).

Not yet implemented or established:
- a repeated-sampling H0 runner and raw artifact schema;
- an independent end-to-end H0 artifact auditor;
- multi-replicate prior-scale sensitivity;
- matched-model inclusion SBC;
- fixed-truth H0 coverage and H0 grid-convergence confirmatory results;
- held-out predictive validation of pre-posterior calibration diagnostics.

**Decision:** the computational plumbing has passed its current CI gates. The H0 scientific-calibration stage remains open. Do not describe the one-replicate smoke test as coverage, SBC, or evidence of estimator superiority.
