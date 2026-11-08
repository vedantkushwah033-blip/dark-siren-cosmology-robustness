# Design: propagate finite-validation host-density error into dark-siren H0 inference

**Date:** 2026-11-08  
**Target branch:** `research/independent-posterior-audit-claim-fix`  
**Status:** design only; implementation is gated on the full density-recovery confirmatory artifact passing its independent audit.  
**Scope:** synthetic low-redshift Euclidean, fixed-detected-count likelihood. No LVK/GWTC data, real survey validation, detector-response claim, or observational H0 estimate.

## 1. Why a separate H0 propagation stage is required

The repository now has a finite-validation density-recovery runner and independent auditor. That stage estimates the merger-weighted parent-host redshift density from a catalogue subject to a separately generated inclusion process. It does not calculate H0 posteriors.

The existing `host_density_benchmark.py` is **not** a safe drop-in for this next stage:

- It generates a different parent/catalogue DGP through `generate_synthetic_population` and `catalogue_inclusion_probability`, rather than the finite-validation benchmark's declared two-overdensity parent law and logistic inclusion law.
- Its arms are oracle, catalogue-only, and known-completeness; it does not fit inclusion probabilities on a disjoint validation sample.
- It uses the legacy `bandwidth_factor=0.8` interface. That is an absolute SciPy covariance factor, not (0.8\times) Scott's rule, and it is not the frozen `scott_multiplier` protocol used by the density-recovery experiment.
- It therefore cannot establish propagation of *this project's* finite-validation errors into H0 inference.

The H0 stage must consume the same parent-frame, validation, retention and fitted-probability design as the density benchmark, or independently regenerate it from the same frozen DGP and RNG protocol. It must not silently combine the old and new simulations.

## 2. Scientific estimand

Let (g_0(z)) be the normalized analytic merger-weighted parent-host density from the finite-validation density protocol:

[
g_0(z)=\frac{\lambda(z)m(z)}
{\int_{0.01}^{0.40}\lambda(u)m(u)\,du},
\qquad m(z)=1+1.5x(z).
]

For each replicate, the catalogue-only, known-inclusion IPW, estimated degree-1 IPW, estimated degree-2 IPW, and parent-sample KDE arms produce candidate densities (widehat g_a(z)). The analytic (g_0) is an additional exact reference arm, not an estimated method.

The propagation question is:

> How do finite-validation inclusion estimation and its induced host-density reconstruction errors affect bias, interval width, frequentist coverage, and calibration diagnostics of a fixed-count dark-siren H0 posterior under the declared synthetic selection model?

The initial H0 stage is a **plug-in** study: each estimated logistic fit supplies a point estimate (widehat\pi(z)), which determines an IPW density. It does not propagate logistic coefficient uncertainty. Results must be labelled accordingly. The observed-information covariance from the current fitter is an asymptotic model-based covariance, not a Bayesian posterior, and must not be sampled as though it were one.

## 3. Data-generating process: separate catalogue inclusion from GW detection

For each Monte Carlo replicate:

1. Generate the parent frame, validation permutation and catalogue-retention uniforms using the finite-validation RNG stream convention:
   - stream 1: parent redshifts;
   - stream 2: validation permutation;
   - stream 3: catalogue-retention uniforms.
2. Apply the same declared (pi(z)) law, nested validation sizes and disjoint validation/analysis partition as the density-recovery benchmark.
3. Generate one set of GW event observations from the **analytic** (g_0(z)), independently of catalogue-retention uniforms and validation assignment.
4. Reuse that exact event set for every density arm in that replicate. This paired design isolates changes in the assumed host density from event-noise variation.

For a fixed true value (H_{0,*}), latent event redshifts are drawn from (g_0(z)). Distances obey

[
d\mid z,H_{0,*}\sim
\mathcal N\!\left(\frac{c z}{H_{0,*}},\sigma_d^2\right),
\qquad D=\{d\ge d_{\rm th}\}.
]

The observed event sample is conditioned on detection. A transparent generator may use rejection sampling from (z\sim g_0) and the Gaussian distance model until (N_{\rm det}) events satisfy the threshold. An equivalent implementation may draw (z) directly from the detection-weighted density proportional to

[
g_0(z)\,
\Phi\!\left(\frac{cz/H_{0,*}-d_{\rm th}}{\sigma_d}\right)
]

and then draw (d\) from the corresponding lower-truncated normal. The two implementations must be shown to agree in a generator-level test. Save the latent redshifts for audit/scoring only; never pass them to the inference likelihood.

Catalogue inclusion and GW detection are distinct processes. The inclusion model is used only to construct a candidate host density. The GW threshold is used only in the event generator and selection-conditioned likelihood. Do not multiply the catalogue inclusion probability into the GW detection efficiency, and do not add a second missing-host mixture correction after IPW without a separately derived generative model.

## 4. Fixed-count selection-conditioned likelihood

For each candidate density (widehat g_a(z)), compute

[
q_k(H_0;\widehat g_a)=
\int_{z_{\min}}^{z_{\max}}
\mathcal N\!\left(d_k\mid\frac{cz}{H_0},\sigma_d^2\right)
\widehat g_a(z)\,dz
]

and

[
\beta(H_0;\widehat g_a)=
\int_{z_{\min}}^{z_{\max}}
\Phi\!\left(\frac{cz/H_0-d_{\rm th}}{\sigma_d}\right)
\widehat g_a(z)\,dz.
]

The conditional likelihood for a fixed number (N_{\rm det}) of detected events is

[
\log L_a(H_0)=
\sum_{k=1}^{N_{\rm det}}\log q_k(H_0;\widehat g_a)
-N_{\rm det}\log\beta(H_0;\widehat g_a).
]

This is the same mathematical form as `selection_aware_toy_likelihood`. It is a fixed-count conditional likelihood: it does not infer an event rate or use a Poisson count term. The density used in each event numerator must also be used in that arm's (eta) denominator. The analytic-reference arm uses (g_0) in both places.

Use stable log-sum-exp or equivalent quadrature for the event mixtures. Independently reconstruct at least one replicate/arm's event numerators and (eta(H_0)) using a second implementation that does not import the production likelihood function.

## 5. Frozen arm definitions

For each replicate, validation size, curvature, and predeclared bandwidth multiplier, evaluate:

| Arm | Host density input | Interpretation |
|---|---|---|
| Analytic truth | (g_0(z)) | Exact model reference; not an estimator |
| Catalogue-only | KDE from retained non-validation catalogue | Uncorrected operational baseline |
| Known-inclusion IPW | KDE using true (pi(z_i)) | Oracle diagnostic only |
| Estimated linear IPW | KDE using validation-fitted degree-1 (widehat\pi(z_i)) | Finite-validation plug-in method |
| Estimated quadratic IPW | KDE using validation-fitted degree-2 (widehat\pi(z_i)) | Finite-validation plug-in method |
| Parent-sample KDE | KDE from the full parent frame | Oracle finite-sample diagnostic, not an operational method |

True inclusion probabilities and latent event redshifts are available to the simulator/auditor only. They must not leak into either estimated-inclusion arm or the catalogue-only arm.

Use `scott_multiplier` consistently. The primary candidate is 0.15, with 0.40 and 1.00 as predeclared sensitivity values, matching the density-recovery protocol. Do not tune the multiplier on H0 outcomes. If the density confirmatory stage changes or rejects this bandwidth protocol, freeze a revised rule before H0 calibration experiments begin.

## 6. Posterior and calibration estimands

For a declared prior (p(H_0)), compute

[
p_a(H_0\mid d,D)\propto
p(H_0)\exp\{\log L_a(H_0)\}
]

on a sufficiently fine H0 grid. A possible initial fixed-truth design is (H_{0,*}=70) km s(^{-1}) Mpc(^{-1}), a uniform prior on [55,85], and a grid-convergence-tested grid. These are proposed design values, not frozen until the implementation protocol is reviewed.

For each arm report, at minimum:

- posterior mean and median;
- equal-tail 68% and 95% intervals;
- interval widths;
- bias and RMSE across paired replicates;
- empirical 68% and 95% fixed-truth coverage;
- posterior-integral and boundary-mass diagnostics;
- event-generation acceptance rate and minimum/maximum detection efficiency;
- fit/density attempt status and failure rate;
- paired replicate-level contrasts against catalogue-only and analytic truth.

Fixed-truth coverage is a frequentist repeated-sampling property. A posterior PIT computed at a fixed (H_{0,*}) is a repeated-sampling calibration diagnostic, not simulation-based calibration (SBC). If SBC is later added, draw (H_0) from the declared prior for each replicate and audit the rank distribution separately. Do not label fixed-truth coverage or a fixed-truth PIT histogram as SBC.

Coverage estimates require Monte Carlo uncertainty intervals. For example, with 100 replicates the standard error of a nominal 68% coverage estimate is about 4.7 percentage points; 100 replicates are a pilot/precision check, not a definitive calibration result. Choose the final replicate count from a predeclared Monte Carlo precision target, and report the uncertainty rather than treating a point estimate near 68% as proof of calibration.

## 7. Staged implementation gates

### Gate A — finish density-stage validation first

- Run the exact 100-replicate, (N=4,000), 801-point density confirmatory design.
- Download the artifact and inspect its manifest, raw partitions, fit ledger, density vectors, and independent audit report.
- Confirm that all planned cells are present and that failure/non-evaluable states are retained.
- Reconstruct the primary paired contrasts from replicate-level records.
- Do not proceed to a scientific H0 claim if this gate fails.

The exact-head PR CI currently verifies the density-recovery **smoke** benchmark (one replicate, (N=1,200), 201-point grid) and its independent audit. That is useful software validation but is not the 100-replicate confirmatory experiment.

### Gate B — numerical likelihood validation

- Test density normalization, nonnegative support and grid convergence.
- Verify (eta(H_0;g)\in(0,1]) and independently integrate it.
- Independently reconstruct event-mixture numerators and selected log likelihoods.
- Test the likelihood under (g=g_0), where the simulation and inference host laws match.
- Verify that the selected likelihood is maximized in the expected region under high-information synthetic settings, without requiring exact recovery in every noisy replicate.
- Compare two event-generation algorithms (rejection and detection-weighted/truncated-normal) on their implied distributions.

### Gate C — paired fixed-truth H0 pilot

- Begin with a small deterministic pilot to test the full event-to-posterior pipeline.
- Use common event data across all host-density arms within each replicate.
- Use the exact finite-validation split and fitted probabilities from the same replicate.
- Preserve all attempt statuses and posterior-grid arrays/checksums.
- Run independent artifact reconstruction before interpreting any arm differences.

### Gate D — calibration experiment

- Freeze the prior, (H_{0,*}), event count, distance noise, detection threshold, redshift grid, H0 grid and replicate count before execution.
- Establish H0-grid and redshift-quadrature convergence independently of host-density comparisons.
- Estimate coverage Monte Carlo uncertainty and predeclare the required precision.
- Separate the primary contrast family from secondary curvature/bandwidth analyses.
- Report both conditional-on-success posterior performance and attempt-level failure rates.

### Gate E — coefficient uncertainty (separate method)

The plug-in stage does not propagate finite-validation coefficient uncertainty. A later uncertainty-aware stage must define a genuine probability model for logistic coefficients, for example

[
p(\theta\mid V)\propto
L_{\rm validation}(\theta)p(\theta),
]

and integrate over that posterior:

[
L_{\rm mod}(H_0)=
\int L(H_0;\theta)\,p(\theta\mid V)\,d\theta.
]

Use posterior draws or a validated integration scheme with the correct weighting. Do not average log likelihoods when the target is an average likelihood; use a stable log-mean-exp of likelihood values. Do not add the coefficient prior twice when using posterior draws. Do not substitute Gaussian draws from the inverse observed information and call them a Bayesian posterior. State whether the cut/factorized posterior is an approximation or derive the full joint model.

## 8. Independent H0 artifact auditor

The H0 auditor must be a separate implementation and must not import the production likelihood, posterior summarizer, event generator, density estimator, or benchmark runner. It should independently verify:

1. Configuration, commit, software versions, seed and manifest hashes.
2. Parent-frame, validation, analysis-catalogue and holdout partition identities.
3. Truth inclusion probabilities and retention uniforms against the declared DGP.
4. Validation-fitted coefficients and fitted probabilities from stored validation records.
5. That estimated arms never use oracle probabilities and catalogue-only uses unit weights.
6. Host-density vectors, normalization, support, bandwidth metadata and checksums.
7. Event-generation distribution, detection rule and latent-truth separation.
8. Event-mixture numerators, (eta(H_0)), selected likelihood and posterior normalization.
9. Posterior summaries, interval coverage, PIT values and paired contrasts.
10. Attempt-level status accounting and missing/failed outputs.
11. Grid convergence and Monte Carlo uncertainty calculations.

Adversarial tests should alter an event distance, threshold, host-density value, fitted probability, arm label, posterior density, H0 grid, (eta) vector, summary statistic or checksum and require the auditor to reject the artifact.

The audit report must state explicitly that numerical reconstruction does not validate the scientific adequacy of the toy DGP or any real gravitational-wave detector model.

## 9. Explicit non-claims

Even after all stages pass, the defensible conclusion is limited to the declared synthetic low-redshift model and tested parameter regimes. Do not claim:

- a real-data H0 measurement;
- validation against LVK/GWTC event catalogues or real galaxy surveys;
- a realistic detector response or event-rate model;
- universal superiority or failure of IPW;
- that the plug-in stage propagates coefficient uncertainty;
- that fixed-truth coverage is SBC;
- that a numerically audited posterior is automatically scientifically calibrated.

The candidate contribution remains a controlled, auditable study of how finite-validation catalogue-inclusion estimation and its induced host-density errors propagate into a simplified dark-siren H0 posterior.
