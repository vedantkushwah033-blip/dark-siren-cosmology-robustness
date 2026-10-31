# Estimated-Inclusion Study: Confirmatory Design Draft v1

**Status:** proposed protocol; not preregistered; no confirmatory run has been executed.  
**Date:** 2026-10-31  
**Scope:** controlled synthetic experiment only. This is not a real galaxy-survey or LVK analysis and does not estimate the real Hubble constant.

## 1. Research question

Under a controlled low-redshift dark-siren toy model, how do finite external validation samples and catalogue-inclusion-model misspecification affect inverse-probability-weighted host-density reconstruction and downstream fixed-truth $H_0$ inference?

The narrow target is the reliability of estimated catalogue-inclusion correction—not a claim that inverse-probability weighting solves all catalogue incompleteness problems.

## 2. Generative scenarios

Use two prespecified data-generating scenarios:

1. **Primary, correctly specified:** logistic catalogue inclusion with intercept $0.25$, standardized-redshift slope $-1.2$, and quadratic coefficient $0$.
2. **Nonlinear stress:** same intercept and slope, with quadratic coefficient $1.2$ in the existing centered basis $x^2-1/3$.

The centering convention is part of the model definition. It is centered under a uniform reference measure on $[-1,1]$, not necessarily under the simulated parent-redshift distribution.

Use the existing toy simulation defaults unless the planning run exposes a documented numerical or positivity problem:

- Parent population: 4,000 galaxies per replicate.
- Detected toy GW events: 30 per replicate.
- External validation sample sizes: 100, 300, and 1,000, nested within replicate.
- Fixed truth: $H_0=70$.
- Toy distance noise: 140 Mpc.
- Toy observed-distance threshold: 250 Mpc.
- $H_0$ grid: 50–90 with 81 points.
- Redshift grid: 0.01–0.40 with 401 points.

The scalar observed-distance threshold is a simplified toy selection model, not injection-calibrated detector selection.

## 3. Estimator arms

For each replicate and each scenario, compare:

- **Catalogue-only:** unweighted observed-catalogue host-density baseline.
- **Oracle:** inverse-probability correction using the true inclusion probabilities; benchmark only, not an implementable estimator.
- **Estimated-linear-$m$:** fitted linear-logistic inclusion model from the external validation sample of size $m$.
- **Estimated-quadratic-$m$:** fitted quadratic-logistic inclusion model from the same nested validation sample of size $m$.

The estimated arms must fit only the external validation redshifts and retention labels. True inclusion probabilities may be used for the oracle arm and for scoring inclusion-probability error, never as fitted-model inputs.

Keep catalogue retention and GW event detection as separate stages. Do not add a hidden-host mixture term to the IPW likelihood unless a coherent joint generative model is derived; otherwise it risks double-counting incompleteness.

## 4. Pairing and randomization

Within each replicate, reuse the same parent population, catalogue-retention realization, detected-event realization, and nested validation samples across estimator arms. Preserve the existing deterministic component-seed schedule.

Use independent root seeds for the two DGP scenarios and record the exact seed, software versions, and configuration in each artifact. Do not drop failed fits; retain every attempted arm and report failure counts.

## 5. Estimands

### Primary point-estimation contrast

For replicate $r$, validation size $m$, and estimated arm $a$:

$$
D_{r,m,a}=|\widehat H_{0,r}^{\,catalogue}-H_0^{true}|
          -|\widehat H_{0,r,m,a}-H_0^{true}|.
$$

Positive values favor arm $a$ **only in absolute posterior-mean error**. They do not imply improved interval calibration or coverage.

Report the paired mean, standard error, 95% interval, number of attempted pairs, number of successful pairs, and failures by arm.

### Secondary model-family contrast

$$
Q_{r,m}=|\widehat H_{0,r,m}^{linear}-H_0^{true}|
       -|\widehat H_{0,r,m}^{quadratic}-H_0^{true}|.
$$

Positive values favor the quadratic model in absolute posterior-mean error only.

### Other prespecified outcomes

- Fixed-truth coverage of the nominal 68% interval at $H_0=70$.
- Mean interval width.
- Host-density integrated squared error.
- Inclusion-probability mean absolute error, scored on the full parent population but not used in fitting.
- Effective sample size, maximum inverse weight, weight concentration, support overlap, and positivity diagnostics.
- Fit failure and boundary-limited rates.

Coverage at this single fixed truth is not prior-averaged calibration and is not simulation-based calibration (SBC). Do not conflate these quantities.

## 6. Monte Carlo precision plan

A minimum of **2,000 replicates per DGP scenario** is proposed for the final study. For a coverage probability near 0.68, the binomial Monte Carlo standard error at $N=2{,}000$ is approximately 0.0104, corresponding to a normal 95% Monte Carlo half-width of about 0.0204 (2.04 percentage points). This is a precision calculation, not a power claim.

For the continuous paired primary contrasts, estimate the replicate-level paired-difference standard deviation in a separate planning pilot. Use a provisional Monte Carlo half-width target of 0.5 km/s/Mpc:

$$
N_D=\left\lceil\left(\frac{1.96s_D}{0.5}\right)^2\right\rceil.
$$

Set the final replicate count per scenario to the maximum of 2,000 and the largest $N_D$ across prespecified primary contrasts. Freeze that count before the confirmatory run. If it exceeds available compute, report the limitation and revise the precision target before looking at confirmatory outcomes; do not silently lower $N$.

## 7. Planning run before confirmation

First run a **non-confirmatory planning pilot of 100–250 replicates per scenario** at the full default configuration. Its purpose is to measure runtime, optimizer failure rates, positivity behavior, and paired-difference variability—not to make scientific claims.

Use planning root seeds 20261031 (linear DGP) and 20261032 (quadratic stress). Keep these separate from the proposed confirmatory root seeds 20261101 and 20261102.

After the planning pilot:
1. Independently audit both generated artifacts and their seed-level provenance.
2. Check that the intended validation-size × model-family matrix is complete.
3. Inspect failed fits, boundary-limited fits, positivity, ESS, and weight concentration without silently excluding difficult replicates.
4. Compute the replicate-level paired-difference SDs and finalize $N$ using Section 6.
5. Freeze the protocol, seeds, code commit, and run configuration before confirmatory execution.

## 8. Analysis and reporting safeguards

- Report both DGP scenarios separately; do not pool them into one headline estimate.
- Report marginal coverage with binomial/Wilson intervals and paired coverage differences separately.
- Report paired primary contrasts with their paired uncertainty; do not treat arms as independent samples.
- State explicitly that successful-pair contrasts condition on both fits succeeding.
- Report all attempted, successful, and failed denominators.
- Do not apply a universal ESS cutoff; show the full diagnostics and interpret them in context.
- Keep fixed-truth coverage distinct from SBC.
- Do not claim improved coverage from a positive $D$ or $Q$.
- Do not call the scalar threshold model a realistic GW detector selection function.
- Do not claim CHIMERA validation, real LVK validation, or a new observational $H_0$ measurement.

## 9. Execution status

The implementation and both smoke pipelines have passed three-version CI and independent artifact/provenance audits at the documented commit. Those smoke runs used only 2 and 3 replicates, respectively, and are execution checks—not scientific evidence.

**No planning pilot at the full default configuration and no confirmatory campaign have yet been run.** This draft must be reviewed and frozen before confirmatory execution.
