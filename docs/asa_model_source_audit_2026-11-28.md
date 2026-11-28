# ASA Paper: Source-Level Mathematical and Design Audit
**Audit date:** 28 November 2026  
**Branch:** `research/asa-model-audit-20261128`  
**Parent checkpoint:** `research/asa-paper-checkpoint-20261128`  
**Purpose:** source-level audit of the hidden-host confirmatory candidate model and the two main factorial runners. This is a methods audit, not a new experiment, independent scientific validation, or manuscript.

## 1. Files inspected

- `src/dark_siren_robustness/hardened_inference.py`
- `scripts/run_host_misspecification_sensitivity_surface.py`
- `scripts/run_hidden_host_assumption_sensitivity.py`
- `docs/asa_paper_checkpoint_2026-11-28.md`
- `docs/host_misspecification_results_2026-10-26.md`
- `docs/known_host_calibration_run_20261025.md`
- `docs/selection_sbc_audit_results.md`
- `docs/fixed_truth_coverage_map_results.md`
- `docs/matched_fixed_h0_coverage_results.md`

The inspection was against the named research checkpoint branch. No source code or results were changed.

## 2. Mathematical model reconstructed from implementation

Let the listed-host redshifts be (z_1=0.02,z_2=0.04), with listed relative weights (q_1=0.7,q_2=0.3). Let the hidden component have true mass (f) and redshift (z_m). The host label is drawn first:
[
P(J=1)=(1-f)q_1,quad P(J=2)=(1-f)q_2,quad P(J=m)=f.
]
At fixed injected (H_0), the distance mean for host (j) is (mu_j(H_0)=cz_j/H_0). The generator draws (D^*\sim N(\mu_j(H_0),\sigma^2)) repeatedly until (D^*>0). Thus
[
p(d\mid H_0,J=j,D^*>0)
=\frac{\phi((d-\mu_j(H_0))/\sigma)}{\sigma\,\Phi(\mu_j(H_0)/\sigma)},\quad d>0.
]
Because rejection is performed **after the host has been selected and within that fixed host**, the host-label probabilities remain the stated probabilities. They are not reweighted by host-specific detection probabilities in this generator.

For an analyst-assumed hidden fraction \(\hat f\) and hidden redshift \(\hat z_m\), the implemented likelihood is
[
L(H_0\mid d)=
(1-\hat f)\sum_{j=1}^2q_j
\frac{\phi((d-cz_j/H_0)/\sigma)}{\sigma\Phi(cz_j/(H_0\sigma))}
+\hat f
\frac{\phi((d-c\hat z_m/H_0)/\sigma)}{\sigma\Phi(c\hat z_m/(H_0\sigma))}.
]
The posterior is this likelihood multiplied by a uniform prior on the specified interval and normalized numerically on the H0 grid.

### Audit conclusion on this normalization
For the *positive-truncated, host-first* generator above, the component-specific denominator \(\Phi(\mu_j/\sigma)\) is mathematically appropriate and the host mixture weights remain unchanged. This is not the same observation-selection mechanism as conditioning a complete host-mixture draw on a global detection event. Under global event selection, the denominator generally involves the mixture-level detection probability, and the host probabilities among detected events change. The two mechanisms must not be blended.

## 3. Source-level findings

### F1 — Likelihood and generator agree on positive truncation (PASS, scoped)
Both factorial runners redraw a non-positive distance for the already selected host. The inference implementation subtracts \(\log\Phi(\mu_j/\sigma)\) separately for each component when `positive_truncated=True`. This is internally coherent for the stated host-first rejection design.

**Scope:** this is an algebra/code consistency finding. It does not validate a physical detector-selection model.

### F2 — “Modelled hidden host” in the surface runner is an oracle arm (PASS with mandatory naming correction)
In `run_host_misspecification_sensitivity_surface.py`, the `hidden_host_modelled` arm passes the observation's **true** hidden fraction and true hidden redshift directly into `posterior_h0`. It is therefore an oracle / correctly specified hidden-population benchmark, not an estimated or practically deployable correction.

**Required manuscript terminology:** call it “oracle hidden-population specification” (or “oracle correctly specified hidden component”). Do not call it an estimated hidden-host method.

### F3 — Assumption factorial correctly separates truth from analyst assumptions (PASS, with interpretation limit)
In `run_hidden_host_assumption_sensitivity.py`, the data generator stores true hidden fraction/redshift. The `oracle_hidden` arm uses those values; `base_omitted` sets hidden fraction to zero; all other arms use fixed analyst-assumed fractions and redshifts from the prespecified grid. This is a valid controlled misspecification comparison for the specified synthetic family.

**Interpretation limit:** the assumed values are fixed sensitivity arms, not estimates learned from a validation sample. The experiment does not assess uncertainty in estimating hidden-population parameters.

### F4 — Pairing is within observation and within design cell (PASS)
Each generated observation is passed through every inference arm before moving to the next replicate. The bootstrap resamples replicate indices and applies the same indices to the compared arms. This supports paired within-cell contrasts, provided the raw rows preserve replicate identity and the independent auditor verifies alignment.

### F5 — The factorial's noise scale is fixed conditional on the injected truth (MUST STATE)
Both runners set
[
\sigma=\eta\,c\,(0.03)/H_{0,\mathrm{true}},
]
where \(\eta\) is the selected noise fraction. The same resulting \(\sigma\) is supplied to every inference arm for that observation. It is not re-estimated and is not candidate-H0-dependent within a posterior evaluation. The paper must state this exactly; it must not describe the experiment as fitting an unknown distance-noise model.

### F6 — The experiment is repeated single-event inference, not a multi-event hierarchical analysis (MUST STATE)
Each replicate generates one host label and one distance observation, then computes a one-event H0 posterior. Repetition across Monte Carlo replicates estimates operating characteristics; it does not create a joint posterior from many gravitational-wave events. Do not describe the present implementation as a full hierarchical population analysis.

### F7 — Grid and prior integration are aligned for the currently declared supports (PASS, bounded)
The runners use a 4,001-point grid from 25 to 125. The declared prior endpoints [50,90], [40,100], and [30,120] align with grid points. The inference code normalizes by trapezoidal integration over in-prior grid points and computes summaries on the positive-density support, avoiding artificial half-bin mass from exterior zero points at the hard prior boundary.

**Remaining numerical requirement:** the existing 4,001-vs-8,001 check is only a small five-replicate-per-cell screen. It does not certify quantile stability for all 278,400 posterior fits or coverage classifications. Full confirmatory grid refinement remains required.

### F8 — The fixed-truth and SBC modules target different probability statements (PASS, separation required)
The older interval-selected calibration module draws H0 from its prior and conditions observations on a toy detection interval. The hidden-host factorials use positive-truncated host-specific Gaussian distances. Prior-predictive posterior-rank calibration and fixed-truth repeated-sampling coverage are different estimands. The paper must present the old SBC module separately, not pool its records with hidden-host factorial outcomes or use its global rank result to certify pointwise coverage.

## 4. Important unresolved methodological issues

1. **No estimation uncertainty:** the hidden fraction and hidden redshift are fixed analyst inputs in the factorial arms. The experiments compare oracle, omission, and deliberately misspecified fixed assumptions; they do not yet compare an estimated hidden-population model with an uncertainty-integrated model.
2. **No realistic selection function:** positive-distance truncation is a mathematical support condition, not a gravitational-wave detection model. The hidden-host factorial does not include detector sensitivity, duty cycle, waveform dependence, or a realistic event-selection function.
3. **One-event toy likelihood:** no multi-event hierarchical combination, population-level H0 inference, peculiar velocities, redshift uncertainty, or full cosmological parameter covariance is included.
4. **Discrete host support:** the two listed redshifts and one hidden redshift are deliberately stylized. The hidden component is a point mass in redshift, not a continuous latent galaxy population.
5. **Noise model known to analyst:** \(\sigma\) is fixed and passed to inference. Robustness to misspecified or estimated noise is not established.
6. **Prior sensitivity is substantive:** changing support changes the inferential problem, not merely numerical stability. Report bias, RMSE, coverage, and width together; do not pick a support because it gives the most attractive coverage.
7. **Multiplicity / selection of contrasts:** the factorials contain many cells and contrasts. A confirmatory subset and primary contrasts must be prespecified; exploratory cellwise intervals must be labelled exploratory. Do not interpret unadjusted cellwise 95% intervals as a familywise error-controlled discovery procedure.
8. **Artifact verification boundary:** the checkpoint reports successful independent artifact reconstruction for the three main factorials. This source inspection did not download and re-execute every artifact. Cite the checkpoint's audit scope accurately and retain its raw manifests and hashes.

## 5. What the current evidence can support

Subject to reproducing the checkpoint's archived artifact audits, the evidence can support a **synthetic operating-characteristics study** of how omission and fixed misspecification of a discrete hidden-host component affect one-event H0 posterior summaries under a known-noise, positive-truncated, low-redshift distance model, and how those results interact with prior support.

It cannot support claims that:
- the method is validated for real LVK/GWTC data;
- the oracle arm is a usable estimator;
- hidden-population parameters have been successfully estimated;
- detector selection or galaxy-catalogue incompleteness has been realistically corrected;
- global SBC proves pointwise coverage;
- the model yields a real H0 measurement or a complete hierarchical dark-siren analysis.

## 6. Required corrections to the planned ASA paper

1. Use “oracle” consistently for the true-parameter hidden-host arm.
2. State the host-first positive-truncation mechanism and its component-specific normalization explicitly.
3. Separate the older interval-selected SBC/coverage module from the positive-truncated hidden-host factorials.
4. Label the current work as repeated one-event synthetic inference.
5. Describe \(\sigma\) as known and fixed within each truth/noise cell.
6. Present prior-support changes as a sensitivity trade-off, never as a selected optimum.
7. Distinguish fixed-parameter frequentist coverage from prior-predictive SBC.
8. Do not draft claims of estimator estimation, uncertainty propagation, real-data validity, or detector-selection correction.

## 7. Next highest-value research step

Before launching larger Monte Carlo runs, conduct a **known-answer posterior validation suite** for the exact positive-truncated mixture likelihood:
- one-host reduction: a one-component unknown-host calculation must equal the known-host posterior;
- zero-hidden-fraction reduction: \(\hat f=0\) must equal the catalogue-only likelihood exactly;
- duplicate-component invariance: splitting one host's weight across identical redshift entries must not change the posterior;
- hidden-component relabelling invariance: permuting mixture components and their paired weights must not change the likelihood;
- direct-density versus log-density agreement on a safe numerical range;
- independent high-accuracy quadrature for posterior normalization and selected quantiles;
- limiting behaviour as \(\mu/\sigma\) becomes large, where positive-truncation normalization approaches one.

These should be deterministic tests with analytic or independently computed expected values, not Monte Carlo-only checks. After they pass, perform the full grid-refinement audit on a frozen confirmatory subset, then choose replication counts from explicit Monte Carlo precision targets.

---
**Audit disposition:** source-level coherence is encouraging for the declared toy model. The most important contribution is not a new likelihood theorem; it is a careful, reproducible characterization of when hidden-population omission and misspecification degrade posterior operating characteristics, with honest separation of prior-predictive calibration and pointwise coverage. The statistical novelty review and external domain review remain open.
