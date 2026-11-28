# ASA Dark-Siren Paper — Research Checkpoint
**Checkpoint date:** 28 November 2026  
**Repository:** `vedantkushwah033-blip/dark-siren-cosmology-robustness`  
**Checkpoint branch:** `research/asa-paper-checkpoint-20261128`  
**Base commit:** `9cb4adab7d9c8a26d9d7a04e7d27ab38b7f1f1a8` on `research/d001-lower-tail-oracle-theory-20261127`

> This is a research-status record, not a manuscript and not a claim that the project is publication-ready. All results described here are synthetic. They do not validate real LVK/GWTC inference or constitute a measurement of the Hubble constant.

## 1. Current recommended ASA paper

### Recommended working title
**Calibration and Robustness of Bayesian Dark-Siren Hubble-Constant Inference under Host Ambiguity and Hidden-Population Misspecification: A Synthetic Simulation Study**

The earlier calibration-only title, “Posterior Calibration and Statistical Reliability of Host-Ambiguous Dark-Siren Cosmological Inference: A Simulation Study,” is now considered too narrow for the strongest completed work. The project contains a substantial hidden-host misspecification and prior-support programme, so the ASA paper should connect calibration diagnostics with robustness while keeping their distinct generative models separate.

### Central research question
How do omission or misspecification of an unobserved host population, together with prior support, affect point estimation and uncertainty coverage for H0 in controlled dark-siren simulations, and what do prior-predictive calibration diagnostics reveal—or fail to reveal—about pointwise fixed-truth performance?

### Intended contribution
A carefully scoped statistical simulation study demonstrating:
- how hidden-host population omission and incorrect hidden-host assumptions affect posterior bias, RMSE, interval width and coverage;
- how matching the hidden fraction but misspecifying its redshift can still degrade inference;
- how prior support trades off pointwise coverage, RMSE and interval width;
- why global prior-predictive rank calibration is not a substitute for fixed-truth coverage analysis.

Do not claim a new general Bayesian calibration theorem, real-data validation, or universal superiority of an estimator.

## 2. Scope and scientific boundaries

The main hidden-host experiments use a simplified low-redshift Euclidean model with discrete candidate hosts and one hidden-host redshift component. The separate older calibration experiments use different data-generating mechanisms.

Explicitly state that the current study does **not** include:
- real LVK/GWTC observations or a real H0 measurement;
- a realistic detector-selection function in the hidden-host factorials;
- galaxy-catalogue inclusion-probability estimation or inverse-probability weighting in those factorials;
- peculiar velocities, redshift uncertainty, full cosmological parameter covariance, or a complete hierarchical multi-event population model.

Do not merge numerical results from the positive-truncated hidden-host experiments with results from the older interval-selected calibration experiment as if they came from one model.

## 3. Completed hidden-host simulation programme

The following three corrected-revision experiments were audited from their archived GitHub Actions artifacts. The recorded common source revision is `ecb260886716c54fdd2d5e47390701949f5a18d5`; the corresponding workflows completed successfully.

| Experiment | Design | Unique generated observations | Posterior fits |
|---|---|---:|---:|
| Host-misspecification sensitivity surface | 108 cells; H0={60,70,80}; noise={5%,10%,20%}; hidden fraction={0,.15,.30,.45}; hidden redshift={.025,.030,.035}; 200 replicates/cell | 21,600 | 86,400 |
| Hidden-host assumption sensitivity | 48 cells; H0={55,70,85}; noise={5%,20%}; true hidden fraction={0,.15,.30,.45}; true hidden redshift={.025,.035}; 200 replicates/cell; 11 inference arms | 9,600 | 105,600 |
| Prior-support × hidden-host factorial | 36 cells; H0={55,70,85}; noise={5%,20%}; hidden fraction={0,.30,.45}; hidden redshift={.025,.035}; 3 prior supports; 200 replicates/cell | 7,200 | 86,400 |
| **Total** | | **38,400** | **278,400** |

### Artifact provenance
| Experiment | Workflow run ID | Artifact ID | ZIP SHA-256 prefix |
|---|---:|---:|---|
| Host-misspecification surface | 40417072115 | 12696905126 | `c2c2c485d6b5…` |
| Hidden-host assumptions | 40417072147 | 12696226361 | `ee684f9ee3b5…` |
| Prior-support factorial | 40417072086 | 12695694723 | `0e82705ebc2c…` |

Related repository records:
- `docs/host_misspecification_sensitivity_surface_results_2026-10-26.md`
- `docs/hidden_host_assumption_sensitivity_results_2026-10-26.md`
- `docs/prior_support_hidden_host_interaction_results_2026-10-26.md`
- `docs/corrected_hidden_host_rerun_provenance_2026-10-27.md`

### Independent audit checks recorded on 28 Nov 2026
- Downloaded ZIP SHA-256 digests matched the digests reported by GitHub artifact metadata.
- All 38,400 unique observations were regenerated from declared NumPy seeds and generator rules. Host indices, redshifts and sigma values matched exactly; observed distances differed by at most (5.68\times10^{-14}), consistent with CSV floating-point serialization.
- All 278,400 per-run posterior rows were independently grouped and the saved cell-level mean error, RMSE, 68%/95% coverage and mean interval-width summaries were reconstructed. Maximum discrepancy was about (1.4\times10^{-14}).
- At true hidden fraction zero, omitted-host and oracle/modelled-host posterior summaries were identical in the checked negative-control rows: 5,400 surface rows, 7,200 prior-support rows and 2,400 assumption-factorial rows.
- Numerical integration of the positive-truncated component density over (d>0) returned 1, with maximum absolute error (2.22\times10^{-16}) over checked combinations.
- A 5-replicate-per-cell 4,001-versus-8,001-point grid screen over the three experiment families found maximum absolute quantile difference 0.00144 H0 units. This is a screening check only, not a full grid-refinement audit of all raw replicates.

These checks establish internal provenance and arithmetic consistency for the archived synthetic experiments. They do not establish astrophysical realism or external validity.

## 4. Main hidden-host results recorded

### A. Omission of the true hidden population
Equal-cell averages over 12 cells per true hidden fraction in the hidden-host assumption experiment:

| True hidden fraction | Inference arm | RMSE | 68% coverage | 95% coverage |
|---:|---|---:|---:|---:|
| 0.00 | Omitted/oracle negative control | 7.378 | 67.0% | 97.1% |
| 0.15 | Omitted | 9.059 | 63.1% | 92.9% |
| 0.15 | Oracle | 8.536 | 66.4% | 98.2% |
| 0.30 | Omitted | 10.549 | 57.3% | 87.4% |
| 0.30 | Oracle | 9.147 | 68.4% | 97.8% |
| 0.45 | Omitted | 11.355 | 53.5% | 82.4% |
| 0.45 | Oracle | 9.070 | 71.2% | 98.4% |

At true hidden fraction 0.45, omission minus oracle increased RMSE by 2.286 (paired 95% interval 1.938 to 2.670), reduced 68% coverage by 17.7 percentage points (−19.7 to −15.7), and reduced 95% coverage by 16.0 percentage points (−17.4 to −14.6).

The oracle arm is supplied the true hidden fraction and redshift. It is a simulation-only benchmark, not a deployable estimator.

### B. Incorrect hidden-host assumptions
Two recorded contrasts:
- True hidden fraction 0, but analyst assumes 45% hidden at z=0.030: RMSE increases by 4.409 relative to the oracle/omitted model (paired 95% interval 4.160 to 4.662).
- True hidden fraction 0.30, analyst assumes fraction 0.30 but redshift 0.030: RMSE increases by 0.888 (95% interval 0.660 to 1.128), while 68% coverage falls by 6.2 percentage points (95% interval −7.8 to −4.7 points).

This supports the important interpretation that matching the hidden fraction alone does not remove redshift misspecification. An incorrectly specified extra component can perform worse than omitting it.

### C. Prior-support trade-off
For true hidden fraction 0.45 and the catalogue-only base-weight arm:

| Uniform prior support | RMSE | 68% coverage | 95% coverage | Mean 68% interval width |
|---|---:|---:|---:|---:|
| [50,90] | 11.868 | 50.5% | 81.0% | 16.39 |
| [40,100] | 12.957 | 75.3% | 89.1% | 26.87 |
| [30,120] | 16.783 | 82.3% | 96.0% | 40.32 |

Wider support improves empirical coverage in this scenario while increasing RMSE and interval width. Do not claim that any one support is universally preferable. Prior-support effects require careful interpretation because changing support changes both prior information and posterior summaries.

## 5. Main hidden-host model to describe accurately

The listed-host redshifts are (z_1=0.02) and (z_2=0.04), with relative listed weights (q_1=0.7) and (q_2=0.3). The hidden component has fraction (f) and redshift (z_m).

The generating host probabilities are:
[
(P(J=1),P(J=2),P(J=m))=((1-f)0.7,(1-f)0.3,f).
]

For injected H0 and selected host (J=j), the distance is generated from a Gaussian with mean (mu_j(H_0)=cz_j/H_0) and standard deviation (sigma), redrawing for the same host until (d>0). Thus the component density is:
[
p(d\mid H_0,J=j,d>0)=\frac{\phi(d;\mu_j(H_0),\sigma)}{\Phi(\mu_j(H_0)/\sigma)},\quad d>0.
]

Because truncation is performed conditional on the already-selected host, the host-mixture probabilities remain unchanged. The analyst's mixture likelihood is:
[
L(H_0\mid d)=(1-\hat f)\sum_jq_jp(d\mid H_0,z_j,d>0)+\hat f p(d\mid H_0,\hat z_m,d>0).
]

The posterior is proportional to this likelihood times a uniform prior over the chosen support. The absolute sigma is set per fixed-truth cell as noise_fraction × (c\times0.03/H_{0,true}) and is supplied as a fixed noise input to inference. Describe this exactly; do not imply an unknown or candidate-dependent noise law.

## 6. Separate older calibration module — keep separate

The older two-host calibration experiments are not the same data-generating model as the hidden-host factorials:
- One prior-predictive SBC experiment has no detection selection.
- Another conditions on a toy interval detection event (0<d<180) Mpc with (sigma=20) Mpc.
- The fixed-truth coverage map uses the latter setup.
- The hidden-host factorials instead use positive-truncated Gaussian components.

For the 2,000-replicate selection-normalized rank experiment, recorded global posterior-CDF ranks had mean 0.5009 and KS p-value 0.962. Rank-derived prior-predictive interval coverage was 67.85% at nominal 68% and 94.6% at nominal 95%.

However, exploratory 5-unit H0 bins showed strong pointwise heterogeneity: 68% coverage was 27.2% in [50,55), 98.3% in [70,75), and 13.2% in [85,90). This illustrates why global prior-predictive calibration does not guarantee good fixed-truth coverage at every H0.

The rank test is prior-predictive conditional on fixed host-mixture assumptions. It is not joint SBC over hidden-host parameters. A non-significant KS test is not proof of calibration.

## 7. What remains — ordered work plan

### P0 — before any manuscript drafting
1. **Literature and novelty review (not complete).** Compare directly against dark-siren catalogue incompleteness, missing-host/host-population reconstruction, Bayesian model misspecification, SBC, prior sensitivity and frequentist coverage literature. Build a claim-by-claim overlap matrix. Do not assert novelty until this is done.
2. **Freeze the confirmatory estimands and model-family separation.** Specify primary outcomes, primary contrasts, scenario cells, inference arms, exclusions and which existing results are exploratory. Keep prior-predictive SBC and fixed-truth robustness as separate targets.
3. **Get domain review.** Ask a knowledgeable gravitational-wave/cosmology researcher to review the toy data-generating process, meaning of mixture weights, and limits of the positive-truncation model.
4. **Resolve numerical validation.** The current 4,001-vs-8,001 comparison is only a five-replicate-per-cell screen. Run a prespecified full grid-refinement audit on the selected primary cells/artifacts and check posterior quantiles, means, coverage classification and tail behaviour.
5. **Reconcile Monte Carlo uncertainty.** At 200 replications per cell, a nominal 68% coverage estimate has approximate 95% MC half-width of 6.5 percentage points per cell. At 1,000 replications it is about 2.9 points. Decide replication targets from a prespecified precision requirement, accounting for multiple cells and paired contrasts.

### P1 — confirmatory evidence
6. Increase replications for the selected primary hidden-host and prior-support cells; do not blindly rerun every exploratory cell.
7. Recompute all primary summaries and paired confidence intervals from raw run-level records using a frozen analysis script.
8. Add prior sensitivity only where scientifically defensible and clearly label it as sensitivity analysis, not a search for a preferred result.
9. Validate likelihood/posterior calculations independently with known-answer cases and a second numerical implementation where feasible.
10. Re-run the relevant automated tests and resolve or explain failures. A separate latest full-suite record in project context reported 1,093 passed, 3 failed, 1 skipped; those failures must not be represented as a clean full-suite pass. Check whether the failures affect this paper's analysis path.

### P2 — manuscript
11. Freeze the final tables and figures only after confirmatory analysis.
12. Write methods before results interpretation; document random seeds, software versions, priors, numerical grids, data-generation rules, exclusions and code/archived artifacts.
13. Include a limitations section that states the synthetic, low-redshift, discrete-host and simplified selection assumptions.
14. Have an independent reader reproduce at least one primary result from raw data and scripts before submission.

## 8. Immediate next actions

Recommended next sequence:
1. Audit the exact source scripts and analysis code for the three hidden-host experiments against the written model above.
2. Construct a focused literature map and novelty/overlap table.
3. Propose and freeze a small confirmatory subset with explicit primary estimands and MCSE targets.
4. Only then run additional confirmatory simulations.

No new confirmatory simulation workflow was launched during the 28 November checkpoint. No repository source code was modified as part of the audit. This checkpoint is documentation only; no changes were made to `main`, no pull request was opened, and no existing experiment artifacts were deleted.

## 9. Important cautions
- The three archived hidden-host experiments and their independent reconstruction are strong internal evidence, but do not establish real-world validity.
- The older SBC result is a prior-predictive calibration diagnostic, not a blanket proof of posterior reliability.
- The oracle arm is not an implementable method.
- The prior-support table shows a trade-off, not a universally optimal prior.
- The paper's scientific value depends on a defensible literature gap and a clear statistical contribution, not merely the number of simulations.
- Do not report the manuscript as complete or ready for submission until the remaining work above is completed.
