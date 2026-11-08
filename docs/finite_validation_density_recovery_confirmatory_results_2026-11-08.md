# Confirmatory finite-validation density-recovery results (2026-11-08)

## Scope and status

This report records two completed, repository-native **synthetic host-redshift density-recovery** experiments. It does not report an (H_0) posterior result.

- Disjoint finite-parent validation: 100 replicates, 13,500 estimator attempts.
- Fixed-analysis-catalogue validation: 100 replicates, 4,500 estimator attempts.
- Both full raw artifacts passed their independent numerical reconstruction auditors.
- All planned attempts in both artifacts were SUCCESS; there were no FAILED, NON_EVALUABLE, or NOT_RUN attempts.
- The experiments use different sampling designs and answer related but distinct questions. Their estimates must not be pooled as if they were one design.

The DGP is synthetic: a two-overdensity parent redshift intensity on (z\in[0.01,0.40]), with logistic catalogue inclusion
\[
\pi(z)=\operatorname{expit}\{\alpha-x+\gamma(x^2-c)\},
\quad x=(z-0.01)/0.39,
\]
where (c=E_{\rm parent}[x^2]\), and the intercept is quadrature-calibrated to mean retention 0.25. The curvature values are \(\gamma\in\{-1.2,0,+1.2\}\).

## Run and artifact provenance

### Disjoint finite-parent design

- GitHub Actions run: [#42476772981](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/42476772981)
- Executed commit: `b7322d1c6e14c3f000109c0920bceda41a6cbf4a`
- Artifact: `finite-validation-density-recovery-confirmatory-py3.12` (artifact ID `13616821671`)
- Design: 100 replicates; parent frame size 4,000; validation sizes 100, 300, 1,000; gamma -1.2, 0, +1.2; degrees 1 and 2; Scott multipliers 0.15, 0.40, 1.00; 801-point grid.
- Auditor: **PASS**; 13,500/13,500 density estimates reconstructed; 900 parent partitions checked; 100 RNG stream sets checked; 4 manifest files verified; 13,500 SUCCESS attempts.
- Audit limitations: numerical reconstruction is not scientific validation of the DGP; no (H_0) likelihood/posterior/coverage is audited; production and auditor share SciPy's `gaussian_kde` primitive.

### Fixed-analysis-catalogue design

- GitHub Actions run: [#42476989338](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/42476989338)
- Executed commit: `c18de7b013cdf583a710dba552daf796fad0e994`
- Artifact: `fixed-analysis-catalogue-density-confirmatory-py3.12` (artifact ID `13616760895`)
- Design: 100 replicates; fixed analysis parent frame size 4,000; independent validation frame size 1,000; nested validation sizes 100, 300, 1,000; gamma -1.2, 0, +1.2; degrees 1 and 2; primary Scott multiplier 0.15; 801-point grid.
- Auditor: **PASS**; 4,500/4,500 density estimates reconstructed; 300 frame/scenario partitions checked; 100 five-stream RNG sets checked; 4 manifest files verified; 4,500 SUCCESS attempts.
- Confirmatory analysis: predeclared primary interaction and 18 Holm-adjusted secondary contrasts were generated. The primary interaction was independently recomputed from the raw attempt CSV and matched the published analysis JSON exactly.
- Audit limitations: synthetic density recovery only; production and auditor share SciPy's `gaussian_kde` primitive; validation assumes complete binary retention-label ascertainment.

Both workflow artifacts have 30-day retention at upload. Preserve local copies if longer-term reproducibility is required.

## Primary fixed-analysis-catalogue result

The analysis catalogue is generated once per replicate/scenario and held fixed across validation sizes. The validation frame is independent. Thus this design estimates the effect of additional validation information conditional on a fixed analysis catalogue.

The predeclared primary estimand at validation size 300, degree-1 fitted inclusion model, and Scott multiplier 0.15 is
\[
\Psi=E[\Delta_{+1.2,r}-\Delta_{-1.2,r}],
\qquad
\Delta_{\gamma,r}=\operatorname{IAE}_{\rm estimated\ IPW}
-\operatorname{IAE}_{\rm catalogue-only}.
\]

Results:

- Mean interaction: **+0.277893**
- Monte Carlo SE: **0.005418**
- Replicate-level percentile-bootstrap 95% interval: **[0.267450, 0.288651]**
- Paired-(t) 95% interval: **[0.267142, 0.288644]**
- Complete replicate pairs: **100/100**
- Failures in all four contributing estimator/scenario cells: **0**

Positive \(\Psi\) means the estimated-IPW-versus-catalogue-only density-recovery contrast is less favourable under positive than negative curvature. This supports the predeclared **synthetic, sign-asymmetric curvature interaction** in this particular DGP and design. It does not establish that IPW is generally harmful under positive curvature or generally beneficial under negative curvature.

### Fixed-catalogue degree-1 paired contrasts

All values are estimated-IPW IAE minus catalogue-only IAE. Negative values favour estimated IPW for this density metric; positive values favour catalogue-only. Each cell has 100 complete replicate pairs. Intervals below are paired-(t) 95% intervals.

| Curvature | Validation n | Mean paired difference | MC SE | 95% interval |
|---:|---:|---:|---:|---:|
| -1.2 | 100 | -0.221862 | 0.007778 | [-0.237295, -0.206429] |
| -1.2 | 300 | -0.259266 | 0.004270 | [-0.267739, -0.250793] |
| -1.2 | 1000 | -0.278635 | 0.003574 | [-0.285727, -0.271543] |
| 0.0 | 100 | -0.027729 | 0.008780 | [-0.045151, -0.010306] |
| 0.0 | 300 | -0.073940 | 0.004663 | [-0.083192, -0.064688] |
| 0.0 | 1000 | -0.093501 | 0.003771 | [-0.100983, -0.086020] |
| +1.2 | 100 | +0.066190 | 0.008497 | [+0.049330, +0.083050] |
| +1.2 | 300 | +0.018627 | 0.003656 | [+0.011373, +0.025882] |
| +1.2 | 1000 | +0.002452 | 0.002232 | [-0.001976, +0.006879] |

At positive curvature, the degree-1 estimated-IPW penalty diminishes toward zero as validation size increases. At negative curvature, its advantage increases in magnitude in this design. The zero-curvature case is correctly specified by the linear model, and estimated IPW improves over catalogue-only in these simulated conditions.

The 18 secondary contrasts were tested with Holm adjustment as specified in the analysis plan. The table above is descriptive and shows pointwise intervals; adjusted p-values are in the workflow artifact's `confirmatory_analysis.json`.

## Disjoint finite-parent design

Here validation units are selected from the same finite parent frame and excluded from the analysis catalogue. Consequently the analysis catalogue changes with validation size. This represents an operational finite-parent trade-off rather than a pure information-only effect.

### Degree-1 paired contrasts at primary bandwidth 0.15

| Curvature | Validation n | Mean paired difference | MC SE | 95% paired-t interval |
|---:|---:|---:|---:|---:|
| -1.2 | 100 | -0.222821 | 0.008271 | [-0.239233, -0.206410] |
| -1.2 | 300 | -0.257871 | 0.005015 | [-0.267823, -0.247920] |
| -1.2 | 1000 | -0.272966 | 0.003965 | [-0.280833, -0.265099] |
| 0.0 | 100 | -0.036827 | 0.007906 | [-0.052513, -0.021140] |
| 0.0 | 300 | -0.072161 | 0.005021 | [-0.082124, -0.062198] |
| 0.0 | 1000 | -0.089254 | 0.003648 | [-0.096494, -0.082015] |
| +1.2 | 100 | +0.054600 | 0.008923 | [+0.036894, +0.072306] |
| +1.2 | 300 | +0.017184 | 0.003206 | [+0.010822, +0.023546] |
| +1.2 | 1000 | +0.000569 | 0.001902 | [-0.003205, +0.004343] |

Mean analysis-catalogue sizes were 979.35, 929.12, and 752.72 for validation sizes 100, 300, and 1,000 respectively. This loss of analysis units is a defining part of the disjoint design, not an incidental nuisance.

The replicate-derived curvature interaction at validation size 300 and multiplier 0.15 was +0.275055 (MC SE 0.005561) for degree 1. This is a descriptive secondary calculation from the raw paired outputs, not the separately predeclared fixed-catalogue primary estimand.

### Bandwidth sensitivity

For the disjoint design at validation size 300, degree 1:

| Scott multiplier | Gamma -1.2 | Gamma 0 | Gamma +1.2 |
|---:|---:|---:|---:|
| 0.15 | -0.2579 | -0.0722 | +0.0172 |
| 0.40 | -0.2181 | -0.0456 | +0.0069 |
| 1.00 | -0.0736 | -0.0053 | +0.0019 |

The sign pattern persists in these three predeclared bandwidth settings, but the magnitude is strongly bandwidth-sensitive. This is a material qualification: these findings are conditional on the specified weighted/reflected KDE protocol and bandwidth grid, not a bandwidth-invariant statement about IPW.

## Cross-design comparison

At validation size 300, degree 1, Scott multiplier 0.15, the disjoint and fixed-catalogue contrasts are numerically close:

| Curvature | Disjoint design | Fixed catalogue |
|---:|---:|---:|
| -1.2 | -0.257871 | -0.259266 |
| 0.0 | -0.072161 | -0.073940 |
| +1.2 | +0.017184 | +0.018627 |

This close agreement is encouraging as a robustness check across these two simulated sampling arrangements. It is not a formal equivalence test, and the two designs have distinct random frames and estimands. In particular, the disjoint design loses analysis units as validation size grows, whereas the fixed-catalogue design does not.

## What the experiments establish — and what they do not

### Supported, narrowly

1. The two corrected-RNG simulation pipelines ran to completion for 100 replicates and their independent auditors reconstructed every recorded density estimate and declared raw design component.
2. Under this synthetic DGP, estimated IPW's density-IAE contrast relative to catalogue-only depends materially on curvature sign, validation sample size, fitted model degree, and KDE bandwidth.
3. The predeclared fixed-catalogue curvature interaction is positive and estimated precisely in this simulation.
4. The corresponding disjoint-design paired contrasts show a similar qualitative pattern, despite the validation-induced reduction in analysis-catalogue size.

### Not established

- No (H_0) likelihood or posterior has been evaluated by these experiments.
- No posterior calibration, simulation-based calibration (SBC), or repeated-sampling coverage result follows from density IAE.
- No real galaxy survey, LVK/GWTC event, or observational (H_0) analysis has been validated.
- No universal ranking or superiority of estimated IPW is established.
- The known-inclusion IPW and parent-sample KDE arms are simulation diagnostics, not operational methods available when true inclusion probabilities or the full parent frame are unknown.
- The observed curvature asymmetry is a result for this DGP, retention level, redshift interval, parent intensity, validation design, logistic model family, merger weighting, and KDE protocol. It is not yet a general theorem or causal decomposition of the error mechanism.
- All attempts succeeded in these 100-replicate runs. This does not demonstrate that separation or non-evaluability is rare under more severe finite-validation regimes.

## Next scientific work

1. Preserve these raw artifacts and analysis JSONs as the density-stage evidence.
2. Add stress regimes that vary retention, curvature magnitude, validation size, and parent-frame shape, with predeclared interaction families and independent seed blocks. Do not tune on the same replicates used for confirmatory claims.
3. Examine whether probability calibration, fitted-coefficient uncertainty, and weight-tail diagnostics predict the density-recovery failures or sign changes. Treat these as mechanistic analyses, not post-hoc proof.
4. Implement the modular/cut coefficient-uncertainty integration stage. For validation data (V), form a coefficient distribution (p(\theta\mid V)), then for each draw reconstruct the host density and evaluate the full selection-conditioned likelihood:
   \[
   L_{\rm marg}(H_0)=\int L(H_0;\theta)p(\theta\mid V)\,d\theta.
   \]
   The GW detection selection correction must be evaluated within each coefficient draw before marginalization.
5. Validate grid convergence separately from posterior calibration. Then run fixed-truth SBC and repeated-sampling coverage experiments with explicit failure accounting.

**Decision:** density-recovery stage = **completed for these two declared synthetic designs, with both raw artifacts independently audited**. (H_0) propagation and posterior-calibration stage = **not yet implemented/validated**.


## Additional diagnostics and Bayesian-prior parameterization issue

### Weight stability is not a sufficient diagnostic

The fixed-catalogue raw attempt artifact shows a concrete case where weights look well behaved but density reconstruction still worsens. At gamma=+1.2, validation n=300, degree 1 and Scott multiplier 0.15:

- Mean Kish ESS is approximately 991 out of approximately 1,002 analysis-catalogue objects.
- Mean maximum normalized weight is approximately 0.0012.
- Mean top-1% weight mass is approximately 0.0122.
- Nevertheless, estimated-IPW minus catalogue-only IAE is +0.018627 (100 paired replicates).

Thus high ESS and low weight concentration do not establish that the fitted inclusion function is correctly specified. Probability-model calibration and density recovery must remain separate audit dimensions.

For degree 2 at gamma=-1.2 and validation n=100, the fixed-catalogue arm has mean Kish ESS about 683, minimum ESS about 40 across replicates, mean maximum normalized weight about 0.0072, and mean top-1% weight mass about 0.061. The degree-2 IAE contrast is -0.102249, less favourable than the degree-1 contrast -0.221862 in the same cell. This is consistent with greater finite-validation variability in the more flexible fit, but the present summaries do not by themselves identify a causal mechanism.

### Prior basis must not inherit validation-sample scaling

Code inspection of the current frequentist fitter confirms that its numerical design matrix uses the validation sample's own redshift mean and standard deviation. That is acceptable as an MLE parameterization. It is not appropriate to place the same fixed Gaussian prior directly on those sample-standardized coefficients across repeated validation samples: the induced prior on the redshift-domain inclusion curve would then vary with the realized validation sample.

For Bayesian coefficient sampling, define a fixed domain coordinate
[
u=(z-z_{\min})/(z_{\max}-z_{\min})\in[0,1],
]
and place the prior on coefficients in the fixed basis `[1,u]` or `[1,u,u^2]`. An equivalent implementation may transform coefficients from the MLE's sample-standardized basis into this fixed basis before applying the prior. The prior and its scales must be declared before the H0 confirmatory experiment.

### Independent arithmetic audit added

A separate analysis auditor has now been added at `src/dark_siren_robustness/fixed_analysis_catalogue_analysis_audit.py`, with a CLI wrapper and adversarial tests. It independently reconstructs the primary interaction, bootstrap and paired-t intervals, all 18 secondary contrasts, and Holm adjustment from the raw attempt ledger. It does not import the production analysis function.

The downloaded confirmatory analysis JSON was also independently recomputed locally from the raw CSV: the primary result and all secondary means, p-values and Holm-adjusted p-values matched to floating-point precision (maximum absolute discrepancy below 2.3e-16). The new automated auditor is wired into the fixed-catalogue workflow; its GitHub Actions result remains the final verification gate for the new code.


## Downstream implementation status update (2026-11-08)

The historical status sentence above reflects the state when the density-recovery results were first written. Since then, finite-validation coefficient-to-H0 propagation has been implemented as **focused, test-covered research software**:

- A fixed-basis proper-prior logistic posterior is available through both a multi-chain MCMC sampler and a separate Student-t importance sampler.
- Weighted coefficient proposals now flow through IPW host-density reconstruction, event-mixture likelihood evaluation, draw-specific GW detection efficiency, per-draw selection correction, and self-normalized marginalization.
- The vectorized likelihood implementation is tested against the existing single-density selection-likelihood function. Tests explicitly reject the average-beta shortcut.
- H0 posterior normalization requires an explicit prior density and reports 68%/95% intervals and a boundary-mass diagnostic.

Focused workflow [#42479890034](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/42479890034) passed the inclusion-posterior and H0 propagation test suite on Python 3.10, 3.11 and 3.12 at PR head `c47ab95e26adaede691f0d71306eaf3b517fa892`. The tests include a real weighted importance-posterior object passed through the density and selected-likelihood pipeline. The general numerical-diagnostics workflow for that head was cancelled, so no full-repository CI pass is claimed.

**Revised status:** H0 propagation code exists and its focused computational tests pass. There is still no confirmatory finite-validation H0 experiment, no SBC, no fixed-truth H0 coverage study, and no evidence that estimated-inclusion IPW improves H0 calibration. The prior scale remains provisional. The project is therefore at **implementation/smoke-validation**, not scientific calibration or coverage.


### End-to-end DGP smoke and H0 marginalization Monte Carlo diagnostics

A further end-to-end smoke test now constructs a finite parent frame, calibrates the synthetic quadratic inclusion law, draws an outcome-independent validation sample, excludes validation units from the analysis catalogue, fits the proper-prior inclusion posterior, simulates threshold-detected GW distances from the known synthetic host law, propagates 1,000 weighted coefficient proposals through IPW KDE and the draw-specific GW selection correction, and normalizes the resulting H0 posterior under an explicit flat-on-H0 grid prior. It asserts partition disjointness, non-evaluable lifecycle checks, finite efficiencies, normalized posterior density and sensible interval ordering. It does **not** use the single smoke posterior as a coverage result.

The H0 propagation layer now also reports disjoint-batch Monte Carlo diagnostics for the marginalized selected log likelihood: 10 batch estimates, per-batch importance ESS, log-scale batch MCSE, and explicit `PASS` / `LOW_BATCH_ESS` / `INSUFFICIENT_DRAWS` status. The default minimum batch ESS is 20. This is separate from the importance sampler's validation-evidence MCSE.

Focused workflow [#42480779870](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/42480779870) passed on Python 3.10, 3.11 and 3.12 at PR head `54ac5a99eba3e8cef3eaffee9a492c6dcf91628e`, including the end-to-end smoke, vectorized/reference likelihood equivalence, per-draw selection-order test, posterior normalization, batch MCSE, and memory-guard checks.

**Still not a scientific H0 result:** this is a single synthetic computational smoke, not a confirmatory repeated-sampling experiment. No SBC, PIT calibration, fixed-truth H0 coverage, prior-scale sensitivity study, or real-data inference has been completed.
