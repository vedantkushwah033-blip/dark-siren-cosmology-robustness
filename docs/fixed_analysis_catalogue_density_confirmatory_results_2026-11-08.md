# Confirmatory fixed-analysis-catalogue density-recovery results

**Date:** 2026-11-08  
**Status:** completed and independently audited  
**Scientific scope:** synthetic density recovery only; no H0 posterior-calibration result.

## Executive result

The 100-replicate fixed-analysis-catalogue (FAC) experiment completed successfully. The independent artifact auditor returned **PASS** for all 4,500 estimator attempts and all 300 replicate-by-curvature frame scenarios. All 4,500 attempts were SUCCESS; there were no FAILED, NON_EVALUABLE, or NOT_RUN attempts.

The predeclared pilot-informed symmetric-curvature interaction was

\[
\Psi=\Delta_{+1.2}-\Delta_{-1.2}=0.277893,
\]

with replicate-level Monte Carlo SE 0.005418, paired-t 95% interval [0.267142, 0.288644], and 10,000-resample replicate-ID percentile-bootstrap 95% interval [0.267450, 0.288651]. Positive values mean the estimated-IPW-versus-catalogue-only density-recovery contrast was less favourable under positive than negative curvature.

The result is a synthetic density-recovery finding. It does **not** establish H0 posterior calibration, coverage, observational validity, or universal superiority or inferiority of IPW.

## Reproducibility record

- Repository: `vedantkushwah033-blip/dark-siren-cosmology-robustness`
- Workflow run: [42476290113](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/42476290113)
- Workflow commit: `77461eaad4d794c1f9c9cc0eaa1cad33267a1836`
- Artifact ID: `13615594133`
- Artifact name: `fixed-analysis-catalogue-density-confirmatory-py3.12`
- Artifact SHA-256 (ZIP): `b71f68a21a8f5867f8de094c0fe3de0a3507e9e5797a3bf22ef0f7426ac78fcf`
- Runtime: Python 3.12.16, NumPy 2.5.4, SciPy 1.18.1, pandas 3.1.1.
- Raw files: `attempt_metrics.csv` (4,500 rows), `frame_units.csv` (1,500,000 rows), `summary.json`, `config.json`, `manifest.json`, `audit_report.json`, and `confirmatory_analysis.json`.

The downloaded ZIP SHA-256 matched GitHub's artifact digest. All four manifest-covered raw/config/summary files matched their declared byte counts and SHA-256 checksums.

## Design

- 100 simulation replicates.
- Per replicate and curvature: 4,000 analysis-parent units and an independent validation frame of 1,000 units.
- Nested validation sizes: 100, 300, and 1,000.
- Curvature: \(\gamma\in\{-1.2,0,+1.2\}\).
- Target retention: 0.25.
- Logistic inclusion fit degrees: 1 and 2.
- Scott multiplier: 0.15.
- Redshift grid: 801 points on [0.01, 0.40].
- Five arms: catalogue-only, known-inclusion IPW (oracle diagnostic), parent-frame KDE (oracle diagnostic), estimated degree-1 IPW, and estimated degree-2 IPW.

The analysis catalogue is generated once per replicate/curvature and remains identical across validation sizes. Thus validation-size contrasts condition on a fixed analysis catalogue and isolate added validation information. This is deliberately distinct from the disjoint finite-parent design, where validation size also reduces the available analysis catalogue.

The inclusion law is
\[
\pi(z)=\operatorname{expit}\{\alpha-x+\gamma(x^2-c)\},
\quad x=(z-0.01)/(0.40-0.01),
\]
where \(\alpha\) is calibrated to target retention 0.25 and \(c=0.31149625211621723\) is the parent-law expectation of \(x^2\). The degree-1 logistic fit is misspecified when curvature is nonzero; degree 2 matches the quadratic logit specification.

Five independent random streams were generated with `SeedSequence([seed, replicate, stream_id])`: analysis redshifts, validation redshifts, validation permutation, analysis retention uniforms, and validation retention uniforms.

## Primary estimand and result

For replicate \(r\) and curvature \(\gamma\), define
\[
\Delta_{r,\gamma}
=\mathrm{IAE}_{\text{estimated degree-1 IPW}}
-\mathrm{IAE}_{\text{catalogue-only}},
\]
at validation size 300 and Scott multiplier 0.15. The primary interaction is
\[
\Psi_r=\Delta_{r,+1.2}-\Delta_{r,-1.2}.
\]

| Quantity | Result |
|---|---:|
| Planned replicate IDs | 100 |
| Complete paired replicate IDs | 100 |
| Incomplete replicate IDs | 0 |
| Mean interaction | +0.277893 |
| Replicate Monte Carlo SE | 0.005418 |
| Paired-t 95% interval | [0.267142, 0.288644] |
| Percentile-bootstrap 95% interval | [0.267450, 0.288651] |
| Bootstrap resamples | 10,000 |

The primary interaction was selected after exploratory pilot evidence and is explicitly **pilot-informed**, not an independent hypothesis-generation exercise.

## Secondary contrasts

Each entry is estimated-IPW IAE minus catalogue-only IAE. Negative values favour estimated IPW for this density-recovery metric; positive values favour catalogue-only. Intervals are pointwise paired-t 95% intervals. P-values are Holm-adjusted across the full predeclared family of 18 contrasts.

| Curvature | Validation n | Degree | Mean paired IAE difference | MC SE | Pointwise 95% interval | Holm-adjusted p |
|---:|---:|---:|---:|---:|---:|---:|
| -1.2 | 100 | 1 | -0.221862 | 0.007778 | [-0.237295, -0.206429] | 2.25e-48 |
| -1.2 | 100 | 2 | -0.102249 | 0.021993 | [-0.145887, -0.058612] | 5.15e-05 |
| -1.2 | 300 | 1 | -0.259266 | 0.004270 | [-0.267739, -0.250793] | 6.04e-79 |
| -1.2 | 300 | 2 | -0.214221 | 0.008564 | [-0.231213, -0.197229] | 1.93e-43 |
| -1.2 | 1,000 | 1 | -0.278635 | 0.003574 | [-0.285727, -0.271543] | 2.01e-89 |
| -1.2 | 1,000 | 2 | -0.269323 | 0.004408 | [-0.278070, -0.260575] | 3.53e-79 |
| 0.0 | 100 | 1 | -0.027729 | 0.008780 | [-0.045151, -0.010306] | 0.00632 |
| 0.0 | 100 | 2 | +0.083879 | 0.021062 | [+0.042087, +0.125671] | 0.000521 |
| 0.0 | 300 | 1 | -0.073940 | 0.004663 | [-0.083192, -0.064688] | 6.92e-28 |
| 0.0 | 300 | 2 | -0.040933 | 0.007875 | [-0.056559, -0.025307] | 8.68e-06 |
| 0.0 | 1,000 | 1 | -0.093501 | 0.003771 | [-0.100983, -0.086020] | 3.77e-43 |
| 0.0 | 1,000 | 2 | -0.081516 | 0.004187 | [-0.089824, -0.073207] | 1.53e-34 |
| +1.2 | 100 | 1 | +0.066190 | 0.008497 | [+0.049330, +0.083050] | 6.75e-11 |
| +1.2 | 100 | 2 | +0.132089 | 0.018948 | [+0.094492, +0.169685] | 3.20e-09 |
| +1.2 | 300 | 1 | +0.018627 | 0.003656 | [+0.011373, +0.025882] | 1.17e-05 |
| +1.2 | 300 | 2 | +0.029817 | 0.005964 | [+0.017982, +0.041651] | 1.49e-05 |
| +1.2 | 1,000 | 1 | +0.002452 | 0.002232 | [-0.001976, +0.006880] | 0.549 |
| +1.2 | 1,000 | 2 | -0.001579 | 0.002849 | [-0.007231, +0.004073] | 0.581 |

All 18 secondary contrasts were evaluable with 100 paired replicates each. I independently recomputed all 18 paired means, raw paired-t p-values, and Holm-adjusted p-values from `attempt_metrics.csv`; maximum absolute discrepancies versus the saved analysis were 0 for means, 2.12e-22 for raw p-values, and 1.70e-21 for Holm-adjusted p-values. Because all 18 tests were evaluable, the separate Holm-family hardening change for unavailable tests does not alter this run's result.

## Interpretation

1. **Negative curvature (\(\gamma=-1.2\)):** estimated IPW substantially reduced IAE relative to catalogue-only at all three validation sizes. The mean paired difference became more negative as validation size increased.
2. **Zero curvature:** estimated degree-1 IPW improved on catalogue-only across all validation sizes. Degree-2 IPW was worse at n=100 but improved at n=300 and n=1,000.
3. **Positive curvature (\(\gamma=+1.2\)):** estimated IPW was worse than catalogue-only at n=100 and n=300 for both fitted degrees. At n=1,000, the paired intervals include zero for both degrees; the observed differences are near zero.
4. The sign reversal is not explained by a simple collapse in weight ESS: in the positive-curvature cells, median Kish ESS is high (approximately 959–994 across the displayed degree/validation-size cells), while the finite-validation density contrast can still be positive. This suggests that weight-concentration diagnostics alone do not diagnose the misspecification-related density error in this experiment. This is a descriptive observation, not a validated predictive diagnostic.
5. The known-inclusion IPW and parent-frame KDE arms are oracle diagnostics only. They are not deployable estimators when the full parent frame and true inclusion probabilities are unavailable.

These findings are conditional on this synthetic parent law, logistic inclusion mechanism, retention target, curvature range, validation design, KDE bandwidth, and finite grid. They do not establish a general ranking of IPW versus alternatives.

## Audit and limitations

The independent auditor reported:

- `status: PASS`
- 4,500/4,500 density estimates reconstructed.
- 300/300 replicate-by-curvature frame scenarios checked.
- 4,500 SUCCESS; zero failures/non-evaluable/not-run attempts.
- Exact reconstruction of the five declared RNG streams, frame redshifts, retention uniforms, validation permutation/ranks, DGP inclusion probabilities, Bernoulli labels, fitted logistic score/objective, KDE densities, IAE values, attempt factorial, and status/contrast summaries.

The auditor independently reconstructs data and calculations but shares SciPy's `gaussian_kde` primitive with production; this is not an independent KDE algorithm implementation. The validation frame assumes complete binary retention-label ascertainment. The design is synthetic and has no real-survey validation.

Most importantly, the experiment uses plug-in fitted inclusion probabilities for density recovery. It does **not** propagate inclusion-coefficient uncertainty into the H0 likelihood. The inverse-information covariance from the logistic fit is not itself a Bayesian posterior. No H0 posterior or coverage claim should be inferred from this report.

## Next research gate

The next stage should be a separate, explicitly modular uncertainty-propagation experiment:

1. Construct a proper validation-data posterior for inclusion coefficients, with a declared prior and checks for separation/weak identification; do not relabel the inverse-information covariance as a posterior.
2. For each coefficient draw, construct the IPW host density from the fixed observed analysis catalogue, preserving the same KDE/bandwidth conventions.
3. Evaluate the event likelihood and GW detection-selection factor on the same draw. Catalogue inclusion and GW detection selection must remain separate mechanisms.
4. Marginalize the likelihood over coefficient draws using log-sum-exp / log-mean-exp:
   \[
   L_{\mathrm{marg}}(H_0)
   =\int L(H_0;\theta)\,p(\theta\mid V)\,d\theta.
   \]
   This is a modular/cut construction unless a fully joint model is explicitly derived; do not average bootstrap likelihoods and call that Bayesian marginalization.
5. Validate numerical integration and Monte Carlo convergence, then conduct fixed-truth posterior calibration/coverage with enough independent replicates. Keep SBC, frequentist coverage, PIT, and grid convergence as separate diagnostics.

The FAC result is a necessary density-recovery gate, not the final scientific endpoint.
