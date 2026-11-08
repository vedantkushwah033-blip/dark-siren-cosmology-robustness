# Predeclared analysis plan: 100-replicate fixed-analysis-catalogue benchmark

**Frozen before the confirmatory run.** This protocol is informed by earlier exploratory disjoint-design pilots and the five-replicate FAC pilot; it is not an independent hypothesis-generation exercise.

## Confirmatory design

- 100 independent simulation replicates.
- Analysis parent frame: 4,000 objects per replicate/scenario.
- Independent validation frame: 1,000 objects.
- Nested validation sizes: 100, 300, 1,000.
- Curvature values: -1.2, 0, +1.2.
- Target retention: 0.25.
- Estimated inclusion models: degree 1 and degree 2.
- Primary Scott multiplier: 0.15.
- Redshift grid: 801 points over [0.01, 0.40].
- Arms: catalogue-only, known-inclusion IPW diagnostic, analysis-parent KDE diagnostic, estimated degree-1 IPW, estimated degree-2 IPW.
- Distinct SeedSequence streams for analysis redshifts, validation redshifts, validation ordering, analysis retention uniforms and validation retention uniforms.

The analysis catalogue is generated once per replicate/scenario and remains identical at all validation sizes. The validation frame is independent of the analysis frame. This design estimates the effect of validation information conditional on a fixed analysis catalogue; it does not estimate the operational trade-off in the disjoint finite-parent design.

## Primary estimand

For replicate (r), curvature (gamma), and fixed validation size (n=300), degree (d=1), and Scott multiplier (h=0.15), define

[
Delta_{r,gamma}
=mathrm{IAE}_{r,gamma,mathrm{estimated degree 1 IPW}}
-mathrm{IAE}_{r,gamma,mathrm{catalogue only}}.
]

The primary symmetric-curvature interaction is

[
Psi_r=Delta_{r,+1.2}-Delta_{r,-1.2},
qquad
Psi=E(Psi_r).
]

Positive (Psi) means the estimated-IPW-versus-catalogue-only density-recovery contrast is less favourable under positive than negative curvature, at the same absolute curvature magnitude, validation size, model degree, bandwidth and analysis-frame size.

This interaction is chosen because the (pm1.2) scenarios are symmetric curvature stresses and the same replicate-level random-number streams pair the scenarios. It is pilot-informed: earlier exploratory work indicated a possible sign-asymmetric misspecification regime. The confirmatory interpretation must therefore be explicitly limited to this predeclared synthetic contrast.

## Primary uncertainty calculation

- Construct (Psi_r) only when catalogue-only and estimated degree-1 IPW IAE are available for both (gamma=-1.2) and (+1.2) in replicate r.
- Report planned replicate count, complete paired replicate count, and all failure rates in the four contributing cells.
- Report the mean (Psi_r), its replicate-level standard error, a two-sided 95% percentile bootstrap interval with 10,000 replicate-ID resamples, and a paired t interval.
- Bootstrap whole replicate IDs, carrying both curvature scenarios and both estimator arms together. Do not resample galaxies, validation units, grid points, or fit coefficients as if they were independent replicates.
- Do not impute a numeric IAE for failed attempts. If any required arm fails, that replicate is absent from the complete primary interaction pair and the failure is separately reported.

## Secondary analysis

Report all 18 estimated-IPW-minus-catalogue-only contrasts at the primary bandwidth (3 curvatures × 3 validation sizes × 2 degrees), with paired replicate-level estimates, pointwise 95% intervals and complete-pair counts. Use Holm adjustment for the family of 18 secondary tests; distinguish adjusted p-values from pointwise intervals. Report sensitivity bandwidths only as secondary analyses if run.

Report conditional-success mean IAE separately from per-cell success probability and the full FAILED / NON_EVALUABLE / NOT_RUN ledger. Do not rank estimators solely by conditional-success means.

## Interpretation gates

The 100-replicate FAC experiment can support a synthetic density-recovery statement only. It cannot establish H0 posterior calibration, posterior coverage, real-survey validity, detector selection validity, or universal superiority of IPW. The known-inclusion arm is an oracle simulation diagnostic. The parent-frame KDE is a diagnostic unavailable when the full analysis parent frame is unknown.

Before any H0 propagation, separately verify numerical convergence, density recovery, coefficient-uncertainty treatment, selection-conditioned likelihood marginalization, and fixed-truth coverage/SBC. Catalogue inclusion and GW detection selection remain distinct processes.


## Independent RNG-stream audit gate (2026-11-08)

The fixed-analysis-catalogue auditor now reconstructs all five declared SeedSequence([seed, replicate, stream_id]) streams independently:

1. analysis-frame redshifts;
2. independent validation-frame redshifts;
3. validation-frame permutation/ranks;
4. analysis-frame retention uniforms;
5. validation-frame retention uniforms.

It checks the reconstructed values against the raw frame-unit artifact, in addition to checking that Bernoulli labels follow the stored uniforms and DGP probabilities. This guards against a common failure mode in which a raw uniform or frame coordinate is altered while the corresponding Bernoulli label and manifest checksum are also updated.

The regression suite includes an adversarial retention-uniform tampering test that refreshes the manifest hash and still requires the independent auditor to reject the artifact. The confirmatory artifact is not considered fully audited until this new check has passed on the full raw output.
