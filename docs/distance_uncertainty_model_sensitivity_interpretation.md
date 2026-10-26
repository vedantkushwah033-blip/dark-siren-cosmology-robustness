# Interpretation limits for the distance-uncertainty sensitivity experiment

## Scope

This is a synthetic, one-known-host, low-redshift experiment with
\(d=c z/H_0\), a uniform prior on \(H_0\in[50,90]\), and Gaussian distance
noise. It is a controlled model comparison, not an LVK analysis, not a
realistic gravitational-wave distance posterior, and not a measurement of
\(H_0\).

## The three arms are not interchangeable

- **A — fixed realized sigma.** The generator first uses the *true* \(H_0\)
  to calculate the true distance, then sets \(\sigma_A=f d_{\rm true}\).
  The likelihood holds this realized sigma fixed while evaluating candidate
  \(H_0\). This is a conditional benchmark that is useful for isolating the
  likelihood-shape question. Because the generated sigma depends on the
  simulation truth, however, this arm is an oracle-style benchmark: it must
  not be presented as an uncertainty prescription that an analyst could
  necessarily construct when the true \(H_0\) is unknown.
- **B — candidate-dependent fractional sigma.** The observation is generated
  using the same realization as A, but the likelihood evaluates
  \(\sigma(H_0)=f c z/H_0\) for each candidate and includes the Gaussian
  normalization term \(-\log\sigma(H_0)\). A-versus-B is therefore a paired
  comparison of two likelihood treatments of the same observation, under the
  chosen fractional-noise generator. It is not a universal ranking of
  uncertainty models.
- **C — fixed reference sigma.** The noise is generated using
  \(\sigma_C=f c z/70\) and then held fixed in the likelihood. Except at
  true \(H_0=70\), its data-generating noise scale differs from A/B. Comparisons
  involving C therefore combine a likelihood assumption with a change in the
  generator; they are not a pure likelihood-only contrast.

## Reading the summaries

1. Coverage is the repeated-simulation fraction of central posterior intervals
   containing the fixed generating value. With 300 replicates per cell, it is
   noisy; always report the Wilson interval and avoid ranking models on small
   differences without uncertainty estimates.
2. Bias and RMSE are descriptive properties of this specified simulation,
   prior, and finite grid. They do not establish performance under a realistic
   population of gravitational-wave events.
3. Truths at 50 and 90 lie exactly on the prior boundaries. Their coverage,
   posterior summaries, and errors are boundary-sensitive and must be shown
   separately from interior truths.
4. The experiment uses a uniform prior and a low-redshift approximation
   \(d=c z/H_0\). It omits realistic GW likelihoods, redshift uncertainty,
   selection effects, catalogue incompleteness, host-population uncertainty,
   and cosmological-model complexity.
5. The CI experiment is only a smoke test with two replicates per cell. Its
   output verifies execution and auditability, not scientific conclusions.
   The default 300-replicate experiment must be run and its full artifacts
   independently audited before reporting numerical findings.

## Appropriate conclusion

The defensible question is: *How do these explicitly defined uncertainty
assumptions alter posterior summaries and repeated-simulation coverage in this
particular controlled model?* The experiment cannot establish that one arm is
generally correct for dark-siren cosmology or that the observed differences
would carry over to real gravitational-wave data.
