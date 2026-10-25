# Prior-measure sensitivity: protocol

**Question.** In a two-host low-redshift mixture, how does the prior measure on H₀ change the integrated posterior mass assigned to each host, and how does that interact with reversing catalogue weights?

This is a synthetic component-evidence diagnostic, not a real gravitational-wave analysis.

## Prespecified design

- Injected H₀: 70 km s⁻¹ Mpc⁻¹.
- Host redshifts: z₁ = 0.02 and z₂ ∈ {0.025, 0.04, 0.06}.
- Generative host weights: (0.7, 0.3); reversed inference weights: (0.3, 0.7).
- Fractional distance uncertainty: {3%, 10%, 20%}.
- Prior support: narrow [50,90] and broad [20,220] km s⁻¹ Mpc⁻¹.
- Prior measures: uniform in H₀ and log-uniform in H₀ (density proportional to 1/H₀).
- Grid: 10–250 km s⁻¹ Mpc⁻¹ at 0.025 spacing.
- Default: 1,000 observations per redshift/noise cell; base seed 20261029; 9,000 raw observations, 36 summary rows, 8 host-mixture posterior component-mass calculations per observation.
- Each synthetic observation, host identity, and sigma is shared across all four combinations of prior support and both prior measures, and across correct/reversed weights.

## Estimand

For each raw observation, calculate each host component's integrated contribution

$$
Z_i = w_i \int_{H_{\min}}^{H_{\max}}
p(d_{\rm obs}\mid H_0,z_i)\,\pi(H_0)\,dH_0,
$$

then normalize the two contributions to component posterior masses $q_i=Z_i/(Z_1+Z_2)$. The main contrast is the paired replicate-level difference in host-2 posterior mass under reversed versus correct weights. Report its mean and deterministic percentile-bootstrap 95% interval, resampling replicate observations within each design cell.

The uniform prior has constant density within the support. The log-uniform prior has density proportional to 1/H₀; the prior's normalizing constant is common to both host components in a given arm and cancels in their normalized masses.

## Validation gates

- Every raw record must match the fixed truth, design grid, host-index/redshift mapping, and sigma formula.
- The independent auditor reconstructs each host component mass directly from raw observation, redshifts, sigma, weights, prior support, and prior measure.
- Component masses must be finite, lie in [0,1], and sum to one.
- Exactly 9,000 raw rows and 36 summaries are expected for the default design.
- All paired means and bootstrap intervals must be independently reconstructed; raw and summary hashes must match metadata.
- Results must distinguish prior-measure sensitivity from prior-support truncation and must not claim that either prior is universally preferred.

## Interpretation limits

This is a fixed-truth, two-host, low-redshift Gaussian-distance toy model. It does not model galaxy-catalogue selection, realistic redshift errors, gravitational-wave detector selection, a population of events, or the full cosmological distance-redshift relation. A prior-measure contrast here is a diagnostic of this model, not a general result for real dark-siren H₀ inference.
