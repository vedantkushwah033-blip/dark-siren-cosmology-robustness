# Prior-measure sensitivity: independently audited results (2026-10-25)

**Status: PASS.** The 1,000-replicate experiment completed and its independent audit reconstructed all host component masses, paired summary estimates, and bootstrap confidence intervals from the raw observations.

- **Workflow:** [successful run 40171791086](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40171791086)
- **Artifact:** [raw observations, summaries, metadata, and audit report](https://api.github.com/repos/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/artifacts/12591150145/zip)
- **Protocol:** [prior-measure sensitivity protocol](prior_measure_sensitivity_protocol.md)
- **Audit:** PASS; 9,000 raw observations, 36 summary rows, raw and summary hashes verified.
- **Seed:** 20261029. Each of the nine redshift/noise cells has 1,000 generated observations, each evaluated under correct and reversed host weights, two prior supports, and two prior measures.

## Main estimand

For each observation, the integrated evidence for host component $i$ is

$$
Z_i = w_i \int_{H_{\min}}^{H_{\max}}
p(d_{\rm obs}\mid H_0,z_i)\,\pi(H_0)\,dH_0,
$$

and the normalized host-component posterior mass is $q_i=Z_i/(Z_1+Z_2)$. The table below reports the mean paired difference in host-2 mass, reversed weights minus correct weights. Intervals are deterministic percentile-bootstrap 95% intervals based on 4,000 paired resamples.

## All design-cell contrasts

| $z_2$ | Distance noise | Broad log-uniform | Broad uniform | Narrow log-uniform | Narrow uniform |
|---:|---:|---:|---:|---:|---:|
| 0.025 | 3% | +0.4000 [0.4000, 0.4000] | +0.3958 [0.3958, 0.3958] | +0.3719 [0.3681, 0.3755] | +0.3764 [0.3729, 0.3797] |
| 0.025 | 10% | +0.4000 [0.4000, 0.4000] | +0.3959 [0.3959, 0.3960] | +0.3561 [0.3521, 0.3600] | +0.3590 [0.3551, 0.3626] |
| 0.025 | 20% | +0.3994 [0.3993, 0.3995] | +0.3975 [0.3974, 0.3976] | +0.3682 [0.3658, 0.3706] | +0.3700 [0.3676, 0.3722] |
| 0.040 | 3% | +0.4000 [0.4000, 0.4000] | +0.3620 [0.3620, 0.3620] | +0.0000 [0.0000, 0.0000] | +0.0000 [0.0000, 0.0000] |
| 0.040 | 10% | +0.3984 [0.3978, 0.3989] | +0.3664 [0.3659, 0.3670] | +0.0056 [0.0043, 0.0069] | +0.0069 [0.0054, 0.0086] |
| 0.040 | 20% | +0.3874 [0.3854, 0.3892] | +0.3734 [0.3719, 0.3748] | +0.1116 [0.1043, 0.1189] | +0.1149 [0.1074, 0.1229] |
| 0.060 | 3% | +0.3821 [0.3792, 0.3847] | +0.3351 [0.3331, 0.3370] | +0.0000 [0.0000, 0.0000] | +0.0000 [0.0000, 0.0000] |
| 0.060 | 10% | +0.3502 [0.3446, 0.3553] | +0.3338 [0.3299, 0.3375] | +0.0000 [0.0000, 0.0000] | +0.0000 [0.0000, 0.0000] |
| 0.060 | 20% | +0.3332 [0.3273, 0.3387] | +0.3205 [0.3157, 0.3251] | +0.0025 [0.0015, 0.0036] | +0.0032 [0.0020, 0.0046] |

Each cell shows mean difference [95% bootstrap CI]. “Broad” denotes prior support [20,220]; “narrow” denotes [50,90] km s⁻¹ Mpc⁻¹. Uniform means constant prior density in H₀; log-uniform means density proportional to 1/H₀. The intervals that collapse to a point reflect effectively deterministic contrasts at the displayed precision in this toy design, not uncertainty-free inference in general.

## Findings

### 1. Prior measure materially changes host-component posterior mass

At $z_2=0.06$ and 10% distance noise under the broad prior, the mean host-2 posterior mass under correct weights was 0.4174 for a uniform-$H_0$ prior and 0.2166 for a log-uniform prior. Under reversed weights, those means were 0.7512 and 0.5668, respectively. Thus the choice of prior measure changes the component-mass allocation substantially even with the same observations, redshifts, support, and catalogue weights.

This is consistent with a local narrow-likelihood argument: for a component with mode well inside support, the likelihood's width in $H_0$ is approximately proportional to $z_i$. Integrating against a uniform-$H_0$ density therefore gives component evidence approximately proportional to $w_i z_i$. A log-uniform density contributes a factor $1/H_0$, approximately cancelling that width factor near the mode, so component evidence is approximately proportional to $w_i$. These are approximations; finite noise, support boundaries, and overlapping components matter.

### 2. Prior support can mask weight sensitivity

For $z_2=0.04$ or $0.06$, the narrow prior excludes one or both alternative-host modes for many observations. In those cells, the paired host-2-mass contrast is zero or small, even though broad-support contrasts remain large. A near-zero effect under the narrow prior must not be interpreted as general evidence that catalogue weights do not matter: it can be caused by truncating a component from the allowed parameter space.

### 3. Reversing weights changes host-2 mass, but mass is not H₀ calibration

The experiment establishes how host-component mass shifts in this specified mixture. It does **not** establish the direction or magnitude of bias, RMSE, or coverage in a full cosmological posterior. A component-mass diagnostic is useful for understanding the mechanism but cannot substitute for posterior-level calibration across a prespecified multi-truth design.

## Independent audit and reproducibility

The workflow's independent audit:
- validated design keys, host-index/redshift mapping, injected truth, sigma, and complete replicate blocks;
- independently reintegrated each host component likelihood under each prior support and prior measure;
- checked finite component masses and normalization;
- reconstructed all 36 paired summary means and 4,000-resample bootstrap intervals;
- verified SHA-256 hashes for the raw and summary CSV files.

The recorded audit report says `PASS`. The raw CSV and `run_metadata.json` are in the linked artifact.

## Limitations and next step

This is a fixed-$H_0$, two-host, low-redshift Gaussian-distance toy model. It omits catalogue selection, realistic redshift uncertainties, detector selection, multiple-event population inference, and the full cosmological distance-redshift relation. The next defensible extension is a prespecified multi-truth experiment that reports posterior-level bias, RMSE, interval width, and coverage alongside these component-mass diagnostics. Any comparison of prior measures should retain the same synthetic data and clearly distinguish sensitivity to prior measure from sensitivity to prior support.
