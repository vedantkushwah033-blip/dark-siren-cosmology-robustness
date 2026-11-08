# Follow-up protocol: distance-uncertainty model sensitivity

**Status:** prespecified follow-up protocol; no new numerical results are claimed here.  
**Parent experiment:** known-host calibration surface, pinned source simulator commit `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`.  
**Purpose:** test whether the calibration pattern is sensitive to the way distance uncertainty is generated and used in inference. This remains a synthetic, one-host diagnostic, not an observational gravitational-wave analysis.

## Why this follow-up is needed

The baseline generator draws

\[
d_{\mathrm{obs}}\sim\mathcal N(d_\star,\sigma_\star^2),
\qquad d_\star=cz/H_{0,\star},\quad \sigma_\star=f d_\star,
\]

where the fractional uncertainty \(f\) is fixed for each experimental arm. The baseline inference supplies the realized \(\sigma_\star\) to a Gaussian likelihood and holds it constant over the candidate values of \(H_0\):

\[
\mathcal L_{\mathrm{fixed}\,\sigma}(H_0)\propto
\exp\!\left[-\frac{(d_{\mathrm{obs}}-cz/H_0)^2}{2\sigma_\star^2}\right].
\]

That is a conditional likelihood for a known, fixed uncertainty. But in this simulator, the uncertainty is itself constructed as a deterministic function of the true distance. If the intended measurement model says the standard deviation varies with the candidate distance, the corresponding likelihood instead has \(\sigma(H_0)=fcz/H_0\), including its normalization term:

\[
\mathcal L_{\mathrm{varying}\,\sigma}(H_0)=
\frac{1}{\sqrt{2\pi}\,\sigma(H_0)}
\exp\!\left[-\frac{(d_{\mathrm{obs}}-cz/H_0)^2}{2\sigma(H_0)^2}\right].
\]

These are different statistical models. Neither should be labelled universally correct without stating what the uncertainty represents. This follow-up compares them under explicit synthetic assumptions; it does not retroactively invalidate or overwrite the baseline.

## Prespecified experimental arms

Use the same truths \(H_{0,\star}\in\{50,55,60,65,70,75,80,85,90\}\), fractional-noise settings \(f\in\{0.05,0.10,0.15,0.20\}\), uniform prior \([50,90]\), and 1,601-point posterior grid. Retain the original 300-replicate run for a direct audit, then use at least 1,000 replicates per cell for the confirmatory comparison if computationally feasible. Treat boundary truths 50 and 90 separately.

### Arm A — baseline conditional fixed uncertainty

Reproduce the current implementation exactly: generate \(d_{\mathrm{obs}}\) with \(\sigma_\star=f d_\star\), then hold the realized \(\sigma_\star\) fixed over the inference grid. This is the compatibility/reference arm.

### Arm B — parameter-dependent uncertainty likelihood

Generate with the same data-generating rule as Arm A. For each candidate \(H_0\), use \(\mu(H_0)=cz/H_0\) and \(\sigma(H_0)=f\mu(H_0)\) in the normalized Gaussian likelihood. Include the \(-\log \sigma(H_0)\) term; omitting it would define a third, different likelihood. Record any numerical stabilization and posterior-grid normalization.

Interpretation: Arm B asks what changes when the analysis explicitly treats the scale of the measurement uncertainty as dependent on the candidate model distance. It is a sensitivity model, not a claim that all real distance posteriors are Gaussian or have this exact scaling.

### Arm C — fixed-absolute-uncertainty control

Define a reference distance scale before simulation, using \(H_{0,\mathrm{ref}}=70\) and the host's redshift: \(\sigma_{\mathrm{ref}}=f c z/H_{0,\mathrm{ref}}\). Generate \(d_{\mathrm{obs}}\sim\mathcal N(cz/H_{0,\star},\sigma_{\mathrm{ref}}^2)\) and use that same \(\sigma_{\mathrm{ref}}\) in inference. This control removes the truth-dependent scale of the noise while preserving a stated noise level relative to a fixed reference model. Because its generator differs from Arms A/B away from the reference truth, comparisons across arms are sensitivity comparisons, not a pure within-dataset likelihood swap.

## Randomness and pairing

- Keep host redshift and standardized Gaussian deviate paired within each truth/replicate across arms wherever the data-generating assumptions allow it.
- For Arms A and B, use the exact same observed distance values and uncertainties so their contrast isolates the likelihood treatment.
- Arm C uses the same standardized deviate but a separately defined absolute uncertainty; document this difference.
- Do not compare individual posterior-grid points as if they were independent observations. Replicate is the resampling unit.
- Save the standard normal deviate explicitly in raw records. Do not rely only on rounded distances to verify pairing.

## Primary outcomes

For each truth/noise/arm cell, report:
1. Nominal 68% and 95% equal-tailed interval coverage, success counts, and Wilson intervals.
2. Bias and RMSE of posterior median and posterior mean.
3. Mean 68% and 95% interval widths.
4. Posterior-median pull summaries, with a clear statement that posterior SD-based pulls are descriptive rather than a universal calibration test.
5. Prior-boundary diagnostics, separately from interior-truth results.
6. Paired replicate-level differences between Arms A and B for coverage indicator, absolute error, squared error, and interval width, with replicate-bootstrap uncertainty intervals.

Do not average boundary and interior truths into one headline calibration number. Do not choose the best arm after seeing results without labelling the comparison exploratory or validating it on a separate simulation seed.

## Required checks before interpretation

- Unit tests verify the normalized Arm B likelihood, including the \(-\log\sigma(H_0)\) term, and confirm it reduces to the intended Gaussian form when uncertainty is constant.
- A recovery test uses a very small fractional error and an interior truth to check that posterior concentration is sensible; this is a diagnostic, not proof of calibration.
- An analytic or high-resolution numerical reference test compares grid summaries against an independently implemented calculation.
- Check all values are finite, all standard deviations positive, interval bounds nested and within prior support, and the correct number of unique truth/noise/replicate/arm records exists.
- Audit the stored standardized noise: \((d_{\mathrm{obs}}-d_\star)/\sigma_\star\) for Arms A/B and \((d_{\mathrm{obs}}-d_\star)/\sigma_{\mathrm{ref}}\) for Arm C should match the stored deviate up to numerical precision.
- Run tests and the full supported Python matrix. Store commit SHA, dependency versions, artifact SHA-256, seeds, failed-run counts, and the exact likelihood equations in the generated report.
- Report differences and Monte Carlo uncertainty, not only point estimates. A result is not evidence of a real-world Hubble-constant bias.

## Decision rules

- If Arms A and B differ materially, conclude only that fixed- versus parameter-dependent uncertainty treatment matters under this specified synthetic model.
- If Arms A and B agree, conclude only that this tested contrast had limited impact over this grid and noise range; do not infer general robustness.
- If coverage remains heterogeneous across truths in both arms, investigate nonlinear distance-to-parameter mapping, bounded-prior effects, grid/discretization, and other model assumptions in separate controlled tests.
- If results change materially with grid resolution or seed, treat the original 300-replicate surface as insufficient for a stable quantitative claim.
- Keep the original baseline artifact immutable; publish the sensitivity experiment as a separate versioned result.
