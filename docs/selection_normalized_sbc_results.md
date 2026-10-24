# Toy selection-normalized posterior and SBC audit

**Status:** initial selected-data SBC run completed; the independent adaptive-quadrature regression test was added afterward and its CI result must be checked separately.  
**SBC workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40068645358  
**Scope:** toy noisy-distance threshold and fixed two-host mixture only; not a realistic GW detector selection model.

## Question

Does the selection-normalized posterior behave consistently with the specified toy data-generating process when the data are generated conditional on detection?

## Model

The candidate host redshifts are 0.02 and 0.04, with host weights 0.70 and 0.30. At each fixed injected H0, a host is drawn from those weights and a Gaussian noisy distance is drawn; the host-distance draw is repeated until the observed distance satisfies d_obs < 180 Mpc. This samples observations conditional on detection **without rejecting or reweighting the injected H0 prior draws**.

The detection efficiency is

$$
\beta(H_0)=\sum_j w_j\Phi\!\left(\frac{180\,\mathrm{Mpc}-cz_j/H_0}{20\,\mathrm{Mpc}}\right).
$$

Inference uses the host-mixture likelihood divided by beta(H0) and a uniform prior on H0 in [50,90] km s^-1 Mpc^-1. The threshold is an intentionally simplified toy rule, not a model of matched-filter SNR or the LVK search pipeline.

## Design

- 2,000 prior-predictive replicates.
- Fixed random seed: 20261027.
- 4,001-point computational grid.
- Gaussian distance uncertainty: 20 Mpc.
- Selection threshold: 180 Mpc.
- Diagnostics: posterior CDF-rank mean/variance, Kolmogorov–Smirnov test against Uniform(0,1), and a ten-decile chi-square diagnostic.

## Results

| Diagnostic | Result |
|---|---:|
| Replicates | 2,000 |
| Mean posterior CDF rank | 0.50093 |
| Rank variance | 0.08390 |
| KS statistic | 0.01117 |
| KS p-value | 0.96189 |
| Decile chi-square p-value | 0.94681 |
| Minimum beta(H0) on prior grid | 0.69948 |
| Maximum beta(H0) on prior grid | 0.99709 |

At the conventional 5% level, neither uniformity diagnostic rejects uniform ranks in this seeded run. The rank mean and variance are close to the Uniform(0,1) reference values of 0.5 and 1/12. The efficiency range shows that the selection normalization is not constant over the prior, so it is materially represented in this toy setup.

These diagnostics are evidence consistent with correct computation under the specified model, not proof of exact calibration. The KS and decile diagnostics use the same ranks and are not independent validations.

## Numerical validation

The analytic selection tests compare beta(H0) to the weighted sum of normal CDF values, verify that efficiency approaches one for a very high threshold, and check the no-selection limit of the posterior. An additional test compares selected-posterior quantiles against an independently written adaptive-quadrature calculation using scipy.integrate.quad and scipy.optimize.brentq.

The focused selection-normalization unit tests passed before the adaptive-quadrature regression test was added: [workflow run 40068498540](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40068498540). The SBC workflow passed: [run 40068645358](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40068645358). The new adaptive-quadrature test and the full test suite must pass on the latest commit before treating the implementation checkpoint as fully validated.

## Limitations and next steps

1. The threshold is applied to noisy distance, not a physical GW detection statistic.
2. The redshifts and host weights are fixed and known; there is no catalogue incompleteness in this experiment.
3. There are no realistic detector sensitivity, duty cycle, waveform, antenna-pattern, or observing-run effects.
4. This prior-predictive SBC run is not a fixed-truth coverage study and does not prove frequentist coverage at each H0.
5. Next compare selected-data inference with and without the beta(H0) correction under paired mock observations, then cross detection selection with catalogue incompleteness only after both branches are independently checked.

## Reproducibility

- Implementation: src/dark_siren_robustness/selection.py
- Simulation: scripts/run_selection_sbc_audit.py
- Tests: tests/test_selection_normalization.py, tests/test_selection_sbc_audit.py
- Workflow: .github/workflows/selection-sbc-audit.yml
- Artifact: selection-sbc-audit-results (per-run CSV and summary CSV; retention is time-limited).
