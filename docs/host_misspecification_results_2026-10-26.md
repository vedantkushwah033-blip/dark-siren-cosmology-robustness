> **Revalidation note — 2026-10-27:** This experiment has now been rerun from corrected revision `0ac70cba1708be82c7b994145f4b0507101afa05`; its independent reconstruction audit passed. The new 4,000-row raw artifact agrees with the historical artifact to floating-point precision (maximum absolute difference in posterior summaries about `1.5e-13`; all scenario-level metrics agree at the displayed precision). The historical run records source revision `5caad04955d735583f77d2c56e455aa28ad43b40`, which predates the later transient log-space double-scaling regression, so this particular report's numbers are reproduced rather than superseded. This does **not** validate the other hidden-host experiments; see the [impact audit](missing_host_mixture_likelihood_impact_audit_2026-10-27.md).

# Host-population misspecification stress-test results

**Run date:** 2026-10-26  
**Source revision:** `5caad04955d735583f77d2c56e455aa28ad43b40`  
**Workflow:** [Host-weight and omitted-host misspecification audit](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40278173997)  
**Artifact:** [Raw per-replicate records and scenario summaries](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40278173997/artifacts/12631380406)  
**Artifact SHA-256:** `c627c12bcad2c41db75769f3566dd2a5cc1b1bd3e9988332760ac901462f094c`

## Corrected-revision rerun and artifact provenance

- **Corrected source revision:** `0ac70cba1708be82c7b994145f4b0507101afa05`
- **Workflow run:** [40416457719 — passed](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40416457719)
- **Artifact:** [corrected raw records and summaries](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40416457719/artifacts/12694878568)
- **Artifact ZIP SHA-256:** `e56cbeebe4b353eab78bfb533117142fc138f163e0719575264ee416ba910c20`
- **Raw `per_run.csv` SHA-256:** `829fc5cb6f61a33d8ec78eb88aae000b5f034131f7b4aa034bfa956dc985cb5b`
- **Raw `summary.csv` SHA-256:** `ff1148ef0938564810249237ecf41ca7dd6deacc63b200c8f6e6f5a8be2dd8ce`
- **Independent reconstruction:** passed for all 4,000 seeded records and four scenario summaries.

The historical artifact SHA-256 was `c627c12bcad2c41db75769f3566dd2a5cc1b1bd3e9988332760ac901462f094c`. Comparing both downloaded artifacts shows the maximum absolute difference across per-run posterior interval summaries is below `1.5e-13`; differences are floating-point rounding only. The table below is therefore reproduced by the corrected revision at the precision reported.

## Design and verification

The fixed synthetic truth is (H_0=70) km s\(^{-1}\) Mpc\(^{-1}\), with a uniform prior on [50, 90], two catalogue redshifts 0.02 and 0.04, and a hidden component at redshift 0.03 in the omitted-host experiments. The noise scale is fixed at (0.10 c(0.03)/70). Non-positive simulated distances are redrawn, and the inference likelihood conditions each Gaussian component on positive distance.

The experiment used seed `20261027` and 1,000 replicates per scenario: 4,000 posterior fits total. The independent auditor regenerated the seeded data, reconstructed all 4,000 posterior summaries, recalculated the coverage indicators and aggregate metrics, and verified the Wilson interval endpoints. The workflow passed all steps.

This is an independent *implementation* audit of the declared toy model, not an independent validation of its scientific assumptions or real gravitational-wave/catalogue data.

## Results

Errors and widths are in km s\(^{-1}\) Mpc\(^{-1}\). Coverage intervals are pointwise two-sided 95% Wilson score intervals for each empirical coverage proportion.

| Scenario | Mean median error | RMSE of median | 68% coverage (95% Wilson interval) | 95% coverage (95% Wilson interval) | Mean 68% width | Mean 95% width |
|---|---:|---:|---:|---:|---:|---:|
| Correct catalogue weights | +1.257 | 6.648 | 0.730 (0.702–0.757) | 0.987 (0.978–0.992) | 15.361 | 27.838 |
| Reversed catalogue weights | +1.251 | 6.654 | 0.743 (0.715–0.769) | 0.990 (0.982–0.995) | 15.451 | 28.126 |
| Hidden host omitted | +1.284 | 8.562 | 0.664 (0.634–0.693) | 0.993 (0.986–0.997) | 16.993 | 29.238 |
| Hidden host modelled | +1.527 | 4.718 | 0.885 (0.864–0.903) | 0.988 (0.979–0.993) | 22.175 | 33.051 |

## Interpretation, with limits

1. **Omitting the hidden host is associated with worse median-estimator RMSE in this design.** RMSE rises from 4.718 when the hidden component is modelled to 8.562 when it is omitted. This is a descriptive contrast for this fixed truth, prior, host support, noise scale, seed, and model; it is not a universal estimate of the cost of catalogue incompleteness.
2. **Modelling the hidden component substantially broadens the intervals.** Mean 68% interval width rises from 16.993 to 22.175, and empirical 68% coverage rises from 0.664 to 0.885. That pattern is consistent with a more conservative posterior under this setup, but it does not alone establish generally calibrated uncertainty.
3. **Reversing the two catalogue weights has little apparent effect on these aggregate metrics here.** The two redshifts and likelihood structure may make this setup weakly informative about the weight reversal. This should not be interpreted as evidence that host weights are generally unimportant.
4. **The empirical 95% coverage is high in all four arms (0.987–0.993).** The pointwise Wilson intervals lie above the nominal 0.95 level for this simulation. Fixed-parameter frequentist coverage is not guaranteed to equal nominal Bayesian credible level, so this is a property to investigate—not by itself a proof of model failure.
5. **The 68% results differ across scenarios**, but the Wilson intervals are marginal intervals, not paired intervals for scenario differences. The scenarios were generated sequentially from one seeded RNG stream rather than as matched replicate-level contrasts. Avoid treating interval overlap or non-overlap as a formal test of differences.

The estimates summarize one predeclared simulation design with one fixed (H_0) truth and one seed. They do not cover uncertainty in the assumed host population, selection process, redshift-distance relation, noise model, prior choice, or numerical approximation. Next validation should vary the injected (H_0), hidden-host fraction/redshift, noise scale, and prior support, with a predeclared paired design where comparisons are intended, and should report uncertainty for scenario contrasts.

## Reproduction

```bash
python scripts/run_host_misspecification_audit.py --replicates 1000 --seed 20261027 --output-dir results/host-misspecification-audit
python scripts/audit_host_misspecification_audit.py --replicates 1000 --seed 20261027 --input-dir results/host-misspecification-audit
```

See the [independent audit protocol](host_misspecification_independent_audit_protocol_2026-10-26.md) for the audit's precise scope and limitations.
