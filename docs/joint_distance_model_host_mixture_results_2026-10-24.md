# Joint distance-model × host-mixture experiment: results and independent audit

**Status:** completed synthetic diagnostic; raw outputs independently checked on 2026-10-24.  
**Workflow:** [GitHub Actions run 40116696053](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40116696053)  
**Artifact:** `joint-distance-host-mixture-results` (artifact ID `12569675504`; SHA-256 `8848c80fc9448466edfddc6ba7d310f6b3c7796317fab436563abb225605b0ae`).  
**Scope warning:** synthetic model-comparison diagnostic only; not an analysis of real gravitational-wave events or a measurement of (H_0).

## Design

The simulation crosses three true values ((H_0=60,70,80)), three redshift pairs ((0.02,0.04)), ((0.02,0.10)), and ((0.10,0.40)), and three intrinsic host-weight pairs ((0.5,0.5)), ((0.7,0.3)), and ((0.9,0.1)). Each of the 27 cells has 300 simulated observations. Each observation is analyzed using four arms:

1. flat-(Lambda)CDM distance relation with two hosts (`lcdm_two_host`);
2. low-redshift (d_L\approx cz/H_0) relation with two hosts (`lowz_two_host`);
3. flat-(Lambda)CDM relation with a collapsed effective-redshift host (`lcdm_collapsed`);
4. low-redshift relation with a collapsed host (`lowz_collapsed`).

The truth generator uses flat (Lambda)CDM with (Omega_m=0.3), selects one host according to the configured weights for each observation, and adds Gaussian distance noise with (sigma=max(20,mathrm{Mpc},0.05\times) the weighted mean distance). Inference uses a uniform (H_0) prior on [40, 100] and a 4,001-point grid on [25, 125]. The baseline deliberately omits detection truncation and selection correction to isolate the distance-model and host-representation effects.

## Output and independent integrity audit

- Raw inference rows: **32,400** (8,100 observations × four arms).
- Summary rows: **108** (27 cells × four arms).
- Paired-contrast rows: **135** (27 cells × five contrasts).
- Every cell/arm contains exactly 300 replicates.
- All numeric raw fields are finite; all posterior medians fall inside the prior; 68% and 95% interval ordering checks passed.
- Independently recomputed every reported summary from `per_run.csv`: maximum absolute discrepancy was approximately (1.1\times10^{-14}), consistent with floating-point rounding.
- Independently recomputed all 135 paired-contrast means from matched replicate-level absolute errors: maximum discrepancy was approximately (3.6\times10^{-15}).
- Independently regenerated all 5,000-draw bootstrap intervals per contrast using the documented deterministic seed rule: maximum endpoint discrepancies were below (8\times10^{-15}).
- These checks validate the arithmetic and internal consistency of the exported artifact; they do not validate the scientific adequacy of the simplified generative model.

## Aggregate results

Across all 8,100 simulated observations (averaged over cells):

| Inference arm | Mean absolute error in (H_0) | 68% interval coverage | 95% interval coverage |
|---|---:|---:|---:|
| `lcdm_two_host` | 6.238 | 75.1% | 96.9% |
| `lowz_two_host` | 7.946 | 62.4% | 83.0% |
| `lcdm_collapsed` | 21.675 | 13.5% | 28.3% |
| `lowz_collapsed` | 21.044 | 14.5% | 30.0% |

The error units are km s⁻¹ Mpc⁻¹ for (H_0). Coverage fractions are empirical repeated-sampling proportions, not probabilities that any particular interval contains the truth. The global aggregates combine distinct truths and configurations, so they should not be read as universal performance guarantees.

## Main findings

### 1. Collapsing a two-host mixture into one effective-redshift host was the dominant failure in this design

The collapsed arms had overall MAE around 21 km s⁻¹ Mpc⁻¹ and very poor interval coverage. By contrast, the matched two-host flat-(Lambda)CDM arm had MAE about 6.24 and 95% coverage about 96.9%. Averaged across the 27 configurations, the paired increase in absolute error from host collapse was 15.44 for the flat-(Lambda)CDM relation and 13.10 for the low-redshift relation. This is evidence for the importance of representing host ambiguity in this particular synthetic setup—not proof that any specific real-catalogue weighting method is correct.

### 2. Distance-model mismatch mattered most for the widely separated redshift pair

For hosts at (z=(0.02,0.04)), average MAE was 7.76 for the two-host (Lambda)CDM arm and 7.79 for the two-host low-redshift arm. For (z=(0.02,0.10)), the corresponding values were 6.76 and 7.39. For the more widely separated (z=(0.10,0.40)) pair, they were 4.20 and 8.66; 95% coverage fell from 95.6% under the matched (Lambda)CDM arm to 63.5% under the low-redshift approximation. This pattern is consistent with the low-redshift approximation becoming inadequate at higher redshift in this setup.

### 3. The effects interact; no universal correction ranking should be claimed

The paired difference-in-differences estimand was defined as
[
(\mathrm{MAE}_{\text{lowz,collapsed}}-\mathrm{MAE}_{\text{LCDM,collapsed}})
-
(\mathrm{MAE}_{\text{lowz,two-host}}-\mathrm{MAE}_{\text{LCDM,two-host}}).
]
Its average across the 27 configurations was -2.34, but individual configurations ranged from about -8.31 to +0.40. Thus the joint effect is configuration-dependent. A single pooled interaction average should not replace per-cell contrasts and their uncertainty intervals.

## Limitations and next checks

- One simulated observation is generated from one selected host at a time; the likelihood marginalizes over the configured host mixture. This is a simplified host-ambiguity model.
- There is no gravitational-wave detector likelihood, event-level selection normalization, galaxy-catalogue incompleteness process, peculiar-velocity model, or real survey data in this experiment.
- (Omega_m) is fixed at 0.3; cosmological parameter uncertainty is not marginalized.
- The grid and noise model are controlled simplifications. The result cannot establish observational bias or quantify a real (H_0) systematic.
- The next useful checks are posterior-grid refinement, an independent likelihood/distance implementation, and an extension that introduces selection and catalogue incompleteness separately after the baseline interaction has been preserved.

**Reproducibility:** workflow artifact contains `per_run.csv`, `summary.csv`, `paired_contrasts.csv`, and `experiment_manifest.json`. The workflow succeeded, and the exported statistics and bootstrap intervals passed the independent checks described above.
