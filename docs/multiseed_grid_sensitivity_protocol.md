# Multi-seed paired grid-sensitivity protocol

## Objective

Extend the single-seed grid-sensitivity screening study to check whether numerical PIT differences and aggregate calibration summaries are stable across independently seeded synthetic datasets. This is a numerical robustness exercise under the repository's declared synthetic model, not validation against astrophysical data.

## Design and execution

The runner `scripts/run_multiseed_grid_sensitivity.py` defaults to three seeds (20261025, 20261026, 20261027), 300 replicates per seed and grid, and posterior grids of 1,201, 2,401, and 4,801 points. The default design therefore generates 2,700 replicate-grid evaluations across 9 seed/grid combinations, with 900 distinct simulated replicate records across seeds and paired observations reused within each seed across grid resolutions. A quick local/CI smoke design can use two grids, two seeds, and 100 replicates per seed.

Example:

```bash
python scripts/run_multiseed_grid_sensitivity.py \
  --replicates 300 \
  --seeds 20261025,20261026,20261027 \
  --grid-points 1201,2401,4801 \
  --output-dir artifacts/multiseed_grid_sensitivity
```

For every seed, the runner delegates to the single-seed experiment. Each grid's raw records and summary are independently replayed and audited, and the paired comparison table is audited before the seed's results are aggregated. A failed audit raises an error and prevents writing the new aggregate success tables.

## Outputs

- `seed_<seed>/grid_<points>/`: raw replicate CSV, summary JSON, and independent audit report for each grid.
- `seed_<seed>/grid_sensitivity.csv`: within-seed paired contrasts against the finest requested grid.
- `seed_<seed>/grid_sensitivity_audit.json`: independent reconstruction audit for that seed-level paired comparison table.
- `per_seed_paired_contrasts.csv`: one row per seed and grid, retaining seed-level PIT summaries and within-seed paired differences.
- `multiseed_grid_sensitivity.csv`: one row per grid, aggregating across seeds.
- `multiseed_audit_report.json`: independent reconstruction of the across-seed aggregate statistics and completeness checks.
- `multiseed_metadata.json`: design, audit state, and interpretation metadata.

The aggregate table reports the mean and sample standard deviation of per-seed PIT means, plus the standard error across seeds. It separately reports the mean and standard error across seeds of each seed's mean paired PIT difference from the finest grid. This keeps seed-to-seed variability distinct from the within-seed standard error of paired replicate-level differences. With only two or three seeds, these uncertainty estimates are themselves imprecise and should be treated as descriptive.

## Interpretation and limits

1. Paired differences within each seed compare the same simulated observations, so they isolate posterior-grid discretization conditional on that synthetic dataset.
2. Cross-seed variation probes sensitivity to independent random realizations but does not vary the model assumptions, selection regime, host mixture, or noise distribution.
3. The finest grid is a numerical reference, not a proven exact posterior. Predefine acceptable numerical tolerances before interpreting any convergence claim.
4. A PIT mean near 0.5 or a non-significant uniformity test does not establish calibration by itself.
5. Synthetic matched-model checks do not establish astrophysical validity for real gravitational-wave likelihoods, detection selection, or galaxy catalogues.
6. Publication-facing claims should include more seeds, justified Monte Carlo precision, alternative noise/selection regimes, and an independent review of the statistical model.

## Selection conditioning and estimand (important distinction)

The generator draws truth from the declared uniform prior first, then rejection-samples host and distance conditional on detection at that fixed truth. The inference uses a detection-conditioned likelihood proportional to p(distance | H0) / P(detection | H0). These choices can form a coherent **conditional-on-truth detection experiment**: at each fixed truth, the accepted observations follow the selected likelihood, while the experiment deliberately keeps the simulated truth distribution equal to the declared prior.

That is not the same sampling design as drawing events from a natural detected-event population. In a population of events selected by detection, the distribution of truth among detections is generally reweighted by P(detection | H0); it is not automatically the original prior. The current experiment therefore should not be described as a full prior-predictive simulation of the detected population. The estimand must be stated explicitly, and the conditional likelihood, data-generation process, and target distribution must be derived together before interpreting PIT uniformity. Independent replay verifies implementation fidelity, not the appropriateness of the estimand for a real astrophysical analysis.
