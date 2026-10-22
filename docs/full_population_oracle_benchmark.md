# Full-population oracle benchmark

## Purpose

The existing missing-host control includes a hybrid oracle scenario: true host-generation weights are used for the hidden component, but the observed component retains uniform inference weights. This additional benchmark provides a distinct, internally consistent reference by using the full simulated galaxy population and the true luminosity host-generation weights in one normalized posterior.

## Scenarios and outputs

The script `scripts/run_full_population_oracle.py` generates 500 mock galaxies per replicate, selects the true host using luminosity weights, simulates a noisy gravitational-wave distance, and creates a 70%-complete random or redshift-dependent mask. For this oracle posterior, inference deliberately uses all 500 galaxies and the true luminosity weights; the mask is recorded for selection context but does not remove galaxies from the oracle likelihood.

It writes:
- `scenario_summary.csv` and `scenario_summary.json`: scenario-level bias, RMSE, interval width, coverage, Monte Carlo standard errors, and metadata.
- `per_run_results.csv`: one row per mock universe.

## Interpretation safeguards

This is a **simulation-only upper-information benchmark**, not a method available to a real analysis with an incomplete galaxy catalogue. Because the full galaxy population and true host weights are supplied to the oracle likelihood, it must not be presented as a practical missing-host correction or as evidence that a real survey can recover hidden galaxies. It is useful as a reference point for quantifying the information lost by catalogue incompleteness.

The two mask scenarios are deliberately **paired**: replicate index (r) uses the same catalogue and GW distance draw in both scenarios, while only the mask mechanism differs. Since the oracle likelihood ignores the mask and uses the same full catalogue and weights, posterior summaries must match exactly within each pair. Realized completeness and whether the host falls in the observed subset may differ. The workflow checks this invariance as a regression guard.

The benchmark is included in the manual `factorial-sensitivity` workflow. A passing unit-test workflow alone does not mean the full experiment has executed; check the manual workflow run and inspect its uploaded artifact before reporting any numerical result.

## Reproducibility

The original simulator is pinned to `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`. The workflow passes the configured base seed and replicate count to this benchmark. Both mask scenarios use the same seed block, making the oracle posterior a paired, deterministic invariance check.
