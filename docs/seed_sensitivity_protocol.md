# Seed-sensitivity workflow

The workflow `.github/workflows/original-coverage-seed-sensitivity.yml` runs
the pinned original-model reproduction independently for seeds 2026 and 31027.
Both runs use the same replicate count and the same source commit.

## Run it

1. Open the repository's Actions tab.
2. Select **original-coverage-seed-sensitivity**.
3. Choose **Run workflow** and set replicates to 1000 initially.
4. Download both artifacts:
   - `original-coverage-seed-2026`
   - `original-coverage-seed-31027`
5. Compare the `summary.json` and `coverage_curve.csv` files.

Workflow page:
https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/workflows/original-coverage-seed-sensitivity.yml

## Interpretation

The two seeds are independent Monte Carlo replications of the same simulation
design, not independent observational datasets. Differences between their
coverage estimates should be judged against Monte Carlo uncertainty. For a
coverage proportion near 0.94, the approximate binomial Monte Carlo standard
error is 0.75 percentage points at 1000 repetitions and 0.34 percentage
points at 5000 repetitions.

Do not treat a small difference between seeds as evidence of a physical
effect. First check that the coverage estimates, interval widths, bias, and
RMSE are stable within expected Monte Carlo variation. If the runs differ
more than expected, inspect random-number handling and the data-generation
pipeline before interpreting the science.

This workflow has been added to the repository, but its runs have not been
triggered or independently inspected yet. Do not cite its output as a result
until the artifacts are available and reviewed.
