# Joint prior-performance scorecard protocol

## Purpose

The prior-support and interval-width experiments should be interpreted together. This scorecard joins empirical coverage, Wilson uncertainty bounds, bias, RMSE, and posterior interval-width summaries into one row per truth × prior × method.

It is a descriptive diagnostic and does not select a prior automatically.

## Design

The script runs the same paired synthetic experiment with truths 52, 70, and 88; prior supports [50,90], [40,100], and [30,110]; 300 replicates per truth; seed 20261024; and a 4,001-point posterior grid. It joins each coverage/bias/RMSE summary with the interval-width summary for the identical experiment design.

The output contains 18 rows (3 truths × 3 priors × 2 methods) and includes absolute coverage error relative to nominal 68% and 95% levels.

## Reproduction

```bash
python scripts/run_joint_prior_scorecard.py --replicates 300 --seed 20261024 --grid-size 4001 --output results/joint-prior-scorecard.csv
python -m pytest tests/test_joint_prior_scorecard.py -q
```

The GitHub Actions workflow tests Python 3.10, 3.11, and 3.12 and uploads the CSV artifact.

## Important implementation note

The scorecard computes coverage/bias/RMSE and interval-width summaries using the same seed and data-generating code. However, the current script executes the data-generating audit twice to obtain the summary and width inputs. Because both runs use the same seed, they should be deterministic and paired; future optimization can reduce this to a single simulation pass without changing results. The scorecard does not replace independent validation of the likelihood, selection function, or prior choice.

## Interpretation

Compare coverage with width, bias, and RMSE together. Coverage close to nominal is not sufficient if intervals are uninformatively broad, and narrow intervals are not desirable if they systematically miss the truth. No prior is declared preferred merely because it minimizes one score.
