# Configuration-driven nonlinear inclusion screen

The reproducible entry point is:

```bash
python scripts/run_estimated_inclusion_curvature_screen_from_config.py \
  --config configs/estimated_inclusion_curvature_screen_v1.json
```

The JSON manifest is schema version 1. It explicitly records replicate count, parent population size, detected-event count, validation sizes, root seed, numerical grid sizes, curvature levels and output directory. Unknown/missing keys, duplicate sizes/curvatures, invalid types and out-of-range validation sizes are rejected before the experiment is launched.

The wrapper computes SHA-256 over the canonical JSON representation (sorted keys, compact separators). After the underlying runner completes, it records the exact configuration and digest in `screen_report.json`, and writes `experiment_config.canonical.json` alongside the raw condition artifacts. The underlying runner retains responsibility for simulation execution, per-condition independent auditing, and refusing to overwrite a nonempty output directory.

The supplied manifest is the planned 20-replicate screening design:
- parent population: 4,000;
- detected events: 30;
- validation sizes: 100, 300, 1,000;
- curvature: -2, -1, 0, +1, +2;
- H0/redshift grids: 41/201 points;
- root seed: 20261031.

This configuration is not evidence that the screen has run. The screen remains non-confirmatory. Do not increase replication count or proceed to a confirmatory campaign until the output from this exact manifest is complete, independently audited, and reviewed for positivity, failure rates and paired denominators.

The manifest hash identifies the declared configuration, not the source-code environment. For a fully reproducible archival run, also record the source commit, Python/dependency versions, platform, and artifact hashes. These are separate provenance fields and must not be inferred from the configuration digest.
