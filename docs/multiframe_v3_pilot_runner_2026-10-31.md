# Multi-frame v3 pilot runner

**Status:** executable pilot runner implemented; no confirmatory experiment or scientific conclusion is claimed.

Run from the repository root after installing the project:

```bash
python scripts/run_multiframe_v3_pilot.py \
  --manifest configs/multiframe_v3_smoke_manifest.json \
  --output-dir artifacts/multiframe_v3_smoke
```

The manifest declares explicit scenario IDs, frames per scenario, outer repetitions, sample sizes, bootstrap replicates, nominal confidence level, and the three supported method IDs:
- `wald_srswor`
- `ordinary_percentile_bootstrap`
- `studentized_pseudopop_srswor`

The runner generates all finite frames, then one method-agnostic SRSWOR sample per scenario/frame/sample-size/replicate key. All three methods therefore receive identical sampled units. Stochastic methods use a stable SHA-256-derived method coordinate with the v3 method seed domain, so changing method order does not change a method's random stream.

Outputs:
- `frames.csv`: all realized finite-frame units and provenance.
- `replicates.csv`: one row per method and outer-sample key, including paired indices, seed lineage, interval endpoints, success/failure, and studentized inner-attempt accounting.
- `frame_summaries.csv`: per-frame conditional coverage, interval success, unconditional successful-and-covering rate, and width.
- `equal_frame_summaries.csv`: equal-frame-weighted descriptive summaries; frames, not pooled outer attempts, are the unit of equal weighting.
- `manifest.json`: declared design, software metadata, row counts, pilot-only status, and SHA-256 checksums.

The equal-frame summaries are descriptive. With only two frames in the smoke manifest they are not precision estimates and must not be interpreted as confirmatory. A separate independent v3 auditor and a precision-based planning run remain release gates.
