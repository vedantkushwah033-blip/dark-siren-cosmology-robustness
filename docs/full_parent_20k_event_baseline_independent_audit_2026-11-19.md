# Independent raw-seed audit of the 20,000-draw full-parent baseline

This audit regenerates all 100 archived parent populations from the immutable 600-cell source artifact, then independently draws 200 detected four-event samples per parent using the archived fixed-truth design and deterministic seed schedule. It independently reconstructs all 20,000 selected-event posterior curves, quantiles, coverage indicators, interval widths, pooled metrics, and parent-cluster summaries.

The audit targets the completed baseline artifact from workflow 44436218047 (artifact 14490011968). It is intentionally separate from the baseline runner and uses a separately coded sampler and selected-event posterior calculation.

## Scope

This is a fixed-truth four-event full-parent diagnostic. It is not joint SBC, not a catalogue-inclusion comparison, and not observational H0 inference. Passing the audit verifies raw-seed reproducibility and arithmetic; it does not establish nominal calibration.
