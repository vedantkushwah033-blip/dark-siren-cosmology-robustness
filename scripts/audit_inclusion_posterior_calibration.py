"""Independently audit multi-replicate inclusion-posterior calibration artifacts.

This checks orchestration, paired/nested design, PIT/coverage summaries and
per-replicate audit attestations. It does not replace the raw-data auditor for
each replicate or establish scientific calibration by itself.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _pit(frame: pd.DataFrame, truth: float, column: str) -> float:
    grid = frame["H0"].to_numpy(dtype=float)
    density = frame[column].to_numpy(dtype=float)
    if (grid.ndim != 1 or density.shape != grid.shape or grid.size < 2
        or np.any(~np.isfinite(grid)) or np.any(~np.isfinite(density))
        or np.any(np.diff(grid) <= 0) or np.any(density < 0)
        or not np.isfinite(truth) or truth < grid[0] or truth > grid[-1]):
        raise AssertionError("invalid posterior grid, density, or generating H0")
    increments = 0.5 * (density[1:] + density[:-1]) * np.diff(grid)
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise AssertionError("posterior density has invalid integrated mass")
    return float(np.interp(truth, grid, cdf / cdf[-1]))


def _assert_nested_validation(artifact_dirs: list[Path]) -> dict[str, int]:
    usable = [path for path in artifact_dirs if path.exists()]
    if len(usable) < 2:
        return {"paired_sizes_checked": 0, "analysis_files_checked": 0}
    configs = [json.loads((path / "config.json").read_text()) for path in usable]
    analysis_hashes = []
    analysis_files = ("analysis_catalogue.csv", "gw_events_observed.csv")
    for path in usable:
        for name in analysis_files:
            analysis_hashes.append((name, _sha256(path / "observable" / name)))
    for name in analysis_files:
        hashes = {value for file_name, value in analysis_hashes if file_name == name}
        if len(hashes) != 1:
            raise AssertionError(f"analysis {name} changed across validation sizes")
    reference_cfg = configs[0]["analysis_generator_config"]
    for cfg in configs[1:]:
        if cfg["analysis_generator_config"] != reference_cfg:
            raise AssertionError("analysis generator configuration changed across validation sizes")

    frame_reference = None
    selected_prefix = []
    for path, cfg in zip(usable, configs):
        truth = pd.read_csv(path / "truth" / "validation_frame_truth.csv")
        observable = pd.read_csv(path / "observable" / "validation_observed.csv")
        required = ["frame_id", "redshift", "inclusion_probability", "retained"]
        current = truth[required].copy()
        if frame_reference is None:
            frame_reference = current
        elif not frame_reference.equals(current):
            raise AssertionError("validation parent frame changed across nested sample sizes")
        ids = observable["validation_id"].to_numpy(dtype=int)
        if len(ids) != int(cfg["validation_generator_config"]["n_validation"]):
            raise AssertionError("validation observable row count differs from config")
        selected_prefix.append(ids)
        truth_by_id = truth.set_index("frame_id")
        matched = truth_by_id.loc[ids]
        if not np.allclose(matched["redshift"], observable["redshift"], rtol=0, atol=1e-14):
            raise AssertionError("validation redshifts do not match selected frame IDs")
        if not np.array_equal(matched["retained"].to_numpy(), observable["retained"].to_numpy()):
            raise AssertionError("validation labels do not match selected frame IDs")
    order = np.argsort([len(ids) for ids in selected_prefix])
    ordered = [selected_prefix[index] for index in order]
    for smaller, larger in zip(ordered, ordered[1:]):
        if not np.array_equal(smaller, larger[:len(smaller)]):
            raise AssertionError("validation samples are not literal nested permutation prefixes")
    return {"paired_sizes_checked": len(usable), "analysis_files_checked": len(analysis_files)}


def _wilson(successes: int, trials: int) -> tuple[float, float]:
    if trials <= 0:
        return (float("nan"), float("nan"))
    z = 1.959963984540054
    p = successes / trials
    denominator = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / denominator
    half = z * np.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / denominator
    return (max(0.0, center - half), min(1.0, center + half))


def audit(root: Path) -> dict:
    root = root.resolve()
    planned = pd.read_csv(root / "planned_cells.csv")
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    summary = json.loads((root / "calibration_summary.json").read_text(encoding="utf-8"))
    if planned["cell_id"].duplicated().any() or ledger["cell_id"].duplicated().any():
        raise AssertionError("planned or ledger cell IDs are duplicated")
    planned_ids, ledger_ids = set(planned["cell_id"]), set(ledger["cell_id"])
    if planned_ids != ledger_ids:
        raise AssertionError("planned cells and replicate ledger do not match exactly")
    allowed = {"SUCCESS", "NON_EVALUABLE", "FAILED", "NOT_RUN"}
    if not set(ledger["status"].dropna()).issubset(allowed):
        raise AssertionError("ledger contains an unknown status")
    if int(summary["planned_cells"]) != len(planned):
        raise AssertionError("summary planned-cell count is incorrect")
    attempted = int((ledger["status"] != "NOT_RUN").sum())
    if int(summary["attempted"]) != attempted:
        raise AssertionError("summary attempted count is incorrect")
    status_counts = ledger["status"].value_counts().to_dict()
    if {str(k): int(v) for k, v in summary["status_counts"].items()} != {
        str(k): int(v) for k, v in status_counts.items()
    }:
        raise AssertionError("summary status counts do not match replicate ledger")

    artifact_paths: dict[str, Path] = {}
    per_replicate_audits = 0
    for row in ledger.itertuples(index=False):
        status = str(row.status)
        if status in {"SUCCESS", "NON_EVALUABLE"}:
            path = Path(row.artifact_dir)
            if not path.is_dir():
                raise AssertionError(f"missing artifact directory for {row.cell_id}")
            report_path = path / "audit_report.json"
            report = json.loads(report_path.read_text(encoding="utf-8"))
            if report.get("status") != "PASS" or getattr(row, "audit_status", None) != "PASS":
                raise AssertionError(f"per-replicate audit did not pass for {row.cell_id}")
            config_path = path / "config.json"
            if getattr(row, "config_sha256", "") and _sha256(config_path) != row.config_sha256:
                raise AssertionError(f"config checksum mismatch for {row.cell_id}")
            artifact_paths[str(row.cell_id)] = path
            per_replicate_audits += 1

            if status == "SUCCESS":
                convergence = json.loads((path / "quadrature_convergence.json").read_text())
                density = pd.read_csv(path / "h0_posterior_density_by_order.csv")
                final_order = int(convergence["orders_evaluated"][-1])
                final = convergence["posterior_by_order"][-1]
                if convergence["status"] != "CONVERGED":
                    raise AssertionError(f"SUCCESS row lacks quadrature convergence: {row.cell_id}")
                for key in ("h0_q025", "h0_q05", "h0_q16", "h0_median", "h0_q84", "h0_q95", "h0_q975"):
                    if not np.isclose(float(getattr(row, key)), float(final[key]), rtol=1e-10, atol=1e-10):
                        raise AssertionError(f"ledger {key} differs from posterior artifact: {row.cell_id}")
                experiment = str(row.experiment)
                if experiment == "H0_PIT_FIXED_INCLUSION_TRUTH":
                    pit_value = _pit(density, float(row.h0_true), f"density_order_{final_order}")
                    if not np.isclose(pit_value, float(row.pit), rtol=1e-9, atol=1e-9):
                        raise AssertionError(f"ledger PIT differs from posterior reconstruction: {row.cell_id}")
                elif experiment == "FIXED_TRUTH":
                    for level, lower, upper in ((68, "h0_q16", "h0_q84"), (90, "h0_q05", "h0_q95"), (95, "h0_q025", "h0_q975")):
                        expected = int(float(getattr(row, lower)) <= float(row.h0_true) <= float(getattr(row, upper)))
                        if expected != int(getattr(row, f"covered_{level}")):
                            raise AssertionError(f"ledger coverage flag mismatch: {row.cell_id}")
                        width = float(getattr(row, upper)) - float(getattr(row, lower))
                        if not np.isclose(width, float(getattr(row, f"width_{level}")), rtol=1e-10, atol=1e-10):
                            raise AssertionError(f"ledger interval width mismatch: {row.cell_id}")

    # Verify paired analysis data and nested validation prefixes within each
    # experiment/floor/replicate block, without combining different designs.
    paired_blocks = 0
    group_cols = ["experiment", "catalogue_floor", "replicate"]
    for _, group in ledger.groupby(group_cols, dropna=False):
        paths = [artifact_paths[str(cell)] for cell in group["cell_id"] if str(cell) in artifact_paths]
        if len(paths) >= 2:
            _assert_nested_validation(paths)
            paired_blocks += 1

    # Recompute headline summary values from successful ledger rows.
    for key, item in summary["summary_by_cell"].items():
        subset = ledger[
            (ledger["experiment"] == item["experiment"])
            & np.isclose(ledger["catalogue_floor"].astype(float), float(item["catalogue_floor"]))
            & (ledger["n_validation"] == int(item["n_validation"]))
        ]
        successful = subset[subset["status"] == "SUCCESS"]
        if int(item["attempted"]) != len(subset) or int(item["success"]) != len(successful):
            raise AssertionError(f"summary attempt/success counts mismatch for {key}")
        if item["experiment"] == "H0_PIT_FIXED_INCLUSION_TRUTH" and len(successful):
            pit_values = successful["pit"].to_numpy(dtype=float)
            if not np.isclose(float(item["pit_mean"]), float(np.mean(pit_values)), rtol=1e-10, atol=1e-10):
                raise AssertionError(f"summary PIT mean mismatch for {key}")
            if len(pit_values) > 1 and not np.isclose(float(item["pit_variance"]), float(np.var(pit_values, ddof=1)), rtol=1e-10, atol=1e-10):
                raise AssertionError(f"summary PIT variance mismatch for {key}")
            bins = len(item["pit_bin_counts"])
            histogram_counts, _ = np.histogram(pit_values, bins=np.linspace(0.0, 1.0, bins + 1))
            if not np.array_equal(histogram_counts, np.asarray(item["pit_bin_counts"], dtype=int)):
                raise AssertionError(f"summary PIT histogram mismatch for {key}")
            for level in (68, 90, 95):
                successes = int(successful[f"h0_pit_interval_included_{level}"].sum())
                rate = successes / len(successful)
                if not np.isclose(rate, float(item[f"prior_predictive_interval_inclusion_{level}"]), rtol=1e-10, atol=1e-10):
                    raise AssertionError(f"summary prior-predictive interval inclusion mismatch for {key}, level={level}")
                interval = _wilson(successes, len(successful))
                if not np.allclose(interval, item[f"prior_predictive_interval_inclusion_{level}_wilson95"], rtol=1e-10, atol=1e-10):
                    raise AssertionError(f"summary PIT interval Wilson interval mismatch for {key}, level={level}")
        if item["experiment"] == "FIXED_TRUTH" and len(successful):
            bias = (successful["h0_median"] - successful["h0_true"]).to_numpy(dtype=float)
            if not np.isclose(float(item["mean_posterior_median_bias"]), float(np.mean(bias)), rtol=1e-10, atol=1e-10):
                raise AssertionError(f"summary posterior-median bias mismatch for {key}")
            if not np.isclose(float(item["median_absolute_posterior_median_bias"]), float(np.median(np.abs(bias))), rtol=1e-10, atol=1e-10):
                raise AssertionError(f"summary absolute median bias mismatch for {key}")
            for level in (68, 90, 95):
                successes = int(successful[f"covered_{level}"].sum())
                rate = successes / len(successful)
                if not np.isclose(rate, float(item[f"coverage_{level}"]), rtol=1e-10, atol=1e-10):
                    raise AssertionError(f"summary coverage mismatch for {key}, level={level}")
                interval = _wilson(successes, len(successful))
                if not np.allclose(interval, item[f"coverage_{level}_wilson95"], rtol=1e-10, atol=1e-10):
                    raise AssertionError(f"summary coverage Wilson interval mismatch for {key}, level={level}")
                if not np.isclose(float(successful[f"width_{level}"].median()), float(item[f"median_interval_width_{level}"]), rtol=1e-10, atol=1e-10):
                    raise AssertionError(f"summary interval width mismatch for {key}, level={level}")

    return {
        "status": "PASS",
        "planned_cells": int(len(planned)),
        "attempted_cells": attempted,
        "status_counts": {str(k): int(v) for k, v in status_counts.items()},
        "per_replicate_audits_passed": per_replicate_audits,
        "paired_design_blocks_checked": paired_blocks,
        "scientific_limit": "Orchestration and arithmetic integrity only; conditional H0 PIT and fixed-truth coverage are not joint SBC or real-data validation.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = audit(args.artifact_dir)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
