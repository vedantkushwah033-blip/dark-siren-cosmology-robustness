"""Independent integrity audit for preparation-only synthetic shared-input bundles.

This auditor checks frozen hashes, catalogue inclusion truth, validation-frame
joins/nesting, and the GW detection trace. It imports no production generator,
inference, KDE, posterior, or lifecycle code.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit

C_KM_S = 299_792.458
REQUIRED = {
    "observable/catalogue_observed.csv",
    "observable/gw_events_observed.csv",
    "truth/parent_truth.csv",
    "truth/gw_event_truth.csv",
    "truth/gw_generation_trace_truth.csv",
    "truth/validation_frame_truth.csv",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _same(actual, expected, label: str, *, atol: float = 1e-12) -> None:
    a, e = np.asarray(actual), np.asarray(expected)
    if a.shape != e.shape:
        raise AssertionError(f"{label}: shape mismatch {a.shape} != {e.shape}")
    if not np.allclose(a.astype(float), e.astype(float), rtol=1e-11, atol=atol):
        delta = float(np.max(np.abs(a.astype(float) - e.astype(float))))
        raise AssertionError(f"{label}: values differ (max abs={delta})")


def audit(bundle_dir: Path) -> dict:
    root = Path(bundle_dir).resolve()
    outer_path = root / "manifest.json"
    shared_path = root / "shared_input_manifest.json"
    config_path = root / "config.json"
    if not outer_path.is_file() or not shared_path.is_file() or not config_path.is_file():
        raise AssertionError("bundle manifest/config files are incomplete")
    outer = json.loads(outer_path.read_text(encoding="utf-8"))
    shared = json.loads(shared_path.read_text(encoding="utf-8"))
    config = json.loads(config_path.read_text(encoding="utf-8"))

    if outer.get("schema_version") != 1 or outer.get("status") != "PASS":
        raise AssertionError("outer manifest schema/status is invalid")
    if outer.get("preparation_only") is not True or outer.get("inference_executed") is not False:
        raise AssertionError("bundle is not explicitly marked preparation-only")
    if config.get("status") != "FROZEN_INPUTS_ONLY_NO_INFERENCE_EXECUTED":
        raise AssertionError("config does not identify a preparation-only bundle")
    if (root / "results").exists() or list(root.rglob("*.npz")):
        raise AssertionError("preparation-only bundle contains inference outputs")

    entries = outer.get("files")
    if not isinstance(entries, dict) or set(entries) != {
        str(p.relative_to(root)) for p in root.rglob("*")
        if p.is_file() and p.name != "manifest.json"
    }:
        raise AssertionError("outer manifest file set differs from bundle contents")
    for relative, meta in entries.items():
        path = (root / relative).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise AssertionError(f"manifest path is missing or escapes bundle: {relative}")
        if path.stat().st_size != meta.get("bytes") or _sha256(path) != meta.get("sha256"):
            raise AssertionError(f"outer manifest hash/size mismatch: {relative}")
    if len(entries) != outer.get("file_count"):
        raise AssertionError("outer manifest file_count mismatch")
    if outer.get("shared_input_manifest_sha256") != _sha256(shared_path):
        raise AssertionError("outer manifest does not bind the shared-input manifest")

    sizes = config.get("validation_sizes")
    if (not isinstance(sizes, list) or not sizes
            or any(type(v) is not int or v < 5 for v in sizes)
            or sizes != sorted(set(sizes))):
        raise AssertionError("config validation-size plan is invalid")
    if shared.get("schema_version") != 1 or shared.get("status") != "FROZEN":
        raise AssertionError("shared-input manifest schema/status is invalid")
    if shared.get("validation_sizes") != sizes:
        raise AssertionError("shared-input manifest validation sizes differ from config")
    expected_files = REQUIRED | {
        f"observable/validation_n{n}" + ".csv" for n in sizes
    }
    if set(shared.get("files", {})) != expected_files:
        raise AssertionError("shared-input manifest file set is invalid")
    for relative, meta in shared["files"].items():
        path = (root / relative).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise AssertionError(f"shared-input path missing/escaping: {relative}")
        if path.stat().st_size != meta.get("bytes") or _sha256(path) != meta.get("sha256"):
            raise AssertionError(f"shared-input hash/size mismatch: {relative}")

    provenance = shared.get("provenance", {})
    for key, value in {
        "seed": config["seed"], "event_seed": config["event_seed"],
        "validation_seed": config["validation_seed"],
        "validation_sample_seed": config["validation_sample_seed"],
        "n_parent": config["n_parent"], "n_detected": config["n_detected"],
        "h0_true": config["h0_true"], "sigma_mpc": config["sigma_mpc"],
        "detection_threshold_mpc": config["detection_threshold_mpc"],
        "n_validation_frame": config["validation_frame"],
        "curvature_gamma": config["curvature_gamma"],
        "curvature_beta": config["curvature_beta"],
        "inclusion_model_degree": config["inclusion_model_degree"],
    }.items():
        if provenance.get(key) != value:
            raise AssertionError(f"provenance/config mismatch: {key}")

    parent = pd.read_csv(root / "truth/parent_truth.csv")
    catalogue = pd.read_csv(root / "observable/catalogue_observed.csv")
    frame = pd.read_csv(root / "truth/validation_frame_truth.csv")
    if len(parent) != int(config["n_parent"]) or not np.array_equal(
        parent["object_id"].to_numpy(int), np.arange(len(parent))
    ):
        raise AssertionError("parent population size or object IDs are invalid")
    zmin, zmax = 0.01, 0.40
    if not np.all(np.isfinite(parent["redshift"])) or np.any(
        (parent["redshift"] < zmin) | (parent["redshift"] > zmax)
    ):
        raise AssertionError("parent redshifts are outside the declared domain")
    if not np.all(np.isfinite(parent["inclusion_probability"])) or np.any(
        (parent["inclusion_probability"] <= 0) | (parent["inclusion_probability"] >= 1)
    ):
        raise AssertionError("parent inclusion probabilities are invalid")

    calibration = config["quadratic_calibration"]
    gamma, beta = float(config["curvature_gamma"]), float(config["curvature_beta"])
    center, scale = 0.205, 0.195
    x = (parent["redshift"].to_numpy(float) - center) / scale
    expected_p = expit(
        float(calibration["alpha"]) + beta * x
        + gamma * (x**2 - float(calibration["x2_center"]))
    )
    _same(parent["inclusion_probability"], expected_p, "parent inclusion probabilities")
    expected_catalogue = parent.loc[parent["retained"].astype(int) == 1, ["object_id", "redshift"]]
    expected_catalogue = expected_catalogue.reset_index(drop=True)
    if not np.array_equal(catalogue["object_id"].to_numpy(int), expected_catalogue["object_id"].to_numpy(int)):
        raise AssertionError("observed catalogue IDs do not equal retained parent IDs")
    _same(catalogue["redshift"], expected_catalogue["redshift"], "catalogue redshifts")

    if len(frame) != int(config["validation_frame"]) or not np.array_equal(
        frame["frame_id"].to_numpy(int), np.arange(len(frame))
    ):
        raise AssertionError("validation frame size or frame IDs are invalid")
    xf = (frame["redshift"].to_numpy(float) - center) / scale
    expected_pf = expit(
        float(calibration["alpha"]) + beta * xf
        + gamma * (xf**2 - float(calibration["x2_center"]))
    )
    _same(frame["inclusion_probability"], expected_pf, "validation-frame inclusion probabilities")

    prior_min, prior_max = map(float, config["h0_prior_support"])
    if not (np.isfinite(prior_min) and np.isfinite(prior_max)
            and prior_min < float(config["h0_true"]) < prior_max):
        raise AssertionError("H0 truth/prior support is invalid")
    event_obs = pd.read_csv(root / "observable/gw_events_observed.csv")
    event_truth = pd.read_csv(root / "truth/gw_event_truth.csv")
    trace = pd.read_csv(root / "truth/gw_generation_trace_truth.csv")
    n_events = int(config["n_detected"])
    if len(event_obs) != n_events or len(event_truth) != n_events:
        raise AssertionError("detected GW event count differs from config")
    if not np.array_equal(event_obs["event_id"].to_numpy(int), np.arange(n_events)):
        raise AssertionError("GW event IDs are invalid")
    if np.any(event_obs["distance_mpc"].to_numpy(float) < float(config["detection_threshold_mpc"])):
        raise AssertionError("observed GW event violates the detection threshold")
    if not np.array_equal(event_truth["event_id"].to_numpy(int), np.arange(n_events)):
        raise AssertionError("GW truth event IDs are invalid")
    hosts = parent.set_index("object_id").loc[event_truth["host_object_id"].to_numpy(int)]
    _same(event_truth["host_redshift"], hosts["redshift"], "GW host redshifts")
    _same(
        event_truth["latent_mean_distance_mpc"],
        C_KM_S * hosts["redshift"].to_numpy(float) / float(config["h0_true"]),
        "GW latent mean distances",
    )
    selected = trace.loc[trace["selected_as_event"].astype(int) == 1]
    if len(selected) != n_events or not np.all(selected["detected"].to_numpy(int) == 1):
        raise AssertionError("GW generation trace selected-event accounting is invalid")
    _same(event_obs["distance_mpc"], selected["distance_mpc"], "GW selected distances")
    if not np.array_equal(
        event_truth["host_object_id"].to_numpy(int), selected["host_object_id"].to_numpy(int)
    ):
        raise AssertionError("GW generation trace host IDs disagree with event truth")

    previous_ids = None
    for n in sizes:
        validation = pd.read_csv(root / f"observable/validation_n{n}.csv")
        if len(validation) != n or validation["validation_id"].duplicated().any():
            raise AssertionError(f"validation sample size/IDs are invalid for n={n}")
        ids = validation["validation_id"].to_numpy(int)
        if np.any((ids < 0) | (ids >= len(frame))):
            raise AssertionError(f"validation IDs out of range for n={n}")
        joined = validation.merge(
            frame[["frame_id", "redshift", "retained", "selected_for_validation"]],
            left_on="validation_id", right_on="frame_id", validate="one_to_one",
            suffixes=("_observed", "_truth"),
        )
        _same(joined["redshift_observed"], joined["redshift_truth"], f"validation redshifts n={n}")
        _same(joined["retained_observed"], joined["retained_truth"], f"validation outcomes n={n}")
        if not np.all(joined["selected_for_validation"].to_numpy(int) == 1):
            raise AssertionError(f"validation rows not marked selected in truth frame for n={n}")
        if previous_ids is not None and not np.array_equal(previous_ids, ids[:len(previous_ids)]):
            raise AssertionError("validation slices are not nested prefixes")
        previous_ids = ids

    return {
        "status": "PASS", "preparation_only": True, "inference_executed": False,
        "outer_manifest_files_verified": len(entries),
        "shared_input_files_verified": len(shared["files"]),
        "parent_rows_verified": len(parent), "catalogue_rows_verified": len(catalogue),
        "validation_frame_rows_verified": len(frame),
        "validation_sizes_verified": sizes, "detected_events_verified": n_events,
        "scope": "Raw synthetic input integrity only; no likelihood or posterior calibration was evaluated.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle_dir", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = audit(args.bundle_dir)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
