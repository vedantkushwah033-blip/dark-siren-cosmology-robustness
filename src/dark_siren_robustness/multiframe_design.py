"""Validated scenario manifests and deterministic finite-frame generation.

This is the first v3 data-generation layer. It deliberately stops before
interval construction: the existing v2 statistical methods are not silently
reused under a new schema until their contracts are adapted and audited.
"""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from numbers import Integral, Real
from typing import Any

import numpy as np

from dark_siren_robustness.multiframe_seed_schedule import frame_seed

MANIFEST_SCHEMA = "linear_finite_population_multiframe_manifest_v3"
FRAME_SCHEMA = "linear_finite_population_frame_v3"
SUPPORTED_FAMILIES = {"normal", "gamma", "lognormal"}


def _integer(name: str, value: Any, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be an integer, not a boolean")
    value = int(value)
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return value


def _finite_real(name: str, value: Any, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number, not a boolean")
    value = float(value)
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return value


def validate_multiframe_manifest(manifest: Mapping[str, Any]) -> dict[str, Any]:
    """Validate and normalize the explicit scenario/frame experiment design.

    Expected fields: schema, root_seed, population_size, frames_per_scenario,
    outer_repetitions, bootstrap_replicates, nominal_level, sample_sizes,
    methods, scenarios. Each scenario has scenario_id, name, family, params.
    """
    if not isinstance(manifest, Mapping):
        raise TypeError("manifest must be a mapping")
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError(f"schema must equal {MANIFEST_SCHEMA}")
    normalized = dict(manifest)
    normalized["root_seed"] = _integer("root_seed", manifest.get("root_seed"))
    normalized["population_size"] = _integer(
        "population_size", manifest.get("population_size"), 2
    )
    normalized["frames_per_scenario"] = _integer(
        "frames_per_scenario", manifest.get("frames_per_scenario"), 2
    )
    normalized["outer_repetitions"] = _integer(
        "outer_repetitions", manifest.get("outer_repetitions"), 1
    )
    normalized["bootstrap_replicates"] = _integer(
        "bootstrap_replicates", manifest.get("bootstrap_replicates"), 2
    )
    level = _finite_real("nominal_level", manifest.get("nominal_level"))
    if not 0 < level < 1:
        raise ValueError("nominal_level must lie strictly between zero and one")
    normalized["nominal_level"] = level

    sizes = manifest.get("sample_sizes")
    if not isinstance(sizes, Sequence) or isinstance(sizes, (str, bytes)) or not sizes:
        raise TypeError("sample_sizes must be a non-empty sequence")
    clean_sizes = [_integer("sample_size", n, 2) for n in sizes]
    if len(set(clean_sizes)) != len(clean_sizes):
        raise ValueError("sample_sizes must be unique")
    if any(n >= normalized["population_size"] for n in clean_sizes):
        raise ValueError("every sample_size must be smaller than population_size")
    normalized["sample_sizes"] = clean_sizes

    methods = manifest.get("methods")
    if not isinstance(methods, Sequence) or isinstance(methods, (str, bytes)) or len(methods) < 2:
        raise TypeError("methods must contain at least two method records")
    method_ids = []
    for method in methods:
        if not isinstance(method, Mapping) or not isinstance(method.get("method_id"), str) or not method["method_id"]:
            raise ValueError("each method requires a non-empty string method_id")
        method_ids.append(method["method_id"])
    if len(set(method_ids)) != len(method_ids):
        raise ValueError("method_id values must be unique")
    normalized["methods"] = [dict(x) for x in methods]

    scenarios = manifest.get("scenarios")
    if not isinstance(scenarios, Sequence) or isinstance(scenarios, (str, bytes)) or not scenarios:
        raise TypeError("scenarios must be a non-empty sequence")
    scenario_ids = []
    clean_scenarios = []
    for item in scenarios:
        if not isinstance(item, Mapping):
            raise TypeError("each scenario must be a mapping")
        sid = _integer("scenario_id", item.get("scenario_id"))
        if sid in scenario_ids:
            raise ValueError("scenario_id values must be unique")
        scenario_ids.append(sid)
        family = item.get("family")
        if family not in SUPPORTED_FAMILIES:
            raise ValueError(f"unsupported distribution family: {family!r}")
        name = item.get("name")
        if not isinstance(name, str) or not name:
            raise ValueError("each scenario requires a non-empty name")
        params = item.get("params")
        if not isinstance(params, Mapping):
            raise TypeError("scenario params must be a mapping")
        p = dict(params)
        if family == "normal":
            p["loc"] = _finite_real("normal loc", p.get("loc"))
            p["scale"] = _finite_real("normal scale", p.get("scale"), positive=True)
        elif family == "gamma":
            p["shape"] = _finite_real("gamma shape", p.get("shape"), positive=True)
            p["scale"] = _finite_real("gamma scale", p.get("scale"), positive=True)
        else:
            p["meanlog"] = _finite_real("lognormal meanlog", p.get("meanlog"))
            p["sigmalog"] = _finite_real("lognormal sigmalog", p.get("sigmalog"), positive=True)
        clean_scenarios.append({
            "scenario_id": sid, "name": name, "family": family, "params": p
        })
    normalized["scenarios"] = clean_scenarios
    normalized["schema"] = MANIFEST_SCHEMA
    return normalized


def generate_finite_frame(
    manifest: Mapping[str, Any], scenario_id: int, frame_id: int
) -> dict[str, Any]:
    """Generate one finite frame and auditable metadata from manifest coordinates."""
    design = validate_multiframe_manifest(manifest)
    sid = _integer("scenario_id", scenario_id)
    fid = _integer("frame_id", frame_id)
    matches = [s for s in design["scenarios"] if s["scenario_id"] == sid]
    if not matches:
        raise ValueError(f"unknown scenario_id: {sid}")
    scenario = matches[0]
    seed = frame_seed(design["root_seed"], sid, fid)
    rng = np.random.default_rng(seed)
    n = design["population_size"]
    params = scenario["params"]
    if scenario["family"] == "normal":
        values = rng.normal(params["loc"], params["scale"], size=n)
    elif scenario["family"] == "gamma":
        values = rng.gamma(params["shape"], params["scale"], size=n)
    else:
        values = rng.lognormal(params["meanlog"], params["sigmalog"], size=n)
    values = np.asarray(values, dtype=np.float64)
    if values.shape != (n,) or not np.all(np.isfinite(values)):
        raise FloatingPointError("generated finite frame has invalid shape or non-finite values")
    values.setflags(write=False)
    return {
        "schema": FRAME_SCHEMA,
        "scenario_id": sid,
        "scenario_name": scenario["name"],
        "frame_id": fid,
        "family": scenario["family"],
        "seed": seed,
        "population_size": n,
        "population_mean": float(np.mean(values)),
        "population_variance": float(np.var(values, ddof=0)),
        "values": values,
    }


def generate_frame_grid(manifest: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Generate the complete scenario × frame grid in explicit manifest order."""
    design = validate_multiframe_manifest(manifest)
    return [
        generate_finite_frame(design, scenario["scenario_id"], frame_id)
        for scenario in design["scenarios"]
        for frame_id in range(design["frames_per_scenario"])
    ]
