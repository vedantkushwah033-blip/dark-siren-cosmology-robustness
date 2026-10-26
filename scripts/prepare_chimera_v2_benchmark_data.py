"""Download and checksum a small, opt-in CHIMERA v2.0 benchmark fixture.

This script is intentionally not run in ordinary CI: the default fixture is about
482 MB and contains public scientific data that should remain outside Git history.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

RECORD_URL = "https://zenodo.org/records/17143346"
FILES = {
    "injections": {
        "filename": "injections_Ninj_2e7_O5Like_snr20.h5",
        "url": f"{RECORD_URL}/files/injections_Ninj_2e7_O5Like_snr20.h5",
        "md5": "213ab2655b9fac5fb8cdc022fe8bd555",
        "size_bytes": 156_100_000,
        "role": "selection-effect injections",
    },
    "parent_catalog": {
        "filename": "parent_galaxy_catalog_1.6M_dVdz.h5",
        "url": f"{RECORD_URL}/files/parent_galaxy_catalog_1.6M_dVdz.h5",
        "md5": "c59d9feece4e0deab2b2b544cdb86303",
        "size_bytes": 133_800_000,
        "role": "parent galaxy sample",
    },
    "pixelated_pe": {
        "filename": "PE_O5Like_snr20_pixelated_nsidelist512-256-128-64-32-16-8_meanpixels15_skyconf0.9_nestFalse.h5",
        "url": f"{RECORD_URL}/files/PE_O5Like_snr20_pixelated_nsidelist512-256-128-64-32-16-8_meanpixels15_skyconf0.9_nestFalse.h5",
        "md5": "165914665a87ebb971fad2ed3bc66f15",
        "size_bytes": 192_200_000,
        "role": "pixelated GW posterior samples",
    },
    "nonpixelated_pe": {
        "filename": "PE_O5Like_snr20.h5",
        "url": f"{RECORD_URL}/files/PE_O5Like_snr20.h5",
        "md5": "f887599e7e11b19ce21039ecb51616c7",
        "size_bytes": 180_100_000,
        "role": "non-pixelated GW posterior samples",
    },
}


def hash_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_verified(url: str, destination: Path, expected_md5: str, timeout: int = 60) -> dict:
    """Stream to a temporary file; atomically publish only after MD5 verification."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=destination.name + ".", suffix=".part", dir=destination.parent)
    os.close(fd)
    temp_path = Path(temp_name)
    md5 = hashlib.md5()
    size = 0
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response, temp_path.open("wb") as output:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                output.write(chunk)
                md5.update(chunk)
                size += len(chunk)
        actual_md5 = md5.hexdigest()
        if actual_md5.lower() != expected_md5.lower():
            raise ValueError(f"MD5 mismatch for {destination.name}: expected {expected_md5}, got {actual_md5}")
        temp_path.replace(destination)
        return {"size_bytes": size, "md5": actual_md5, "sha256": hash_file(destination, "sha256")}
    finally:
        if temp_path.exists():
            temp_path.unlink()


def prepare(output_dir: Path, selected=("injections", "parent_catalog", "pixelated_pe"), timeout=60):
    unknown = sorted(set(selected) - set(FILES))
    if unknown:
        raise ValueError(f"Unknown CHIMERA data keys: {unknown}")
    manifest = {
        "source_record": RECORD_URL,
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
        "license_and_redistribution_reviewed": False,
        "note": "Do not commit these large data files to Git. Review upstream licence and benchmark terms before redistribution.",
        "files": {},
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    for key in selected:
        spec = FILES[key]
        destination = output_dir / spec["filename"]
        if destination.exists():
            actual_md5 = hash_file(destination, "md5")
            if actual_md5 != spec["md5"]:
                raise ValueError(f"Existing file has incorrect MD5: {destination}; move it aside and retry")
            result = {"size_bytes": destination.stat().st_size, "md5": actual_md5,
                      "sha256": hash_file(destination, "sha256"), "status": "already_verified"}
        else:
            result = download_verified(spec["url"], destination, spec["md5"], timeout)
            result["status"] = "downloaded_and_verified"
        manifest["files"][key] = {
            "filename": spec["filename"], "url": spec["url"], "role": spec["role"],
            "expected_md5": spec["md5"], **result,
        }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="external_data/chimera_v2")
    parser.add_argument("--files", nargs="+", choices=tuple(FILES),
                        default=["injections", "parent_catalog", "pixelated_pe"])
    parser.add_argument("--include-nonpixelated-pe", action="store_true")
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    selected = list(args.files)
    if args.include_nonpixelated_pe and "nonpixelated_pe" not in selected:
        selected.append("nonpixelated_pe")
    if args.dry_run:
        total = sum(FILES[key]["size_bytes"] for key in selected)
        print(json.dumps({"record": RECORD_URL, "selected": selected,
                          "approx_total_bytes": total, "approx_total_mib": round(total / 1024**2, 1)}, indent=2))
        return
    manifest = prepare(Path(args.output_dir), selected, args.timeout)
    print(json.dumps({"manifest": str(Path(args.output_dir) / "manifest.json"),
                      "verified_files": list(manifest["files"])}, indent=2))


if __name__ == "__main__":
    main()
