#!/usr/bin/env python3
"""Download only generic network/perturbation resources for blinded Phase 4A."""

from __future__ import annotations

import csv
import hashlib
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "data" / "raw" / "phase4a"
DEST.mkdir(parents=True, exist_ok=True)

RESOURCES = [
    (
        "ligand_target_matrix_nsga2r_final.rds",
        "https://zenodo.org/records/7074291/files/ligand_target_matrix_nsga2r_final.rds",
        "NicheNet-v2:Zenodo:7074291",
    ),
    (
        "omnipath_lr_crosscheck_20260815.tsv",
        "https://omnipathdb.org/interactions?format=tsv&genesymbols=1&datasets=omnipath,pathwayextra,ligrecextra&organisms=9606&fields=sources,references,curation_effort,datasets&license=academic",
        "OmniPath:LR-crosscheck:2026-08-15",
    ),
    (
        "omnipath_signed_activity_flow_20260815.tsv",
        "https://omnipathdb.org/interactions?format=tsv&genesymbols=1&datasets=omnipath&organisms=9606&directed=1&signed=1&fields=sources,references,curation_effort,datasets&license=academic",
        "OmniPath:signed-activity-flow:2026-08-15",
    ),
    (
        "collectri_signed_20260815.tsv",
        "https://omnipathdb.org/interactions?format=tsv&genesymbols=1&datasets=collectri&organisms=9606&directed=1&signed=1&resources=CollecTRI&fields=sources,references,curation_effort,datasets&license=academic",
        "OmniPath:CollecTRI:2026-08-15",
    ),
    (
        "cytosig_signature_centroid_core.tsv",
        "https://raw.githubusercontent.com/data2intelligence/CytoSig/master/CytoSig/signature.centroid",
        "CytoSig:core-signature:GitHub-master:2026-08-15",
    ),
]
KNOWN_MD5 = {"ligand_target_matrix_nsga2r_final.rds": "b09606b04b2d4490418d9028c0e58b9f"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    if destination.is_file() and destination.stat().st_size > 0:
        print(f"CACHED {destination.name} ({destination.stat().st_size} bytes)", flush=True)
        return
    part = destination.with_suffix(destination.suffix + ".part")
    command = ["curl.exe", "-L", "--fail", "--retry", "4", "--retry-delay", "3", "-C", "-", "-o", str(part), url]
    print(f"DOWNLOAD {url}", flush=True)
    completed = subprocess.run(command, cwd=ROOT)
    if completed.returncode:
        raise RuntimeError(f"curl failed ({completed.returncode}): {url}")
    os.replace(part, destination)


def main() -> int:
    manifest: list[dict[str, object]] = []
    for filename, url, source in RESOURCES:
        path = DEST / filename
        download(url, path)
        if filename in KNOWN_MD5 and md5(path) != KNOWN_MD5[filename]:
            path.unlink()
            download(url, path)
            if md5(path) != KNOWN_MD5[filename]:
                raise RuntimeError(f"official MD5 mismatch after clean retry: {filename}")
        manifest.append({
            "filename": filename,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "url": url,
            "source_accession": source,
            "download_date": "2026-08-15",
        })
    manifest_path = DEST / "download_manifest.csv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)
    print(f"Phase 4A resources ready: {len(manifest)} files; {sum(int(r['bytes']) for r in manifest)} bytes", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
