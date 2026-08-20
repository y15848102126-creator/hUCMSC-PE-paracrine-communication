#!/usr/bin/env python3
"""Download only public sender-side Phase 3 inputs; never reads receiver outputs."""

from __future__ import annotations

import csv
import gzip
import hashlib
import os
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlparse


ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "data" / "raw" / "phase3"
DEST.mkdir(parents=True, exist_ok=True)


def read_soft(accession: str) -> str:
    with gzip.open(ROOT / "data" / "raw" / f"{accession}_family.soft.gz", "rt", encoding="utf-8", errors="replace") as handle:
        return handle.read()


def sample_urls(accession: str, keep_title) -> list[tuple[str, str]]:
    text = read_soft(accession)
    out: list[tuple[str, str]] = []
    for block in re.split(r"(?=\^SAMPLE = )", text):
        gsm = re.search(r"^\^SAMPLE = (\S+)", block, re.M)
        title = re.search(r"^!Sample_title = (.+)$", block, re.M)
        if not gsm or not title or not keep_title(title.group(1).strip()):
            continue
        for match in re.finditer(r"^!Sample_supplementary_file_\d+ = (\S+)", block, re.M):
            url = match.group(1).replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
            out.append((url, f"{accession}:{gsm.group(1)}:{title.group(1).strip()}"))
    return out


def series_urls(accession: str) -> list[tuple[str, str]]:
    text = read_soft(accession)
    return [
        (m.group(1).replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov"), f"{accession}:series_supplement")
        for m in re.finditer(r"^!Series_supplementary_file = (\S+)", text, re.M)
    ]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, destination: Path) -> None:
    invalid_omnipath = destination.name.startswith("omnipath_intercell") and destination.is_file() and b"Unknown argument" in destination.read_bytes()[:1000]
    if destination.is_file() and destination.stat().st_size > 0 and not invalid_omnipath:
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
    records: list[tuple[str, str, str]] = []
    for url, source in sample_urls("GSE182158", lambda _: True):
        records.append((url, Path(urlparse(url).path).name, source))
    for url, source in sample_urls("GSE199071", lambda title: title.upper().startswith("HUCMSC")):
        records.append((url, Path(urlparse(url).path).name, source))
    for url, source in series_urls("GSE117837"):
        records.append((url, Path(urlparse(url).path).name, source))
    records.extend([
        ("https://zenodo.org/records/7074291/files/lr_network_human_21122021.rds", "lr_network_human_21122021.rds", "NicheNet:Zenodo:7074291"),
        ("https://omnipathdb.org/intercell?format=tsv&fields=genesymbol&entity_types=protein&license=academic", "omnipath_intercell_20260810.tsv", "OmniPath:intercell_API:2026-08-10"),
    ])
    if len(records) != 33 + 12 + 1 + 2:
        raise AssertionError(f"unexpected resource count: {len(records)}")
    manifest = []
    for url, filename, source in records:
        path = DEST / filename
        download(url, path)
        manifest.append({
            "filename": filename,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "url": url,
            "source_accession": source,
            "download_date": "2026-08-10",
        })
    with (DEST / "download_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)
    print(f"Phase 3 sender inputs ready: {len(manifest)} files; {sum(r['bytes'] for r in manifest)} bytes", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
