#!/usr/bin/env python3
"""Download the frozen public Phase 2A gene-set and regulon resources."""

from __future__ import annotations

import csv
import hashlib
import json
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/raw/phase2a_resources"

RESOURCES = {
    "h.all.v2026.1.Hs.symbols.gmt": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/h.all.v2026.1.Hs.symbols.gmt",
    "c2.cp.reactome.v2026.1.Hs.symbols.gmt": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/c2.cp.reactome.v2026.1.Hs.symbols.gmt",
    "c5.go.bp.v2026.1.Hs.symbols.gmt": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/c5.go.bp.v2026.1.Hs.symbols.gmt",
}


def download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    request = urllib.request.Request(url, headers={"User-Agent": "hUCMSC-PE-Phase2A/1.0"})
    with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as handle:
        while chunk := response.read(1024 * 1024):
            handle.write(chunk)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name, url in RESOURCES.items():
        path = OUT / name
        download(url, path)
        manifest.append({"resource": name, "url": url, "bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "version": "MSigDB 2026.1.Hs", "download_date": "2026-08-09"})
    query = urllib.parse.urlencode({"datasets": "collectri", "organisms": "9606", "fields": "sources,references,curation_effort", "genesymbols": "1", "format": "tsv"})
    url = "https://omnipathdb.org/interactions?" + query
    raw_path = OUT / "collectri_human_genesymbols_omnipath.tsv"
    download(url, raw_path)
    # Keep an exact raw snapshot and record its schema; parsing is done in R.
    with raw_path.open("r", encoding="utf-8", errors="replace") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    manifest.append({"resource": raw_path.name, "url": url, "bytes": raw_path.stat().st_size, "sha256": hashlib.sha256(raw_path.read_bytes()).hexdigest(), "version": "OmniPath CollecTRI snapshot 2026-08-09", "download_date": "2026-08-09", "columns": ";".join(header)})
    (OUT / "resource_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
