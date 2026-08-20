#!/usr/bin/env python3
"""Download small, machine-readable Phase 2A.2 provenance and GEO metadata."""

from __future__ import annotations

import hashlib
import json
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/raw/phase2a2"
OUT.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "hUCMSC-PE-reproducibility-audit/2A2 (metadata only)"}

URLS = {
    "admati_figshare_23264102.json": "https://api.figshare.com/v2/articles/23264102",
    "admati_figshare_23264165.json": "https://api.figshare.com/v2/articles/23264165",
    "admati_github_tree.json": "https://api.github.com/repos/zeiselamit/PE_2023/git/trees/main?recursive=1",
    "GSE282038_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282038/soft/GSE282038_family.soft.gz",
    "GSE267340_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267340/soft/GSE267340_family.soft.gz",
    "GSE298119_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298119/soft/GSE298119_family.soft.gz",
    "GSE282038_suppl_listing.html": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282038/suppl/",
    "GSE267340_suppl_listing.html": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267340/suppl/",
    "GSE298119_suppl_listing.html": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298119/suppl/",
    "zheng_frontiers_full.html": "https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2025.1638603/full",
}

SEARCHES = {
    "admati_gds_esearch.json": ("gds", '37572658[PMID] OR Admati[Author] OR \"10.1016/j.medj.2023.07.005\"'),
    "admati_sra_esearch.json": ("sra", '37572658[All Fields] OR \"10.1016/j.medj.2023.07.005\"[All Fields] OR \"Inbal Admati\"[All Fields]'),
    "admati_bioproject_esearch.json": ("bioproject", '37572658[All Fields] OR \"10.1016/j.medj.2023.07.005\"[All Fields] OR \"Inbal Admati\"[All Fields]'),
    "admati_biosample_esearch.json": ("biosample", '37572658[All Fields] OR \"10.1016/j.medj.2023.07.005\"[All Fields] OR \"Inbal Admati\"[All Fields]'),
}


def fetch(url: str, path: Path) -> dict[str, object]:
    request = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(request, timeout=180) as response:
        payload = response.read()
    path.write_bytes(payload)
    return {
        "filename": path.name,
        "url": url,
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def main() -> int:
    manifest: list[dict[str, object]] = []
    for filename, url in URLS.items():
        manifest.append(fetch(url, OUT / filename))
    for filename, (database, term) in SEARCHES.items():
        query = urllib.parse.urlencode({"db": database, "term": term, "retmode": "json", "retmax": 100})
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + query
        manifest.append(fetch(url, OUT / filename))
    (OUT / "metadata_download_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"downloaded": len(manifest), "bytes": sum(int(x["bytes"]) for x in manifest)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
