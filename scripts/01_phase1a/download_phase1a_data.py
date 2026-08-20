#!/usr/bin/env python3
"""Download Phase 1A metadata and processed expression evidence; never FASTQ/SRA archives."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import time
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
P1RAW = RAW / "phase1a"
CFG = json.loads((ROOT / "config" / "phase1a_sources.json").read_text(encoding="utf-8"))
USER_AGENT = "hUCMSC-PE-Phase1A-freeze/1.0 (processed-data and metadata audit)"


def fetch(url: str, path: Path, force: bool = False) -> None:
    if path.exists() and path.stat().st_size > 0 and not force:
        return
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    error: Exception | None = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("wb") as handle:
                    shutil.copyfileobj(response, handle)
            if path.stat().st_size == 0:
                raise RuntimeError(f"empty response from {url}")
            return
        except Exception as exc:
            error = exc
            if path.exists():
                path.unlink()
            time.sleep(2**attempt)
    raise RuntimeError(f"failed to download {url}: {error}")


def geo_soft_url(accession: str) -> str:
    stem = accession[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{accession}/soft/{accession}_family.soft.gz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    P1RAW.mkdir(parents=True, exist_ok=True)

    downloads: list[tuple[str, Path, str, str]] = []
    for accession in CFG["geo_soft_accessions"]:
        downloads.append((geo_soft_url(accession), RAW / f"{accession}_family.soft.gz", "GEO_SOFT", "PROGRAMMATIC"))
    for pmc in CFG["europe_pmc_full_text"]:
        downloads.append((
            f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML",
            P1RAW / f"{pmc}_fulltext.xml",
            "EUROPE_PMC_XML", "PROGRAMMATIC",
        ))
    for pmc in CFG["europe_pmc_supplements"]:
        downloads.append((
            f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/supplementaryFiles",
            P1RAW / f"{pmc}_supplementary.zip",
            "EUROPE_PMC_SUPPLEMENT", "PROGRAMMATIC",
        ))
    downloads.append((
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({
            "db": "pubmed", "id": ",".join(CFG["pubmed_ids"]), "retmode": "xml",
        }),
        P1RAW / "pubmed_phase1a.xml",
        "PUBMED_XML", "PROGRAMMATIC",
    ))
    for item in CFG["files"]:
        downloads.append((item["url"], P1RAW / item["name"], item["kind"], item.get("retrieval", "PROGRAMMATIC")))

    manifest = []
    for index, (url, path, kind, retrieval) in enumerate(downloads, 1):
        print(f"[{index}/{len(downloads)}] {kind}: {path.name}")
        if retrieval == "BROWSER_OPTIONAL":
            if not path.exists():
                print("  optional browser-acquired supplement absent; continuing because the freeze build does not depend on it")
                continue
            if path.read_bytes()[:2] != b"PK":
                raise RuntimeError(f"optional Office supplement is not a valid ZIP-based Office file: {path}")
        else:
            fetch(url, path, args.force)
        manifest.append({
            "kind": kind,
            "path": str(path.relative_to(ROOT)).replace("\\", "/"),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "source_url": url,
            "retrieval": retrieval,
        })
    (P1RAW / "download_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"Phase 1A cache ready: {P1RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
