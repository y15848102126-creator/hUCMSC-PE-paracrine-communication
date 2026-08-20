#!/usr/bin/env python3
"""Download only metadata and small audit evidence; never expression FASTQ/SRA files."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
CONFIG = ROOT / "config" / "audit_accessions.json"
USER_AGENT = "hUCMSC-PE-phase0-audit/1.0 (metadata-only)"


def geo_prefix(accession: str) -> str:
    return accession[:-3] + "nnn"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fetch(url: str, destination: Path, force: bool = False, retries: int = 4) -> None:
    if destination.exists() and destination.stat().st_size > 0 and not force:
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".partial")
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as output:
                while True:
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    output.write(block)
            partial.replace(destination)
            return
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_error = exc
            if partial.exists():
                partial.unlink()
            time.sleep(2**attempt)
    raise RuntimeError(f"Failed to download {url}: {last_error}")


def head_size(url: str) -> tuple[str, str]:
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return str(response.headers.get("Content-Length", "UNRESOLVED")), str(response.status)
    except Exception as exc:  # recorded, never silently discarded
        return "UNRESOLVED", f"ERROR:{type(exc).__name__}"


def series_supplement_urls(soft_path: Path) -> list[str]:
    urls: list[str] = []
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Series_supplementary_file = "):
                urls.append(line.rstrip().split(" = ", 1)[1].replace("ftp://", "https://"))
    return urls


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="redownload existing metadata")
    args = parser.parse_args()
    cfg = json.loads(CONFIG.read_text(encoding="utf-8"))
    RAW.mkdir(parents=True, exist_ok=True)

    downloads: list[tuple[str, Path, str]] = []
    for accession in cfg["main_accessions"] + cfg["support_accessions"]:
        url = (
            f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_prefix(accession)}/"
            f"{accession}/soft/{accession}_family.soft.gz"
        )
        downloads.append((url, RAW / f"{accession}_family.soft.gz", "GEO_SOFT"))

    for accession in cfg["sra_accessions"]:
        downloads.append((
            f"https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc={accession}",
            RAW / f"{accession}_sra_runinfo.csv",
            "SRA_RUNINFO",
        ))

    pmids = ",".join(cfg["pubmed_ids"])
    downloads.append((
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=pubmed&id={pmids}&retmode=xml",
        RAW / "pubmed_batch.xml",
        "PUBMED_XML",
    ))
    for pmc in cfg["pmc_full_text_ids"]:
        downloads.append((
            f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML",
            RAW / f"{pmc}_fulltext.xml",
            "EUROPE_PMC_XML",
        ))

    downloads.extend([
        (
            "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8715893/supplementaryFiles",
            RAW / "PMC8715893_supplementary.zip",
            "EUROPE_PMC_SUPPLEMENT",
        ),
        (
            "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12092795/supplementaryFiles",
            RAW / "PMC12092795_supplementary.zip",
            "EUROPE_PMC_SUPPLEMENT",
        ),
        (
            "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5519nnn/GSM5519462/suppl/"
            "GSM5519462_U01_barcodes.tsv.gz",
            RAW / "GSM5519462_U01_barcodes.tsv.gz",
            "GEO_BARCODE_EVIDENCE",
        ),
        (
            "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5519nnn/GSM5519463/suppl/"
            "GSM5519463_U02_barcodes.tsv.gz",
            RAW / "GSM5519463_U02_barcodes.tsv.gz",
            "GEO_BARCODE_EVIDENCE",
        ),
    ])

    for index, (url, path, kind) in enumerate(downloads, 1):
        print(f"[{index}/{len(downloads)}] {kind}: {path.name}", flush=True)
        # PubMed records can gain identifiers/corrections and the requested PMID list
        # may change; this small batch is intentionally refreshed on every online run.
        fetch(url, path, force=args.force or kind == "PUBMED_XML")

    manifest_rows = []
    for url, path, kind in downloads:
        manifest_rows.append({
            "kind": kind,
            "local_file": path.name,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "source_url": url,
        })
    with (RAW / "download_manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_rows[0].keys())
        writer.writeheader()
        writer.writerows(manifest_rows)

    size_rows = []
    for accession in cfg["main_accessions"]:
        for url in series_supplement_urls(RAW / f"{accession}_family.soft.gz"):
            length, status = head_size(url)
            size_rows.append({
                "geo_accession": accession,
                "supplement_url": url,
                "content_length_bytes": length,
                "http_status": status,
            })
    with (RAW / "supplement_sizes.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=size_rows[0].keys())
        writer.writeheader()
        writer.writerows(size_rows)

    print(f"Metadata cache ready: {RAW}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
