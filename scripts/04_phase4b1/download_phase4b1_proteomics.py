#!/usr/bin/env python3
"""Download only processed/search-output and small open supplementary files for Phase 4B.1."""

from __future__ import annotations

import hashlib
import json
import os
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "data/raw/phase4b1"
DEST.mkdir(parents=True, exist_ok=True)
CONTACT_EMAIL = os.environ.get("CONTACT_EMAIL", "").strip()
USER_AGENT = "hUCMSC-PE-Phase4B1/1.0" + (f" ({CONTACT_EMAIL})" if CONTACT_EMAIL else "")
CANDIDATES = ["ADAM17", "AGRN", "COL18A1", "DCN", "ENPP1", "FURIN", "GDF11", "GRN", "HSPG2", "MDK", "NAMPT", "NID1", "PSEN1", "SERPINE1", "TIMP1", "TIMP2", "WNT5A"]


def fetch(url: str, out: Path, accept: str = "*/*") -> Path:
    if out.exists() and out.stat().st_size > 0:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": accept})
    error = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=180) as response, out.open("wb") as handle:
                while True:
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    handle.write(block)
            return out
        except Exception as exc:
            error = exc
            out.unlink(missing_ok=True)
            time.sleep(2 ** attempt)
    raise RuntimeError(f"download failed: {url}: {error}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(manifest: list[dict], label: str, url: str, path: Path) -> None:
    manifest.append({"label": label, "url": url, "local_file": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha256(path)})


def main() -> int:
    manifest: list[dict] = []

    for pxd in ["PXD056371", "PXD020948", "PXD022174"]:
        url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{pxd}/files"
        path = fetch(url, DEST / f"{pxd}_pride_files.json", "application/json")
        record(manifest, f"{pxd}_FILE_INVENTORY", url, path)

    for pxd in ["PXD033723", "PXD036694"]:
        url = f"https://www.iprox.cn/proxi/datasets/{pxd}"
        path = fetch(url, DEST / f"{pxd}_iprox_proxi.json", "application/json")
        record(manifest, f"{pxd}_FILE_INVENTORY", url, path)

    selected = {
        "PXD056371_mzTab": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD056371/mzTab.mzTab",
        "PXD020948_MaxQuant_txt": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/04/PXD020948/txt.zip",
        "PXD033723_search_output": "https://download.iprox.org/IPX0004396000/IPX0004396001/1.txt",
    }
    for label, url in selected.items():
        suffix = ".zip" if url.endswith(".zip") else (".mztab" if "mzTab" in url else ".txt")
        if label == "PXD033723_search_output":
            zero = DEST / f"{label}{suffix}"
            manifest.append({"label": f"{label}_REPOSITORY_ZERO_BYTE", "url": url, "local_file": str(zero.relative_to(ROOT)), "bytes": zero.stat().st_size if zero.exists() else 0, "sha256": sha256(zero) if zero.exists() else "NOT_AVAILABLE", "error": "iProX endpoint advertises Content-Length: 0"})
            continue
        try:
            path = fetch(url, DEST / f"{label}{suffix}")
            record(manifest, label, url, path)
        except RuntimeError as exc:
            manifest.append({"label": f"{label}_FAILED", "url": url, "local_file": "NOT_DOWNLOADED", "bytes": 0, "sha256": "NOT_AVAILABLE", "error": str(exc)})

    pmcs = ["PMC11522688", "PMC7694919", "PMC7510303", "PMC13369709", "PMC13204574"]
    for pmc in pmcs:
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML"
        xml_path = fetch(url, DEST / f"{pmc}_fullText.xml", "application/xml")
        record(manifest, f"{pmc}_FULLTEXT_XML", url, xml_path)
        xml = xml_path.read_text(encoding="utf-8", errors="replace")
        hrefs = sorted(set(re.findall(r'(?:xlink:href|href)="([^"]+)"', xml)))
        for href in hrefs:
            name = Path(urllib.parse.urlparse(href).path).name
            lower = name.lower()
            is_supp = any(token in lower for token in ["moesm", "supp", "mmc", "esm", "additional", "data_s", "table_s"])
            allowed = lower.endswith((".xlsx", ".xls", ".csv", ".tsv", ".txt", ".docx", ".zip", ".pdf", ".tif"))
            if not name or not is_supp or not allowed:
                continue
            if href.startswith("http"):
                supp_url = href
            else:
                supp_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/bin/{name}"
            try:
                supp_path = fetch(supp_url, DEST / "supplements" / f"{pmc}_{name}")
                record(manifest, f"{pmc}_SUPPLEMENT", supp_url, supp_path)
            except RuntimeError as exc:
                manifest.append({"label": f"{pmc}_SUPPLEMENT_FAILED", "url": supp_url, "local_file": "NOT_DOWNLOADED", "bytes": 0, "sha256": "NOT_AVAILABLE", "error": str(exc)})

    query = " OR ".join(f"gene_exact:{g}" for g in CANDIDATES)
    fields = "accession,id,gene_names,protein_name"
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({"query": f"({query}) AND organism_id:9606 AND reviewed:true", "format": "tsv", "fields": fields, "size": 100})
    mapping = fetch(url, DEST / "uniprot_candidate_mapping.tsv", "text/tab-separated-values")
    record(manifest, "UNIPROT_REVIEWED_MAPPING", url, mapping)

    (DEST / "download_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    total = sum(int(x.get("bytes", 0)) for x in manifest)
    print(f"PHASE4B1_DOWNLOAD_OK files={sum(1 for x in manifest if int(x.get('bytes', 0)) > 0)} bytes={total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
