#!/usr/bin/env python3
"""Download Phase 0B evidence; processed expression is opt-in, never FASTQ/SRA archives."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "phase0b"
CFG = json.loads((ROOT / "config" / "phase0b_sources.json").read_text(encoding="utf-8"))
USER_AGENT = "hUCMSC-PE-Phase0B-audit/1.0 (metadata feasibility audit)"


def fetch(url: str, path: Path, force: bool = False) -> None:
    if path.exists() and path.stat().st_size > 0 and not force:
        return
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                payload = response.read()
            if not payload:
                raise RuntimeError(f"empty response from {url}")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
            return
        except Exception as exc:  # network retry is bounded and auditable
            last_error = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"failed to download {url}: {last_error}")


def geo_soft_url(accession: str) -> str:
    stem = accession[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{accession}/soft/{accession}_family.soft.gz"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--include-processed",
        action="store_true",
        help="also cache the public Admati processed matrix (~257 MB compressed)",
    )
    args = parser.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    downloads: list[tuple[str, Path, str]] = []

    for pmc in CFG["pmc_full_text_ids"]:
        downloads.append((
            f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML",
            RAW / f"{pmc}_fulltext.xml", "EUROPE_PMC_XML",
        ))
    for pmc in CFG["pmc_supplement_ids"]:
        downloads.append((
            f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/supplementaryFiles",
            RAW / f"{pmc}_supplementary.zip", "EUROPE_PMC_SUPPLEMENT",
        ))
    for accession in CFG["geo_soft_accessions"]:
        downloads.append((geo_soft_url(accession), RAW / f"{accession}_family.soft.gz", "GEO_SOFT"))

    pmid_query = ",".join(CFG["pubmed_ids"])
    downloads.append((
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode({
            "db": "pubmed", "id": pmid_query, "retmode": "xml",
        }), RAW / "pubmed_phase0b.xml", "PUBMED_XML",
    ))
    discovery_terms = {
        "gds_pe_single_cell": '(preeclampsia OR pre-eclampsia) AND (single cell OR single-cell OR single nucleus OR snRNA-seq) AND Homo sapiens[Organism]',
        "gds_78_subject_title": '"Single-cell mapping of maternal-fetal cross-talk in preeclampsia"',
        "gds_hucmsc_single_cell": '(umbilical cord OR Wharton jelly) AND (mesenchymal stem cell OR mesenchymal stromal cell) AND (single cell OR single-cell)',
    }
    for name, term in discovery_terms.items():
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode({
            "db": "gds", "term": term, "retmode": "json", "retmax": 500,
        })
        downloads.append((url, RAW / f"{name}_esearch.json", "NCBI_GDS_ESEARCH"))
    for version in CFG["dbgap_versions"]:
        url = "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?" + urllib.parse.urlencode({"study_id": version})
        downloads.append((url, RAW / f"{version}.html", "DBGAP_PUBLIC_STUDY_PAGE"))
        downloads.extend([
            (f"https://www.ncbi.nlm.nih.gov/gap/sstr/api/v1/study/{version}/summary",
             RAW / f"{version}_sstr_summary.json", "DBGAP_SSTR_API"),
            (f"https://www.ncbi.nlm.nih.gov/gap/sstr/api/v1/study/{version}/subjects?data_format=row&length=1000&start=0",
             RAW / f"{version}_sstr_subjects.json", "DBGAP_SSTR_API"),
        ])
    downloads.extend([
        ("https://api.github.com/repos/JustMoveOnnn/preeclampsia/contents/single_cell_matrix/data",
         RAW / "yang2023_github_contents.json", "GITHUB_API"),
        ("https://api.github.com/repos/hypaik/PePT_vignette/git/trees/main?recursive=1",
         RAW / "gse290578_github_tree.json", "GITHUB_API"),
        ("https://api.github.com/repos/zeiselamit/PE_2023/git/trees/main?recursive=1",
         RAW / "admati2023_github_tree.json", "GITHUB_API"),
        ("https://api.figshare.com/v2/articles/23264102",
         RAW / "admati2023_figshare_article.json", "FIGSHARE_API"),
        ("https://api.figshare.com/v2/articles/23264165",
         RAW / "admati2023_trophoblast_figshare_article.json", "FIGSHARE_API"),
        ("https://raw.githubusercontent.com/JustMoveOnnn/preeclampsia/main/single_cell_matrix/data/readme.md",
         RAW / "yang2023_readme.md", "GITHUB_RAW"),
        ("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=GSE290578",
         RAW / "GSE290578_sra_runinfo.csv", "SRA_RUNINFO"),
        ("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA643879",
         RAW / "PRJNA643879_sra_runinfo.csv", "SRA_RUNINFO"),
        ("https://cell.ucsf.edu/snPlacenta/", RAW / "normal_atlas_portal.html", "PUBLIC_DATA_PORTAL"),
        ("https://assets-eu.researchsquare.com/files/rs-8254581/v1_covered_ac1f2319-439f-4f0e-b096-2c6a7bfeaf46.pdf",
         RAW / "rs-8254581-v1.pdf", "PREPRINT_PDF"),
        ("https://assets-eu.researchsquare.com/files/rs-8254581/v1/838318f119a18d771dcb6b4f.pdf",
         RAW / "rs-8254581-v1_supplementary_information.pdf", "PREPRINT_SUPPLEMENT_PDF"),
        ("https://assets-eu.researchsquare.com/files/rs-8254581/v1/fa4398f676e015653b4e6f75.xlsx",
         RAW / "rs-8254581-v1_supplementary_tables_3_to_8.xlsx", "PREPRINT_SUPPLEMENT_XLSX"),
        ("https://ngdc.cncb.ac.cn/gsa-human/browse/HRA005090",
         RAW / "HRA005090.html", "GSA_HUMAN_STUDY_PAGE"),
        ("https://ngdc.cncb.ac.cn/gsa-human/ajaxb/indinstudy?accession=HRA005090",
         RAW / "HRA005090_individuals.json", "GSA_HUMAN_API"),
        ("https://ngdc.cncb.ac.cn/gsa-human/ajaxb/runinstudy?accession=HRA005090",
         RAW / "HRA005090_runs.json", "GSA_HUMAN_API"),
        ("https://ngdc.cncb.ac.cn/gsa-human/browse/HRA004699",
         RAW / "HRA004699.html", "GSA_HUMAN_STUDY_PAGE"),
        ("https://ngdc.cncb.ac.cn/gsa-human/browse/HRA003297",
         RAW / "HRA003297.html", "GSA_HUMAN_STUDY_PAGE"),
        ("https://db.cngb.org/search/project/CNP0000562/",
         RAW / "CNP0000562.html", "CNGB_PROJECT_PAGE"),
        ("https://api.elsevier.com/content/article/PII:S0161589023001402?httpAccept=text%2Fxml",
         RAW / "Jiao2023_elsevier.xml", "PUBLISHER_ARTICLE_XML"),
    ])
    if args.include_processed:
        downloads.append((
            "https://ndownloader.figshare.com/files/41003240",
            RAW / "sc_PE_allcells_with_metadata_29-May-2023.txt.zip",
            "PUBLIC_PROCESSED_MATRIX",
        ))

    manifest = []
    for index, (url, path, kind) in enumerate(downloads, 1):
        print(f"[{index}/{len(downloads)}] {kind}: {path.name}", flush=True)
        fetch(url, path, force=args.force or kind in {"NCBI_GDS_ESEARCH", "PUBMED_XML", "DBGAP_PUBLIC_STUDY_PAGE", "GITHUB_API", "PUBLIC_DATA_PORTAL", "SRA_RUNINFO"})
        manifest.append({
            "kind": kind, "local_file": path.name, "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "source_url": url,
        })
    with (RAW / "download_manifest.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)
    print(f"Phase 0B metadata cache ready: {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
