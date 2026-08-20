#!/usr/bin/env python3
"""Download small authoritative records for reproducible Phase 4B review."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "data/raw/phase4b"
DEST.mkdir(parents=True, exist_ok=True)
CANDIDATES = ["ADAM17", "AGRN", "COL18A1", "DCN", "ENPP1", "FURIN", "GDF11", "GRN", "HSPG2", "MDK", "NAMPT", "NID1", "PSEN1", "SERPINE1", "TIMP1", "TIMP2", "WNT5A"]
DATE = "2026-08-15"
TOOL = "hUCMSC_PE_Phase4B"
CONTACT_EMAIL = os.environ.get("CONTACT_EMAIL", "").strip()
USER_AGENT = f"{TOOL}/1.0" + (f" ({CONTACT_EMAIL})" if CONTACT_EMAIL else "")


def fetch(url: str, destination: Path, accept: str = "application/json") -> bytes:
    if destination.is_file() and destination.stat().st_size > 0:
        return destination.read_bytes()
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": accept})
    error = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                payload = response.read()
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
            return payload
        except Exception as exc:  # recorded and retried; no silent skipping
            error = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"download failed: {url}: {error}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def text_of(node: ET.Element | None) -> str:
    return "" if node is None else "".join(node.itertext()).strip()


def pubmed_records(xml_payload: bytes, candidate: str, family: str) -> list[dict[str, str]]:
    root = ET.fromstring(xml_payload)
    rows = []
    for article in root.findall(".//PubmedArticle"):
        citation = article.find("MedlineCitation")
        record = citation.find("Article") if citation is not None else None
        if citation is None or record is None:
            continue
        pmid = text_of(citation.find("PMID"))
        title = text_of(record.find("ArticleTitle"))
        abstract = " ".join(text_of(n) for n in record.findall("Abstract/AbstractText"))
        journal = text_of(record.find("Journal/Title"))
        year = text_of(record.find("Journal/JournalIssue/PubDate/Year")) or text_of(record.find("Journal/JournalIssue/PubDate/MedlineDate"))
        doi = ""
        for identifier in article.findall("PubmedData/ArticleIdList/ArticleId"):
            if identifier.attrib.get("IdType") == "doi":
                doi = text_of(identifier)
        publication_types = ";".join(text_of(x) for x in record.findall("PublicationTypeList/PublicationType"))
        rows.append({"candidate": candidate, "query_family": family, "PMID": pmid, "DOI": doi, "year": year, "title": title, "abstract": abstract, "journal": journal, "publication_types": publication_types, "source_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
    return rows


def main() -> int:
    manifest = []
    uniprot_rows = []
    for gene in CANDIDATES:
        query = f"(gene_exact:{gene}) AND (organism_id:9606) AND (reviewed:true)"
        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode({"query": query, "format": "json", "size": 5})
        path = DEST / "uniprot" / f"{gene}.json"
        payload = fetch(url, path)
        data = json.loads(payload)
        results = data.get("results", [])
        exact = [r for r in results if gene in [x.get("geneName", {}).get("value") for x in r.get("genes", [])]]
        chosen = exact[0] if exact else (results[0] if results else {})
        uniprot_rows.append({"candidate": gene, "uniprot_accession": chosen.get("primaryAccession", "NOT_FOUND"), "reviewed": chosen.get("entryType", ""), "record_count": len(results), "json_file": str(path.relative_to(ROOT)), "sha256": sha256(path), "source_url": url})
        manifest.append([str(path.relative_to(ROOT)), path.stat().st_size, sha256(path), url, DATE])
        time.sleep(0.15)
    with (DEST / "uniprot_record_registry.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(uniprot_rows[0]))
        writer.writeheader(); writer.writerows(uniprot_rows)

    query_templates = {
        "HUCMSC_PROTEIN_SOURCE": '("{gene}"[Title/Abstract]) AND ("umbilical cord mesenchymal"[Title/Abstract] OR "Wharton jelly mesenchymal"[Title/Abstract] OR hUCMSC[Title/Abstract] OR UCMSC[Title/Abstract]) AND (secretome[Title/Abstract] OR "conditioned medium"[Title/Abstract] OR proteomic*[Title/Abstract] OR ELISA[Title/Abstract] OR "western blot"[Title/Abstract] OR extracellular[Title/Abstract])',
        "EMPIRICAL_PERTURBATION": '("{gene}"[Title/Abstract]) AND (recombinant[Title/Abstract] OR treatment[Title/Abstract] OR stimulation[Title/Abstract] OR perturbation[Title/Abstract] OR knockout[Title/Abstract]) AND (transcriptom*[Title/Abstract] OR RNA-seq[Title/Abstract] OR microarray[Title/Abstract] OR expression[Title/Abstract]) AND (placent*[Title/Abstract] OR trophoblast*[Title/Abstract] OR macrophage*[Title/Abstract] OR fibroblast*[Title/Abstract] OR stromal[Title/Abstract])',
        "PE_DIRECT": '("{gene}"[Title/Abstract]) AND (preeclampsia[Title/Abstract] OR pre-eclampsia[Title/Abstract])',
        "PLACENTA_CONTEXT": '("{gene}"[Title/Abstract]) AND (placenta[Title/Abstract] OR trophoblast*[Title/Abstract] OR Hofbauer[Title/Abstract] OR decidual[Title/Abstract] OR "pregnancy hypertension"[Title/Abstract])',
        "MSC_PE_NOVELTY": '("{gene}"[Title/Abstract]) AND ("mesenchymal stem"[Title/Abstract] OR "mesenchymal stromal"[Title/Abstract] OR hUCMSC[Title/Abstract] OR UCMSC[Title/Abstract] OR "MSC-derived"[Title/Abstract]) AND (preeclampsia[Title/Abstract] OR pre-eclampsia[Title/Abstract])',
        "MSC_RELATED_MECHANISM": '("{gene}"[Title/Abstract]) AND ("mesenchymal stem"[Title/Abstract] OR "mesenchymal stromal"[Title/Abstract] OR hUCMSC[Title/Abstract] OR UCMSC[Title/Abstract]) AND (secretome[Title/Abstract] OR "conditioned medium"[Title/Abstract] OR extracellular[Title/Abstract] OR exosome*[Title/Abstract] OR vesicle*[Title/Abstract])'
    }
    search_log, all_records = [], []
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    for gene in CANDIDATES:
        for family, template in query_templates.items():
            query = template.format(gene=gene)
            params = {"db": "pubmed", "term": query, "retmode": "json", "retmax": 50, "sort": "relevance", "tool": TOOL}
            if CONTACT_EMAIL:
                params["email"] = CONTACT_EMAIL
            search_url = base + "esearch.fcgi?" + urllib.parse.urlencode(params)
            search_path = DEST / "pubmed" / f"{gene}__{family}__esearch.json"
            search_payload = fetch(search_url, search_path)
            search = json.loads(search_payload)["esearchresult"]
            ids = search.get("idlist", [])
            count = int(search.get("count", 0))
            records = []
            fetch_path = ""
            if ids:
                efetch_params = {"db": "pubmed", "id": ",".join(ids), "retmode": "xml", "tool": TOOL}
                if CONTACT_EMAIL:
                    efetch_params["email"] = CONTACT_EMAIL
                fetch_url = base + "efetch.fcgi?" + urllib.parse.urlencode(efetch_params)
                fetch_file = DEST / "pubmed" / f"{gene}__{family}__records.xml"
                fetch_payload = fetch(fetch_url, fetch_file, "application/xml")
                records = pubmed_records(fetch_payload, gene, family)
                fetch_path = str(fetch_file.relative_to(ROOT))
                manifest.append([fetch_path, fetch_file.stat().st_size, sha256(fetch_file), fetch_url, DATE])
            all_records.extend(records)
            search_log.append({"candidate": gene, "database_source": "PubMed_NCBI_EUTILS", "query_family": family, "query": query, "search_date": DATE, "result_count": count, "screened_records": len(records), "included_records": "PENDING_MANUAL_PRIMARY_STUDY_SCREEN", "inclusion_exclusion_reason": "All retrieved titles/abstracts retained for equal-depth screening; reviews not eligible as sole evidence", "search_url": search_url, "raw_search_file": str(search_path.relative_to(ROOT)), "raw_records_file": fetch_path, "source_url": search_url})
            manifest.append([str(search_path.relative_to(ROOT)), search_path.stat().st_size, sha256(search_path), search_url, DATE])
            time.sleep(0.36)
    with (DEST / "pubmed_search_log_raw.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(search_log[0]))
        writer.writeheader(); writer.writerows(search_log)
    with (DEST / "pubmed_candidate_records.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_records[0]) if all_records else ["candidate", "query_family", "PMID", "DOI", "year", "title", "abstract", "journal", "publication_types", "source_url"])
        writer.writeheader(); writer.writerows(all_records)

    pride_queries = ["umbilical cord mesenchymal secretome", "Wharton jelly mesenchymal secretome", "umbilical cord mesenchymal conditioned medium", "Wharton jelly mesenchymal proteomics"]
    pride_rows = []
    for index, query in enumerate(pride_queries, 1):
        url = "https://www.ebi.ac.uk/pride/ws/archive/v3/search/projects?" + urllib.parse.urlencode({"keyword": query, "page": 0, "pageSize": 100})
        path = DEST / "pride" / f"query_{index}.json"
        status, result_count = "DOWNLOADED", "UNRESOLVED_SCHEMA"
        try:
            payload = fetch(url, path)
            parsed = json.loads(payload)
            embedded = parsed.get("_embedded", {}) if isinstance(parsed, dict) else {}
            projects = embedded.get("projects", []) if isinstance(embedded, dict) else []
            result_count = parsed.get("page", {}).get("totalElements", len(projects)) if isinstance(parsed, dict) else len(parsed)
            manifest.append([str(path.relative_to(ROOT)), path.stat().st_size, sha256(path), url, DATE])
        except Exception as exc:
            status = f"API_ERROR:{type(exc).__name__}"
        pride_rows.append({"database_source": "PRIDE_ARCHIVE_API_V3", "query": query, "search_date": DATE, "result_count": result_count, "status": status, "raw_file": str(path.relative_to(ROOT)) if path.exists() else "", "source_url": url})
    with (DEST / "pride_search_log_raw.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pride_rows[0]))
        writer.writeheader(); writer.writerows(pride_rows)
    with (DEST / "download_manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle); writer.writerow(["file", "bytes", "sha256", "url", "download_date"]); writer.writerows(manifest)
    print(f"PHASE4B_DOWNLOAD_OK uniprot={len(uniprot_rows)} pubmed_searches={len(search_log)} records={len(all_records)} pride_queries={len(pride_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
