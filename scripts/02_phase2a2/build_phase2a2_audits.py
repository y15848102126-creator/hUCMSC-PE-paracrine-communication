#!/usr/bin/env python3
"""Build source-traceable Phase 2A.2 provenance and external-cohort audits."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/phase2a2"
OUT = ROOT / "results/02_phase2a2"
ADM = "https://doi.org/10.6084/m9.figshare.23264102.v1"
CODE = "https://github.com/zeiselamit/PE_2023"
PAPER = "https://doi.org/10.1016/j.medj.2023.07.005"
ZHENG = "https://doi.org/10.3389/fimmu.2025.1638603"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def esearch_count(name: str) -> int:
    payload = json.loads((RAW / name).read_text(encoding="utf-8"))
    return int(payload["esearchresult"]["count"])


def main() -> int:
    config = json.loads((ROOT / "config/phase2a2_analysis.json").read_text(encoding="utf-8"))
    for rel, expected in config["history_policy"]["frozen_sha256"].items():
        observed = sha256(ROOT / rel)
        if observed != expected:
            raise RuntimeError(f"Frozen history changed: {rel}: {observed} != {expected}")

    fig = json.loads((RAW / "admati_figshare_23264102.json").read_text(encoding="utf-8"))
    nuclei = json.loads((RAW / "admati_figshare_23264165.json").read_text(encoding="utf-8"))
    public_zip = ROOT / "data/raw/phase0b/sc_PE_allcells_with_metadata_29-May-2023.txt.zip"
    with zipfile.ZipFile(public_zip) as archive:
        member = archive.infolist()[0]
        if len(archive.infolist()) != 1:
            raise RuntimeError("Unexpected Figshare ZIP member count")
    file_info = fig["files"][0]
    audit = [
        {
            "record_type": "FIGSHARE_ARTICLE_INVENTORY", "repository": "Figshare", "article_or_accession": "23264102.v1",
            "filename": file_info["name"], "file_id": file_info["id"], "bytes": file_info["size"],
            "checksum_type": "MD5 supplied/computed", "checksum": file_info["computed_md5"], "format": file_info["mimetype"],
            "used_in_phase2a": "YES", "expression_state": "DESCRIPTION_CALLS_UMI_COUNTS_BUT_VALUES_ARE_NORMALIZED_AND_CEILED",
            "evidence": re.sub("<[^>]+>", "", fig["description"]), "source_url": ADM,
        },
        {
            "record_type": "LOCAL_DOWNLOADED_FILE", "repository": "Figshare download", "article_or_accession": "file:41003240",
            "filename": public_zip.name, "file_id": 41003240, "bytes": public_zip.stat().st_size,
            "checksum_type": "SHA256", "checksum": sha256(public_zip), "format": "ZIP",
            "used_in_phase2a": "YES", "expression_state": "PER_CELL_10000_NORMALIZED_THEN_CEIL",
            "evidence": "Exact Phase 2A input archive; checksum matches frozen config.", "source_url": ADM,
        },
        {
            "record_type": "ZIP_MEMBER", "repository": "Figshare download", "article_or_accession": "file:41003240",
            "filename": member.filename, "file_id": "", "bytes": member.file_size,
            "checksum_type": "ZIP_CRC32", "checksum": f"{member.CRC:08x}", "format": "tab-delimited genes x cells plus metadata",
            "used_in_phase2a": "YES", "expression_state": "PER_CELL_10000_NORMALIZED_THEN_CEIL",
            "evidence": "Sole member; Phase 2A parser read this table.", "source_url": ADM,
        },
        {
            "record_type": "NUMERICAL_LAYER_AUDIT", "repository": "Local reproducibility audit", "article_or_accession": "23264102.v1",
            "filename": member.filename, "file_id": "", "bytes": member.file_size, "checksum_type": "", "checksum": "", "format": "integer-valued expression",
            "used_in_phase2a": "YES", "expression_state": "PER_CELL_10000_NORMALIZED_THEN_CEIL",
            "evidence": "86,752/86,752 cells matched the ceil-to-10,000 signature; matrix sum 962,152,952 versus published total_molecules sum 646,369,597.",
            "source_url": "results/02_phase2a1/phase2a1_risk_flags.csv|data/interim/phase2a1/count_layer_provenance_summary.csv",
        },
        {
            "record_type": "AUTHOR_INTERNAL_LOADER", "repository": "GitHub", "article_or_accession": "zeiselamit/PE_2023",
            "filename": "load_sc_PE_data_and_save_v1.m", "file_id": "", "bytes": "", "checksum_type": "", "checksum": "", "format": "MATLAB source",
            "used_in_phase2a": "PROVENANCE_ONLY", "expression_state": "INTERNAL_CELL_RANGER_FILTERED_FEATURE_BC_INPUTS",
            "evidence": "Loader references internal matrix.mtx, barcodes.tsv and features.tsv paths; these files are absent from the public GitHub tree.", "source_url": CODE,
        },
        {
            "record_type": "AUTHOR_DOWNSTREAM_TRANSFORM", "repository": "GitHub", "article_or_accession": "zeiselamit/PE_2023",
            "filename": "multiple public MATLAB scripts", "file_id": "", "bytes": "", "checksum_type": "", "checksum": "", "format": "MATLAB source",
            "used_in_phase2a": "PROVENANCE_ONLY", "expression_state": "CEIL(DATA/COLUMN_SUM*10000)",
            "evidence": "Public scripts use data = ceil(data./repmat(sum(data),length(data(:,1)),1)*10e3), matching the deposited value behavior.", "source_url": CODE,
        },
        {
            "record_type": "LINKED_FIGSHARE_ARTICLE", "repository": "Figshare", "article_or_accession": "23264165.v1",
            "filename": nuclei["files"][0]["name"], "file_id": nuclei["files"][0]["id"], "bytes": nuclei["files"][0]["size"],
            "checksum_type": "MD5", "checksum": nuclei["files"][0]["computed_md5"], "format": nuclei["files"][0]["mimetype"],
            "used_in_phase2a": "NO", "expression_state": "SEPARATE_TROPHOBLAST_NUCLEI_DATASET_NOT_RAW_SCRNA_COUNTS",
            "evidence": "The linked article is the separate trophoblast single-nucleus dataset, not the underlying Admati scRNA Cell Ranger files.",
            "source_url": "https://doi.org/10.6084/m9.figshare.23264165.v1",
        },
    ]
    write_csv(OUT / "provenance/admati_expression_layer_audit.csv", audit)

    counts = {k: esearch_count(k) for k in [
        "admati_gds_esearch.json", "admati_sra_esearch.json", "admati_bioproject_esearch.json", "admati_biosample_esearch.json"
    ]}
    search = [
        {"repository": "Figshare 23264102.v1", "query_or_scope": "complete article file inventory", "result": "ONE_NORMALIZED_CEILED_TABLE", "public_raw_umi": "NO", "cell_ranger_matrix": "NO", "fastq": "NO", "evidence": "One ZIP only; article has no related materials or linked files.", "source_url": ADM},
        {"repository": "Figshare linked article 23264165.v1", "query_or_scope": "author-linked nuclei article", "result": "SEPARATE_SNRNA_TABLE", "public_raw_umi": "NO_FOR_ADMATI_SCRNA", "cell_ranger_matrix": "NO", "fastq": "NO", "evidence": "Separate trophoblast nuclei dataset.", "source_url": "https://doi.org/10.6084/m9.figshare.23264165.v1"},
        {"repository": "GitHub", "query_or_scope": "complete public repository tree", "result": "CODE_ONLY_NO_MATRIX_OR_FASTQ", "public_raw_umi": "NO", "cell_ranger_matrix": "NO", "fastq": "NO", "evidence": "Internal filtered_feature_bc_matrix paths are referenced but not committed.", "source_url": CODE},
        {"repository": "NCBI GEO/GDS", "query_or_scope": "PMID, DOI and author ESearch", "result": f"NOT_FOUND_COUNT_{counts['admati_gds_esearch.json']}", "public_raw_umi": "NO_RECORD_FOUND", "cell_ranger_matrix": "NO_RECORD_FOUND", "fastq": "NO_RECORD_FOUND", "evidence": "Programmatic NCBI ESearch on audit date.", "source_url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds"},
        {"repository": "NCBI SRA", "query_or_scope": "PMID, DOI and author ESearch", "result": f"NOT_FOUND_COUNT_{counts['admati_sra_esearch.json']}", "public_raw_umi": "NO_RECORD_FOUND", "cell_ranger_matrix": "NO_RECORD_FOUND", "fastq": "NO_RECORD_FOUND", "evidence": "Programmatic NCBI ESearch on audit date.", "source_url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra"},
        {"repository": "NCBI BioProject", "query_or_scope": "PMID, DOI and author ESearch", "result": f"NOT_FOUND_COUNT_{counts['admati_bioproject_esearch.json']}", "public_raw_umi": "NO_RECORD_FOUND", "cell_ranger_matrix": "NO_RECORD_FOUND", "fastq": "NO_RECORD_FOUND", "evidence": "Programmatic NCBI ESearch on audit date.", "source_url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject"},
        {"repository": "NCBI BioSample", "query_or_scope": "PMID, DOI and author ESearch", "result": f"NOT_FOUND_COUNT_{counts['admati_biosample_esearch.json']}", "public_raw_umi": "NO_RECORD_FOUND", "cell_ranger_matrix": "NO_RECORD_FOUND", "fastq": "NO_RECORD_FOUND", "evidence": "Programmatic NCBI ESearch on audit date.", "source_url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=biosample"},
        {"repository": "ENA/EGA", "query_or_scope": "publication DOI/PMID and author-linked accessions", "result": "NO_AUTHOR_LINKED_ACCESSION_IDENTIFIED", "public_raw_umi": "NO_RECORD_FOUND", "cell_ranger_matrix": "NO_RECORD_FOUND", "fastq": "NO_RECORD_FOUND", "evidence": "No ENA/EGA accession is given by the paper, Figshare record or author repository; no NCBI study alias was found for ENA mirroring.", "source_url": PAPER},
        {"repository": "Published data availability", "query_or_scope": "MedJ article statement", "result": "EXPRESSION_FIGSHARE_CODE_GITHUB_ADDITIONAL_INFO_ON_REQUEST", "public_raw_umi": "NO", "cell_ranger_matrix": "NO", "fastq": "NO", "evidence": "On-request material is not publicly recoverable.", "source_url": PAPER},
        {"repository": "AUDIT_CONCLUSION", "query_or_scope": "all sources above", "result": "PUBLIC_MATRIX_NORMALIZED_RAW_NOT_PUBLIC", "public_raw_umi": "NO", "cell_ranger_matrix": "NO", "fastq": "NO", "evidence": "The description/code contradiction is documented, not guessed away.", "source_url": f"{ADM}|{CODE}|{PAPER}"},
    ]
    write_csv(OUT / "provenance/admati_raw_count_search.csv", search)

    zheng_samples = [
        ("GSE282038", "GSM8634701", "PE001", "EOPE", "31+4", "YES", "UNRESOLVED", "CESAREAN; paper says elective", "PRJNA1186695", "SRX26736654", "GEO labels EOPE; Phase 0 paper supplement identifies FGR and GA"),
        ("GSE298119", "GSM9008678", "PE002", "PE_CASE", "UNRESOLVED", "UNRESOLVED", "UNRESOLVED", "CESAREAN; elective only at paper cohort level", "PRJNA1268385", "SRX28967107", "GEO series calls this a term PE cohort; EOPE attribution is not sample-resolved"),
        ("GSE298119", "GSM9008679", "PE003", "PE_CASE", "UNRESOLVED", "UNRESOLVED", "UNRESOLVED", "CESAREAN; elective only at paper cohort level", "PRJNA1268385", "SRX28967108", "GEO series calls this a term PE cohort; EOPE attribution is not sample-resolved"),
        ("GSE267340", "GSM8264272", "CONTROL", "CONTROL", "38+4", "NO_EVIDENCE", "UNRESOLVED", "CESAREAN; paper says elective", "PRJNA1111104", "SRX24539333", "External healthy control embedded in a term GDM/macrosomia series"),
        ("GSE298119", "GSM9008680", "CTL2", "CONTROL", "UNRESOLVED", "NO_EVIDENCE", "UNRESOLVED", "CESAREAN; elective only at paper cohort level", "PRJNA1268385", "SRX28967109", "GEO series calls this a term normotensive control"),
    ]
    zrows = []
    for gse, gsm, pid, disease, ga, fgr, sex, delivery, bp, srx, note in zheng_samples:
        zrows.append({
            "study": "Zheng_2025", "reported_cohort": "3 EOPE + 2 control", "gse": gse, "gsm": gsm, "patient_id": pid,
            "disease_label_geo": disease, "eoPE_definition_or_status": "CONFIRMED" if pid == "PE001" else ("NOT_APPLICABLE" if disease == "CONTROL" else "UNRESOLVED_CONTRADICTED_BY_TERM_SERIES"),
            "gestational_age_weeks": ga, "FGR": fgr, "fetal_sex": sex, "delivery_mode": delivery,
            "processed_raw_matrix": "PUBLIC_10X_BARCODES_FEATURES_MATRIX", "FASTQ": "PUBLIC_SRA", "raw_count_availability": "PUBLIC_CELL_RANGER_STYLE_COUNTS",
            "cell_annotations_with_barcode_map": "NOT_PUBLICLY_IDENTIFIED", "sequencing_batch_or_series": gse, "bioproject": bp, "sra_experiment": srx,
            "overlap": "SAME_REPORTED_FIVE_SUBJECT_COHORT_ACROSS_THREE_GSE; NOT_INDEPENDENT_GSE_EVIDENCE",
            "cohort_reusability": "FAIL_TARGETED_EOPE_GATE", "reason": note,
            "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}|{ZHENG}",
        })
    write_csv(OUT / "external_validation/zheng_eope_dataset_audit.csv", zrows)

    shared = read_csv(ROOT / "results/02_phase2a/programs/shared_pe_programs.csv")
    not_run = []
    for row in shared:
        not_run.append({
            "validation_study": "Zheng_2025", "celltype": row["celltype"], "collection": row["collection"], "gene_set": row["gene_set"], "pathway": row["pathway"],
            "legacy_direction": "UP_IN_PE_PROGRAM" if float(row["NES_EOPE"]) > 0 else "DOWN_IN_PE_PROGRAM",
            "mapping_status": "NOT_MAPPABLE", "PE_n": 3, "control_n": 2, "patient_level_effect": "", "direction": "", "exact_or_permutation_P": "", "BH_FDR": "",
            "validation_status": "NOT_RUN_EOPE_ESTIMAND_CONTRADICTION", "reason": "PE002/PE003 are described as term PE in GSE298119 but EOPE in the paper; their GA/FGR metadata and public cell-annotation barcode map are absent. Frozen validation gate failed.",
            "source_url": f"{ZHENG}|https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282038|https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267340|https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298119",
        })
    write_csv(OUT / "external_validation/zheng_eope_targeted_validation.csv", not_run)

    legacy_modules = {r["pathway"]: r for r in read_csv(ROOT / "results/02_phase2a1/shared_program_modules.csv") if r["record_type"] == "ORIGINAL_GENE_SET"}
    yang = read_csv(ROOT / "results/02_phase2a1/yang_lope_replication.csv")
    yang_updated = []
    for row in yang:
        module = legacy_modules[row["pathway"]]
        item = dict(row)
        item.update({
            "evidence_version": "PHASE2A1_VALUES_PRESERVED_WITH_PHASE2A2_LABEL", "legacy_hypothesis_status": "LEGACY_HYPOTHESIS_MODULE",
            "program_module": module["program_module"], "module_label": module["module_label"], "rerun_status": "NOT_RERUN_VALUES_PRESERVED",
        })
        yang_updated.append(item)
    write_csv(OUT / "external_validation/yang_lope_updated_evidence.csv", yang_updated)
    print(json.dumps({"admati_provenance_rows": len(audit), "raw_search_rows": len(search), "zheng_samples": len(zrows), "zheng_tests_not_run": len(not_run), "yang_rows_preserved": len(yang_updated)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
