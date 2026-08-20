#!/usr/bin/env python3
"""Build the Phase 0 audit tables from GEO/SRA metadata plus cited judgment config."""

from __future__ import annotations

import csv
import gzip
import io
import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "results" / "00_dataset_audit"
DOCS = ROOT / "docs"
CFG = json.loads((ROOT / "config" / "audit_accessions.json").read_text(encoding="utf-8"))
JUDGMENTS = json.loads((ROOT / "config" / "dataset_judgments.json").read_text(encoding="utf-8"))
UNRESOLVED = "UNRESOLVED"
GSE182158_DONORS = {
    "A01": ("38", "Female"), "A02": ("46", "Female"), "A03": ("32", "Female"),
    "B01": ("33", "Male"), "B02": ("43", "Male"), "B03": ("27", "Female"),
    "D01": ("33", "Male"), "D02": ("28", "Male"), "D03": ("22", "Male"),
    "U01": ("28", "Female"), "U02": ("37", "Female"),
}


def first(mapping: dict[str, list[str]], key: str, default: str = UNRESOLVED) -> str:
    values = mapping.get(key, [])
    return values[0] if values else default


def parse_soft(accession: str) -> dict:
    path = RAW / f"{accession}_family.soft.gz"
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}; run download_metadata.py first")
    series: dict[str, list[str]] = defaultdict(list)
    platforms: dict[str, dict[str, list[str]]] = {}
    samples: list[dict] = []
    current_kind = None
    current = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith("^SERIES = "):
                current_kind, current = "series", series
            elif line.startswith("^PLATFORM = "):
                pid = line.split(" = ", 1)[1]
                platforms[pid] = defaultdict(list)
                current_kind, current = "platform", platforms[pid]
            elif line.startswith("^SAMPLE = "):
                sample = {"gsm": line.split(" = ", 1)[1], "fields": defaultdict(list)}
                samples.append(sample)
                current_kind, current = "sample", sample["fields"]
            elif line.startswith("^"):
                current_kind, current = None, None
            elif line.startswith("!") and " = " in line and current is not None:
                key, value = line[1:].split(" = ", 1)
                if current_kind == "series" and key.startswith("Series_"):
                    current[key].append(value)
                elif current_kind == "platform" and key.startswith("Platform_"):
                    current[key].append(value)
                elif current_kind == "sample" and key.startswith("Sample_"):
                    current[key].append(value)
    return {"accession": accession, "series": dict(series), "platforms": platforms, "samples": samples}


def characteristics(sample: dict) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in sample["fields"].get("Sample_characteristics_ch1", []):
        if ":" in item:
            key, value = item.split(":", 1)
            result[key.strip().lower()] = value.strip()
    return result


def char_get(chars: dict[str, str], *keys: str) -> str:
    for key in keys:
        if key.lower() in chars and chars[key.lower()].strip():
            return chars[key.lower()].strip()
    return UNRESOLVED


def classify_disease(chars: dict[str, str], source: str, title: str) -> str:
    exact = char_get(chars, "diagnosis", "clinical group", "disease", "disease state", "subject status", "condition", "group")
    if exact != UNRESOLVED:
        return exact
    text = f"{source} {title}".lower()
    if "no pe" in text or "non-pe" in text or "control" in text or "healthy" in text:
        return "control/non-PE (derived from exact title)"
    if " pe" in f" {text}" or "preeclamp" in text or "eclamps" in text:
        return "PE (derived from exact title)"
    return UNRESOLVED


def classify_subtype(disease: str) -> str:
    value = disease.lower()
    if "early" in value or "eope" in value:
        return "early onset"
    if "late" in value or "lope" in value:
        return "late onset"
    if "severe" in value or "serious" in value:
        return "severe"
    if "pe" in value or "preeclamp" in value or "eclamp" in value:
        return "unspecified"
    return "NOT_APPLICABLE"


def raw_available(sample: dict) -> str:
    fields = sample["fields"]
    has_sra = any("sra" in value.lower() for value in fields.get("Sample_relation", []))
    has_supp = bool([k for k in fields if k.startswith("Sample_supplementary_file")])
    return "YES" if has_sra or has_supp or first(fields, "Sample_type", "").lower() == "rna" else "UNRESOLVED"


def processed_available(sample: dict) -> str:
    return "YES" if any(k.startswith("Sample_supplementary_file") for k in sample["fields"]) else "SERIES_LEVEL_OR_UNRESOLVED"


def twin_subject_id(accession: str, title: str, gsm: str) -> str:
    if accession != "GSE272342":
        return gsm
    parts = [part.strip() for part in title.split(",")]
    if len(parts) < 2:
        return gsm
    cohort, code = parts[0], parts[1]
    if cohort.startswith(("D ", "M ")):
        number = re.match(r"(\d+)", code)
        return f"GSE272342_TWIN_PREG_{number.group(1) if number else code}"
    return f"GSE272342_SINGLETON_{code}"


def make_sample_row(dataset_accession: str, sample: dict, origin_gse: str | None = None,
                    membership: str = "direct_gsm", reanalysis_of: str = "") -> dict[str, str]:
    fields = sample["fields"]
    chars = characteristics(sample)
    gsm = sample["gsm"]
    title = first(fields, "Sample_title")
    source = first(fields, "Sample_source_name_ch1")
    origin = origin_gse or dataset_accession
    tissue = char_get(chars, "tissue")
    if tissue == UNRESOLVED:
        tissue = source
    disease = classify_disease(chars, source, title)
    ga = char_get(chars, "gestational age (weeks)", "ga (week)", "gestational age_at_delivery_(weeks.days)")
    if "ga (week)" in chars:
        day = chars.get("ga (day)", "0")
        ga = f"{chars['ga (week)']}w+{day}d"
    delivery_ga = char_get(chars, "gestational age_at_delivery_(weeks.days)")
    maternal_age = char_get(chars, "maternal age", "age of_mother_(years)", "age")
    fgr = char_get(chars, "iugr diagnosis", "fgr", "sga")
    labor = char_get(chars, "labor status", "attempted vaginal delivery")
    delivery_mode = char_get(chars, "mode of delivery")
    fetal_sex = char_get(chars, "fetal_sex", "fetal sex", "infant gender", "assigned neonatal_sex", "gender", "sex")
    ancestry = char_get(chars, "maternal_ancestry_self_report", "maternal ethnicity", "maternal race", "ethnicity", "race")
    if "maternal ethnicity" in chars and "maternal race" in chars:
        ancestry = f"{chars['maternal ethnicity']}; {chars['maternal race']}"
    multiplicity = char_get(chars, "singleton or_twin")
    passage = char_get(chars, "passage", "passages")
    treatment = char_get(chars, "treatment")
    donor = char_get(chars, "donor id", "donor")
    donor_age = UNRESOLVED
    donor_sex = UNRESOLVED
    biological_unit = "cell" if dataset_accession == "GSE117837" else (
        "placenta_library" if "placenta" in tissue.lower() or "placenta" in source.lower() else "donor_library"
    )
    cell_id = title if biological_unit == "cell" else "NOT_APPLICABLE"
    treatment_duration = "12 h" if dataset_accession == "GSE117837" and treatment != "naïve" else "NOT_APPLICABLE"
    if dataset_accession == "GSE182158" and title in GSE182158_DONORS:
        donor = title
        donor_age, donor_sex = GSE182158_DONORS[title]
        passage = "P1 or P2; exact donor mapping UNRESOLVED"
    if dataset_accession == "GSE282038" and title == "PE001":
        ga, delivery_ga, maternal_age, fgr, labor, delivery_mode = "31+4 weeks", "31+4 weeks", "30", "FGR", "not in active labor", "cesarean"
    if origin == "GSE267340" and title == "CONTROL" and dataset_accession == "GSE282038":
        disease, ga, delivery_ga, maternal_age = "healthy external control", "38+4 weeks", "38+4 weeks", "36"
    source_fields = [tissue, disease, ga, maternal_age, donor_age, donor_sex, fgr, labor, delivery_mode, fetal_sex, ancestry, multiplicity, passage, treatment]
    known = sum(value not in {UNRESOLVED, "", "NOT_APPLICABLE"} for value in source_fields)
    reanalysis = reanalysis_of or (origin if origin != dataset_accession else "")
    return {
        "dataset_accession": dataset_accession,
        "sample_accession": gsm,
        "sample_title": title,
        "analytical_membership": membership,
        "origin_gse": origin,
        "biological_unit_type": biological_unit,
        "biological_subject_id": twin_subject_id(dataset_accession, title, gsm),
        "donor_id": donor,
        "cell_id": cell_id,
        "tissue_source": tissue,
        "exact_biological_material": source,
        "disease_group": disease,
        "pe_subtype": classify_subtype(disease),
        "gestational_age": ga,
        "delivery_gestational_age": delivery_ga,
        "maternal_age": maternal_age,
        "donor_age": donor_age,
        "donor_sex": donor_sex,
        "fgr_iugr_status": fgr,
        "labor_status": labor,
        "delivery_mode": delivery_mode,
        "fetal_sex": fetal_sex,
        "ancestry_ethnicity": ancestry,
        "singleton_twin_status": multiplicity,
        "placental_sampling_location": char_get(chars, "sampling location", "placental sampling location"),
        "passage": passage,
        "treatment": treatment,
        "treatment_duration": treatment_duration,
        "platform": first(fields, "Sample_platform_id"),
        "raw_available": raw_available(sample),
        "processed_available": processed_available(sample),
        "technical_replicate": "YES: re-sequenced GSE203507 placenta" if dataset_accession == "GSE272342" and title.split(",")[1].strip() in {"S06", "S18", "S20"} else "NO_EVIDENCE",
        "reanalysis_of": reanalysis or "NONE",
        "metadata_completeness": f"{known}/14 key fields present (NOT_APPLICABLE excluded)",
        "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}" + ("|https://pmc.ncbi.nlm.nih.gov/articles/PMC8715893/" if dataset_accession == "GSE182158" else ""),
        "source_accession": gsm + (";PMID34965030;Supplementary Table 1" if dataset_accession == "GSE182158" else ""),
        "source_record_type": "GEO_SAMPLE_SOFT;PUBLICATION_SUPPLEMENT_XLSX" if dataset_accession == "GSE182158" else "GEO_SAMPLE_SOFT",
        "source_evidence": "exact submitted Sample_characteristics; derived labels explicitly marked" + ("; donor age/sex from publication Supplementary Table 1" if dataset_accession == "GSE182158" else ""),
    }


def parse_pubmed() -> dict[str, dict[str, str]]:
    path = RAW / "pubmed_batch.xml"
    if not path.exists():
        return {}
    root = ET.parse(path).getroot()
    result = {}
    for article in root.findall("./PubmedArticle"):
        pmid = article.findtext("./MedlineCitation/PMID") or ""
        title_node = article.find("./MedlineCitation/Article/ArticleTitle")
        title = "".join(title_node.itertext()) if title_node is not None else ""
        year = article.findtext("./MedlineCitation/Article/Journal/JournalIssue/PubDate/Year")
        if not year:
            year = article.findtext("./MedlineCitation/Article/Journal/JournalIssue/PubDate/MedlineDate") or UNRESOLVED
        ids = {node.attrib.get("IdType", ""): node.text or "" for node in article.findall("./PubmedData/ArticleIdList/ArticleId")}
        result[pmid] = {"title": title, "year": year, "doi": ids.get("doi", UNRESOLVED), "pmc": ids.get("pmc", UNRESOLVED)}
    return result


def verify_primary_source_evidence() -> None:
    """Machine-check high-risk facts that originate in paper supplements/full text."""
    msc_supplement = RAW / "PMC8715893_supplementary.zip"
    if not msc_supplement.exists():
        raise FileNotFoundError(f"Missing {msc_supplement}; run download_metadata.py")
    with zipfile.ZipFile(msc_supplement) as outer:
        members = [name for name in outer.namelist() if name.endswith("CTM2-11-e650-s025.xlsx")]
        assert len(members) == 1, "GSE182158 donor-information XLSX not uniquely found"
        xlsx_bytes = outer.read(members[0])
    with zipfile.ZipFile(io.BytesIO(xlsx_bytes)) as workbook:
        spreadsheet_ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        shared_strings = []
        if "xl/sharedStrings.xml" in workbook.namelist():
            shared_root = ET.fromstring(workbook.read("xl/sharedStrings.xml"))
            for item in shared_root.findall("x:si", spreadsheet_ns):
                shared_strings.append("".join(node.text or "" for node in item.findall(".//x:t", spreadsheet_ns)))
        sheet_root = ET.fromstring(workbook.read("xl/worksheets/sheet1.xml"))
        donor_rows = []
        for row in sheet_root.findall(".//x:sheetData/x:row", spreadsheet_ns):
            values = []
            for cell in row.findall("x:c", spreadsheet_ns):
                value_node = cell.find("x:v", spreadsheet_ns)
                value = value_node.text if value_node is not None and value_node.text is not None else ""
                if cell.attrib.get("t") == "s" and value:
                    value = shared_strings[int(value)]
                elif cell.attrib.get("t") == "inlineStr":
                    value = "".join(node.text or "" for node in cell.findall(".//x:t", spreadsheet_ns))
                values.append(value.strip())
            if values:
                donor_rows.append(values)
    observed_donors = {row[0]: (row[2], row[3]) for row in donor_rows[1:] if len(row) >= 4}
    assert observed_donors == GSE182158_DONORS, observed_donors

    supplement = RAW / "PMC12092795_supplementary.zip"
    if not supplement.exists():
        raise FileNotFoundError(f"Missing {supplement}; run download_metadata.py")
    with zipfile.ZipFile(supplement) as outer:
        members = [name for name in outer.namelist() if name.endswith("MOESM6_ESM.docx")]
        assert len(members) == 1, "GSE282038 clinical supplement DOCX not uniquely found"
        docx_bytes = outer.read(members[0])
    with zipfile.ZipFile(io.BytesIO(docx_bytes)) as docx:
        document_xml = docx.read("word/document.xml")
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    root = ET.fromstring(document_xml)
    tables = []
    for table in root.findall(".//w:tbl", ns):
        rows = []
        for row in table.findall("./w:tr", ns):
            cells = []
            for cell in row.findall("./w:tc", ns):
                cells.append("".join(node.text or "" for node in cell.findall(".//w:t", ns)).strip())
            rows.append(cells)
        tables.append(rows)
    clinical = next(table for table in tables if table and table[0][:2] == ["Sample", "Group"] and table[0][2].replace(" ", "") == "Age,years")
    by_sample = {row[0]: row for row in clinical[1:]}
    assert by_sample["1"][1] == "EOPE" and by_sample["1"][4:6] == ["31 + 4", "27 + 4"]
    assert by_sample["9"][1] == "Healthy" and by_sample["9"][4] == "38 + 4"

    licensing = RAW / "PMC6506509_fulltext.xml"
    if not licensing.exists():
        raise FileNotFoundError(f"Missing {licensing}; run download_metadata.py")
    licensing_text = " ".join("".join(ET.parse(licensing).getroot().itertext()).split())
    assert "10 ng/mL" in licensing_text and "12 h" in licensing_text


def sra_stats(accession: str) -> tuple[int, float, str]:
    path = RAW / f"{accession}_sra_runinfo.csv"
    if not path.exists():
        return 0, 0.0, UNRESOLVED
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    size_mb = sum(float(row.get("size_MB") or 0) for row in rows)
    studies = sorted({row.get("SRAStudy", "") for row in rows if row.get("SRAStudy")})
    return len(rows), size_mb / 1024, ";".join(studies) or UNRESOLVED


HEAD_FALLBACK = {
    "GSE117837": 6155545, "GSE131355": 5468160, "GSE173193": 503316480,
    "GSE182158": 1554186240, "GSE190639": 174671, "GSE199071": 167833600,
    "GSE233634": 5597176, "GSE234729": 10209449, "GSE272342": 4874240,
    "GSE275980": 3111693, "GSE282038": 127528960, "GSE329173": 389713920,
    "GSE75010": 697139200 + 34204908,
}


def supplement_stats(accession: str, parsed: dict) -> tuple[list[str], int]:
    files = parsed["series"].get("Series_supplementary_file", [])
    size_path = RAW / "supplement_sizes.csv"
    if size_path.exists():
        with size_path.open(encoding="utf-8-sig", newline="") as handle:
            rows = [row for row in csv.DictReader(handle) if row["geo_accession"] == accession]
        sizes = [int(row["content_length_bytes"]) for row in rows if row["content_length_bytes"].isdigit()]
        if sizes:
            return files, sum(sizes)
    return files, HEAD_FALLBACK.get(accession, 0)


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_overlap(parsed_all: dict[str, dict]) -> list[dict[str, str]]:
    rows = []
    relation_ids = {
        value.rsplit(": ", 1)[-1]
        for value in parsed_all["GSE75010"]["series"].get("Series_relation", [])
        if value.startswith("Reanalysis of:")
    }
    for source in ["GSE30186", "GSE10588", "GSE24129", "GSE25906", "GSE43942", "GSE4707", "GSE44711"]:
        source_ids = {sample["gsm"] for sample in parsed_all[source]["samples"]}
        overlap = sorted(relation_ids & source_ids)
        rows.append({
            "dataset_a": "GSE75010", "dataset_b": source, "overlap_type": "exact GSM reanalysis",
            "overlap_n": str(len(overlap)), "overlap_accessions": ";".join(overlap), "evidence_status": "CONFIRMED",
            "data_leakage_rule": "same biological samples; never use as independent train/validation cohorts",
            "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75010",
            "source_accession": f"GSE75010;{source}", "source_evidence": "GSE75010 Series_relation intersected with source-series GSM",
        })
    rows.extend([
        {"dataset_a": "GSE234729", "dataset_b": "GSE186257", "overlap_type": "exact GSM reanalysis", "overlap_n": "12",
         "overlap_accessions": ";".join(value.rsplit(": ", 1)[-1] for value in parsed_all["GSE234729"]["series"].get("Series_relation", []) if value.startswith("Reanalysis of:")),
         "evidence_status": "CONFIRMED", "data_leakage_rule": "shared samples; count once", "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234729",
         "source_accession": "GSE234729;GSE186257", "source_evidence": "GEO Series_relation and GSE186257 GSM mapping"},
        {"dataset_a": "GSE272342", "dataset_b": "GSE203507", "overlap_type": "same placenta re-extracted/re-sequenced", "overlap_n": "3",
         "overlap_accessions": "S06=GSM6175121;S18=GSM6175133;S20=GSM6175135", "evidence_status": "CONFIRMED",
         "data_leakage_rule": "biological overlap; protocol-comparison libraries are not independent validation", "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272342",
         "source_accession": "GSE272342;GSE203507", "source_evidence": "GEO Overall design"},
        {"dataset_a": "GSE282038", "dataset_b": "GSE267340", "overlap_type": "external control dependency (not an exact duplicate)", "overlap_n": "1 control",
         "overlap_accessions": "GSM8264272 CONTROL", "evidence_status": "CONFIRMED", "data_leakage_rule": "cross-series disease/batch confounding; no internal control",
         "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12092795/", "source_accession": "GSE282038;GSE267340;PMID40394177", "source_evidence": "paper Data Availability and Supplementary Table S4"},
        {"dataset_a": "GSE282038", "dataset_b": "GSE298119", "overlap_type": "later cross-series integration (not an exact duplicate)", "overlap_n": "3 added donors",
         "overlap_accessions": "GSM9008678 PE002;GSM9008679 PE003;GSM9008680 CTL2", "evidence_status": "CONFIRMED", "data_leakage_rule": "report series provenance and batch; pooled total is 3 PE + 2 controls",
         "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12702750/", "source_accession": "GSE282038;GSE298119;PMID41403942", "source_evidence": "paper Data Availability and GEO SOFT"},
    ])
    return rows


RISK_ROWS = [
    ("GSE75010", "OVERLAP_7_SOURCE_GSE", "CRITICAL", "CONFIRMED", "330-sample matrix contains 173 exact GSM from seven source GSEs", "Keep one independence group; prohibit source-GSE train/validation reuse"),
    ("GSE234729", "DIRECT_VS_ANALYTIC_COUNT", "HIGH", "RESOLVED", "111 direct GSM plus 12 linked GSE186257 GSM explain the paper's 123 samples", "Download both series and deduplicate exact GSM"),
    ("GSE234729", "CLINICAL_METADATA_GAP", "HIGH", "UNRESOLVED", "GA, maternal age, delivery and FGR fields are absent for 111 direct GSM", "Obtain author/paper clinical supplement before adjusted analysis"),
    ("GSE190639", "EOPE_GA_CONFOUNDING", "CRITICAL", "CONFIRMED", "EOPE median delivery GA ~30.6 weeks vs healthy ~39.0", "Do not interpret naive EOPE-vs-term-control contrast as disease-only"),
    ("GSE272342", "NONINDEPENDENT_TWINS", "CRITICAL", "CONFIRMED", "32 placentas come from 16 twin pregnancies", "Use pregnancy-level blocking/random effects; patient n is not placenta n"),
    ("GSE272342", "MULTIPLICITY_CONFOUNDING", "HIGH", "CONFIRMED", "All 18 controls are twin placentas while case set includes singleton and twin samples", "Restrict/stratify by multiplicity and chorionicity"),
    ("GSE272342", "GSE203507_RESEQUENCING", "HIGH", "CONFIRMED", "S06/S18/S20 are re-extracted/resequenced placentas from GSE203507", "Do not use GSE203507 as independent validation for these subjects"),
    ("GSE272342", "FGR_STATUS_GAP", "MEDIUM", "UNRESOLVED", "Birthweight interval is present but explicit pregnancy-level FGR status is absent", "Obtain the clinical supplement before FGR-stratified use"),
    ("GSE173193", "SMALL_DONOR_N", "CRITICAL", "CONFIRMED", "PE n=2 and control n=2; other groups are GDM n=2 and AMA n=2", "Auxiliary replication only; cells are not donor replicates"),
    ("GSE173193", "GA_AND_CLINICAL_GAP", "HIGH", "UNRESOLVED", "Gestational age, delivery, FGR and uncomplicated-control status are absent", "Do not use as sole discovery cohort"),
    ("GSE282038", "NO_INTERNAL_CONTROL", "CRITICAL", "CONFIRMED", "GSE contains one EOPE case and zero controls", "External validation only"),
    ("GSE282038", "CROSS_SERIES_CONTROL", "CRITICAL", "CONFIRMED", "Primary control is GSE267340 CONTROL; later work also uses GSE298119", "Model/report series batch; do not claim within-series contrast"),
    ("GSE282038", "GA_FGR_CONFOUNDING", "CRITICAL", "CONFIRMED", "EOPE+FGR case 31+4 weeks vs external healthy control 38+4", "Disease, FGR, gestational age and series are inseparable"),
    ("GSE282038", "GSE298119_CLINICAL_GAP", "HIGH", "UNRESOLVED", "Exact gestational age and FGR status for PE002, PE003 and CTL2 are absent from submitted GEO sample metadata", "Obtain the linked study clinical table before cross-series integration"),
    ("GSE329173", "CASE_ONLY", "CRITICAL", "CONFIRMED", "Three severe PE/eclampsia cases and no controls", "Cell-state validation only; no PE-vs-normal design"),
    ("GSE329173", "PUBLICATION_CLINICAL_GAP", "HIGH", "UNRESOLVED", "No linked publication/PMID/DOI and no gestational-age or delivery metadata were found", "Obtain a citable publication or author-supplied clinical manifest before use"),
    ("GSE182158", "UC_DONOR_N2", "HIGH", "CONFIRMED", "U01/U02 are only two hUC-MSC donors", "Use donor-aware inference and external atlas support"),
    ("GSE182158", "PASSAGE_MAPPING_GAP", "MEDIUM", "UNRESOLVED", "Publication Supplementary Table 1 resolves donor age/sex, but protocol states P1/P2 without exact donor-specific passage mapping", "Retain passage as P1-or-P2; do not impute an exact passage"),
    ("GSE117837", "CELL_LABEL_COUNT_DISCREPANCY", "HIGH", "UNRESOLVED", "GEO labels give 203 naive/158 stimulated; paper says 202/159", "Resolve against original authors/source manifest; do not reassign a cell ad hoc"),
    ("GSE117837", "DONOR_PASSAGE_CONFOUNDING", "HIGH", "CONFIRMED", "Donor1 only P5; donor2 supplies P0/P2/P5", "Licensing comparison must be donor/passage aware"),
    ("GSE131355", "NO_REPLICATION", "HIGH", "UNRESOLVED", "One sample per tissue label and donor identities absent", "Supplementary descriptive use only"),
    ("GSE233634", "DONOR_INDEPENDENCE_UNKNOWN", "HIGH", "UNRESOLVED", "Triplicates lack GEO donor identifiers", "Do not assume biological triplicates"),
]


def build_report(dataset_rows: list[dict], overlap_rows: list[dict], risk_rows: list[dict], availability_rows: list[dict]) -> str:
    lines = [
        "# Phase 0 Data Feasibility Audit Report",
        "",
        f"Audit date: {CFG['audit_date']}",
        "Scope: metadata and data-availability feasibility only; no DEG, WGCNA, machine learning, CellChat, NicheNet, GSEA, or final biological inference was performed.",
        "",
        "## Executive decision",
        "",
        "**GO_WITH_MODIFICATIONS**",
        "",
        "Bulk PE discovery/replication and hUC-MSC atlas/licensing components are feasible, but the currently audited PE single-cell cohorts do not provide a robust, internally controlled primary discovery cohort. Phase 1 should proceed only after adding a better PE scRNA-seq cohort or formally downgrading PE single-cell work to external cell-state validation. GSE75010 must remain one composite dependence group, and all gestational-age/multiplicity constraints below are mandatory.",
        "",
        "## Methods and evidence standard",
        "",
        "GEO family SOFT was parsed programmatically for every GSE/GSM. SRA RunInfo supplied run counts and compressed download estimates. PubMed/Europe PMC supplied publication identifiers and paper methods. The GSE282038 clinical values were checked against Supplementary Tables S4/S5, and all 11 GSE182158 donor age/sex values were checked against Supplementary Table 1. Missing fields remain `UNRESOLVED`; no samples were removed or relabeled to improve results.",
        "",
        "## Dataset-by-dataset evaluation",
        "",
    ]
    for row in dataset_rows:
        status = "recommended conditional inclusion" if row["proposed_role"] not in {"EXCLUDE", "PENDING"} else "pending/exclude"
        lines.extend([
            f"### {row['geo_accession']} - {row['proposed_role']}", "",
            f"**Assessment:** {status}. {row['role_restrictions']}", "",
            f"**Confirmed facts:** {row['analytical_total_samples']}; cases `{row['case_count']}`, controls `{row['control_count']}`. Material: {row['exact_biological_material']}. Assay/platform: {row['assay_type']} / {row['platform']}.", "",
            f"**Clinical comparability:** subtype `{row['pe_subtype']}`; GA `{row['gestational_age']}`; control status `{row['controls_genuine_normotensive_uncomplicated']}`; GA comparability `{row['control_ga_comparable']}`.", "",
            f"**Unresolved/limitations:** {row['sample_metadata_completeness']}. {row['evidence_notes']}", "",
            f"**Sources:** {row['source_url']}", "",
        ])
    lines.extend(["## Duplicate and data-leakage risk", "", "Confirmed dependencies:", ""])
    for row in overlap_rows:
        lines.append(f"- {row['dataset_a']} ↔ {row['dataset_b']}: {row['overlap_type']}, n={row['overlap_n']}. {row['data_leakage_rule']}")
    lines.extend([
        "", "GSE75010 and its seven source GSEs are one dependence group. GSE234729 shares 12 exact samples with GSE186257. GSE272342 re-sequences three GSE203507 placentas. Cross-series controls for GSE282038 are dependencies even though they are not duplicate biospecimens.",
        "", "## Gestational-age matching risk", "",
        "- GSE190639: EOPE (~30.6 weeks median) versus healthy term (~39.0) is strongly mismatched.",
        "- GSE282038: EOPE+FGR case 31+4 weeks versus external healthy control 38+4; disease, FGR, GA and series batch are inseparable.",
        "- GSE75010: cohort-wide control GA is higher than PE and matching varies by source cohort; use study/GA-aware models only.",
        "- GSE234729: GA is not available in GEO sample metadata and must be obtained before adjusted modeling.",
        "- GSE272342: GA is available, but pregnancy multiplicity/chorionicity and non-independent twin placentas dominate design risk.",
        "", "## PE subtype heterogeneity", "",
        "Severe PE (GSE234729), EOPE/LOPE (GSE190639), mixed PE in singleton/twin pregnancies (GSE272342), EOPE+FGR (GSE282038), and case-only severe PE/eclampsia (GSE329173) must not be pooled as a homogeneous disease. Subtype labels must remain explicit in every Phase 1 design matrix.",
        "", "## Recommended inclusion and roles", "",
        "- `PRIMARY_BULK`: GSE234729, conditional on clinical-table acquisition and shared-sample handling.",
        "- `BULK_REPLICATION`: GSE272342, pregnancy-blocked and multiplicity-stratified.",
        "- `SUBTYPE_VALIDATION`: GSE75010 as one composite cohort; GSE190639 as a targeted immune-panel subtype check.",
        "- `SCRNA_REPLICATION`: GSE173193, GSE282038, GSE329173; none qualifies as the sole primary PE single-cell discovery cohort.",
        "- `PRIMARY_HUCMSC_ATLAS`: GSE182158 (U01 age 28 female; U02 age 37 female; exact donor-specific passage remains P1-or-P2 unresolved).",
        "- `HUCMSC_LICENSING`: GSE117837, with donor/passage awareness and the one-cell label discrepancy unresolved.",
        "- `SUPPLEMENTARY`: GSE131355, GSE199071, GSE233634, GSE275980.",
        "", "## Recommended exclusion from specific analyses", "",
        "No entire accession is deleted for convenience. Instead: exclude GSE329173 and GSE282038 from any internally controlled PE-vs-normal discovery; exclude GSE173193 as the sole discovery cohort; exclude GSE75010 constituent series as independent validation sets; exclude cell counts and twin co-placentas from patient-level n; exclude GSE131355 from inferential donor-level analyses.",
        "", "## Download scale estimate", "",
        "SRA estimates below are compressed archive estimates from RunInfo, not FASTQ expansion sizes:", "",
        "| Dataset | SRA runs | SRA GiB | GEO supplements |", "|---|---:|---:|---:|",
    ])
    for row in availability_rows:
        lines.append(f"| {row['geo_accession']} | {row['sra_run_count']} | {row['sra_compressed_size_gib']} | {row['geo_supplement_size_mib']} MiB |")
    total_gib = sum(float(row["sra_compressed_size_gib"]) for row in availability_rows if row["sra_compressed_size_gib"] not in {UNRESOLVED, ""})
    lines.extend([
        "", f"All audited SRA raw archives together are approximately **{total_gib:.1f} GiB compressed**. A sensible next phase should start from processed matrices and selectively retrieve raw reads only for datasets that pass design QC.",
        "", "## Unresolved issues", "",
    ])
    unresolved = [row for row in risk_rows if row["status"] == "UNRESOLVED"]
    for row in unresolved:
        lines.append(f"- {row['geo_accession']} / {row['flag_id']}: {row['issue']} Required: {row['required_mitigation']}")
    lines.extend([
        "", "## Phase 1 feasibility", "",
        "Phase 1 is feasible only with the modifications above. The bulk and hUC-MSC sides are adequate for carefully stratified, donor/pregnancy-aware work. The PE single-cell side is not adequate for a definitive discovery analysis because the audited cohorts are tiny, external-control-dependent, gestational-age/FGR-confounded, or case-only. The appropriate Phase 1 gate is therefore to locate/qualify an additional internally controlled PE scRNA-seq cohort or explicitly limit single-cell work to replication and cell-state localization.",
        "", "## Final conclusion", "", "**GO_WITH_MODIFICATIONS**", "",
    ])
    return "\n".join(lines)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    all_accessions = CFG["main_accessions"] + CFG["support_accessions"]
    parsed_all = {accession: parse_soft(accession) for accession in all_accessions}
    pubmed = parse_pubmed()
    verify_primary_source_evidence()

    dataset_rows = []
    for accession in CFG["main_accessions"]:
        parsed = parsed_all[accession]
        judgment = JUDGMENTS[accession]
        platforms = ";".join(first(fields, "Platform_title") for fields in parsed["platforms"].values()) or UNRESOLVED
        series = parsed["series"]
        dataset_rows.append({
            "geo_accession": accession, "category": judgment["category"], "title": first(series, "Series_title"),
            "pmid": judgment["pmid"], "doi": judgment["doi"], "publication": judgment["publication"], "publication_year": judgment["publication_year"],
            "organism": ";".join(sorted({first(sample["fields"], "Sample_organism_ch1") for sample in parsed["samples"]})),
            "tissue_source": judgment["tissue_source"], "exact_biological_material": judgment["exact_biological_material"],
            "assay_type": ";".join(series.get("Series_type", [UNRESOLVED])), "platform": platforms,
            "geo_gsm_count": str(len(parsed["samples"])), "analytical_total_samples": judgment["analytical_total_samples"],
            "actual_downloadable_samples": judgment["actual_downloadable_samples"], "biological_donor_or_pregnancy_n": judgment["biological_subject_n"],
            "case_count": judgment["case_count"], "control_count": judgment["control_count"], "other_count": judgment["other_count"],
            "pe_subtype": judgment["pe_subtype"], "gestational_age": judgment["gestational_age"], "delivery_gestational_age": judgment["delivery_gestational_age"],
            "maternal_age": judgment["maternal_age"], "fgr_iugr_status": judgment["fgr_iugr_status"], "labor_status": judgment["labor_status"],
            "donor_age": "A01 38;A02 46;A03 32;B01 33;B02 43;B03 27;D01 33;D02 28;D03 22;U01 28;U02 37" if accession == "GSE182158" else "NOT_APPLICABLE_OR_UNRESOLVED",
            "donor_sex": "A01 F;A02 F;A03 F;B01 M;B02 M;B03 F;D01 M;D02 M;D03 M;U01 F;U02 F" if accession == "GSE182158" else "NOT_APPLICABLE_OR_UNRESOLVED",
            "delivery_mode": judgment["delivery_mode"], "fetal_sex": judgment["fetal_sex"], "ancestry_ethnicity": judgment["ancestry_ethnicity"],
            "singleton_twin_status": judgment["singleton_twin_status"], "placental_sampling_location": judgment["placental_sampling_location"],
            "raw_data_availability": "YES" if accession in CFG["sra_accessions"] or accession in {"GSE75010", "GSE190639"} else "UNRESOLVED",
            "processed_matrix_availability": "YES" if series.get("Series_supplementary_file") else "UNRESOLVED",
            "sample_metadata_completeness": judgment["sample_metadata_completeness"],
            "controls_genuine_normotensive_uncomplicated": judgment["controls_genuine"], "control_ga_comparable": judgment["control_ga_comparable"],
            "primary_or_reanalysis": judgment["primary_or_reanalysis"], "overlaps_other_geo": judgment["overlaps_other_geo"],
            "proposed_role": judgment["proposed_role"], "role_restrictions": judgment["role_restrictions"],
            "source_url": judgment["source_urls"], "source_accession": f"{accession};{judgment['pmid']}", "source_pmid": judgment["pmid"],
            "source_doi": judgment["doi"], "evidence_notes": judgment["evidence_notes"], "audit_status": "AUDITED",
        })

    sample_rows = []
    for accession in CFG["main_accessions"]:
        sample_rows.extend(make_sample_row(accession, sample) for sample in parsed_all[accession]["samples"])
    # Explicit analytical dependencies not represented by direct GSM lists.
    source_map = {sample["gsm"]: (source, sample) for source in CFG["support_accessions"] for sample in parsed_all[source]["samples"]}
    for dataset in ["GSE75010", "GSE234729"]:
        for relation in parsed_all[dataset]["series"].get("Series_relation", []):
            if relation.startswith("Reanalysis of:"):
                gsm = relation.rsplit(": ", 1)[-1]
                if gsm in source_map:
                    origin, sample = source_map[gsm]
                    sample_rows.append(make_sample_row(dataset, sample, origin, "reused_external_gsm", origin))
    for source in ["GSE267340", "GSE298119"]:
        for sample in parsed_all[source]["samples"]:
            title = first(sample["fields"], "Sample_title")
            if source == "GSE267340" and title != "CONTROL":
                continue
            sample_rows.append(make_sample_row("GSE282038", sample, source, "external_series_reference", source))

    # Add cell-column evidence for U01/U02 without conflating it with donor n.
    barcode_counts = {}
    for gsm, title in [("GSM5519462", "U01"), ("GSM5519463", "U02")]:
        path = RAW / f"{gsm}_{title}_barcodes.tsv.gz"
        if path.exists():
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                barcode_counts[gsm] = sum(1 for _ in handle)
    for row in sample_rows:
        if row["dataset_accession"] == "GSE182158" and row["sample_accession"] in barcode_counts:
            row["source_evidence"] += f"; processed matrix barcode columns={barcode_counts[row['sample_accession']]} (not donor n)"

    overlap_rows = build_overlap(parsed_all)
    risk_rows = [{
        "geo_accession": acc, "flag_id": flag, "severity": severity, "status": status, "issue": issue,
        "evidence": JUDGMENTS[acc]["evidence_notes"], "required_mitigation": mitigation,
        "source_url": JUDGMENTS[acc]["source_urls"], "source_accession": f"{acc};{JUDGMENTS[acc]['pmid']}",
    } for acc, flag, severity, status, issue, mitigation in RISK_ROWS]

    availability_rows = []
    for accession in CFG["main_accessions"]:
        parsed = parsed_all[accession]
        runs, sra_gib, sra_study = sra_stats(accession)
        files, supplement_bytes = supplement_stats(accession, parsed)
        raw_repo = "SRA" if runs else ("GEO" if accession in {"GSE75010", "GSE190639"} else UNRESOLVED)
        availability_rows.append({
            "geo_accession": accession, "raw_available": "YES" if raw_repo != UNRESOLVED else UNRESOLVED,
            "raw_repository": raw_repo, "raw_accession": sra_study, "sra_run_count": str(runs),
            "sra_compressed_size_gib": f"{sra_gib:.2f}" if runs else "0.00",
            "processed_matrix_available": "YES" if files else UNRESOLVED,
            "processed_or_supplement_files": ";".join(Path(url).name for url in files) if files else UNRESOLVED,
            "geo_supplement_size_mib": f"{supplement_bytes / 1024 / 1024:.2f}" if supplement_bytes else UNRESOLVED,
            "download_recommendation": "START_WITH_PROCESSED" if files else "METADATA_REVIEW_BEFORE_RAW",
            "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession};https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc={accession}",
            "source_accession": f"{accession};{sra_study}", "size_basis": "SRA RunInfo size_MB sum; GEO HTTP Content-Length",
        })

    role_rows = [{
        "geo_accession": row["geo_accession"], "proposed_role": row["proposed_role"],
        "phase1_inclusion": "CONDITIONAL" if row["proposed_role"] not in {"EXCLUDE", "PENDING"} else "NO",
        "rationale": row["evidence_notes"], "mandatory_restrictions": row["role_restrictions"],
        "independence_group": (
            "GSE75010_COMPOSITE_330" if row["geo_accession"] == "GSE75010" else
            "GSE234729_GSE186257_SHARED" if row["geo_accession"] == "GSE234729" else
            "GSE282038_GSE267340_GSE298119_CROSS_SERIES" if row["geo_accession"] == "GSE282038" else row["geo_accession"]
        ),
        "source_url": row["source_url"], "source_accession": row["source_accession"],
    } for row in dataset_rows]

    write_csv(OUT / "dataset_registry.csv", dataset_rows)
    write_csv(OUT / "sample_registry.csv", sample_rows)
    write_csv(OUT / "dataset_overlap_matrix.csv", overlap_rows)
    write_csv(OUT / "dataset_risk_flags.csv", risk_rows)
    write_csv(OUT / "data_availability.csv", availability_rows)
    write_csv(OUT / "proposed_dataset_roles.csv", role_rows)
    report = build_report(dataset_rows, overlap_rows, risk_rows, availability_rows)
    (DOCS / "DATASET_AUDIT_REPORT.md").write_text(report, encoding="utf-8")
    print(f"Wrote {len(dataset_rows)} datasets, {len(sample_rows)} sample/dependency rows, {len(risk_rows)} risk flags")
    return 0


if __name__ == "__main__":
    sys.exit(main())
