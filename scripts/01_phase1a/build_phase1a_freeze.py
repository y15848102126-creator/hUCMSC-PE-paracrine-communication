# HISTORICAL PROVENANCE: the GSE30186 cohort-minimum shift-log branch below was withdrawn before formal outcome analysis. Use the Phase 1A.1 frozen normexp/quantile/log2 matrix.
#!/usr/bin/env python3
"""Build the Phase 1A bulk sample/data freeze without performing differential analysis."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
P1RAW = RAW / "phase1a"
INTERIM = ROOT / "data" / "interim" / "phase1a"
MATRIX_DIR = INTERIM / "processed_matrices"
OUT = ROOT / "results" / "01_phase1a"
DESIGN = json.loads((ROOT / "config" / "phase1a_design.json").read_text(encoding="utf-8"))
UNRESOLVED = "UNRESOLVED"
NA_VALUES = {"", "na", "n/a", "null", "not available", "unknown", "unresolved"}

ACCESSIONS = [
    "GSE234729", "GSE186257", "GSE75010", "GSE30186", "GSE10588",
    "GSE24129", "GSE25906", "GSE43942", "GSE4707", "GSE44711",
    "GSE190639", "GSE272342",
]
EXPRESSION_COHORTS = [
    "GSE234729_ANALYTIC_123", "GSE75010_BIOBANK", "GSE30186", "GSE10588",
    "GSE24129", "GSE25906", "GSE43942", "GSE4707", "GSE44711", "GSE190639",
]


def clean(value: object) -> str:
    text = str(value).strip()
    return UNRESOLVED if text.lower() in NA_VALUES else text


def first(mapping: dict[str, list[str]], key: str, default: str = UNRESOLVED) -> str:
    values = mapping.get(key, [])
    return clean(values[0]) if values else default


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_soft(accession: str) -> dict:
    path = RAW / f"{accession}_family.soft.gz"
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}; run download_phase1a_data.py")
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
            result[key.strip().lower()] = clean(value)
    return result


def cget(chars: dict[str, str], *keys: str) -> str:
    for key in keys:
        value = chars.get(key.lower(), UNRESOLVED)
        if value != UNRESOLVED:
            return value
    return UNRESOLVED


def disease_label(accession: str, sample: dict) -> tuple[str, str]:
    fields, chars = sample["fields"], characteristics(sample)
    title = first(fields, "Sample_title")
    source = first(fields, "Sample_source_name_ch1")
    exact = cget(chars, "disease", "diagnosis", "classification", "disease state", "phenotype", "condition", "group", "clinical group")
    text = f"{exact} {title} {source}".lower()
    if accession == "GSE24129" and "fetal growth restriction" in text:
        return "OTHER", "ISOLATED_FGR"
    if any(token in text for token in ["normal labor", "normal rep", "placenta_nl", "normotensive", "control", "healthy", "non-pe", "no pe", "normal pregnancy", "phenotype normal"]):
        return "CONTROL", "NOT_APPLICABLE"
    if any(token in text for token in ["preeclamp", "pre-eclamp", "eopet", "eope", "lope", "eo pe", "pe(", " pe", "preeclamptic"]):
        if any(token in text for token in ["eopet", "eope", "early-onset", "early onset", "pe(eo)"]):
            subtype = "EOPE"
        elif any(token in text for token in ["lope", "late-onset", "late onset", "pe(lo)"]):
            subtype = "LOPE"
        elif "severe" in text:
            subtype = "SEVERE_PE"
        else:
            subtype = "UNSPECIFIED_PE"
        return "PE", subtype
    return UNRESOLVED, UNRESOLVED


def normalize_ga(chars: dict[str, str]) -> str:
    if chars.get("ga (week)") not in {None, UNRESOLVED}:
        day = chars.get("ga (day)", "0")
        return f"{chars['ga (week)']}w+{day}d"
    return cget(chars, "gestational age", "gestational age (weeks)", "gestational age_at_delivery_(weeks.days)")


def pregnancy_id(accession: str, sample: dict) -> str:
    gsm = sample["gsm"]
    title = first(sample["fields"], "Sample_title")
    if accession != "GSE272342":
        return gsm
    parts = [part.strip() for part in title.split(",")]
    if len(parts) < 2:
        return gsm
    if parts[0].startswith("S "):
        return f"GSE272342_SINGLETON_{parts[1]}"
    match = re.match(r"(\d+)", parts[1])
    return f"GSE272342_TWIN_{match.group(1) if match else parts[1]}"


def sample_row(cohort: str, accession: str, sample: dict, provenance: str = "direct") -> dict:
    fields, chars = sample["fields"], characteristics(sample)
    group, subtype = disease_label(accession, sample)
    title = first(fields, "Sample_title")
    ga = normalize_ga(chars)
    delivery_ga = cget(chars, "gestational age_at_delivery_(weeks.days)")
    if delivery_ga == UNRESOLVED and accession in {"GSE75010", "GSE190639", "GSE272342"}:
        delivery_ga = ga
    ancestry = cget(chars, "maternal_ancestry_self_report", "maternal ethnicity", "ethnic group", "ethnicity")
    if "maternal race" in chars:
        ancestry = f"{ancestry}; {chars['maternal race']}" if ancestry != UNRESOLVED else chars["maternal race"]
    row = {
        "dataset": cohort,
        "GSM/sample ID": sample["gsm"],
        "sample_title": title,
        "patient/pregnancy ID if known": pregnancy_id(accession, sample),
        "PE/control": group,
        "PE subtype": subtype,
        "GA": ga,
        "delivery GA": delivery_ga,
        "maternal_age": cget(chars, "maternal age", "age", "age of_mother_(years)"),
        "BMI": cget(chars, "maternal bmi", "bmi"),
        "FGR": cget(chars, "iugr diagnosis", "sga", "fgr"),
        "singleton/twin": cget(chars, "singleton or_twin"),
        "fetal sex": cget(chars, "infant gender", "fetal sex", "fetal_sex", "gender", "assigned neonatal_sex"),
        "ancestry": ancestry,
        "platform": first(fields, "Sample_platform_id"),
        "batch": cget(chars, "batch"),
        "labor": cget(chars, "induction of labor", "attempted vaginal delivery", "labor status"),
        "delivery_mode": cget(chars, "mode of delivery"),
        "parity": cget(chars, "parity", "previous nulliparity"),
        "include_phase1b": "NO",
        "exclusion_reason": "UNRESOLVED_PENDING_RULE",
        "planned_phase1b_use": DESIGN["cohorts"][cohort]["role"],
        "independence_group": pregnancy_id(accession, sample),
        "source_provenance": provenance,
        "source": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={sample['gsm']}",
        "source_accession": f"{accession};{sample['gsm']}",
    }
    return row


def apply_inclusion(row: dict) -> None:
    cohort = row["dataset"]
    is_main = DESIGN["cohorts"][cohort]["phase1b_main"]
    group = row["PE/control"]
    title = row["sample_title"]
    if not is_main:
        row["include_phase1b"] = "NO"
        row["exclusion_reason"] = DESIGN["cohorts"][cohort]["reason"]
        return
    if cohort == "GSE75010_BIOBANK" and "-CH" in title:
        row["include_phase1b"] = "NO"
        row["exclusion_reason"] = "Chronic-hypertension stratum excluded from the primary PE-vs-normotensive contrast; retained for design-specific sensitivity."
    elif cohort == "GSE24129" and group == "OTHER":
        row["include_phase1b"] = "NO"
        row["exclusion_reason"] = "Isolated FGR without PE is not an eligible PE or uncomplicated normotensive-control sample."
    elif group in {"PE", "CONTROL"}:
        row["include_phase1b"] = "YES"
        row["exclusion_reason"] = "NONE"
    else:
        row["include_phase1b"] = "NO"
        row["exclusion_reason"] = "Disease label is not an eligible PE/control label."


def build_sample_freeze(parsed: dict[str, dict]) -> list[dict]:
    rows: list[dict] = []
    for sample in parsed["GSE234729"]["samples"]:
        rows.append(sample_row("GSE234729_ANALYTIC_123", "GSE234729", sample, "GSE234729 direct GSM"))
    rel_gsms = [item.split(": ", 1)[1] for item in parsed["GSE234729"]["series"].get("Series_relation", []) if item.startswith("Reanalysis of: GSM")]
    linked = {sample["gsm"]: sample for sample in parsed["GSE186257"]["samples"]}
    for gsm in rel_gsms:
        row = sample_row("GSE234729_ANALYTIC_123", "GSE186257", linked[gsm], "GSE186257 reanalysis included once in GSE234729 analytical matrix")
        row["dataset"] = "GSE234729_ANALYTIC_123"
        rows.append(row)
    accession_to_cohort = {
        "GSE75010": "GSE75010_BIOBANK", "GSE30186": "GSE30186", "GSE10588": "GSE10588",
        "GSE24129": "GSE24129", "GSE25906": "GSE25906", "GSE43942": "GSE43942",
        "GSE4707": "GSE4707", "GSE44711": "GSE44711", "GSE190639": "GSE190639",
        "GSE272342": "GSE272342",
    }
    for accession, cohort in accession_to_cohort.items():
        for sample in parsed[accession]["samples"]:
            rows.append(sample_row(cohort, accession, sample))
    for row in rows:
        apply_inclusion(row)
    rows.sort(key=lambda r: (r["dataset"], r["GSM/sample ID"]))
    return rows


def read_soft_expression(accession: str) -> pd.DataFrame:
    path = RAW / f"{accession}_family.soft.gz"
    sample_ids: list[str] = []
    probe_ids: list[str] | None = None
    values: list[np.ndarray] = []
    gsm = ""
    in_table = False
    header_seen = False
    current_ids: list[str] = []
    current_values: list[float] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                gsm = line.split(" = ", 1)[1]
            elif line == "!sample_table_begin":
                in_table, header_seen = True, False
                current_ids, current_values = [], []
            elif line == "!sample_table_end" and in_table:
                in_table = False
                if current_ids:
                    if probe_ids is None:
                        probe_ids = current_ids
                    elif probe_ids != current_ids:
                        raise AssertionError(f"Probe order differs across {accession} samples")
                    sample_ids.append(gsm)
                    values.append(np.asarray(current_values, dtype=np.float64))
            elif in_table:
                parts = line.split("\t")
                if not header_seen:
                    header_seen = True
                    if len(parts) < 2 or parts[0] != "ID_REF" or parts[1] != "VALUE":
                        raise AssertionError(f"Unexpected sample table header in {accession}: {parts[:2]}")
                elif len(parts) >= 2:
                    current_ids.append(parts[0])
                    try:
                        current_values.append(float(parts[1]))
                    except ValueError:
                        current_values.append(np.nan)
    if probe_ids is None or not values:
        raise ValueError(f"No embedded expression table in {accession}")
    return pd.DataFrame(np.column_stack(values), index=pd.Index(probe_ids, name="feature_id"), columns=sample_ids)


def read_platform_annotation(accession: str) -> pd.DataFrame:
    path = RAW / f"{accession}_family.soft.gz"
    in_table = False
    header: list[str] = []
    rows: list[list[str]] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line == "!platform_table_begin":
                in_table = True
                header = []
            elif line == "!platform_table_end" and in_table:
                break
            elif in_table:
                parts = line.split("\t")
                if not header:
                    header = parts
                else:
                    if len(parts) < len(header):
                        parts += [""] * (len(header) - len(parts))
                    rows.append(parts[:len(header)])
    if not header:
        raise ValueError(f"No platform annotation in {accession}")
    return pd.DataFrame(rows, columns=header, dtype=str)


def hgnc_lookup() -> tuple[dict[str, tuple[str, str]], dict[str, tuple[str, str]], dict[str, tuple[str, str]]]:
    hgnc = pd.read_csv(P1RAW / "hgnc_complete_set.txt", sep="\t", dtype=str, keep_default_na=False)
    symbol_map: dict[str, tuple[str, str]] = {}
    accession_map: dict[str, tuple[str, str]] = {}
    ensembl_map: dict[str, tuple[str, str]] = {}
    for _, row in hgnc.iterrows():
        symbol, ens = row["symbol"], row["ensembl_gene_id"]
        target = (symbol, ens or UNRESOLVED)
        for value in [symbol, *row["alias_symbol"].split("|"), *row["prev_symbol"].split("|")]:
            value = value.strip()
            if value:
                symbol_map.setdefault(value.upper(), target)
        for value in [*row["ena"].split("|"), *row["refseq_accession"].split("|")]:
            value = value.strip().split(".")[0]
            if value:
                accession_map.setdefault(value.upper(), target)
        if ens:
            ensembl_map[ens.split(".")[0].upper()] = target
    return symbol_map, accession_map, ensembl_map


def symbol_candidates(platform: str, row: pd.Series) -> tuple[str, str, str]:
    if platform == "GPL6244":
        assignment = row.get("gene_assignment", "")
        symbols, ids = [], []
        for segment in assignment.split(" /// "):
            parts = segment.split(" // ")
            if len(parts) >= 2 and parts[1] not in {"", "---"}:
                symbols.append(parts[1])
                ids.append(parts[0])
        return "|".join(dict.fromkeys(ids)), "|".join(dict.fromkeys(symbols)), row.get("GB_LIST", "")
    if platform == "GPL10558":
        return row.get("Entrez_Gene_ID", ""), row.get("Symbol", ""), row.get("RefSeq_ID", "") or row.get("GB_ACC", "")
    if platform == "GPL2986":
        return row.get("GeneID", ""), row.get("Gene Symbol", ""), row.get("GenBank", "")
    if platform == "GPL6102":
        return row.get("Entrez_Gene_ID", ""), row.get("Symbol", ""), row.get("RefSeq_ID", "") or row.get("GB_ACC", "")
    if platform == "GPL10191":
        return row.get("GB_ACC", ""), "", row.get("GB_ACC", "")
    if platform == "GPL1708":
        return row.get("ENSEMBL_ID", "") or row.get("GENE", ""), row.get("GENE_SYMBOL", ""), row.get("REFSEQ", "") or row.get("GB_ACC", "")
    if platform == "GPL31059":
        return row.get("ORF", ""), row.get("ORF", ""), row.get("GB_ACC", "")
    return "", "", ""


def build_platform_mapping(platform: str, annotation: pd.DataFrame, symbol_map: dict, accession_map: dict, ensembl_map: dict, datasets: str, source_acc: str) -> tuple[list[dict], dict[str, str]]:
    result: list[dict] = []
    mapping: dict[str, str] = {}
    hgnc_hash = sha256(P1RAW / "hgnc_complete_set.txt")
    for _, ann in annotation.iterrows():
        probe = clean(ann.get("ID", ""))
        original_id, original_symbol, accession = symbol_candidates(platform, ann)
        symbols = [value.strip() for value in original_symbol.split("|") if value.strip() and value.strip() != "---"]
        resolved = {symbol_map[value.upper()] for value in symbols if value.upper() in symbol_map}
        if not resolved:
            for token in re.split(r"[|,; ]+", accession):
                key = token.strip().split(".")[0].upper()
                if key in accession_map:
                    resolved.add(accession_map[key])
        if not resolved:
            for token in re.split(r"[|,; ]+", original_id):
                key = token.strip().split(".")[0].upper()
                if key in ensembl_map:
                    resolved.add(ensembl_map[key])
        if len(resolved) == 1:
            mapped_symbol, mapped_ens = next(iter(resolved))
            status = "MAPPED_UNAMBIGUOUS"
            mapping[probe] = mapped_symbol
        elif len(resolved) > 1:
            mapped_symbol, mapped_ens = UNRESOLVED, UNRESOLVED
            status = "AMBIGUOUS_MULTI_GENE_EXCLUDED_FROM_COLLAPSE"
        else:
            mapped_symbol, mapped_ens = UNRESOLVED, UNRESOLVED
            status = "UNMAPPED"
        result.append({
            "dataset": datasets, "platform": platform, "original_probe_id": probe,
            "original_gene_id": clean(original_id), "original_symbol": clean(original_symbol),
            "mapped_symbol": mapped_symbol, "mapped_ensembl_gene_id": mapped_ens,
            "mapping_source": f"{platform}+HGNC_COMPLETE_SET",
            "mapping_version": f"HGNC_2026-08-09_{hgnc_hash[:16]}",
            "mapping_status": status,
            "collapse_rule": "MEDIAN_AFTER_TRANSFORM_NO_OUTCOME",
            "source": f"{source_acc};{platform};HGNC_COMPLETE_SET_2026-08-09",
        })
    return result, mapping


def build_rnaseq_mapping(matrix: pd.DataFrame, symbol_map: dict, ensembl_map: dict) -> tuple[list[dict], dict[str, str]]:
    rows, mapping = [], {}
    hgnc_hash = sha256(P1RAW / "hgnc_complete_set.txt")
    for feature in matrix.index.astype(str):
        parts = feature.split("$", 1)
        original_id = parts[0].split(".")[0]
        original_symbol = parts[1] if len(parts) == 2 else ""
        resolved = symbol_map.get(original_symbol.upper()) or ensembl_map.get(original_id.upper())
        if resolved:
            mapped_symbol, mapped_ens = resolved
            status = "MAPPED_UNAMBIGUOUS"
            mapping[feature] = mapped_symbol
        else:
            mapped_symbol, mapped_ens, status = UNRESOLVED, UNRESOLVED, "UNMAPPED"
        rows.append({
            "dataset": "GSE234729_ANALYTIC_123", "platform": "GPL24676",
            "original_probe_id": feature, "original_gene_id": original_id,
            "original_symbol": clean(original_symbol), "mapped_symbol": mapped_symbol,
            "mapped_ensembl_gene_id": mapped_ens,
            "mapping_source": "GSE234729_SUBMITTED_ID+HGNC_COMPLETE_SET",
            "mapping_version": f"GRCh38.89;HGNC_2026-08-09_{hgnc_hash[:16]}",
            "mapping_status": status,
            "collapse_rule": "GENE_INPUT_MEDIAN_IF_DUPLICATE_NO_OUTCOME",
            "source": "GSE234729;GPL24676;HGNC_COMPLETE_SET_2026-08-09",
        })
    return rows, mapping


def transform_matrix(matrix: pd.DataFrame, rule: str) -> pd.DataFrame:
    values = matrix.astype(float)
    if rule == "none":
        return values
    if rule.startswith("log2(value - cohort_min"):
        minimum = np.nanmin(values.to_numpy())
        return np.log2(values - minimum + 1.0)
    if rule.startswith("log2(value + 1") or rule.startswith("log2(normalized_count + 1"):
        return np.log2(values + 1.0)
    raise ValueError(f"Unsupported transform: {rule}")


def collapse_matrix(matrix: pd.DataFrame, feature_mapping: dict[str, str]) -> pd.DataFrame:
    mapped = pd.Series([feature_mapping.get(str(feature), "") for feature in matrix.index], index=matrix.index)
    keep = mapped != ""
    filtered = matrix.loc[keep].copy()
    filtered.insert(0, "mapped_symbol", mapped.loc[keep].to_numpy())
    return filtered.groupby("mapped_symbol", sort=True).median(numeric_only=True)


def write_matrix(path: Path, matrix: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(path, sep="\t", compression="gzip", float_format="%.8g", lineterminator="\n")


def qc_rows(cohort: str, matrix: pd.DataFrame) -> list[dict]:
    numeric = matrix.replace([np.inf, -np.inf], np.nan).dropna(axis=0, how="any")
    variances = numeric.var(axis=1)
    numeric = numeric.loc[variances[variances > 0].nlargest(min(2000, int((variances > 0).sum()))).index]
    if numeric.shape[0] < 2 or numeric.shape[1] < 3:
        return []
    x = numeric.to_numpy(dtype=float, copy=True).T
    x -= x.mean(axis=0, keepdims=True)
    scale = x.std(axis=0, ddof=1)
    scale[scale == 0] = 1
    x /= scale
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    scores = u[:, :2] * s[:2]
    corr = np.corrcoef(x)
    median_corr = np.array([np.median(np.delete(corr[i], i)) for i in range(corr.shape[0])])
    center = np.median(median_corr)
    mad = np.median(np.abs(median_corr - center))
    threshold = center - 3 * (1.4826 * mad if mad > 0 else 0)
    rows = []
    for i, sample in enumerate(numeric.columns):
        flag = "FLAG_REVIEW" if median_corr[i] < threshold else "PASS"
        rows.append({
            "dataset": cohort, "sample_id": sample, "pc1": f"{scores[i, 0]:.8g}",
            "pc2": f"{scores[i, 1]:.8g}", "median_sample_correlation": f"{median_corr[i]:.8g}",
            "outlier_threshold": f"{threshold:.8g}", "qc_flag": flag,
            "action": "FLAG_ONLY_NO_AUTOMATIC_EXCLUSION" if flag != "PASS" else "NONE",
            "method": "Outcome-blind PCA and robust median-correlation rule on up to 2000 highest-variance mapped genes; no disease labels used",
            "source": f"data/interim/phase1a/processed_matrices/{cohort}_gene_level.tsv.gz",
        })
    return rows


def expression_and_mapping(parsed: dict[str, dict], freeze_rows: list[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    symbol_map, accession_map, ensembl_map = hgnc_lookup()
    platform_specs = {
        "GPL6244": ("GSE75010", "GSE75010_BIOBANK|GSE24129"),
        "GPL10558": ("GSE30186", "GSE30186|GSE44711"),
        "GPL2986": ("GSE10588", "GSE10588"),
        "GPL6102": ("GSE25906", "GSE25906"),
        "GPL10191": ("GSE43942", "GSE43942"),
        "GPL1708": ("GSE4707", "GSE4707"),
        "GPL31059": ("GSE190639", "GSE190639"),
    }
    mapping_rows: list[dict] = []
    platform_maps: dict[str, dict[str, str]] = {}
    for platform, (source_acc, datasets) in platform_specs.items():
        ann = read_platform_annotation(source_acc)
        rows, mapping = build_platform_mapping(platform, ann, symbol_map, accession_map, ensembl_map, datasets, source_acc)
        mapping_rows.extend(rows)
        platform_maps[platform] = mapping

    processing_rows: list[dict] = []
    all_qc: list[dict] = []
    freeze_by_cohort = defaultdict(list)
    for row in freeze_rows:
        freeze_by_cohort[row["dataset"]].append(row)

    for cohort in EXPRESSION_COHORTS:
        spec = DESIGN["cohorts"][cohort]
        accession = spec["accession"]
        if cohort == "GSE234729_ANALYTIC_123":
            matrix = pd.read_csv(P1RAW / "GSE234729_deseq2_normalized_filtered.txt.gz", sep="\t", index_col=0)
            title_to_gsm = {}
            for row in freeze_by_cohort[cohort]:
                title_to_gsm[row["sample_title"]] = row["GSM/sample ID"]
            missing = sorted(set(matrix.columns) - set(title_to_gsm))
            if missing:
                raise AssertionError(f"GSE234729 columns not mapped to GSM: {missing}")
            matrix.columns = [title_to_gsm[column] for column in matrix.columns]
            rnaseq_rows, feature_map = build_rnaseq_mapping(matrix, symbol_map, ensembl_map)
            mapping_rows.extend(rnaseq_rows)
        else:
            matrix = read_soft_expression(accession)
            feature_map = platform_maps[spec["platform"]]
        transformed = transform_matrix(matrix, spec["transform"])
        gene_matrix = collapse_matrix(transformed, feature_map)
        mapped_feature_n = sum(str(feature) in feature_map for feature in matrix.index)
        probe_path = MATRIX_DIR / f"{cohort}_feature_level.tsv.gz"
        gene_path = MATRIX_DIR / f"{cohort}_gene_level.tsv.gz"
        write_matrix(probe_path, transformed)
        write_matrix(gene_path, gene_matrix)
        all_qc.extend(qc_rows(cohort, gene_matrix))
        included = [row for row in freeze_by_cohort[cohort] if row["include_phase1b"] == "YES"]
        processing_rows.append({
            "dataset": cohort, "accession": accession, "assay_type": first(parsed[accession]["platforms"].get(spec["platform"], {}), "Platform_technology"),
            "platform": spec["platform"], "expression_unit": spec["expression_unit"],
            "log_transformed_input": "YES" if "log2" in spec["expression_unit"].lower() or "log-ratio" in spec["expression_unit"].lower() else "NO_OR_NOT_EXPLICIT",
            "normalization_state": first(parsed[accession]["samples"][0]["fields"], "Sample_data_processing"),
            "prespecified_qc_transform": spec["transform"],
            "probe_annotation_version": first(parsed[accession]["platforms"].get(spec["platform"], {}), "Platform_title"),
            "gene_identifier_system": "GEO platform annotation; HGNC current symbols and Ensembl IDs",
            "duplicated_probe_rule": "Outcome-blind median across unambiguously mapped probes; ambiguous multi-gene probes excluded",
            "submitted_feature_n": matrix.shape[0], "unambiguously_mapped_feature_n": mapped_feature_n,
            "unmapped_or_ambiguous_feature_n": matrix.shape[0] - mapped_feature_n,
            "mapped_gene_n": gene_matrix.shape[0],
            "missing_gene_policy": "A cohort contributes no effect for an absent/unmapped gene; missing genes are never assigned zero expression or zero effect",
            "sample_qc_information": "Outcome-blind PCA coordinates and robust sample-correlation flags written to bulk_sample_qc.csv; no flagged sample automatically removed",
            "matrix_scope": f"{matrix.shape[0]} submitted features x {matrix.shape[1]} frozen samples; {gene_matrix.shape[0]} mapped genes; {len(included)} main-signature samples",
            "feature_matrix_path": probe_path.relative_to(ROOT).as_posix(), "feature_matrix_sha256": sha256(probe_path),
            "gene_matrix_path": gene_path.relative_to(ROOT).as_posix(), "gene_matrix_sha256": sha256(gene_path),
            "git_tracking": "NOT_TRACKED_REBUILDABLE_INTERIM", "source": spec["source"],
        })

    processing_rows.append({
        "dataset": "GSE272342", "accession": "GSE272342", "assay_type": "high-throughput sequencing",
        "platform": "GPL16791", "expression_unit": "raw HTSeq counts per placenta library",
        "log_transformed_input": "NO", "normalization_state": "raw feature counts; normalization deferred",
        "prespecified_qc_transform": "No matrix generated in the main freeze because this is pregnancy-aware design-specific sensitivity only",
        "probe_annotation_version": "hg38/iGenomes submitted annotation; TopHat2/HTSeq pipeline",
        "gene_identifier_system": "submitted feature-count identifiers; formal harmonization deferred until sensitivity activation",
        "duplicated_probe_rule": "NOT_APPLICABLE", "sample_qc_information": "45 library records frozen; 29 pregnancy independence groups",
        "submitted_feature_n": "UNRESOLVED_UNTIL_SENSITIVITY_ACTIVATION", "unambiguously_mapped_feature_n": "NOT_APPLICABLE",
        "unmapped_or_ambiguous_feature_n": "NOT_APPLICABLE", "mapped_gene_n": "NOT_APPLICABLE",
        "missing_gene_policy": "Deferred with pregnancy-aware sensitivity activation; missing genes will not be assigned zero effect",
        "matrix_scope": "45 placenta libraries from 29 pregnancies; 27 PE and 18 control libraries; 20 PE and 9 control pregnancies",
        "feature_matrix_path": "NOT_GENERATED_ROLE_RESTRICTED", "feature_matrix_sha256": "NOT_APPLICABLE",
        "gene_matrix_path": "NOT_GENERATED_ROLE_RESTRICTED", "gene_matrix_sha256": "NOT_APPLICABLE",
        "git_tracking": "NOT_APPLICABLE", "source": DESIGN["cohorts"]["GSE272342"]["source"],
    })
    mapping_rows.sort(key=lambda r: (r["platform"], r["original_probe_id"]))
    all_qc.sort(key=lambda r: (r["dataset"], r["sample_id"]))
    return mapping_rows, processing_rows, all_qc


def cohort_registry(freeze_rows: list[dict]) -> list[dict]:
    grouped = defaultdict(list)
    for row in freeze_rows:
        grouped[row["dataset"]].append(row)
    result = []
    for cohort, spec in DESIGN["cohorts"].items():
        rows = grouped[cohort]
        included = [row for row in rows if row["include_phase1b"] == "YES"]
        total_group = defaultdict(set)
        included_group = defaultdict(set)
        for row in rows:
            total_group[row["PE/control"]].add(row["independence_group"])
        for row in included:
            included_group[row["PE/control"]].add(row["independence_group"])
        clinical_fields = ["GA", "maternal_age", "BMI", "FGR", "fetal sex", "ancestry", "batch", "labor", "delivery_mode", "singleton/twin"]
        available = []
        for field in clinical_fields:
            n = sum(row[field] != UNRESOLVED for row in rows)
            if n:
                available.append(f"{field}:{n}/{len(rows)}")
        result.append({
            "dataset": cohort, "accession": spec["accession"], "role_phase1b": spec["role"],
            "strategy": "STRATEGY_A_SOURCE_COHORT" if cohort.startswith("GSE") and cohort in {"GSE30186", "GSE10588", "GSE24129", "GSE25906", "GSE43942", "GSE4707", "GSE44711"} else ("STRATEGY_A_DIRECT_BIOBANK" if cohort == "GSE75010_BIOBANK" else "NOT_APPLICABLE"),
            "frozen_library_or_sample_n": len(rows), "frozen_independent_unit_n": len({row["independence_group"] for row in rows}),
            "frozen_pe_independent_n": len(total_group["PE"]), "frozen_control_independent_n": len(total_group["CONTROL"]),
            "phase1b_main_sample_n": len(included), "phase1b_main_pe_n": len(included_group["PE"]),
            "phase1b_main_control_n": len(included_group["CONTROL"]),
            "tissue": "human placenta/chorionic villi; cohort-specific submitted material retained",
            "genome_wide": "NO_TARGETED_PANEL" if cohort == "GSE190639" else "YES",
            "sample_independence": "PREGNANCY_CLUSTERED_45_LIBRARIES_29_PREGNANCIES" if cohort == "GSE272342" else "ONE_PLACENTA_SAMPLE_PER_SUBMITTED_PREGNANCY_NO_DUPLICATES_IDENTIFIED",
            "recoverable_clinical_covariates": ";".join(available) or "NONE_SAMPLE_LEVEL",
            "decision_reason": spec["reason"], "source": spec["source"],
        })
    return result


def overlap_rows(parsed: dict[str, dict]) -> list[dict]:
    source_counts = {
        "GSE30186": 12, "GSE10588": 43, "GSE24129": 16, "GSE25906": 60,
        "GSE43942": 12, "GSE4707": 14, "GSE44711": 16,
    }
    rows = [{
        "dataset_a": "GSE75010_COMPOSITE_330", "dataset_b": "GSE75010_BIOBANK",
        "overlap_n": 157, "overlap_unit": "exact direct GSM/biological samples",
        "overlap_identifiers": "GSM1940492-GSM1940648 (157 direct Series samples)",
        "independence_decision": "SAME_SAMPLES_NOT_INDEPENDENT; Strategy A uses BioBank and never adds the composite",
        "source": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75010",
    }]
    for accession, count in source_counts.items():
        rows.append({
            "dataset_a": "GSE75010_COMPOSITE_330", "dataset_b": accession,
            "overlap_n": count, "overlap_unit": "exact reanalysis-of GSM",
            "overlap_identifiers": f"All {count} eligible PE/control GSM from {accession} included in the GSE75010 composite",
            "independence_decision": "SAME_SAMPLES_NOT_INDEPENDENT; Strategy A uses source cohort and never adds the composite",
            "source": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75010|https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
        })
    rel_gsms = [item.split(": ", 1)[1] for item in parsed["GSE234729"]["series"].get("Series_relation", []) if item.startswith("Reanalysis of: GSM")]
    rows.append({
        "dataset_a": "GSE234729_ANALYTIC_123", "dataset_b": "GSE186257",
        "overlap_n": 12, "overlap_unit": "exact GSM reanalysis",
        "overlap_identifiers": "|".join(rel_gsms),
        "independence_decision": "COUNT_ONCE_INSIDE_GSE234729_ANALYTIC_123; never add GSE186257 as independent evidence",
        "source": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234729|https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186257",
    })
    rows.append({
        "dataset_a": "GSE272342", "dataset_b": "GSE203507", "overlap_n": 3,
        "overlap_unit": "same biological placentas re-extracted/resequenced",
        "overlap_identifiers": "S06|S18|S20",
        "independence_decision": "NOT_INDEPENDENT; do not combine as separate validation samples",
        "source": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272342|https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203507",
    })
    return rows


def exclusion_rows(freeze_rows: list[dict]) -> list[dict]:
    rows = []
    for row in freeze_rows:
        if row["include_phase1b"] == "NO":
            rows.append({
                "dataset": row["dataset"], "sample_id": row["GSM/sample ID"],
                "patient_or_pregnancy_id": row["patient/pregnancy ID if known"],
                "exclusion_scope": "PHASE1B_MAIN_STABLE_DISEASE_SIGNATURE_ONLY",
                "exclusion_reason": row["exclusion_reason"],
                "retained_elsewhere": row["planned_phase1b_use"],
                "outcome_separation_used": "NO", "source": row["source"],
            })
    return rows


def risk_rows(freeze_rows: list[dict], qc: list[dict]) -> list[dict]:
    risks = [
        ("P1A-R001", "CRITICAL", "GSE234729_ANALYTIC_123", "Sample-level gestational age and sequencing batch are absent from the public crosswalk although the paper adjusted for both.", "Do not use an unadjusted de novo contrast. Recover the crosswalk or use published adjusted evidence only.", "OPEN", DESIGN["cohorts"]["GSE234729_ANALYTIC_123"]["source"]),
        ("P1A-R002", "HIGH", "GSE234729_ANALYTIC_123", "GEO and the publication supplement report different STAR/R/DESeq2 versions; aggregate ancestry/control counts also differ across text versions.", "Freeze exact files and hashes; report the discrepancy and do not silently reconcile it.", "OPEN", DESIGN["cohorts"]["GSE234729_ANALYTIC_123"]["source"]),
        ("P1A-R003", "CRITICAL", "GSE75010_COMPOSITE_330", "The composite reuses 157 BioBank samples and 173 samples from seven historical GSE cohorts.", "Use Strategy A or Strategy B, never both. This freeze locks Strategy A.", "MITIGATED_BY_DESIGN", DESIGN["strategy"]["source"]),
        ("P1A-R004", "HIGH", "GSE75010_BIOBANK", "The non-PE group includes chronic-hypertension and preterm-control strata; processing batches of 20 are not exposed sample-by-sample.", "Exclude CH from the main contrast, model GA, and retain preterm-control and latent-batch sensitivity analyses.", "PARTIALLY_MITIGATED", DESIGN["cohorts"]["GSE75010_BIOBANK"]["source"]),
        ("P1A-R005", "HIGH", "GSE30186|GSE10588|GSE43942", "Sample-level gestational age is unavailable.", "Never impute GA; use cohort-wise estimates, heterogeneity diagnostics and leave-one-cohort-out analysis.", "OPEN", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30186|https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10588|https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43942"),
        ("P1A-R006", "HIGH", "GSE44711", "EOPE cases and controls have strong gestational-age separation.", "Subtype validation only; do not include in the core stable-gene definition.", "MITIGATED_BY_ROLE", DESIGN["cohorts"]["GSE44711"]["source"]),
        ("P1A-R007", "HIGH", "GSE272342", "All controls are twin pregnancies while cases include singleton and twin pregnancies; 45 placentas arise from 29 pregnancies.", "Pregnancy-aware design-specific sensitivity only; never count co-twin placentas as independent mothers.", "MITIGATED_BY_ROLE", DESIGN["cohorts"]["GSE272342"]["source"]),
        ("P1A-R008", "MODERATE", "GSE24129", "Isolated FGR is a third disease group and cannot be mislabeled as control or PE.", "Exclude the eight isolated-FGR samples from the main PE contrast but retain them in the freeze.", "MITIGATED_BY_EXCLUSION", DESIGN["cohorts"]["GSE24129"]["source"]),
        ("P1A-R009", "MODERATE", "GSE4707", "Two-color arrays use a pooled reference made from three normal placentas and have only four individual normal arrays.", "Subtype direction validation only.", "MITIGATED_BY_ROLE", DESIGN["cohorts"]["GSE4707"]["source"]),
        ("P1A-R010", "MODERATE", "GSE190639", "NanoString immune panel is targeted and cannot support genome-wide discovery.", "Subtype validation only.", "MITIGATED_BY_ROLE", DESIGN["cohorts"]["GSE190639"]["source"]),
        ("P1A-R011", "HIGH", "ALL_COHORTS", "Microarray/RNA-seq platforms and clinical covariate availability are heterogeneous.", "Analyze each cohort separately and combine effects/evidence only; no cross-study ComBat or pooled discovery matrix.", "MITIGATED_BY_PLAN", "docs/PHASE1B_STATISTICAL_ANALYSIS_PLAN.md"),
    ]
    rows = [{"risk_id": a, "severity": b, "dataset": c, "risk": d, "required_action": e, "status": f, "source": g} for a, b, c, d, e, f, g in risks]
    flagged = [row for row in qc if row["qc_flag"] != "PASS"]
    if flagged:
        rows.append({
            "risk_id": "P1A-R012", "severity": "MODERATE", "dataset": "MULTIPLE",
            "risk": f"Outcome-blind QC flagged {len(flagged)} sample(s) for review; no sample was removed.",
            "required_action": "Inspect raw/processed provenance and sensitivity influence before any exclusion; outcome separation is not a valid deletion criterion.",
            "status": "OPEN_FLAG_ONLY", "source": "results/01_phase1a/bulk_sample_qc.csv",
        })
    return rows


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    MATRIX_DIR.mkdir(parents=True, exist_ok=True)
    parsed = {accession: parse_soft(accession) for accession in ACCESSIONS}
    freeze = build_sample_freeze(parsed)
    mapping, processing, qc = expression_and_mapping(parsed, freeze)
    cohorts = cohort_registry(freeze)
    overlaps = overlap_rows(parsed)
    exclusions = exclusion_rows(freeze)
    risks = risk_rows(freeze, qc)

    write_csv(OUT / "bulk_sample_freeze.csv", freeze)
    write_csv(OUT / "bulk_cohort_registry.csv", cohorts)
    write_csv(OUT / "bulk_overlap_freeze.csv", overlaps)
    write_csv(OUT / "gene_mapping_registry.csv", mapping)
    write_csv(OUT / "bulk_processing_registry.csv", processing)
    write_csv(OUT / "bulk_exclusion_log.csv", exclusions)
    write_csv(OUT / "phase1a_risk_flags.csv", risks)
    write_csv(OUT / "bulk_sample_qc.csv", qc)
    print(f"Wrote Phase 1A freeze: {len(freeze)} sample rows, {len(mapping)} mapping rows, {len(qc)} QC rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
