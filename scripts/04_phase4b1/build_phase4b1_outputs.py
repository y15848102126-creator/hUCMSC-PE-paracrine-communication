#!/usr/bin/env python3
"""Build Phase 4B.1 protein-source correction tables from public protein files.

This script changes only HUCMSC_PROTEIN_SOURCE and mechanically re-applies the
Phase 4B classification predicates. All other evidence dimensions are copied
from the frozen Phase 4B outputs.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw/phase4b1"
OUT = ROOT / "results/04_phase4b1"
CFG = json.loads((ROOT / "config/phase4b1_analysis.json").read_text(encoding="utf-8"))
CANDIDATES = CFG["frozen_candidates"]

UNIPROT = {
    "ADAM17": "P78536", "AGRN": "O00468", "COL18A1": "P39060", "DCN": "P07585",
    "ENPP1": "P22413", "FURIN": "P09958", "GDF11": "O95390", "GRN": "P28799",
    "HSPG2": "P98160", "MDK": "P21741", "NAMPT": "P43490", "NID1": "P14543",
    "PSEN1": "P49768", "SERPINE1": "P05121", "TIMP1": "P01033",
    "TIMP2": "P16035", "WNT5A": "P41221",
}
ALIASES = {
    "SERPINE1": "PAI1;PAI-1;plasminogen activator inhibitor 1",
    "DCN": "PGS2;decorin", "HSPG2": "PGBM;perlecan", "MDK": "MK;midkine",
    "NID1": "nidogen-1;entactin", "TIMP1": "metalloproteinase inhibitor 1",
    "TIMP2": "metalloproteinase inhibitor 2", "WNT5A": "Wnt-5a",
}

URLS = {
    "PXD056371": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD056371",
    "PXD020948": "https://www.ebi.ac.uk/pride/archive/projects/PXD020948",
    "YU2024_SUPP": "https://doi.org/10.3390/ijms25094758",
    "FIGUEROA2025_SUPP": "https://doi.org/10.1186/s12951-024-03088-x",
    "PXD033723": "https://www.iprox.cn/page/project.html?id=IPX0004396000",
    "PMID32967723": "https://pubmed.ncbi.nlm.nih.gov/32967723/",
    "PXD036694": "https://www.iprox.cn/page/project.html?id=IPX0005018000",
    "PXD044276": "https://www.iprox.cn/page/project.html?id=IPX0006849000",
    "PXD022174": "https://www.ebi.ac.uk/pride/archive/projects/PXD022174",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def pxd056371() -> dict[str, dict]:
    """Count direct-accession PSMs in the three hUCMSC-EV runs (ms_run 4-6)."""
    path = RAW / "PXD056371_mzTab.mztab"
    accession_to_gene = {v: k for k, v in UNIPROT.items()}
    records: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    header = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("PSH\t"):
                header = line.rstrip("\n").split("\t")
            elif line.startswith("PSM\t") and header:
                values = line.rstrip("\n").split("\t")
                row = dict(zip(header, values))
                gene = accession_to_gene.get(row.get("accession", ""))
                match = re.search(r"ms_run\[(\d+)\]", row.get("spectra_ref", ""))
                if gene and match and match.group(1) in {"4", "5", "6"}:
                    records[gene].append((match.group(1), row.get("sequence", ""), row.get("accession", "")))
    out = {}
    for gene in CANDIDATES:
        recs = records.get(gene, [])
        psm = Counter(x[0] for x in recs)
        peptides = {run: len({x[1] for x in recs if x[0] == run}) for run in ("4", "5", "6")}
        positive_runs = sum(peptides[r] > 0 for r in peptides)
        out[gene] = {
            "detected": positive_runs > 0,
            "identifier": UNIPROT[gene],
            "unique_peptides": ";".join(f"Huexo{int(r)-3}:{peptides[r]}" for r in ("4", "5", "6")),
            "quantitative": ";".join(f"Huexo{int(r)-3}_PSM:{psm.get(r, 0)}" for r in ("4", "5", "6")),
            "positive_runs": positive_runs,
            "fdr": "protein-group combined FDR score in mzTab; exact q-value not exported",
            "limitations": "Three hUCMSC-EV MS runs; donor independence not reported. Direct accession PSMs required; ambiguity-only WNT5A does not count.",
        }
    return out


def pxd020948() -> dict[str, dict]:
    path = RAW / "PXD020948_MaxQuant_txt.zip"
    with zipfile.ZipFile(path) as archive:
        table = pd.read_csv(archive.open("txt/proteinGroups.txt"), sep="\t", low_memory=False)
    out = {}
    for gene in CANDIDATES:
        symbol_lists = table["Gene names"].fillna("").astype(str).str.split(";")
        mask = symbol_lists.apply(lambda xs: gene in xs)
        rows = table.loc[mask]
        if rows.empty:
            out[gene] = {"detected": False, "identifier": UNIPROT[gene], "unique_peptides": "0;0;0", "quantitative": "0;0;0", "positive_runs": 0, "fdr": "NOT_APPLICABLE", "limitations": "Complete MaxQuant table searched; no unambiguous candidate group in UC-MSC runs."}
            continue
        row = rows.iloc[0]
        peptides = [int(row.get(f"Unique peptides U_MSC_{i}", 0) or 0) for i in (1, 2, 3)]
        intensity = [float(row.get(f"Intensity U_MSC_{i}", 0) or 0) for i in (1, 2, 3)]
        positive_runs = sum(p > 0 and x > 0 for p, x in zip(peptides, intensity))
        out[gene] = {
            "detected": positive_runs > 0,
            "identifier": str(row.get("Majority protein IDs", UNIPROT[gene])),
            "unique_peptides": ";".join(f"U_MSC_{i}:{p}" for i, p in enumerate(peptides, 1)),
            "quantitative": ";".join(f"U_MSC_{i}_intensity:{x:g}" for i, x in enumerate(intensity, 1)),
            "positive_runs": positive_runs,
            "fdr": f"MaxQuant protein-group q={row.get('Q-value', 'NOT_REPORTED')}",
            "limitations": "Three UC-MSC exosome MS runs; cells were commercial and independent donor n is unresolved. ENPP1 is present in one run with one unique peptide only.",
        }
    return out


def yu2024() -> dict[str, dict]:
    path = RAW / "yu2024_supplement/supplementary table 1.xlsx"
    table = pd.read_excel(path, sheet_name="Proteins", header=1)
    descriptions = table["Description"].fillna("").astype(str)
    out = {}
    for gene in CANDIDATES:
        mask = descriptions.str.contains(rf"GN={re.escape(gene)}(?: |$)", regex=True)
        rows = table.loc[mask]
        if rows.empty:
            out[gene] = {"detected": False, "identifier": UNIPROT[gene], "unique_peptides": "0", "quantitative": "WJ-1:0;WJ-2:0", "positive_runs": 0, "fdr": "NOT_APPLICABLE", "limitations": "Complete 1,695-protein Supplementary Table S1 searched; no unambiguous candidate row."}
            continue
        row = rows.iloc[0]
        a1 = float(row["Abundances (Grouped): WJ-1"] or 0)
        a2 = float(row["Abundances (Grouped): WJ-2"] or 0)
        out[gene] = {
            "detected": a1 > 0 or a2 > 0,
            "identifier": str(row["Accession"]),
            "unique_peptides": str(int(row["# Unique Peptides"])),
            "quantitative": f"WJ-1:{a1:g};WJ-2:{a2:g}",
            "positive_runs": int(a1 > 0) + int(a2 > 0),
            "fdr": f"combined experimental q={row['Exp. q-value: Combined']}",
            "limitations": "Two WJ-MSC exosome protein samples; donor independence is not explicit. Study inclusion required FDR<0.01 and >=2 unique peptides.",
        }
    return out


def figueroa2025() -> dict[str, dict]:
    path = RAW / "12951_2024_3088_MOESM5_ESM.xlsx"
    table = pd.read_excel(path, sheet_name="Identified proteins")
    descriptions = table["Fasta headers"].fillna("").astype(str)
    out = {}
    for gene in CANDIDATES:
        mask = descriptions.str.contains(rf"GN={re.escape(gene)}(?: |$)", regex=True)
        rows = table.loc[mask]
        if rows.empty:
            out[gene] = {"detected": False, "identifier": UNIPROT[gene], "unique_peptides": "0", "quantitative": "D1:0;D2:0;D3:0", "positive_runs": 0, "fdr": "NOT_APPLICABLE", "limitations": "Complete 628-row identified-protein table searched; no unambiguous candidate group."}
            continue
        row = rows.iloc[0]
        header = str(row["Fasta headers"])
        # Protein groups spanning distinct frozen candidates are ambiguous.
        present_candidates = [g for g in CANDIDATES if re.search(rf"GN={re.escape(g)}(?: |$)", header)]
        values = [float(row.get(f"LFQ intensity D{i}_CTL", 0) or 0) for i in (1, 2, 3)]
        ambiguous = len(present_candidates) > 1 or (gene == "WNT5A" and re.search(r"GN=WNT5B(?: |$)", header) is not None)
        positive_runs = 0 if ambiguous else sum(v > 0 for v in values)
        acc = re.findall(r"(?:sp|tr)\|([^|]+)\|", header)
        out[gene] = {
            "detected": positive_runs > 0,
            "identifier": ";".join(acc) if acc else UNIPROT[gene],
            "unique_peptides": str(int(row["Peptides"])),
            "quantitative": ";".join(f"D{i}:{v:g}" for i, v in enumerate(values, 1)),
            "positive_runs": positive_runs,
            "fdr": "protein FDR 2.5% (study method)",
            "limitations": "Three independent UC-MSC donors. WNT5A/WNT5B group is ambiguous and not counted; MDK row has zero LFQ in all donors.",
        }
    return out


def dataset_registry() -> pd.DataFrame:
    rows = [
        ["PXD056371", "39472581", "10.1038/s41392-024-01993-z", 2024, "human umbilical-cord MSC", "passage 5", "EV_EXOSOME", "3 HU EV MS runs; donor n unresolved", "MaxQuant/mzTab; LC-MS/MS", "YES", "COMPLETE_PUBLIC_MZTAB", "TECHNICALLY_EVALUABLE", "PXD056371_mzTab.mztab", "Direct accession PSMs allow HU-run-specific calls; NID1 Western blot is reported in the supplement."],
        ["PXD020948", "33246507", "10.1186/s13287-020-02032-8", 2020, "commercial human UC-MSC", "passage 3; 48h serum-free", "EV_EXOSOME", "3 UC-MSC MS runs; donor n unresolved", "MaxQuant/LC-MS-MS", "YES", "COMPLETE_PUBLIC_MAXQUANT", "TECHNICALLY_EVALUABLE", "PXD020948_MaxQuant_txt.zip", "Run-specific unique peptide and intensity columns; source cells may not represent independent donors."],
        ["YU2024_SUPP", "38731977", "10.3390/ijms25094758", 2024, "Wharton's-jelly MSC", "NOT_REPORTED", "EV_EXOSOME", "2 WJ-MSC exosome protein samples; donor independence unresolved", "iTRAQ LC-Q-Exactive", "NO_REPOSITORY_ACCESSION", "COMPLETE_SUPPLEMENTARY_TABLE", "TECHNICALLY_EVALUABLE", "supplementary table 1.xls", "1,695 proteins; FDR<0.01 and >=2 unique peptides."],
        ["FIGUEROA2025_SUPP", "39806427", "10.1186/s12951-024-03088-x", 2025, "human UC-MSC", "passage 5; 48h serum-free", "EV_EXOSOME", "3 independent UC donors", "MaxQuant LFQ; nanoLC-QTOF", "NO_REPOSITORY_ACCESSION", "COMPLETE_SUPPLEMENTARY_TABLE", "TECHNICALLY_EVALUABLE", "12951_2024_3088_MOESM5_ESM.xlsx", "628 identified protein groups; protein FDR 2.5%; donor-specific LFQ."],
        ["PXD033723", "NOT_REPORTED", "10.1002/pmic.202200204", 2023, "human UC-MSC", "NOT_REPORTED", "EV_EXOSOME", "NOT_REPORTED", "LC-MS/MS", "YES", "RAW_PUBLIC_SEARCH_OUTPUT_ZERO_BYTE", "NOT_CURRENTLY_EVALUABLE", "iProX 1.txt", "Publication reports ~4,200 proteins, but public candidate-level search output currently resolves to a zero-byte object; absence cannot be called."],
        ["PMID32967723", "32967723", "10.1186/s13287-020-01931-0", 2020, "umbilical-cord MSC-related material", "preterm and term; see paper", "EV_EXOSOME_OR_RELATED_PEPTIDOME", "3 preterm + 3 term donors", "LC-MS peptidome", "NO", "SELECTED_SUPPLEMENTARY_TABLES_ONLY", "NOT_ASSESSABLE_FOR_ABSENCE", "four supplementary DOCX tables", "Identifier mapping was checked; tables are selected/differential peptide lists, not a complete protein universe."],
        ["PXD036694", "NOT_REPORTED", "NOT_REPORTED", 2022, "BM/AD/UC-MSC conditions", "high-glucose conditioned medium", "SOLUBLE_CONDITIONED_MEDIUM", "NOT_REPORTED", "raw MS plus proprietary .msf", "YES", "RAW_AND_PROPRIETARY_OUTPUT_ONLY", "NOT_CURRENTLY_EVALUABLE", "IPX0005018000", "Relevant UC condition exists, but no flat public identification table supports candidate-level calls."],
        ["PXD044276", "NOT_REPORTED", "10.1038/s41598-024-79063-1", 2024, "human WJ-MSC", "NOT_REPORTED", "CELL_LYSATE", "NOT_REPORTED", "LC-MS/MS", "YES", "CELLULAR_PROTEOME", "CONTEXT_ONLY_NOT_EXTRACELLULAR", "IPX0006849000", "Biological compartment is cell lysate, so it cannot establish extracellular source evidence."],
        ["PXD022174", "NOT_REPORTED", "NOT_REPORTED", 2021, "human bone-marrow clonal MSC", "NOT_REPORTED", "EV_EXOSOME", "NOT_REPORTED", "LC-MS/MS", "YES", "COMPLETE_REPOSITORY", "EXCLUDE_WRONG_MSC_SOURCE", "PXD022174", "Not UC/WJ-MSC."],
    ]
    columns = ["dataset_accession", "PMID", "DOI", "year", "biological_source", "passage_or_condition", "extracellular_compartment", "donor_or_replicate_n", "identification_method", "public_repository", "public_data_state", "technical_evaluability", "file_or_record", "limitations"]
    df = pd.DataFrame(rows, columns=columns)
    df["local_file_size_bytes"] = df["file_or_record"].map({
        "PXD056371_mzTab.mztab": (RAW / "PXD056371_mzTab.mztab").stat().st_size,
        "PXD020948_MaxQuant_txt.zip": (RAW / "PXD020948_MaxQuant_txt.zip").stat().st_size,
        "supplementary table 1.xls": (RAW / "yu2024_supplement/supplementary table 1.xls").stat().st_size,
        "12951_2024_3088_MOESM5_ESM.xlsx": (RAW / "12951_2024_3088_MOESM5_ESM.xlsx").stat().st_size,
    }).fillna("NOT_APPLICABLE")
    df["local_file_sha256"] = df["file_or_record"].map({
        "PXD056371_mzTab.mztab": sha256(RAW / "PXD056371_mzTab.mztab"),
        "PXD020948_MaxQuant_txt.zip": sha256(RAW / "PXD020948_MaxQuant_txt.zip"),
        "supplementary table 1.xls": sha256(RAW / "yu2024_supplement/supplementary table 1.xls"),
        "12951_2024_3088_MOESM5_ESM.xlsx": sha256(RAW / "12951_2024_3088_MOESM5_ESM.xlsx"),
    }).fillna("NOT_APPLICABLE")
    df["source_url"] = df["dataset_accession"].map(URLS)
    return df


def make_detection_matrix(evaluable: dict[str, dict[str, dict]]) -> pd.DataFrame:
    meta = {
        "PXD056371": ["human umbilical-cord MSC", "passage 5", "EV_EXOSOME", "3 HU EV MS runs; donor n unresolved", "LC-MS/MS; MaxQuant mzTab", "baseline"],
        "PXD020948": ["commercial human UC-MSC", "passage 3; 48h serum-free", "EV_EXOSOME", "3 MS runs; donor n unresolved", "LC-MS/MS; MaxQuant", "baseline"],
        "YU2024_SUPP": ["Wharton's-jelly MSC", "NOT_REPORTED", "EV_EXOSOME", "2 protein samples; donor independence unresolved", "iTRAQ LC-Q-Exactive", "baseline"],
        "FIGUEROA2025_SUPP": ["human UC-MSC", "passage 5; 48h serum-free", "EV_EXOSOME", "3 independent donors", "nanoLC-QTOF; MaxQuant LFQ", "baseline"],
    }
    rows = []
    for dataset, calls in evaluable.items():
        for gene in CANDIDATES:
            c = calls[gene]
            rows.append({
                "candidate": gene, "dataset_accession": dataset, "UC_WJ_MSC_source": meta[dataset][0],
                "passage_condition": meta[dataset][1], "extracellular_compartment": meta[dataset][2],
                "donor_or_replicate_n": meta[dataset][3], "protein_detected": "YES" if c["detected"] else "NO",
                "identification_method": meta[dataset][4], "uniprot_or_protein_identifier": c["identifier"],
                "protein_FDR_or_confidence": c["fdr"], "unique_peptide_evidence": c["unique_peptides"],
                "quantitative_abundance": c["quantitative"], "positive_sample_or_run_n": c["positive_runs"],
                "western_or_ELISA_confirmation": "NID1_WESTERN_REPORTED" if dataset == "PXD056371" and gene == "NID1" else "NO_OR_NOT_REPORTED",
                "baseline_vs_engineered": meta[dataset][5], "limitations": c["limitations"], "source_url": URLS[dataset],
            })
    for dataset in ("PXD033723", "PMID32967723", "PXD036694"):
        for gene in CANDIDATES:
            rows.append({
                "candidate": gene, "dataset_accession": dataset, "UC_WJ_MSC_source": "see dataset registry",
                "passage_condition": "NOT_REPORTED", "extracellular_compartment": "EV_EXOSOME" if dataset != "PXD036694" else "SOLUBLE_CONDITIONED_MEDIUM",
                "donor_or_replicate_n": "see dataset registry", "protein_detected": "NOT_ASSESSABLE",
                "identification_method": "public table insufficient for candidate-level calls", "uniprot_or_protein_identifier": UNIPROT[gene],
                "protein_FDR_or_confidence": "NOT_ASSESSABLE", "unique_peptide_evidence": "NOT_ASSESSABLE",
                "quantitative_abundance": "NOT_ASSESSABLE", "positive_sample_or_run_n": "NOT_ASSESSABLE",
                "western_or_ELISA_confirmation": "NOT_ASSESSABLE", "baseline_vs_engineered": "baseline",
                "limitations": "Incomplete, empty, or proprietary public output cannot support a non-detection.", "source_url": URLS[dataset],
            })
    return pd.DataFrame(rows)


def source_summary(detections: pd.DataFrame) -> pd.DataFrame:
    old = pd.read_csv(ROOT / "results/04_phase4b/integration/phase4b_candidate_evidence_matrix.csv").set_index("candidate")
    rows = []
    for gene in CANDIDATES:
        d = detections[(detections.candidate == gene) & (detections.protein_detected == "YES")]
        datasets = sorted(d.dataset_accession.unique())
        direct = bool(datasets)
        if direct:
            cls = "HUCMSC_EV_PROTEIN_DIRECT"
            strength = "MULTI_DATASET" if len(datasets) >= 2 else "SINGLE_DATASET"
            evidence = f"Direct unambiguous EV protein detection in {len(datasets)} technically evaluable dataset(s): {';'.join(datasets)}"
        else:
            cls = "NO_PROTEIN_SOURCE_EVIDENCE"
            strength = "EVIDENCE_ABSENT_NOT_AGAINST"
            evidence = "No unambiguous detection in four technically evaluable UC/WJ-MSC extracellular protein tables"
        previous = old.loc[gene, "HUCMSC_PROTEIN_SOURCE"]
        rows.append({
            "candidate": gene, "corrected_protein_source_classification": cls,
            "extracellular_compartment": "EV_EXOSOME" if direct else "NOT_DEMONSTRATED",
            "direct_dataset_n": len(datasets), "direct_datasets": ";".join(datasets) if datasets else "NONE",
            "evidence_strength": strength, "soluble_CM_direct": "NO", "EV_direct": "YES" if direct else "NO",
            "ECM_direct": "NO", "cell_lysate_only": "NO", "protein_source_evidence": evidence,
            "previous_phase4b_classification": previous,
            "phase4b_false_negative": "YES" if previous == "NO_PROTEIN_SOURCE_EVIDENCE" and direct else "NO",
            "limitations": "EV cargo is not equivalent to a freely soluble paracrine ligand. Absence is not evidence against secretion." if direct else "FURIN/ADAM17 may be below detection, condition-specific, membrane-associated, or absent; no biological negative is inferred.",
            "source_url": "|".join(URLS[x] for x in datasets) if datasets else "results/04_phase4b1/candidate_protein_detection_matrix.csv",
        })
    return pd.DataFrame(rows)


def identifier_mapping_audit(detections: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in CANDIDATES:
        ids = sorted({x for v in detections[(detections.candidate == gene) & (detections.protein_detected == "YES")].uniprot_or_protein_identifier.astype(str) for x in v.split(";") if x and x != "nan"})
        rows.append({
            "candidate": gene, "HGNC_symbol": gene, "reviewed_UniProt_accession": UNIPROT[gene],
            "aliases_checked": ALIASES.get(gene, "protein name and HGNC symbol"),
            "observed_identifiers": ";".join(ids) if ids else "NONE_UNAMBIGUOUS",
            "mapping_status": "UNAMBIGUOUS" if ids else "NO_UNAMBIGUOUS_DETECTION",
            "ambiguity_note": "WNT5A/WNT5B shared protein groups were excluded; direct WNT5A rows in PXD020948 and Yu 2024 remain valid." if gene == "WNT5A" else "No cross-candidate protein group used as positive evidence.",
            "mapping_source": "reviewed human UniProt accession; GN field; documented aliases",
            "source_url": f"https://rest.uniprot.org/uniprotkb/{UNIPROT[gene]}.txt|results/04_phase4b1/candidate_protein_detection_matrix.csv",
        })
    return pd.DataFrame(rows)


def mechanical_reclassification(summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    old_e = pd.read_csv(ROOT / "results/04_phase4b/integration/phase4b_candidate_evidence_matrix.csv")
    old_c = pd.read_csv(ROOT / "results/04_phase4b/integration/phase4b_candidate_classification.csv").set_index("candidate")
    source = summary.set_index("candidate").corrected_protein_source_classification.to_dict()
    direct_classes = {"HUCMSC_EV_PROTEIN_DIRECT", "HUCMSC_SOLUBLE_CM_PROTEIN_DIRECT", "HUCMSC_ECM_PROTEIN_DIRECT"}
    plausible = {"CANONICAL_SOLUBLE_SECRETED", "ECM_ASSOCIATED", "SHED_SOLUBLE_FORM", "EXTRACELLULAR_ENZYME"}
    class_rows = []
    corrected = old_e.copy()
    for i, row in corrected.iterrows():
        gene = row["candidate"]
        corrected.at[i, "HUCMSC_PROTEIN_SOURCE"] = source[gene]
        rev = int(re.search(r"REVERSAL_AXES=(\d+)", row["EMPIRICAL_SIGNED_PERTURBATION"]).group(1))
        conc = int(re.search(r"CONCORDANT_AXES=(\d+)", row["EMPIRICAL_SIGNED_PERTURBATION"]).group(1))
        topo_ok = row["PARACRINE_TOPOLOGY"] in plausible
        protein_direct = source[gene] in direct_classes
        direct_prior = row["DIRECT_MSC_PE_PRECEDENT"] in {"DIRECT_HUCMSC_PE_MECHANISM_ALREADY_SHOWN", "DIRECT_MSC_PE_MECHANISM_ALREADY_SHOWN"}
        high = topo_ok and protein_direct and rev > 0 and conc == 0 and not direct_prior
        context = topo_ok and (protein_direct or source[gene] == "OTHER_MSC_EXTRACELLULAR_PROTEIN") and rev > 0 and row["MIXED_DIRECTION_RISK"] in {"MIXED_DIRECTION_CONTEXT_DEPENDENT", "PREDOMINANTLY_DISEASE_CONCORDANT"}
        computational = source[gene] in {"NO_PROTEIN_SOURCE_EVIDENCE", "HUCMSC_PROTEIN_CELL_ONLY"} and rev == 0
        weak = row["PARACRINE_TOPOLOGY"] in {"MEMBRANE_ASSOCIATED", "INTRACELLULAR_OR_QUESTIONABLE_PARACRINE", "UNCERTAIN"}
        low_novelty = direct_prior
        flags = []
        if high: flags.append("TRIANGULATED_HIGH_PRIORITY")
        if context: flags.append("TRIANGULATED_CONTEXT_DEPENDENT")
        if computational: flags.append("COMPUTATIONAL_ONLY")
        if weak: flags.append("BIOPHYSICALLY_WEAK_PARACRINE")
        if low_novelty: flags.append("KNOWN_MECHANISM_LOW_NOVELTY")
        corrected.at[i, "classification_flags"] = ";".join(flags) if flags else "NO_SUGGESTED_CATEGORY_THRESHOLD_MET"
        if high: primary = "TRIANGULATED_HIGH_PRIORITY"
        elif context: primary = "TRIANGULATED_CONTEXT_DEPENDENT"
        elif weak: primary = "BIOPHYSICALLY_WEAK_PARACRINE"
        elif computational: primary = "COMPUTATIONAL_ONLY"
        elif protein_direct and rev == 0: primary = "PROTEIN_SUPPORTED_BUT_DIRECTION_UNRESOLVED"
        elif rev > 0 and not protein_direct: primary = "PERTURBATION_SUPPORTED_BUT_SOURCE_UNCONFIRMED"
        else: primary = "PARTIAL_EXTERNAL_EVIDENCE"
        class_rows.append({
            "candidate": gene, "previous_primary_classification": old_c.loc[gene, "primary_classification"],
            "corrected_primary_classification": primary, "classification_changed": "YES" if primary != old_c.loc[gene, "primary_classification"] else "NO",
            "TRIANGULATED_HIGH_PRIORITY": "YES" if high else "NO", "TRIANGULATED_CONTEXT_DEPENDENT": "YES" if context else "NO",
            "COMPUTATIONAL_ONLY": "YES" if computational else "NO", "BIOPHYSICALLY_WEAK_PARACRINE": "YES" if weak else "NO",
            "KNOWN_MECHANISM_LOW_NOVELTY": "YES" if low_novelty else "NO", "topology_pass": "YES" if topo_ok else "NO",
            "direct_hucmsc_extracellular_protein": "YES" if protein_direct else "NO",
            "corrected_protein_source_classification": source[gene], "empirical_reversal_axis_n": rev,
            "empirical_concordant_or_conflict_axis_n": conc,
            "final_interpretation": "NO_CANDIDATE_IS_THERAPEUTICALLY_VALIDATED",
            "classification_rule_source": "config/phase4b_analysis.json; unchanged predicates",
            "source_url": "results/04_phase4b/integration/phase4b_candidate_classification.csv|results/04_phase4b1/corrected_hucmsc_protein_source_evidence.csv",
        })
    corrected["source_url"] = corrected["source_url"].astype(str) + "|results/04_phase4b1/corrected_hucmsc_protein_source_evidence.csv"
    return corrected, pd.DataFrame(class_rows)


def risk_flags(summary: pd.DataFrame, classifications: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame([
        ["P4B1_RISK_001", "HIGH", "COMPARTMENT", "All corrected direct evidence is EV/exosome cargo; no technically evaluable complete soluble-CM protein table was recovered.", "Do not translate EV detection into soluble-ligand availability.", "OPEN"],
        ["P4B1_RISK_002", "HIGH", "DONOR_INDEPENDENCE", "PXD056371 and PXD020948 contain three MS runs but independent donor identity is unresolved; Yu 2024 has two WJ samples with donor independence unstated.", "Count datasets and explicit donors separately; do not call run replication donor replication.", "OPEN"],
        ["P4B1_RISK_003", "HIGH", "ENPP1", "ENPP1 is strong in both Yu WJ samples (16 unique peptides) but only one peptide in one of three PXD020948 UC runs; its empirical reversal source remains non-placental.", "Treat mechanical high-priority status as synthesis eligibility, not PE therapeutic proof.", "OPEN"],
        ["P4B1_RISK_004", "MEDIUM", "LOW_DEPTH_CALLS", "MDK is supported by one peptide in one PXD056371 HU run; GDF11, GRN and PSEN1 are near the two-peptide floor in Yu 2024.", "Retain direct calls with explicit strength and require orthogonal protein validation before experiments.", "OPEN"],
        ["P4B1_RISK_005", "HIGH", "PXD033723", "The iProX search-output object for the ~4,200-protein Xu study is currently zero bytes; candidate calls cannot be reproduced.", "Do not infer non-detection; retry repository retrieval in final experimental planning if necessary.", "OPEN"],
        ["P4B1_RISK_006", "MEDIUM", "INCOMPLETE_TABLES", "PMID32967723 supplementary tables are selected peptidome results, and PXD036694 provides raw/proprietary output without a flat identification table.", "Mark all candidate cells NOT_ASSESSABLE rather than NO.", "OPEN"],
        ["P4B1_RISK_007", "MEDIUM", "AMBIGUOUS_GROUPS", "WNT5A appears in WNT5A/WNT5B ambiguity groups in PXD056371/Figueroa; those groups were excluded from positive evidence.", "Use only direct WNT5A evidence from PXD020948 and Yu 2024.", "MITIGATED"],
        ["P4B1_RISK_008", "HIGH", "SCOPE", "Protein-source correction does not update empirical perturbation, PE context, novelty, or mixed-direction evidence.", "Preserve upstream hashes and avoid interpreting source support as direction support.", "MITIGATED"],
    ], columns=["risk_id", "severity", "domain", "risk", "required_mitigation", "status"]).assign(source_url="results/04_phase4b1/candidate_protein_detection_matrix.csv|docs/PHASE4B1_PROTEIN_SOURCE_COMPLETENESS_REPORT.md")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    evaluable = {
        "PXD056371": pxd056371(), "PXD020948": pxd020948(),
        "YU2024_SUPP": yu2024(), "FIGUEROA2025_SUPP": figueroa2025(),
    }
    registry = dataset_registry()
    detections = make_detection_matrix(evaluable)
    summary = source_summary(detections)
    mappings = identifier_mapping_audit(detections)
    evidence, classifications = mechanical_reclassification(summary)
    risks = risk_flags(summary, classifications)

    tables = {
        "hucmsc_proteomics_dataset_registry.csv": registry,
        "candidate_protein_detection_matrix.csv": detections,
        "corrected_hucmsc_protein_source_evidence.csv": summary,
        "protein_identifier_mapping_audit.csv": mappings,
        "corrected_phase4b_candidate_evidence_matrix.csv": evidence,
        "corrected_phase4b_candidate_classification.csv": classifications,
        "phase4b1_risk_flags.csv": risks,
    }
    for name, table in tables.items():
        table.to_csv(OUT / name, index=False, quoting=csv.QUOTE_MINIMAL)
    print(
        "PHASE4B1_OUTPUTS_OK "
        f"evaluable_extracellular_datasets={registry.technical_evaluability.eq('TECHNICALLY_EVALUABLE').sum()} "
        f"direct_candidates={summary.EV_direct.eq('YES').sum()} "
        f"high_priority={classifications.TRIANGULATED_HIGH_PRIORITY.eq('YES').sum()}"
    )


if __name__ == "__main__":
    main()
