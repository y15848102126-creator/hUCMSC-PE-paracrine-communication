#!/usr/bin/env python3
"""Build auditable Phase 4B tables from frozen inputs and screened evidence."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/04_phase4b"
DATE = "2026-08-15"
CANDIDATES = ["ADAM17", "AGRN", "COL18A1", "DCN", "ENPP1", "FURIN", "GDF11", "GRN", "HSPG2", "MDK", "NAMPT", "NID1", "PSEN1", "SERPINE1", "TIMP1", "TIMP2", "WNT5A"]


def u(gene: str, accession: str) -> str:
    return f"https://rest.uniprot.org/uniprotkb/{accession}.json"


TOPOLOGY = {
    "ADAM17": ["YES", 1, "Cell membrane; single-pass type I", "NO_CANONICAL_SOLUBLE_LIGAND", "NO", "ADAM17 itself is membrane anchored; it sheds other ectodomains", "YES", "YES", "Cytoplasmic tail", "No reviewed soluble isoform", "MEMBRANE_ASSOCIATED"],
    "AGRN": ["YES", 0, "Secreted extracellular matrix (isoform 1); membrane isoform 2", "YES_SECRETED_ISOFORM", "YES", "Proteolytic fragments; isoform-dependent membrane attachment", "NO", "YES", "None for secreted isoform", "Secreted agrin and transmembrane agrin", "ECM_ASSOCIATED"],
    "COL18A1": ["YES", 0, "Secreted basement-membrane ECM", "YES_CLEAVED_ENDOSTATIN", "YES", "Cathepsin/elastase processing generates endostatin", "NO", "YES", "NO", "Endostatin/NC1 fragments", "ECM_ASSOCIATED"],
    "DCN": ["YES", 0, "Secreted extracellular matrix", "YES_ECM_BOUND", "YES", "No membrane shedding required", "NO", "YES", "NO", "Multiple splice isoforms; all evidence interpreted conservatively", "ECM_ASSOCIATED"],
    "ENPP1": ["NO", 1, "Cell membrane type II; extracellular catalytic domain", "YES_AFTER_PROTEOLYTIC_CLEAVAGE", "NO", "Secreted form produced by cleavage near residue 103", "YES", "YES", "Cytoplasmic N-terminus", "Reviewed secreted chain 103-925", "EXTRACELLULAR_ENZYME"],
    "FURIN": ["YES", 1, "Trans-Golgi/cell/endosome membrane; secreted form annotated", "YES_SHED_FORM", "NO", "Soluble form depends on processing/shedding", "YES", "YES", "Cytoplasmic tail; predominant trafficking is intracellular/membrane", "No separate reviewed soluble isoform; secreted localization annotated", "SHED_SOLUBLE_FORM"],
    "GDF11": ["YES", 0, "Secreted", "YES_MATURE_LIGAND", "NO", "Furin-like and tolloid-like processing required for activation", "NO", "YES", "NO", "Mature C-terminal ligand", "CANONICAL_SOLUBLE_SECRETED"],
    "GRN": ["YES", 0, "Secreted and lysosomal after endocytosis", "YES", "NO", "Extracellular/lysosomal proteolysis generates granulins", "NO", "YES", "Lysosomal pool also present", "Progranulin and granulin peptides", "CANONICAL_SOLUBLE_SECRETED"],
    "HSPG2": ["YES", 0, "Secreted basement-membrane ECM", "YES_ECM_BOUND_FRAGMENTS", "YES", "Proteolysis generates endorepellin/LG3", "NO", "YES", "NO", "Endorepellin and LG3", "ECM_ASSOCIATED"],
    "MDK": ["YES", 0, "Secreted", "YES", "NO", "No shedding required", "NO", "YES", "NO", "Two reviewed splice isoforms", "CANONICAL_SOLUBLE_SECRETED"],
    "NAMPT": ["NO", 0, "Predominantly nuclear/cytoplasmic; non-classical extracellular pool reported", "NONCLASSICAL_CONTEXT_DEPENDENT", "NO", "No signal peptide; secretion mechanism/context unresolved", "YES_INTRACELLULAR_ENZYME", "MIXED", "YES", "No reviewed secreted isoform", "INTRACELLULAR_OR_QUESTIONABLE_PARACRINE"],
    "NID1": ["YES", 0, "Secreted basement-membrane ECM", "YES_ECM_BOUND", "YES", "No shedding required", "NO", "YES", "NO", "Two reviewed splice isoforms", "ECM_ASSOCIATED"],
    "PSEN1": ["NO", 9, "ER/Golgi/endosome/cell membranes; multi-pass", "NO", "NO", "Endoproteolysis forms membrane gamma-secretase subunits, not a soluble ligand", "YES_INTRAMEMBRANE_PROTEASE", "NO", "YES", "No soluble isoform", "INTRACELLULAR_OR_QUESTIONABLE_PARACRINE"],
    "SERPINE1": ["YES", 0, "Secreted", "YES", "MATRIX_BINDING_POSSIBLE", "No shedding required", "NO", "YES", "NO", "Two splice isoforms", "CANONICAL_SOLUBLE_SECRETED"],
    "TIMP1": ["YES", 0, "Secreted", "YES", "ECM_INTERACTION", "No shedding required", "NO", "YES", "NO", "Canonical secreted protein", "CANONICAL_SOLUBLE_SECRETED"],
    "TIMP2": ["YES", 0, "Secreted", "YES", "ECM_INTERACTION", "No shedding required", "NO", "YES", "NO", "Canonical secreted protein", "CANONICAL_SOLUBLE_SECRETED"],
    "WNT5A": ["YES", 0, "Secreted extracellular space/ECM", "YES_BUT_LIPIDATED", "YES", "Secretion requires glycosylation; lipidation affects carrier/receptor binding", "NO", "YES", "NO", "Two splice isoforms", "CANONICAL_SOLUBLE_SECRETED"],
}


PROTEIN_EVIDENCE = {
    "ADAM17": [("37099959", "EXCLUDED_RECIPIENT_TARGET", "hUCMSC EV miR-26a-5p targets recipient ADAM17; no ADAM17 cargo measurement", "EV_ONLY", "NO_PROTEIN_SOURCE_EVIDENCE")],
    "DCN": [("40386129", "INCLUDED_OTHER_MSC_SECRETOME", "ADSC supernatant proteomics and functional trapping experiments identify decorin as an extracellular antifibrotic component", "SOLUBLE_OR_TOTAL_SECRETOME", "OTHER_MSC_SECRETOME"), ("40438789", "INCLUDED_HUCMSC_CELL_ENGINEERING", "DCN-overexpressing hUC-MSCs; the abstract establishes engineered cellular expression but not a quantitative baseline soluble secretome measurement", "ENGINEERED_CELL_PRODUCT_UNRESOLVED", "OTHER_MSC_SECRETOME")],
    "GDF11": [("40620065", "EXCLUDED_NOT_EV_CARGO", "UC-MSC EV miR-31-5p suppresses macrophage PDGFB and activates endogenous GDF11 downstream; GDF11 is not shown as EV cargo", "EV_DOWNSTREAM_ONLY", "NO_PROTEIN_SOURCE_EVIDENCE")],
    "SERPINE1": [("35218720", "EXCLUDED_RECIPIENT_TARGET", "hUCMSC exosomal miR-148a-3p targets recipient SERPINE1; no SERPINE1 cargo measurement", "EV_ONLY", "NO_PROTEIN_SOURCE_EVIDENCE")],
    "TIMP1": [("39472581", "INCLUDED_HUCMSC_EV_PROTEIN", "Proteomics found TIMP1 highly expressed in both adipose- and umbilical-cord-MSC EVs; TIMP1 reproduced EV effects", "EV_ONLY", "HUCMSC_SECRETOME_DIRECT")],
    "TIMP2": [("31182988", "INCLUDED_ENGINEERED_HUCMSC_EV", "Exosomes were derived from TIMP2-overexpressing hUC-MSCs, but the abstract attributes part of the effect to exosomal SFRP2 and does not quantify endogenous TIMP2 cargo", "ENGINEERED_EV", "HUCMSC_PROTEIN_CELL_ONLY")],
    "WNT5A": [("39565502", "INCLUDED_HUCMSC_EV_PROTEIN", "Proteomics and western blot identify enrichment of WNT5A in term hUCMSC exosomes", "EV_ONLY", "HUCMSC_SECRETOME_DIRECT")],
}


PE_STUDIES = {
    "ADAM17": [("22018416", "PE_MECHANISTIC_AGGRAVATING", "Human placental trophoblast ADAM17 regulates TNF-alpha production; mechanistic placental inflammation evidence, not MSC treatment"), ("37374349", "PE_ASSOCIATION_ONLY", "Placental ADAM17 immunostaining association in PE")],
    "COL18A1": [("17616861", "PE_ASSOCIATION_ONLY", "Placental angiogenesis-gene expression study includes COL18A1/endostatin; observational")],
    "DCN": [("26554635", "PE_MECHANISTIC_AGGRAVATING", "Decorin restrains human trophoblast invasion and is linked mechanistically to PE"), ("21659473", "PE_MECHANISTIC_AGGRAVATING", "Decorin is a VEGFR2-binding antagonist in human extravillous trophoblast")],
    "FURIN": [("23598405", "PE_MECHANISTIC_PROTECTIVE", "Furin is required for human trophoblast syncytialization; loss would be adverse"), ("34289413", "PE_ASSOCIATION_ONLY", "Soluble furin measured in PE/FGR serum; observational")],
    "GDF11": [("41078302", "PE_MECHANISTIC_PROTECTIVE", "GDF11 treatment promotes primary human EVT/HTR8 invasion; RNA-seq identifies ANGPTL4 downstream; serum GDF11 lower in PE"), ("35339026", "PE_MECHANISTIC_AGGRAVATING", "FST deficiency raises trophoblast GDF11 and impairs trophoblast functions through SMAD2/3; compartment-specific results conflict")],
    "MDK": [("42100812", "PE_MECHANISTIC_PROTECTIVE", "A 2026 primary study reports impaired Midkine/Treg signaling contributing to PE; single-study evidence")],
    "NAMPT": [("31610400", "PE_ASSOCIATION_ONLY", "NAMPT inversely associates with NO and positively with sFlt-1 in PE; not an intervention"), ("33944612", "PE_ASSOCIATION_ONLY", "NAMPT genotype/levels associate with NO, sFlt-1 and treatment response")],
    "SERPINE1": [("35703881", "PE_ASSOCIATION_ONLY", "Hypoxia-associated placental fibrosis-gene regulation includes SERPINE1; no ligand intervention"), ("40424676", "PE_MECHANISTIC_PROTECTIVE", "Inflammation-driven miRNA destabilizes SERPINE1 and impairs trophoblast proliferation/invasion, suggesting a context-specific protective cellular role")],
    "TIMP1": [("27256632", "PE_ASSOCIATION_ONLY", "MMP/TIMP imbalance in PE and gestational trophoblastic disease; observational")],
    "TIMP2": [("27256632", "PE_ASSOCIATION_ONLY", "MMP/TIMP imbalance in PE and gestational trophoblastic disease; observational")],
    "WNT5A": [("26865089", "PE_MECHANISTIC_AGGRAVATING", "Recombinant WNT5A inhibits HTR8 trophoblast invasion and is increased in PE placenta"), ("33790866", "PE_MECHANISTIC_AGGRAVATING", "Aspirin benefit in an LPS PE mouse model accompanies inhibition of WNT5A/NF-kB"), ("30177057", "PE_MECHANISTIC_PROTECTIVE", "Impaired EVT WNT5A signaling is linked to poor placentation, creating context conflict")],
}


PE_CLASS = {
    "ADAM17": "PE_MECHANISTIC_AGGRAVATING", "AGRN": "NO_DIRECT_PE_EVIDENCE", "COL18A1": "PE_ASSOCIATION_ONLY",
    "DCN": "PE_MECHANISTIC_AGGRAVATING", "ENPP1": "NO_DIRECT_PE_EVIDENCE", "FURIN": "PE_MECHANISTIC_PROTECTIVE",
    "GDF11": "PE_CONTEXT_CONFLICTING", "GRN": "NO_DIRECT_PE_EVIDENCE", "HSPG2": "NO_DIRECT_PE_EVIDENCE",
    "MDK": "PE_MECHANISTIC_PROTECTIVE", "NAMPT": "PE_ASSOCIATION_ONLY", "NID1": "NO_DIRECT_PE_EVIDENCE",
    "PSEN1": "NO_DIRECT_PE_EVIDENCE", "SERPINE1": "PE_CONTEXT_CONFLICTING", "TIMP1": "PE_ASSOCIATION_ONLY",
    "TIMP2": "PE_ASSOCIATION_ONLY", "WNT5A": "PE_CONTEXT_CONFLICTING",
}


RELATED_MSC = {
    "ADAM17": ("RELATED_MSC_MECHANISM", "hUCMSC EV miRNA targets ADAM17 in non-PE fibrosis; not ADAM17 cargo and not PE"),
    "DCN": ("RELATED_MSC_MECHANISM", "ADSC secretome DCN and engineered hUCMSC-DCN have been studied in fibrosis, not PE"),
    "GDF11": ("RELATED_MSC_MECHANISM", "UC-MSC EVs activate endogenous GDF11 downstream in a metabolic model; GDF11 is not demonstrated cargo"),
    "HSPG2": ("RELATED_MSC_MECHANISM", "MSC exosome studies implicate recipient HSPG2 signaling in scar biology; not hUCMSC HSPG2 cargo or PE"),
    "SERPINE1": ("RELATED_MSC_MECHANISM", "hUCMSC EV miRNA targets recipient SERPINE1 in vascular remodeling; not cargo or PE"),
    "TIMP1": ("RELATED_MSC_MECHANISM", "hUCMSC EV TIMP1/Notch1 mechanism shown in photoaging, not PE"),
    "TIMP2": ("RELATED_MSC_MECHANISM", "TIMP2-engineered hUCMSC exosomes studied in myocardial infarction, not PE"),
    "WNT5A": ("RELATED_MSC_MECHANISM", "hUCMSC exosomal WNT5A mechanism shown in fetal lung development, not PE"),
}


PERTURBATION = {
    ("DCN", "PROGRAM_MODULE_01"): ("EMPIRICAL_DISEASE_CONCORDANT", "23460644", "Decorin increases IFN-gamma stability, STAT1 activation and CXCL10 in endothelial cells/fibroblasts", "OTHER_HUMAN_CELLS_IFNG_NOT_TYPE_I"),
    ("DCN", "PROGRAM_MODULE_04"): ("EMPIRICAL_DISEASE_CONCORDANT", "23460644", "Decorin potentiates IFN-gamma/TNF downstream effector cytokines", "OTHER_HUMAN_CELLS_IFNG_NOT_TYPE_I"),
    ("ENPP1", "PROGRAM_MODULE_01"): ("EMPIRICAL_REVERSAL_SUPPORTED", "37333273", "ENPP1 overexpression/hydrolysis suppresses extracellular cGAMP-STING innate immune signaling with scRNA evidence", "TUMOR_STROMAL_AND_IMMUNE_LOW_PLACENTAL_CONTEXT"),
    ("ENPP1", "PROGRAM_MODULE_04"): ("EMPIRICAL_REVERSAL_SUPPORTED", "37333273", "ENPP1 suppresses paracrine cGAMP-STING immune responses", "TUMOR_STROMAL_AND_IMMUNE_LOW_PLACENTAL_CONTEXT"),
}

CONTEXT_HINTS = {
    "ADAM17": ("22018416", "DISEASE_CONCORDANT_NON_MODULE_SPECIFIC", "Trophoblast TNF-alpha shedding/inflammation"),
    "GDF11": ("30407878", "REVERSAL_HINT_NON_MODULE_SPECIFIC", "GDF11 antagonizes TNF/NF-kB inflammation in mouse macrophages; type-I-IFN module not measured"),
    "MDK": ("40943439", "DISEASE_CONCORDANT_NON_MODULE_SPECIFIC", "MDK deficiency/knockdown reduces LPS-induced TNF/CXCL8; IFN module not measured"),
    "NAMPT": ("30833708", "DISEASE_CONCORDANT_NON_MODULE_SPECIFIC", "NAMPT aggravates inflammation/atherosclerosis; IFN module not measured"),
    "TIMP1": ("39472581", "REVERSAL_HINT_NON_MODULE_SPECIFIC", "hUCMSC EV/TIMP1 reduces UVB inflammation in fibroblast/keratinocyte models; IFN module not measured"),
    "WNT5A": ("26865089", "DISEASE_CONCORDANT_NON_MODULE_SPECIFIC", "WNT5A impairs trophoblast invasion; macrophage IFN-alpha/beta module not measured"),
}


def pubmed_lookup() -> dict[str, dict]:
    frame = pd.read_csv(ROOT / "data/raw/phase4b/pubmed_candidate_records.csv", dtype={"PMID": str})
    frame["PMID"] = frame.PMID.astype(str).str.replace(".0", "", regex=False)
    frame = frame.sort_values(["PMID", "candidate"]).drop_duplicates("PMID")
    return frame.set_index("PMID").to_dict("index")


def citation(pmid: str, lookup: dict[str, dict]) -> dict:
    row = lookup.get(str(pmid), {})
    return {"PMID": str(pmid), "DOI": row.get("DOI", "UNRESOLVED") if pd.notna(row.get("DOI", "")) else "UNRESOLVED", "year": row.get("year", "UNRESOLVED"), "publication": row.get("title", "UNRESOLVED"), "source_url": row.get("source_url", f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/")}


def main() -> int:
    cfg = json.loads((ROOT / "config/phase4b_analysis.json").read_text(encoding="utf-8"))
    frozen = pd.read_csv(ROOT / "results/04_phase4b/freeze/phase4b_frozen_candidates.csv")
    hierarchy = pd.read_csv(ROOT / "results/04_phase4a/integration/sender_receiver_evidence_matrix.csv")
    phase4a_a = hierarchy[hierarchy.axis_tier.eq("TIER_A_DIRECTIONAL_RESCUE_CANDIDATE")].copy()
    lookup = pubmed_lookup()
    for directory in ["topology", "protein", "perturbation", "disease", "novelty", "integration", "qc", "figures"]:
        (OUT / directory).mkdir(parents=True, exist_ok=True)

    # Protein topology.
    reg = pd.read_csv(ROOT / "data/raw/phase4b/uniprot_record_registry.csv")
    accession = reg.set_index("candidate").uniprot_accession.to_dict()
    top_rows = []
    for gene in CANDIDATES:
        vals = TOPOLOGY[gene]
        acc = accession[gene]
        top_rows.append({"candidate": gene, "uniprot_accession": acc, "reviewed_human_record": "YES", "signal_peptide": vals[0], "transmembrane_domain_n": vals[1], "canonical_extracellular_localization": vals[2], "secreted_soluble_form": vals[3], "ECM_association": vals[4], "shedding_or_cleavage_requirement": vals[5], "extracellular_enzymatic_activity": vals[6], "extracellular_annotation": vals[7], "intracellular_localization": vals[8], "known_soluble_isoforms_or_chains": vals[9], "topology_classification": vals[10], "omnipath_phase4a_P1": "YES", "HPA_crosscheck_status": "URL_RECORDED_NOT_CLASSIFICATION_PRIMARY", "evidence_source": "UniProtKB reviewed human record; frozen Phase4A OmniPath intercell annotation", "source_url": f"{u(gene, acc)}|https://www.proteinatlas.org/search/{gene}|results/04_phase4a/freeze/phase4_sender_scopes.csv"})
    topology = pd.DataFrame(top_rows)
    topology.to_csv(OUT / "topology/candidate_secretion_topology.csv", index=False)

    # Protein-source evidence; every candidate gets at least one row.
    protein_rows = []
    protein_class = {}
    for gene in CANDIDATES:
        entries = PROTEIN_EVIDENCE.get(gene, [])
        if not entries:
            protein_class[gene] = "NO_PROTEIN_SOURCE_EVIDENCE"
            protein_rows.append({"candidate": gene, "protein_source_classification": "NO_PROTEIN_SOURCE_EVIDENCE", "evidence_disposition": "EVIDENCE_ABSENT", "extracellular_compartment": "NOT_DEMONSTRATED", "MSC_source": "human UC-MSC searched separately from other MSC", "candidate_intervention_or_exposure": "NONE", "protein_method": "NONE_FOUND", "protein_secretion_evidence": "NO_DIRECT_CANDIDATE_PROTEIN_MEASUREMENT_FOUND", "study_type": "SYSTEMATIC_SEARCH_SENTINEL", "sample_size": "NA", "directness": "ABSENT_NOT_AGAINST", "limitations": "Absence may reflect under-study or incomplete secretome tables", "PMID": "NOT_FOUND", "DOI": "NOT_FOUND", "year": "NA", "publication": "NO_INCLUDED_PRIMARY_STUDY", "source_url": "results/04_phase4b/novelty/literature_search_log.csv"})
            continue
        classes = [e[4] for e in entries]
        protein_class[gene] = "HUCMSC_SECRETOME_DIRECT" if "HUCMSC_SECRETOME_DIRECT" in classes else ("OTHER_MSC_SECRETOME" if "OTHER_MSC_SECRETOME" in classes else ("HUCMSC_PROTEIN_CELL_ONLY" if "HUCMSC_PROTEIN_CELL_ONLY" in classes else "NO_PROTEIN_SOURCE_EVIDENCE"))
        for pmid, disposition, statement, compartment, cls in entries:
            row = {"candidate": gene, "protein_source_classification": cls, "evidence_disposition": disposition, "extracellular_compartment": compartment, "MSC_source": "hUCMSC" if "HUCMSC" in disposition else ("adipose-derived MSC" if "OTHER_MSC" in disposition else "hUCMSC"), "candidate_intervention_or_exposure": "baseline unless ENGINEERED is stated", "protein_method": "proteomics/Western blot where stated in evidence summary", "protein_secretion_evidence": statement, "study_type": "PRIMARY_EXPERIMENTAL", "sample_size": "SEE_PRIMARY_STUDY;often not recoverable from abstract", "directness": "DIRECT_CARGO" if disposition in ["INCLUDED_HUCMSC_EV_PROTEIN", "INCLUDED_OTHER_MSC_SECRETOME"] else "INDIRECT_OR_EXCLUSION", "limitations": "EV evidence is not soluble-CM evidence; engineered expression is not baseline secretion"}
            row.update(citation(pmid, lookup)); protein_rows.append(row)
    protein = pd.DataFrame(protein_rows)
    protein.to_csv(OUT / "protein/hucmsc_protein_source_evidence.csv", index=False)

    # Public proteomics registry: candidate-level detectability without overcalling.
    prot_rows = []
    datasets = [
        ["PMID32967723_CM_PEPTIDOME", "human umbilical-cord MSC conditioned medium; preterm versus term", "3 preterm + 3 term donors", "baseline culture; gestational-source comparison", "CONDITIONED_MEDIUM", "Available DOCX derived tables were structurally extracted; none contains exact candidate symbol/accession/name", "study-reported peptidomics confidence; complete flat identification matrix not exposed", "10.1186/s13287-020-01931-0", "https://pmc.ncbi.nlm.nih.gov/articles/PMC7510303/"],
        ["PXD036694", "umbilical-cord, adipose and bone-marrow MSC conditioned media", "UNRESOLVED biological donor n; proteins found in three replications", "baseline high-glucose HUVEC repair comparison", "CONDITIONED_MEDIUM", "Public PROXI record exposes raw files and proprietary .msf search output; no flat processed identification table", "PSM FDR 1%; TMT-10 plex", "PXD036694", "https://www.iprox.cn/proxi/datasets/PXD036694"],
        ["PMID37671699_TNF_SECRETOME", "TNF-alpha-induced human UC-MSC secretome", "UNRESOLVED_FROM_ABSTRACT", "TNF-alpha licensing", "CONDITIONED_MEDIUM", "Primary LC-MS/MS secretome paper found; reusable processed identification table or PX accession not found", "study-specific; unavailable for candidate-level reanalysis", "10.2217/rme-2023-0085", "https://pubmed.ncbi.nlm.nih.gov/37671699/"],
    ]
    for ds in datasets:
        for gene in CANDIDATES:
            status = "NOT_ASSESSABLE_NO_REUSABLE_FLAT_IDENTIFICATION_TABLE"
            if ds[0] == "PMID32967723_CM_PEPTIDOME":
                status = "NOT_DETECTED_IN_AVAILABLE_DERIVED_SUPPLEMENT_TABLES"
            prot_rows.append({"dataset": ds[0], "candidate": gene, "biological_source": ds[1], "donor_n": ds[2], "culture_condition": ds[3], "fraction": ds[4], "candidate_detected": status, "protein_FDR_or_confidence": ds[6], "quantitative_abundance": "NOT_RECOVERABLE", "processed_table_availability": ds[5], "repository_or_accession": ds[7], "interpretation": "NOT_DETECTED_IN_LIMITED_TABLES_IS_NOT_EVIDENCE_OF_ABSENCE" if "NOT_DETECTED" in status else "NOT_ASSESSABLE", "source_url": ds[8]})
    pd.DataFrame(prot_rows).to_csv(OUT / "protein/public_secretome_proteomics_registry.csv", index=False)

    # Perturbation: exactly the 38 Tier-A axes; generic Phase4A prior is not reused.
    pert_rows = []
    for axis in phase4a_a.sort_values(["ligand", "program_module"]).itertuples(index=False):
        key = (axis.ligand, axis.program_module)
        cls, pmid, statement, context = PERTURBATION.get(key, ("NO_EMPIRICAL_SIGNED_EVIDENCE", "NOT_FOUND", "No independent empirical study measured the frozen receiver program under candidate perturbation", "NOT_EVALUABLE"))
        hint = CONTEXT_HINTS.get(axis.ligand, ("", "NONE", ""))
        row = {"candidate": axis.ligand, "program_module": axis.program_module, "receiver_celltype": axis.receiver_celltype, "receiver_level": axis.receiver_level, "phase4a_internal_label": "TIER_A_DIRECTIONAL_RESCUE_CANDIDATE", "manuscript_facing_label": cfg["manuscript_facing_phase4a_label"], "phase4a_generic_signed_prior": axis.signed_reversal_class, "empirical_signed_classification": cls, "candidate_intervention_or_exposure": statement, "direction_of_perturbation": "OPPOSES_FROZEN_PE_PROGRAM" if cls == "EMPIRICAL_REVERSAL_SUPPORTED" else ("REINFORCES_FROZEN_PE_PROGRAM" if cls == "EMPIRICAL_DISEASE_CONCORDANT" else "NOT_ESTIMABLE"), "receiver_context": context, "molecular_readout": "program-relevant interferon/innate-immune readout" if cls != "NO_EMPIRICAL_SIGNED_EVIDENCE" else "FROZEN_MODULE_NOT_DIRECTLY_MEASURED", "receiver_phenotype": "context-specific; see evidence statement", "species": "human plus model systems" if pmid != "NOT_FOUND" else "NA", "study_type": "PRIMARY_PERTURBATION" if pmid != "NOT_FOUND" else "SYSTEMATIC_SEARCH_SENTINEL", "sample_size": "SEE_PRIMARY_STUDY", "directness": "INDEPENDENT_FROM_PHASE4A" if pmid != "NOT_FOUND" else "ABSENT_NOT_AGAINST", "contextual_nonmodule_PMID": hint[0], "contextual_direction_hint": hint[1], "contextual_hint_summary": hint[2], "limitations": "Non-placental and non-matching interferon contexts have low transferability; hints cannot upgrade the axis"}
        row.update(citation(pmid, lookup) if pmid != "NOT_FOUND" else {"PMID": "NOT_FOUND", "DOI": "NOT_FOUND", "year": "NA", "publication": "NO_INCLUDED_MODULE_MATCHED_PERTURBATION", "source_url": "results/04_phase4b/novelty/literature_search_log.csv"})
        pert_rows.append(row)
    perturb = pd.DataFrame(pert_rows)
    perturb.to_csv(OUT / "perturbation/empirical_signed_perturbation_evidence.csv", index=False)

    # PE disease-context evidence.
    disease_rows = []
    for gene in CANDIDATES:
        studies = PE_STUDIES.get(gene, [])
        if not studies:
            disease_rows.append({"candidate": gene, "candidate_PE_classification": "NO_DIRECT_PE_EVIDENCE", "study_evidence_class": "NO_DIRECT_PE_EVIDENCE", "species": "NA", "cell_or_tissue": "NA", "candidate_intervention_or_exposure": "NONE", "direction_of_perturbation": "NOT_ESTIMABLE", "receiver_phenotype": "NA", "molecular_readout": "NA", "PE_relevance": "NO_INCLUDED_DIRECT_RECORD", "study_type": "SYSTEMATIC_SEARCH_SENTINEL", "sample_size": "NA", "directness": "EVIDENCE_ABSENT_NOT_AGAINST", "limitations": "Absence of direct literature is not evidence against the candidate", "PMID": "NOT_FOUND", "DOI": "NOT_FOUND", "year": "NA", "publication": "NO_INCLUDED_PRIMARY_STUDY", "source_url": "results/04_phase4b/novelty/literature_search_log.csv"})
        for pmid, evidence_class, summary in studies:
            row = {"candidate": gene, "candidate_PE_classification": PE_CLASS[gene], "study_evidence_class": evidence_class, "species": "human" if pmid not in ["33790866"] else "mouse and human cells", "cell_or_tissue": "placenta/trophoblast/clinical pregnancy as stated", "candidate_intervention_or_exposure": summary, "direction_of_perturbation": "PROTECTIVE_CONTEXT" if evidence_class == "PE_MECHANISTIC_PROTECTIVE" else ("AGGRAVATING_CONTEXT" if evidence_class == "PE_MECHANISTIC_AGGRAVATING" else "OBSERVATIONAL_ONLY"), "receiver_phenotype": "trophoblast invasion/inflammation/clinical biomarkers as stated", "molecular_readout": "study-specific", "PE_relevance": "DIRECT_PE_OR_PLACENTAL_MECHANISM", "study_type": "PRIMARY_EXPERIMENTAL" if "MECHANISTIC" in evidence_class else "PRIMARY_OBSERVATIONAL", "sample_size": "SEE_PRIMARY_STUDY", "directness": "MECHANISTIC" if "MECHANISTIC" in evidence_class else "ASSOCIATION_NOT_CAUSATION", "limitations": "Not an hUCMSC treatment experiment; compartment and disease-subtype transferability limited"}
            row.update(citation(pmid, lookup)); disease_rows.append(row)
    disease = pd.DataFrame(disease_rows)
    disease.to_csv(OUT / "disease/pe_candidate_context_evidence.csv", index=False)

    # Novelty audit, one row per frozen candidate.
    novelty_rows = []
    for gene in CANDIDATES:
        cls, summary = RELATED_MSC.get(gene, ("NO_CLOSE_PRIOR_MECHANISM_FOUND", "No direct or near-exact candidate-MSC-PE mechanism found in the equal-depth searches"))
        novelty_rows.append({"candidate": gene, "novelty_classification": cls, "direct_hucmsc_PE_mechanism": "NO", "direct_MSC_PE_mechanism": "NO", "related_MSC_mechanism": "YES" if cls == "RELATED_MSC_MECHANISM" else "NO", "mechanism_summary": summary, "PE_publication_overlap": "NONE_FOUND", "publication_count_used_for_ranking": "NO", "interpretation": "NOVELTY_CLASS_NOT_EFFICACY_SCORE", "limitations": "Search limited to indexed public records through audit date; absence is not proof of novelty", "source_url": "results/04_phase4b/novelty/literature_search_log.csv"})
    novelty = pd.DataFrame(novelty_rows)
    novelty.to_csv(OUT / "novelty/direct_msc_pe_overlap.csv", index=False)

    # Reproducible search log with post-screen included counts.
    search = pd.read_csv(ROOT / "data/raw/phase4b/pubmed_search_log_raw.csv")
    included_map = {}
    for gene, entries in PROTEIN_EVIDENCE.items(): included_map[(gene, "HUCMSC_PROTEIN_SOURCE")] = len(entries)
    for gene, entries in PE_STUDIES.items(): included_map[(gene, "PE_DIRECT")] = len(entries)
    for gene in CANDIDATES:
        included_map[(gene, "MSC_PE_NOVELTY")] = 0
        included_map[(gene, "EMPIRICAL_PERTURBATION")] = len({x[1] for k, x in PERTURBATION.items() if k[0] == gene})
        included_map[(gene, "MSC_RELATED_MECHANISM")] = 1 if gene in RELATED_MSC else 0
        included_map[(gene, "PLACENTA_CONTEXT")] = len(PE_STUDIES.get(gene, []))
    search["included_records"] = [included_map.get((r.candidate, r.query_family), 0) for r in search.itertuples(index=False)]
    search["inclusion_exclusion_reason"] = "Primary experimental/observational records included only when candidate-specific; reviews excluded as sole evidence; false-positive target/cargo records retained in evidence tables"
    pride = pd.read_csv(ROOT / "data/raw/phase4b/pride_search_log_raw.csv")
    extra = []
    for row in pride.itertuples(index=False):
        extra.append({"candidate": "ALL_17_EQUAL_DEPTH", "database_source": row.database_source, "query_family": "PUBLIC_PROTEOMICS_DISCOVERY", "query": row.query, "search_date": row.search_date, "result_count": row.result_count, "screened_records": row.result_count, "included_records": 0, "inclusion_exclusion_reason": "No PRIDE keyword hits; a separate ProteomeXchange/iProX record was found from primary-paper data availability", "search_url": row.source_url, "raw_search_file": row.raw_file, "raw_records_file": "", "source_url": row.source_url})
    extra.extend([
        {"candidate": "ALL_17_EQUAL_DEPTH", "database_source": "PubMed_NCBI_EUTILS", "query_family": "GENERAL_HUCMSC_SECRETOME_PROTEOMICS", "query": "(umbilical cord mesenchymal OR Wharton jelly mesenchymal) AND (secretome OR conditioned medium) AND (proteomics OR mass spectrometry)", "search_date": DATE, "result_count": 10, "screened_records": 10, "included_records": 6, "inclusion_exclusion_reason": "Primary hUCMSC/UCMSC CM or secretome proteomics/peptidomics studies retained; EV-only and non-proteomic records separated", "search_url": "https://pubmed.ncbi.nlm.nih.gov/?term=%28umbilical+cord+mesenchymal+OR+Wharton+jelly+mesenchymal%29+AND+%28secretome+OR+conditioned+medium%29+AND+%28proteomics+OR+mass+spectrometry%29", "raw_search_file": "NOT_SAVED_AD_HOC_QUERY;exact query recorded", "raw_records_file": "data/raw/phase4b/pubmed_candidate_records.csv", "source_url": "https://pubmed.ncbi.nlm.nih.gov/"},
        {"candidate": "ALL_17_EQUAL_DEPTH", "database_source": "ProteomeXchange_iProX_PROXI", "query_family": "PUBLIC_PROTEOMICS_EXACT_ACCESSION", "query": "PXD036694", "search_date": DATE, "result_count": 1, "screened_records": 1, "included_records": 1, "inclusion_exclusion_reason": "Relevant UCMSC/ADSC/BMMSC conditioned-medium TMT dataset; no public flat processed identification table", "search_url": "https://www.iprox.cn/proxi/datasets/PXD036694", "raw_search_file": "NOT_SAVED_API_RESPONSE;URL versioned in registry", "raw_records_file": "", "source_url": "https://www.iprox.cn/proxi/datasets/PXD036694"},
        {"candidate": "ALL_17_EQUAL_DEPTH", "database_source": "Google_Scholar_like", "query_family": "AVAILABILITY_NOTE", "query": "Candidate-identical query families", "search_date": DATE, "result_count": "NOT_AVAILABLE_PROGRAMMATICALLY", "screened_records": 0, "included_records": 0, "inclusion_exclusion_reason": "No reproducible programmatic Google Scholar interface available; PubMed plus repository APIs used", "search_url": "NOT_USED", "raw_search_file": "", "raw_records_file": "", "source_url": "config/phase4b_analysis.json"},
    ])
    search = pd.concat([search, pd.DataFrame(extra)], ignore_index=True)
    search.to_csv(OUT / "novelty/literature_search_log.csv", index=False)

    # Mixed-direction stress test.
    frozen_map = frozen.set_index("ligand").to_dict("index")
    mixed_rows = []
    for gene in CANDIDATES:
        f = frozen_map[gene]
        empirical_rev = perturb[(perturb.candidate.eq(gene)) & perturb.empirical_signed_classification.eq("EMPIRICAL_REVERSAL_SUPPORTED")].program_module.tolist()
        empirical_conc = perturb[(perturb.candidate.eq(gene)) & perturb.empirical_signed_classification.eq("EMPIRICAL_DISEASE_CONCORDANT")].program_module.tolist()
        if gene in ["ADAM17", "NAMPT", "WNT5A"]:
            mixed_class = "PREDOMINANTLY_DISEASE_CONCORDANT"
        elif f["mixed_signed_direction_across_modules"] == "YES" or empirical_conc:
            mixed_class = "MIXED_DIRECTION_CONTEXT_DEPENDENT"
        else:
            mixed_class = "SIGNED_EVIDENCE_UNRESOLVED"
        mixed_rows.append({"candidate": gene, "phase4a_reversal_supported_modules": f["tier_A_modules"], "phase4a_disease_concordant_modules": f["disease_concordant_modules"], "phase4a_mixed_signed_direction": f["mixed_signed_direction_across_modules"], "receiver_evidence_strength": "R1/R2A frozen modules; see phase4_receiver_hierarchy.csv", "empirical_reversal_modules": ";".join(empirical_rev), "empirical_disease_concordant_modules": ";".join(empirical_conc), "PE_context_classification": PE_CLASS[gene], "apparent_benefit_receiver_context_dependent": "YES" if mixed_class in ["MIXED_DIRECTION_CONTEXT_DEPENDENT", "PREDOMINANTLY_DISEASE_CONCORDANT"] else "UNRESOLVED", "mixed_direction_classification": mixed_class, "interpretation": "GENERIC_SIGNED_PRIOR_AND_EXTERNAL_EVIDENCE_ARE_CONTEXT_DEPENDENT_NOT_THERAPEUTIC_PROOF", "source_url": "results/04_phase4a/integration/phase4a_candidate_hierarchy.csv|results/04_phase4b/perturbation/empirical_signed_perturbation_evidence.csv|results/04_phase4b/disease/pe_candidate_context_evidence.csv"})
    mixed = pd.DataFrame(mixed_rows)
    mixed.to_csv(OUT / "integration/mixed_direction_stress_test.csv", index=False)

    # Separate evidence dimensions and non-exclusive deterministic flags.
    evidence_rows, class_rows = [], []
    topo_map = topology.set_index("candidate").topology_classification.to_dict()
    novel_map = novelty.set_index("candidate").novelty_classification.to_dict()
    mixed_map = mixed.set_index("candidate").mixed_direction_classification.to_dict()
    plausible = set(cfg["topology_plausible_for_primary_paracrine"])
    for gene in CANDIDATES:
        f = frozen_map[gene]
        p = perturb[perturb.candidate.eq(gene)]
        empirical_rev_n = int(p.empirical_signed_classification.eq("EMPIRICAL_REVERSAL_SUPPORTED").sum())
        empirical_conc_n = int(p.empirical_signed_classification.isin(["EMPIRICAL_DISEASE_CONCORDANT", "EMPIRICAL_SIGN_CONFLICT"]).sum())
        topo_ok = topo_map[gene] in plausible
        protein_direct = protein_class[gene] == "HUCMSC_SECRETOME_DIRECT"
        direct_prior = novel_map[gene] in ["DIRECT_HUCMSC_PE_MECHANISM_ALREADY_SHOWN", "DIRECT_MSC_PE_MECHANISM_ALREADY_SHOWN"]
        high = topo_ok and protein_direct and empirical_rev_n > 0 and empirical_conc_n == 0 and not direct_prior
        context = topo_ok and protein_class[gene] in ["HUCMSC_SECRETOME_DIRECT", "OTHER_MSC_SECRETOME"] and empirical_rev_n > 0 and mixed_map[gene] in ["MIXED_DIRECTION_CONTEXT_DEPENDENT", "PREDOMINANTLY_DISEASE_CONCORDANT"]
        computational = protein_class[gene] in ["NO_PROTEIN_SOURCE_EVIDENCE", "HUCMSC_PROTEIN_CELL_ONLY"] and empirical_rev_n == 0
        weak = topo_map[gene] in ["MEMBRANE_ASSOCIATED", "INTRACELLULAR_OR_QUESTIONABLE_PARACRINE", "UNCERTAIN"]
        low_novelty = direct_prior
        flags = []
        if high: flags.append("TRIANGULATED_HIGH_PRIORITY")
        if context: flags.append("TRIANGULATED_CONTEXT_DEPENDENT")
        if computational: flags.append("COMPUTATIONAL_ONLY")
        if weak: flags.append("BIOPHYSICALLY_WEAK_PARACRINE")
        if low_novelty: flags.append("KNOWN_MECHANISM_LOW_NOVELTY")
        evidence_rows.append({"candidate": gene, "phase4a_internal_label": f["phase4a_internal_label"], "manuscript_facing_label": f["manuscript_facing_label"], "S1_S2": f["sender_evidence_level"], "licensing_class": f["licensing_context"], "P1_annotation": f["paracrine_scope"], "tier_A_modules": f["tier_A_modules"], "tier_B_modules": f["tier_B_modules"], "disease_concordant_modules": f["disease_concordant_modules"], "mixed_signed_direction_across_modules": f["mixed_signed_direction_across_modules"], "PHASE4A_COMPUTATIONAL_COMPATIBILITY": "TIER_A_FROZEN", "GENERIC_SIGNED_PRIOR": "ONE_GENERIC_OMNIPATH_COLLECTRI_LAYER", "PARACRINE_TOPOLOGY": topo_map[gene], "HUCMSC_PROTEIN_SOURCE": protein_class[gene], "EMPIRICAL_SIGNED_PERTURBATION": "REVERSAL_AXES=" + str(empirical_rev_n) + ";CONCORDANT_AXES=" + str(empirical_conc_n), "PE_DISEASE_CONTEXT": PE_CLASS[gene], "DIRECT_MSC_PE_PRECEDENT": novel_map[gene], "MIXED_DIRECTION_RISK": mixed_map[gene], "RECEIVER_EVIDENCE_STRENGTH": "R1/R2A;cell-localization limitations inherited", "NOVELTY": "EXACT_MSC_PE_NOT_FOUND" if not direct_prior else "LOW", "classification_flags": ";".join(flags) if flags else "NO_SUGGESTED_CATEGORY_THRESHOLD_MET", "source_url": "results/04_phase4b/freeze/phase4b_frozen_candidates.csv|results/04_phase4b/domain_tables"})
        if high:
            primary = "TRIANGULATED_HIGH_PRIORITY"
        elif context:
            primary = "TRIANGULATED_CONTEXT_DEPENDENT"
        elif weak:
            primary = "BIOPHYSICALLY_WEAK_PARACRINE"
        elif computational:
            primary = "COMPUTATIONAL_ONLY"
        elif protein_direct and empirical_rev_n == 0:
            primary = "PROTEIN_SUPPORTED_BUT_DIRECTION_UNRESOLVED"
        elif empirical_rev_n > 0 and not protein_direct:
            primary = "PERTURBATION_SUPPORTED_BUT_SOURCE_UNCONFIRMED"
        else:
            primary = "PARTIAL_EXTERNAL_EVIDENCE"
        class_rows.append({"candidate": gene, "primary_classification": primary, "TRIANGULATED_HIGH_PRIORITY": "YES" if high else "NO", "TRIANGULATED_CONTEXT_DEPENDENT": "YES" if context else "NO", "COMPUTATIONAL_ONLY": "YES" if computational else "NO", "BIOPHYSICALLY_WEAK_PARACRINE": "YES" if weak else "NO", "KNOWN_MECHANISM_LOW_NOVELTY": "YES" if low_novelty else "NO", "topology_pass": "YES" if topo_ok else "NO", "direct_hucmsc_extracellular_protein": "YES" if protein_direct else "NO", "empirical_reversal_axis_n": empirical_rev_n, "empirical_concordant_or_conflict_axis_n": empirical_conc_n, "final_interpretation": "NO_CANDIDATE_IS_THERAPEUTICALLY_VALIDATED", "classification_rule_source": "config/phase4b_analysis.json", "source_url": "results/04_phase4b/integration/phase4b_candidate_evidence_matrix.csv"})
    evidence_matrix = pd.DataFrame(evidence_rows)
    classifications = pd.DataFrame(class_rows)
    evidence_matrix.to_csv(OUT / "integration/phase4b_candidate_evidence_matrix.csv", index=False)
    classifications.to_csv(OUT / "integration/phase4b_candidate_classification.csv", index=False)

    risks = pd.DataFrame([
        ["P4B_RISK_001", "CRITICAL", "EMPIRICAL_DIRECTION", "Only ENPP1 has independent module-matched reversal evidence, in a non-placental tumor context; transferability is low.", "Do not call any candidate empirically validated for PE receiver reversal", "OPEN"],
        ["P4B_RISK_002", "HIGH", "PROTEIN_SOURCE", "Direct hUCMSC candidate protein evidence is EV-only for TIMP1 and WNT5A; soluble conditioned-medium availability is unproven.", "Keep EV and soluble evidence separate", "OPEN"],
        ["P4B_RISK_003", "HIGH", "MIXED_DIRECTION", "Phase4A and external evidence show substantial receiver/context sign conflict; ADAM17, NAMPT and WNT5A lean disease-concordant.", "Retain negative/concordant evidence and receiver-specific labels", "OPEN"],
        ["P4B_RISK_004", "HIGH", "TOPOLOGY", "ADAM17 and PSEN1 are membrane proteins; NAMPT lacks a signal peptide and has context-dependent non-classical secretion.", "Flag biophysically weak paracrine interpretations", "CONTROLLED"],
        ["P4B_RISK_005", "MEDIUM", "PUBLIC_PROTEOMICS", "PXD036694 lacks a public flat identification table and exposes proprietary MSF output; candidate detection is not assessable.", "Do not treat not-assessable as non-detection", "OPEN"],
        ["P4B_RISK_006", "MEDIUM", "PUBLICATION_BIAS", "Well-studied candidates generate more records; sparse evidence may indicate under-study.", "No publication-count ranking; evidence absent separated from evidence against", "CONTROLLED"],
        ["P4B_RISK_007", "MEDIUM", "NOVELTY", "No exact MSC-candidate-PE mechanism was found, but indexing and full-text terminology can miss near-exact precedent.", "Treat novelty as search-bounded, not proven", "OPEN"],
        ["P4B_RISK_008", "CONTROLLED", "MODULE10", "Module10 had zero target-compatible Phase4A axes.", "External mitochondrial literature cannot create an axis", "CONTROLLED"],
        ["P4B_RISK_009", "CONTROLLED", "PHASE4A_FREEZE", "All Phase4A hashes and 17 identities verified before review.", "No hierarchy modification", "CONTROLLED"],
    ], columns=["risk_id", "severity", "domain", "risk", "mitigation", "status"])
    risks["source_url"] = "config/phase4b_analysis.json|results/04_phase4b/domain_tables"
    risks.to_csv(OUT / "qc/phase4b_risk_flags.csv", index=False)
    print(f"PHASE4B_OUTPUTS_OK candidates={len(classifications)} high={classifications.TRIANGULATED_HIGH_PRIORITY.eq('YES').sum()} context={classifications.TRIANGULATED_CONTEXT_DEPENDENT.eq('YES').sum()} computational={classifications.COMPUTATIONAL_ONLY.eq('YES').sum()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
