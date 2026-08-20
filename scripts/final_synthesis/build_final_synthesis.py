#!/usr/bin/env python3
"""Build the frozen final-synthesis tables from existing phase outputs only."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


PHASE4B17 = [
    "ADAM17", "AGRN", "COL18A1", "DCN", "ENPP1", "FURIN", "GDF11", "GRN",
    "HSPG2", "MDK", "NAMPT", "NID1", "PSEN1", "SERPINE1", "TIMP1", "TIMP2", "WNT5A",
]


def read_csv(root: Path, rel: str) -> pd.DataFrame:
    return pd.read_csv(root / rel, dtype=str, keep_default_na=False)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, encoding="utf-8-sig")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def join_nonempty(values) -> str:
    return ";".join(sorted({str(v).strip() for v in values if str(v).strip()}))


def build_phase_registry() -> pd.DataFrame:
    rows = [
        ("0A", "Initial public-data feasibility and overlap audit", "13 GEO candidates including PE bulk/scRNA and hUC-MSC datasets", "dataset/biological subject", "Programmatic GEO/SRA/publication metadata audit", "Traceable evidence; no guessed metadata; no overlap counted as independent", "GO_WITH_MODIFICATIONS", "Repository and auditable registries established; GSE75010 overlap and GSE234729 discrepancies exposed", "Candidate PE scRNA search was incomplete and was superseded by Phase 0B", "FROZEN_HISTORICAL"),
        ("0B", "Rescue PE scRNA availability and audit sender redundancy", "Admati resources, GSE290578, Yang/Jiao/Zhang candidates, normal atlas, GSE182158/GSE199071/GSE117837", "independent pregnancy/donor", "Systematic repository and donor-identity audit", "PRIMARY_PE_SCRNA requires accessible expression plus recoverable donor identity and valid comparison", "GO_TO_PHASE1_WITH_RESTRICTIONS", "Admati public matrix with 26 recoverable donors selected; sender redundancy documented", "The nominal 78-subject preprint cohort was not reusable as the primary public analysis object", "ACTIVE_FOUNDATION"),
        ("1A", "Freeze independent placental bulk cohorts and statistical design", "GSE75010_BIOBANK, GSE30186, GSE10588, GSE24129, GSE25906, GSE43942", "pregnancy/sample", "Cohort-specific preprocessing/QC and overlap freeze", "Independent cohorts analyzed separately; no cross-study pooling or ComBat", "GO_TO_PHASE1B_WITH_RESTRICTIONS", "Six independent core cohorts frozen for cohort-wise synthesis", "Clinical covariates and gestational-age matching were incomplete and heterogeneous", "ACTIVE_FOUNDATION"),
        ("1A.1", "Correct preprocessing and freeze comparable primary estimands", "Six Phase 1A core cohorts", "pregnancy/sample", "Raw-versus-processed audit; platform-appropriate reconstruction; formula freeze", "No universal abs(meta log2FC)>=0.25 cutoff; retain all 18 QC flags", "GO_TO_PHASE1B_WITH_RESTRICTIONS", "GSE30186 reconstructed without arbitrary shift-log; defensible matrices/formulas frozen", "Several historical platforms required submitted processed matrices", "ACTIVE_FOUNDATION"),
        ("1B", "Test for a stable gene-level PE signature across cohorts", "Six Phase 1A.1 core cohorts", "pregnancy/sample", "Cohort-wise limma plus REML random-effects meta-analysis and leave-one-cohort-out", "FDR<0.05, >=75% direction, I2<=60%, LOO direction retained and >=80% LOO FDR<0.10", "NO_GO", "17,731 genes estimable in >=4 cohorts; zero meta-FDR significant and zero STABLE genes", "Cross-platform heterogeneity and limited covariate comparability precluded a universal gene signature", "FROZEN_NEGATIVE"),
        ("2A", "Discover patient-level cell-state PE programs", "Admati public single-cell expression layer", "pregnancy within cell type", "Legacy edgeR pseudobulk and cameraPR", "Subtype-stratified EOPE/LOPE; frozen donor eligibility", "GO_TO_PHASE2B_WITH_RESTRICTIONS (historical)", "20 SHARED_PE hypotheses and subtype-specific programs were generated", "Public layer was normalized-to-about-10,000 and ceiled, not unmodified raw UMI counts", "LEGACY_COUNT_MODEL_DISCOVERY"),
        ("2A.1", "Stress-test legacy receiver hypotheses and audit confounding", "Admati plus Yang LOPE placenta", "pregnancy within cell type", "Clinical overlap tests, frozen-program scores, restricted sensitivities, Yang targeted replication", "Do not refit non-estimable confounders or rediscover pathways", "REVISE_RECEIVER_FRAMEWORK", "Delivery/IUGR positivity failures documented; partial Yang directional evidence; expression-layer contradiction confirmed", "EOPE delivery mode and induction remained non-identifiable", "FROZEN_DIAGNOSTIC"),
        ("2A.2", "Correct receiver framework and retest/rediscover programs", "Admati normalized/ceiled public matrix; Zheng audit; preserved Yang evidence", "pregnancy within cell type", "Donor mean normalized expression plus log transform and limma; corrected ranked gene-set testing", "Tier 1 frozen-20 retest separated from exploratory corrected rediscovery", "GO_TO_PHASE2B_WITH_RESTRICTIONS", "All 20 frozen hypotheses retained direction under the corrected model; 11 modules frozen; Zheng 3+2 cohort not reconstructable", "Clinical confounding remained; external localization was sparse", "ACTIVE_FROZEN"),
        ("2B", "Independently validate frozen receiver programs in bulk placenta", "Six Phase 1A.1 core cohorts", "pregnancy/sample", "Within-cohort rank scores, cohort disease models, random-effects program meta-analysis", "Robust support required FDR<0.05, direction agreement, >=75% cohort agreement, I2<=60%, stable LOO direction", "GO_TO_PHASE3_WITH_RESTRICTIONS", "19/19 sets estimable; zero robust constituent sets/modules; nine directional sets and five directional modules", "Bulk supports tissue-level program direction, not cell-type localization", "FROZEN_LIMITED_DIRECTIONAL"),
        ("3", "Define an independent reproducible hUC-MSC sender ligand universe", "GSE182158, GSE199071, GSE117837", "independent donor and donor×passage stratum", "Donor-aware baseline expression plus within-stratum licensing summaries", "Robust baseline support across datasets/donors; licensing kept descriptive", "GO_TO_PHASE4_WITH_RESTRICTIONS", "214 robust S1/S2 sender candidates frozen (38 S1; 176 S2)", "Only two licensing donors and incomplete protein-level secretion evidence", "ACTIVE_FROZEN"),
        ("4A", "Blindly integrate frozen sender and receiver arms", "214 senders; frozen R1/R2A/R2B receiver hierarchy", "ligand-receiver-module axis with pregnancy-level receptor competence", "NicheNet LR/targets plus generic signed OmniPath-CollecTRI layer", "Separate compatibility from signed reversal; deterministic Tier A/B/C", "GO_TO_PHASE4B_WITH_RESTRICTIONS", "1,284 axes; 791 target-compatible, 63 reversal-supported, 246 disease-concordant; 17 generic-signed-prior Tier A candidates", "Signed prior was generic; 12/17 Tier A candidates had mixed direction; Module 10 had zero target-compatible axes", "ACTIVE_FROZEN"),
        ("4B", "Triangulate the frozen 17 with external topology, protein, perturbation, PE and novelty evidence", "Frozen 17 Phase 4A Tier A candidates", "candidate and candidate-module axis", "Curated topology, candidate-specific literature/repository and perturbation audit", "Unchanged deterministic multi-domain classification", "REVISE_CANDIDATE_FRAMEWORK", "Initial audit found no TRIANGULATED_HIGH_PRIORITY candidate", "Candidate-name-first retrieval missed proteins present only in proteomics tables; protein dimension superseded by 4B.1", "FROZEN_SUPERSEDED_PROTEIN_DIMENSION"),
        ("4B.1", "Correct completeness of hUC-MSC protein-source evidence", "Frozen 17; four technically evaluable UC/WJ-MSC extracellular proteomic resources", "candidate×proteomic dataset", "Proteomics-first identifier-mapped audit", "Only protein-source dimension updated; Phase 4B rules unchanged", "GO_TO_FINAL_SYNTHESIS_WITH_RESTRICTIONS", "15/17 candidates detected in hUC-MSC EV material; none confirmed in complete soluble-CM tables; ENPP1 became sole TRIANGULATED_HIGH_PRIORITY candidate", "EV evidence does not establish soluble secretion, placental exposure, direction, or efficacy", "ACTIVE_FROZEN"),
    ]
    cols = ["phase", "scientific_question", "datasets", "biological_replication_unit", "primary_method", "prespecified_rule", "gate", "key_result", "limitation", "result_status"]
    df = pd.DataFrame(rows, columns=cols)
    df["source_url"] = [
        "docs/DATASET_AUDIT_REPORT.md", "docs/DATASET_AUDIT_PHASE0B_REPORT.md",
        "docs/PHASE1A_BULK_DATA_FREEZE_REPORT.md", "docs/PHASE1A1_PREPROCESSING_AMENDMENT.md",
        "docs/PHASE1B_PE_DISEASE_SIGNATURE_REPORT.md", "docs/PHASE2A_PE_CELLSTATE_DISEASE_PROGRAM_REPORT.md",
        "docs/PHASE2A1_RECEIVER_PROGRAM_STRESS_TEST_REPORT.md", "docs/PHASE2A2_RECEIVER_FRAMEWORK_CORRECTION_REPORT.md",
        "docs/PHASE2B_BULK_PROGRAM_VALIDATION_REPORT.md", "docs/PHASE3_HUCMSC_SENDER_PROGRAM_REPORT.md",
        "docs/PHASE4A_SENDER_RECEIVER_INTEGRATION_REPORT.md", "docs/PHASE4B_EXTERNAL_TRIANGULATION_REPORT.md",
        "docs/PHASE4B1_PROTEIN_SOURCE_COMPLETENESS_REPORT.md",
    ]
    return df


def build_receiver(root: Path) -> pd.DataFrame:
    r = read_csv(root, "results/04_phase4a/freeze/phase4_receiver_hierarchy.csv")
    bulk = read_csv(root, "results/02_phase2b/meta/program_module_validation.csv")
    bmap = bulk.set_index("program_module")["module_classification"].to_dict()
    conf = (
        "EOPE: delivery mode/induction/IUGR NON_ESTIMABLE; GA/fetal sex/maternal age WEAKLY_ESTIMABLE with early-control n=3. "
        "LOPE: C-section sensitivity possible; induction/IUGR NON_ESTIMABLE; residual GA/cell-abundance limits."
    )
    external = {}
    for module in r["program_module"]:
        if module in {"PROGRAM_MODULE_06", "PROGRAM_MODULE_07"}:
            external[module] = "PARTIAL_APPROXIMATE_YANG_LOPE_LOCALIZATION"
        elif module == "PROGRAM_MODULE_08":
            external[module] = "YANG_LOPE_DIRECTIONALLY_DISCORDANT"
        else:
            external[module] = "PRIMARY_ADMATI_LOCALIZATION_ONLY"
    out = pd.DataFrame({
        "module": r["program_module"],
        "module_label": r["module_label"],
        "cell_type_origin": r["celltype"],
        "frozen_direction": r["frozen_direction"],
        "constituent_gene_set_n": r["constituent_gene_set_n"],
        "admati_corrected_evidence": r["CORRECTED_ADMATI_SUPPORT"],
        "yang_external_scrna_evidence": r["EXTERNAL_SCRNA_SUPPORT"],
        "bulk_directional_evidence": r["INDEPENDENT_BULK_PROGRAM_SUPPORT"].eq("BULK_MODULE_DIRECTIONAL").map({True:"YES", False:"NO"}),
        "bulk_robust_evidence": "NO",
        "bulk_module_classification": r["program_module"].map(bmap).fillna("NOT_VALIDATED_PHASE2B_HOLD"),
        "clinical_confounding_limitation": conf,
        "external_localization_confidence": r["program_module"].map(external),
        "cell_type_localization_claim": r["celltype_localization_claim"],
        "final_receiver_level": r["receiver_level"],
        "phase4_analysis_scope": r["phase4a_analysis_scope"],
    })
    def interp(row):
        if row.final_receiver_level == "R1":
            return "PRIMARY_RECEIVER; corrected Admati localization plus partial Yang LOPE directional localization and bulk tissue-level directional support"
        if row.final_receiver_level == "R2A":
            return "PRIMARY_RECEIVER_WITH_ADMATI_LOCALIZATION; bulk tissue-level direction where indicated; no independent localization"
        if row.final_receiver_level == "R2B":
            return "SECONDARY_SENSITIVITY_RECEIVER; corrected Admati and partial Yang localization, without bulk support"
        if row.final_receiver_level == "HOLD":
            return "HOLD_EXTERNAL_DISCORDANCE; not promoted"
        return "NON_PRIMARY_RECEIVER; preserved but not promoted"
    out["final_interpretation"] = out.apply(interp, axis=1)
    out["source_url"] = r["source_url"] + "|results/02_phase2b/meta/program_module_validation.csv|results/02_phase2a1/clinical_confounding_audit.csv"
    return out


def build_sender(root: Path) -> pd.DataFrame:
    base = read_csv(root, "results/03_phase3/sender/frozen_phase4_sender_candidates.csv")
    scope = read_csv(root, "results/04_phase4a/freeze/phase4_sender_scopes.csv")
    hier = read_csv(root, "results/04_phase4a/integration/phase4a_candidate_hierarchy.csv")
    top = read_csv(root, "results/04_phase4b/topology/candidate_secretion_topology.csv")
    prot = read_csv(root, "results/04_phase4b1/corrected_hucmsc_protein_source_evidence.csv")
    scope_keep = scope[[
        "gene", "sender_evidence_level", "baseline_classification", "licensing_classification",
        "omnipath_secreted_annotation", "omnipath_extracellular_annotation",
        "P1_PARACRINE_CORE", "P2_EXTRACELLULAR_EXTENDED", "P3_FULL_LR_SENSITIVITY",
        "primary_scope_membership",
    ]]
    x = base.merge(scope_keep, on="gene", suffixes=("_phase3", "_scope"), validate="one_to_one")
    x = x.merge(hier[["ligand", "best_phase4a_tier"]], left_on="gene", right_on="ligand", how="left", validate="one_to_one")
    x = x.merge(top[["candidate", "topology_classification"]], left_on="gene", right_on="candidate", how="left", validate="one_to_one")
    x = x.merge(prot[["candidate", "EV_direct", "soluble_CM_direct", "corrected_protein_source_classification", "direct_datasets"]], left_on="gene", right_on="candidate", how="left", suffixes=("", "_prot"), validate="one_to_one")
    audited = x["gene"].isin(PHASE4B17)
    out = pd.DataFrame({
        "gene": x["gene"],
        "sender_evidence_level": x["sender_evidence_level_phase3"],
        "baseline_datasets": "GSE182158;GSE199071",
        "licensing_dataset": "GSE117837",
        "gse182158_supported_donor_n": x["gse182_supported_donor_n"],
        "gse199071_supported_donor_n": x["gse199_supported_donor_n"],
        "total_supported_donor_n": x["total_supported_donor_n"],
        "baseline_classification": x["baseline_classification_phase3"],
        "licensing_class": x["licensing_classification_phase3"],
        "secreted_annotation": x["omnipath_secreted_annotation_scope"],
        "extracellular_annotation": x["omnipath_extracellular_annotation_scope"],
        "P1_PARACRINE_CORE": x["P1_PARACRINE_CORE"],
        "P2_EXTRACELLULAR_EXTENDED": x["P2_EXTRACELLULAR_EXTENDED"],
        "P3_FULL_LR_SENSITIVITY": x["P3_FULL_LR_SENSITIVITY"],
        "primary_scope": x["primary_scope_membership"],
        "protein_topology": x["topology_classification"].where(audited, "NOT_AUDITED_PHASE4B_FROZEN17_ONLY"),
        "topology_audit_status": audited.map({True:"AUDITED_FROZEN17", False:"NOT_SYSTEMATICALLY_AUDITED"}),
        "huc_wj_msc_ev_protein_evidence": x["EV_direct"].where(audited, "NOT_SYSTEMATICALLY_AUDITED"),
        "soluble_cm_protein_evidence": x["soluble_CM_direct"].where(audited, "NOT_SYSTEMATICALLY_AUDITED"),
        "protein_source_classification": x["corrected_protein_source_classification"].where(audited, "NOT_SYSTEMATICALLY_AUDITED"),
        "protein_evidence_datasets": x["direct_datasets"].where(audited, "NOT_SYSTEMATICALLY_AUDITED"),
        "phase4_eligible": x["phase4_eligible"],
        "phase4a_best_tier": x["best_phase4a_tier"],
    })
    out["manuscript_facing_sender_interpretation"] = audited.map({True:"TRANSCRIPTOMIC_SENDER_WITH_FROZEN_PHASE4B_PROTEIN_AUDIT", False:"TRANSCRIPTOMIC_SENDER; PROTEIN_SECRETION_NOT_SYSTEMATICALLY_TESTED"})
    out["source_url"] = "results/03_phase3/sender/frozen_phase4_sender_candidates.csv|results/04_phase4a/freeze/phase4_sender_scopes.csv|results/04_phase4b1/corrected_hucmsc_protein_source_evidence.csv"
    return out.sort_values(["sender_evidence_level", "gene"]).reset_index(drop=True)


def build_candidates(root: Path) -> pd.DataFrame:
    ev = read_csv(root, "results/04_phase4b1/corrected_phase4b_candidate_evidence_matrix.csv")
    cl = read_csv(root, "results/04_phase4b1/corrected_phase4b_candidate_classification.csv")
    h = read_csv(root, "results/04_phase4a/integration/phase4a_candidate_hierarchy.csv")
    sc = read_csv(root, "results/04_phase4a/freeze/phase4_sender_scopes.csv")
    tp = read_csv(root, "results/04_phase4b/topology/candidate_secretion_topology.csv")
    pr = read_csv(root, "results/04_phase4b1/corrected_hucmsc_protein_source_evidence.csv")
    lr = read_csv(root, "results/04_phase4a/lr/sender_receiver_lr_compatibility.csv")
    pert = read_csv(root, "results/04_phase4b/perturbation/empirical_signed_perturbation_evidence.csv")
    nov = read_csv(root, "results/04_phase4b/novelty/direct_msc_pe_overlap.csv")
    pe = read_csv(root, "results/04_phase4b/disease/pe_candidate_context_evidence.csv")
    x = ev.merge(cl, on="candidate", suffixes=("", "_classification"), validate="one_to_one")
    x = x.merge(h, left_on="candidate", right_on="ligand", suffixes=("_phase4b", "_phase4a"), validate="one_to_one")
    x = x.merge(sc[["gene", "P1_PARACRINE_CORE", "P2_EXTRACELLULAR_EXTENDED", "P3_FULL_LR_SENSITIVITY"]], left_on="candidate", right_on="gene", validate="one_to_one")
    x = x.merge(tp[["candidate", "topology_classification"]], on="candidate", validate="one_to_one")
    x = x.merge(pr[["candidate", "EV_direct", "soluble_CM_direct", "direct_datasets"]], on="candidate", validate="one_to_one")
    x = x.merge(nov[["candidate", "novelty_classification", "mechanism_summary"]], on="candidate", validate="one_to_one")
    pe_agg = pe.groupby("candidate", as_index=False).agg(pe_context=("candidate_PE_classification", join_nonempty))
    x = x.merge(pe_agg, on="candidate", how="left", validate="one_to_one")
    pert_agg = pert.groupby("candidate", as_index=False).agg(
        empirical_perturbation_class=("empirical_signed_classification", join_nonempty),
        perturbation_context=("receiver_context", join_nonempty),
        empirical_modules=("program_module", join_nonempty),
    )
    x = x.merge(pert_agg, on="candidate", how="left", validate="one_to_one")
    lrf = lr[(lr["ligand"].isin(PHASE4B17)) & (lr["nichenet_lr_support"] == "YES")]
    lr_agg = lrf.groupby("ligand", as_index=False).agg(
        lr_supported_edge_n=("receptor", "size"),
        lr_supported_receptor_n=("receptor", "nunique"),
        lr_supported_receptors=("receptor", join_nonempty),
    )
    x = x.merge(lr_agg, left_on="candidate", right_on="ligand", how="left", suffixes=("", "_lragg"), validate="one_to_one")
    x["lr_supported_edge_n"] = x["lr_supported_edge_n"].fillna("0")
    x["lr_supported_receptor_n"] = x["lr_supported_receptor_n"].fillna("0")
    x["lr_supported_receptors"] = x["lr_supported_receptors"].fillna("")
    x["empirical_perturbation_class"] = x["empirical_perturbation_class"].fillna("NO_EMPIRICAL_SIGNED_EVIDENCE")
    x["perturbation_context"] = x["perturbation_context"].fillna("NONE")
    x["empirical_modules"] = x["empirical_modules"].fillna("NONE")
    def interpretation(row):
        if row["candidate"] == "ENPP1":
            return "Leading experimentally testable communication hypothesis from the frozen cross-dataset framework; context-dependent, EV-source and non-placental reversal evidence only; not a validated therapeutic factor"
        if row["corrected_primary_classification"] == "PROTEIN_SUPPORTED_BUT_DIRECTION_UNRESOLVED":
            return "Direct hUC/WJ-MSC EV protein detection supports source plausibility, but independent signed receiver reversal is unresolved"
        if row["corrected_primary_classification"] == "BIOPHYSICALLY_WEAK_PARACRINE":
            return "Computational compatibility retained, but topology weakens a soluble/paracrine interpretation"
        return "Frozen computational candidate without sufficient independent source-and-direction triangulation"
    def experiment(row):
        if row["candidate"] == "ENPP1":
            return "Confirm across independent hUC-MSC donors; separate soluble versus EV-associated ENPP1; quantify licensing response; expose placental stromal/Hofbauer receivers to CM/EV; ENPP1-specific loss/add-back; test IFN-state reversal and only then PE-relevant function"
        return "Verify extracellular compartment in independent hUC-MSC donors, then perform candidate-specific perturbation in the mapped placental receiver with signed program readouts"
    out = pd.DataFrame({
        "candidate": x["candidate"],
        "phase3_sender_robustness": x["baseline_sender_robustness"],
        "sender_evidence_level": x["S1_S2"],
        "licensing": x["licensing_class"],
        "P1_scope": x["P1_PARACRINE_CORE"], "P2_scope": x["P2_EXTRACELLULAR_EXTENDED"], "P3_scope": x["P3_FULL_LR_SENSITIVITY"],
        "receptor_competent_module_n": x["receptor_competent_module_n"],
        "lr_supported_edge_n": x["lr_supported_edge_n"], "lr_supported_receptor_n": x["lr_supported_receptor_n"], "lr_supported_receptors": x["lr_supported_receptors"],
        "target_compatible_module_n": x["target_compatible_module_n"],
        "generic_signed_prior_reversal_module_n": x["reversal_supported_module_n"],
        "generic_signed_prior_tier_A_modules": x["tier_A_modules_phase4a"],
        "disease_concordant_axis_n": x["disease_concordant_module_n"],
        "disease_concordant_modules": x["disease_concordant_modules_phase4a"],
        "mixed_direction_flag": x["mixed_signed_direction_across_modules_phase4a"],
        "protein_topology": x["topology_classification"],
        "huc_wj_msc_ev_protein_evidence": x["EV_direct"],
        "ev_protein_datasets": x["direct_datasets"],
        "soluble_conditioned_medium_evidence": x["soluble_CM_direct"],
        "empirical_signed_perturbation_evidence": x["empirical_perturbation_class"],
        "empirical_perturbation_modules": x["empirical_modules"],
        "perturbation_biological_context": x["perturbation_context"],
        "pe_context_evidence": x["pe_context"],
        "direct_msc_pe_precedent": x["DIRECT_MSC_PE_PRECEDENT"],
        "novelty": x["novelty_classification"],
        "receiver_evidence_level": x["RECEIVER_EVIDENCE_STRENGTH"],
        "phase4a_internal_label": x["phase4a_internal_label"],
        "manuscript_facing_phase4a_label": x["manuscript_facing_label"],
        "final_deterministic_category": x["corrected_primary_classification"],
        "context_dependent": x["MIXED_DIRECTION_RISK"].eq("MIXED_DIRECTION_CONTEXT_DEPENDENT").map({True:"YES", False:"NO"}),
    })
    out["final_frozen_status"] = out["candidate"].map(lambda g: "MECHANICALLY_TRIANGULATED_LEAD_HYPOTHESIS" if g == "ENPP1" else "FROZEN_PHASE4B1_CLASSIFICATION")
    out["ev_only_protein_source"] = ((out["huc_wj_msc_ev_protein_evidence"] == "YES") & (out["soluble_conditioned_medium_evidence"] != "YES")).map({True:"YES", False:"NO"})
    out["non_placental_empirical_reversal"] = out.apply(lambda r: "YES" if r.candidate == "ENPP1" else ("NOT_ESTABLISHED" if "REVERSAL" not in r.empirical_signed_perturbation_evidence else "CONTEXT_REVIEW_REQUIRED"), axis=1)
    out["direct_pe_experimental_validation"] = "NO"
    out["manuscript_facing_interpretation"] = x.apply(interpretation, axis=1)
    out["required_validation_experiment"] = x.apply(experiment, axis=1)
    out["source_url"] = "results/04_phase4a/integration/phase4a_candidate_hierarchy.csv|results/04_phase4b1/corrected_phase4b_candidate_evidence_matrix.csv|results/04_phase4b1/corrected_phase4b_candidate_classification.csv"
    return out.set_index("candidate").loc[PHASE4B17].reset_index()


def build_negative_results() -> pd.DataFrame:
    rows = [
        ("NEG01", "1B", "Universal individual-gene PE signature", "Zero STABLE genes under preregistered six-cohort criteria", "0 STABLE; 0 meta BH FDR<0.05 among 17,731 genes estimable in >=4 cohorts", "The strict gene-level meta-analysis was negative", "This does not imply absence of PE biology", "FROZEN_NEGATIVE", "results/01_phase1b/meta/stable_pe_genes.csv|docs/PHASE1B_PE_DISEASE_SIGNATURE_REPORT.md"),
        ("NEG02", "2B", "Robust bulk validation of receiver programs", "Zero BULK_ROBUST_SUPPORT constituent sets and zero BULK_MODULE_SUPPORTED modules", "0/19 constituent sets; 0/10 tested modules", "Bulk evidence was limited/directional, not robust", "Directional tissue-level support remains reportable", "FROZEN_NEGATIVE", "results/02_phase2b/meta/program_gene_set_meta_analysis.csv|results/02_phase2b/meta/program_module_validation.csv"),
        ("NEG03", "4A", "Sender compatibility with SCT oxidative-phosphorylation receiver", "PROGRAM_MODULE_10 had zero target-compatible sender axes", "0 axes", "No frozen computational sender connection was found", "External mitochondrial literature cannot retroactively create an axis", "FROZEN_NEGATIVE", "results/04_phase4a/targets/nichenet_target_compatibility.csv|docs/PHASE4A_SENDER_RECEIVER_INTEGRATION_REPORT.md"),
        ("NEG04", "2A.2", "Independent Zheng EOPE validation", "Reported 3 EOPE + 2 control cohort could not be reconstructed with donor-linked expression and annotations", "0 targeted validation results", "No direct reconstructed Zheng validation exists", "This is data non-reusability, not biological contradiction", "FROZEN_UNRESOLVED", "results/02_phase2a2/external_validation/zheng_eope_dataset_audit.csv"),
        ("NEG05", "4B.1", "Soluble conditioned-medium proteomic confirmation", "No complete technically evaluable soluble hUC/WJ-MSC conditioned-medium proteome table was available", "0/17 with complete-table soluble-CM confirmation", "EV evidence must remain compartment-specific", "It does not prove candidates are absent from soluble CM", "ACTIVE_LIMITATION", "results/04_phase4b1/hucmsc_proteomics_dataset_registry.csv|results/04_phase4b1/corrected_hucmsc_protein_source_evidence.csv"),
        ("NEG06", "4B.1", "PE-specific experimental validation of ENPP1", "No direct PE experimental validation identified in the frozen audit", "0 direct PE experiments", "ENPP1 remains a hypothesis", "It is not a validated therapeutic factor or mediator", "ACTIVE_LIMITATION", "results/04_phase4b/disease/pe_candidate_context_evidence.csv|results/04_phase4b/novelty/direct_msc_pe_overlap.csv"),
        ("NEG07", "4B", "Independent placental perturbation confirmation of ENPP1 reversal", "No placental-cell ENPP1 perturbation confirmed reversal of the frozen receiver program", "0 placental perturbation studies", "Existing module-matched reversal evidence is non-placental", "Non-placental cGAMP-STING evidence does not establish placental direction", "ACTIVE_LIMITATION", "results/04_phase4b/perturbation/empirical_signed_perturbation_evidence.csv"),
        ("NEG08", "4A", "Sign consistency across receiver modules", "Twelve of 17 Tier A candidates had disease-concordant signed evidence in at least one other module", "12/17 mixed/disease-concordant elsewhere", "Context-dependent or potentially aggravating axes were retained", "Tier A does not imply uniform rescue", "FROZEN_RISK", "results/04_phase4a/integration/phase4a_candidate_hierarchy.csv|results/04_phase4b/integration/mixed_direction_stress_test.csv"),
    ]
    return pd.DataFrame(rows, columns=["negative_result_id","phase","question","frozen_result","quantitative_value","interpretation","prohibited_inference","status","source_url"])


def build_claims() -> pd.DataFrame:
    rows = [
        ("CLAIM01", "The strict six-cohort gene-level meta-analysis yielded no STABLE PE genes under the preregistered rule.", "No stable universal individual-gene PE signature was identified under the frozen criteria.", "PE has no transcriptomic abnormalities.", "results/01_phase1b/meta/pe_gene_meta_analysis.csv|results/01_phase1b/meta/stable_pe_genes.csv", "Cross-platform effect heterogeneity and incomplete covariates"),
        ("CLAIM02", "Corrected pregnancy-level analysis localized PE-associated programs to selected placental cell states in Admati.", "PE-associated interferon and related programs were localized to selected cell states in the corrected primary scRNA analysis.", "These cell states were independently validated in bulk placenta.", "results/02_phase2a2/corrected_analysis/corrected_program_modules.csv", "Primary localization largely depends on one scRNA dataset"),
        ("CLAIM03", "Yang LOPE provided partial directional support for placental-stromal interferon themes.", "The stromal interferon theme had partial independent LOPE scRNA directional support.", "All receiver cell-type localizations replicated independently.", "results/02_phase2a1/yang_lope_replication.csv|results/04_phase4a/freeze/phase4_receiver_hierarchy.csv", "n=3 versus 3 and conservative mapping; SCT/Hofbauer limitations"),
        ("CLAIM04", "Five receiver modules showed bulk tissue-level directional support; none met robust module criteria.", "Single-cell-defined programs showed limited directional support across independent bulk cohorts.", "Hofbauer or SCT localization was validated by bulk.", "results/02_phase2b/meta/program_module_validation.csv", "Bulk cannot assign cellular origin"),
        ("CLAIM05", "Independent hUC-MSC datasets defined 214 robust transcriptomic S1/S2 sender candidates.", "A donor-aware, cross-dataset transcriptomic sender ligand universe was frozen.", "All 214 ligands are secreted proteins or therapeutic mediators.", "results/03_phase3/sender/frozen_phase4_sender_candidates.csv", "Transcript abundance is not secretion or functional exposure"),
        ("CLAIM06", "Licensing responses were described within donor×passage strata.", "Inflammatory licensing altered availability of some sender transcripts in a donor/passage-dependent manner.", "Licensing-up means therapeutically beneficial.", "results/03_phase3/licensing/licensing_effect_by_donor_passage.csv", "Only two donors; licensing is context, not efficacy"),
        ("CLAIM07", "Blinded integration identified receptor and target compatibility between frozen sender and receiver arms.", "Many frozen ligands were mechanistically compatible with PE receiver programs.", "hUC-MSC ligands reverse PE programs.", "results/04_phase4a/lr/sender_receiver_lr_compatibility.csv|results/04_phase4a/targets/nichenet_target_compatibility.csv", "NicheNet compatibility is largely unsigned"),
        ("CLAIM08", "Generic signed priors yielded 17 Tier A candidates, with mixed direction retained.", "Seventeen candidates satisfied a generic-signed-prior reversal rule.", "The 17 are predicted therapeutic rescue factors.", "results/04_phase4a/integration/phase4a_candidate_hierarchy.csv", "Generic network propagation is not empirical placental reversal"),
        ("CLAIM09", "Fifteen of 17 candidates were detected in hUC/WJ-MSC EV protein resources.", "Direct EV protein evidence supports extracellular-source plausibility for 15 candidates.", "The proteins are soluble conditioned-medium factors.", "results/04_phase4b1/candidate_protein_detection_matrix.csv", "EV and soluble compartments are distinct"),
        ("CLAIM10", "ENPP1 alone satisfied the unchanged TRIANGULATED_HIGH_PRIORITY rule after protein-source correction.", "ENPP1 emerged as the only candidate satisfying the prespecified multi-domain triangulation rule.", "ENPP1 mediates hUC-MSC therapy in PE.", "results/04_phase4b1/corrected_phase4b_candidate_classification.csv", "EV-only source, non-placental perturbation, mixed receiver directions, no PE experiment"),
        ("CLAIM11", "ENPP1 is the mechanically triangulated lead hypothesis.", "ENPP1 is the leading experimentally testable communication hypothesis emerging from the frozen framework.", "ENPP1 is a validated therapeutic factor, protective mechanism, or established PE target.", "results/final_synthesis/final_17_candidate_evidence.csv|results/04_phase4b1/corrected_phase4b_candidate_classification.csv", "Causality, placental exposure, direction and efficacy are untested"),
        ("CLAIM12", "PROGRAM_MODULE_10 had zero target-compatible sender axes under the frozen rule.", "No frozen sender connection was identified for the SCT oxidative-phosphorylation receiver.", "External mitochondrial biology establishes a Phase4A Module10 axis.", "results/04_phase4a/targets/nichenet_target_compatibility.csv", "Frozen negative result cannot be rescued by later literature"),
    ]
    return pd.DataFrame(rows, columns=["claim_id","supported_claim","maximum_allowed_wording","prohibited_overclaim","evidence_source","key_limitation"]).assign(source_url=lambda d: d["evidence_source"])


def build_figures() -> pd.DataFrame:
    rows = [
        ("F1A","Figure 1","Study design and evidence-arm separation","results/final_synthesis/phase_registry.csv","phase, question, gate, status","phase","13 phases","descriptive gate sequence","Sender and receiver evidence remained independent until Phase4A","Workflow summary only","MAIN"),
        ("F1B","Figure 1","Frozen negative and correction checkpoints","results/final_synthesis/negative_results_register.csv","phase, frozen_result, status","frozen analytical result","8 registered negatives","counts and categorical gates","Negative results and corrections shaped the analysis path","Not an efficacy endpoint","MAIN"),
        ("F2A","Figure 2","Admati expression-layer provenance correction","results/02_phase2a2/provenance/admati_expression_layer_audit.csv","file, value behavior, raw-count status","public expression file","audited public files","distribution/checksum audit","The public layer required continuous-expression modeling","True UMI counts were not public","MAIN"),
        ("F2B","Figure 2","Frozen-20 corrected retest","results/02_phase2a2/corrected_analysis/frozen20_corrected_retest.csv","module/gene set, EOPE and LOPE effect, FDR, classification","pregnancy within cell type","20 frozen hypotheses","limma effect and within-family BH FDR","Legacy hypotheses were retested under the corrected model","Same cohort; not independent validation","MAIN"),
        ("F2C","Figure 2","Corrected receiver modules","results/02_phase2a2/corrected_analysis/corrected_program_modules.csv","module, cell type, constituents, classification","nonredundant module","11 modules","membership overlap and corrected support","Receiver hypotheses were redundancy-reduced","Module names do not imply causal mechanisms","MAIN"),
        ("F3A","Figure 3","Receiver evidence hierarchy","results/final_synthesis/final_receiver_evidence.csv","module, Admati, Yang, bulk, final level","module","11 frozen modules","evidence-presence matrix","R1/R2A/R2B/HOLD hierarchy integrates separate evidence domains","Bulk does not validate localization","MAIN"),
        ("F3B","Figure 3","Bulk program validation","results/02_phase2b/meta/program_module_validation.csv","robust/directional constituent counts, module class","module across six cohorts","10 tested modules","random-effects constituent-set criteria","Five modules had directional but no robust support","Zero robust modules","MAIN"),
        ("F3C","Figure 3","Yang LOPE targeted evidence","results/02_phase2a1/yang_lope_replication.csv","cell mapping, module, effect direction, adjusted evidence","pregnancy","LOPE 3 vs control 3","patient-level targeted score comparison","Partial external stromal IFN support","Small n and incomplete mapping","SUPPLEMENTARY"),
        ("F4A","Figure 4","Sender universe and donor support","results/final_synthesis/final_sender_evidence.csv","S1/S2, donor counts, baseline class","candidate ligand","214 candidates","frozen donor-support rule","Independent datasets define a reproducible transcriptomic sender universe","Not protein secretion","MAIN"),
        ("F4B","Figure 4","Licensing response classes","results/03_phase3/licensing/licensing_ligand_classification.csv","gene, donor/passages, licensing class","donor×passage","2 donors; four frozen strata","within-stratum effect summaries","Licensing is donor/passage context","Not therapeutic direction","MAIN"),
        ("F5A","Figure 5","Blinded sender-receiver compatibility","results/04_phase4a/integration/sender_receiver_evidence_matrix.csv","ligand, module, receptor, target, signed class","ligand-module axis","1,284 axes","frozen deterministic evidence axes","Compatibility was broad but signed reversal limited","Network priors are not causal proof","MAIN"),
        ("F5B","Figure 5","Generic signed and disease-concordant axes","results/04_phase4a/signed/signed_reversal_analysis.csv","ligand, receiver, reversal class","ligand-module axis","frozen eligible axes","signed propagation class","Mixed directions remain visible","Generic signed prior only","MAIN"),
        ("F5C","Figure 5","Module 10 frozen negative","results/04_phase4a/targets/nichenet_target_compatibility.csv","module, target compatibility","ligand-module axis","PROGRAM_MODULE_10; 214 senders","target-compatibility rule","No sender axis met target compatibility","Cannot be replaced by external literature","MAIN"),
        ("F6A","Figure 6","Protein-source triangulation","results/04_phase4b1/candidate_protein_detection_matrix.csv","candidate, dataset, compartment, detection confidence","candidate×dataset","17 candidates; 4 evaluable extracellular resources","identifier-mapped protein detection","EV protein evidence materially corrected the source layer","No complete soluble-CM proteome","MAIN"),
        ("F6B","Figure 6","Final 17-candidate evidence matrix","results/final_synthesis/final_17_candidate_evidence.csv","separate evidence domains and final category","candidate","17 frozen candidates","unchanged deterministic categories","Only ENPP1 meets the mechanical high-priority rule","Category is not therapeutic validation","MAIN"),
        ("F6C","Figure 6","ENPP1 evidence and validation gaps","results/final_synthesis/final_17_candidate_evidence.csv","ENPP1 source, context, receiver, gaps","one lead hypothesis","ENPP1","evidence-domain presence/absence","ENPP1 is a testable cross-dataset hypothesis","EV-only; non-placental perturbation; no PE experiment","MAIN"),
    ]
    return pd.DataFrame(rows, columns=["panel_id","figure","panel_title","source_file","variables","statistical_unit","n","exact_statistic","interpretation","limitation","placement"]).assign(source_url=lambda d: d.source_file)


def build_supplementary(root: Path) -> pd.DataFrame:
    specs = [
        ("ST01","Dataset feasibility registry","results/00_dataset_audit/dataset_registry.csv","dataset","Phase0A/0B dataset provenance and eligibility"),
        ("ST02","Dataset overlap audit","results/00_dataset_audit/dataset_overlap_matrix.csv","dataset pair","GSE75010 and other overlap/leakage safeguards"),
        ("ST03","Sample and donor registries","results/00_dataset_audit/sample_registry.csv|results/02_phase2a/metadata/patient_registry.csv","sample/pregnancy","Frozen biological identities and labels"),
        ("ST04","Bulk preprocessing decisions","results/01_phase1a1/formal_phase1b_matrix_registry.csv|results/01_phase1a1/preprocessing_comparison_registry.csv","cohort","Formal matrix and reconstruction choices"),
        ("ST05","Phase1B gene-level negative meta-analysis","results/01_phase1b/meta/pe_gene_meta_analysis.csv|results/01_phase1b/meta/stable_pe_genes.csv","gene","Complete gene-level outcome including zero STABLE result"),
        ("ST06","Receiver expression-layer correction audit","results/02_phase2a2/provenance/admati_expression_layer_audit.csv|results/02_phase2a2/provenance/admati_raw_count_search.csv","file/resource","Why the legacy count model was superseded"),
        ("ST07","Corrected receiver gene and program statistics","results/02_phase2a2/corrected_analysis/corrected_gene_statistics.csv|results/02_phase2a2/corrected_analysis/frozen20_corrected_retest.csv","gene/program×cell type","Corrected donor-level receiver results"),
        ("ST08","Phase2B program validation","results/02_phase2b/meta/program_gene_set_meta_analysis.csv|results/02_phase2b/meta/program_module_validation.csv","constituent set/module","Independent bulk directional and negative robust evidence"),
        ("ST09","Frozen sender ligand universe","results/03_phase3/sender/frozen_phase4_sender_candidates.csv","candidate ligand","All 214 S1/S2 candidates"),
        ("ST10","Licensing donor and passage strata","results/03_phase3/licensing/licensing_effect_by_donor_passage.csv|results/03_phase3/licensing/licensing_ligand_classification.csv","donor×passage×gene","Licensing context without cell-level pseudoreplication"),
        ("ST11","Complete Phase4A communication axes","results/04_phase4a/integration/sender_receiver_evidence_matrix.csv|results/04_phase4a/lr/sender_receiver_lr_compatibility.csv","ligand×module/receptor","Compatible, insufficient, concordant and reversal axes retained"),
        ("ST12","Frozen 17 candidate external evidence","results/final_synthesis/final_17_candidate_evidence.csv","candidate","All evidence domains and deterministic category"),
        ("ST13","Complete hUC-MSC protein-source audit","results/04_phase4b1/candidate_protein_detection_matrix.csv|results/04_phase4b1/protein_identifier_mapping_audit.csv","candidate×dataset","EV/CM/ECM compartments and identifier mappings"),
        ("ST14","Negative results register","results/final_synthesis/negative_results_register.csv","negative result","Frozen null, unavailable and discordant outcomes"),
        ("ST15","Literature search log","results/04_phase4b/novelty/literature_search_log.csv","candidate×query","Reproducible frozen Phase4B evidence search"),
        ("ST16","Claim ceiling matrix","results/final_synthesis/claim_ceiling_matrix.csv","claim","Permitted and prohibited manuscript wording"),
    ]
    rows=[]
    for sid,title,sources,unit,purpose in specs:
        paths=sources.split("|")
        counts=[]
        for rel in paths:
            p=root/rel
            counts.append(str(max(sum(1 for _ in p.open(encoding="utf-8-sig"))-1,0)) if p.exists() else "MISSING")
        rows.append((sid,title,sources,unit,";".join(counts),purpose,"REQUIRED","Include full small tables or archived large-table link/checksum",sources))
    return pd.DataFrame(rows,columns=["supplement_id","title","source_files","statistical_unit","source_row_counts","purpose","delivery_status","presentation_note","source_url"])


def build_roadmap() -> pd.DataFrame:
    rows = [
        (1,"Confirm ENPP1 protein across independent hUC-MSC donors","At least three new biological donors; orthogonal immunoassay plus targeted MS","Donor-level reproducible detection with negative controls","Failure indicates source non-reproducibility","NOT_PERFORMED"),
        (2,"Distinguish soluble versus EV-associated ENPP1","Matched EV-depleted CM, purified EV, cell lysate and media blank; nanoparticle/EV QC","Compartment-resolved abundance and recovery","EV restriction determines exposure model","NOT_PERFORMED"),
        (3,"Quantify ENPP1 after inflammatory licensing","Same donor and passage, naive versus IFN-gamma/TNF-alpha; repeated donors","Donor-paired ENPP1 protein effect with uncertainty","Do not infer benefit from up/down licensing alone","NOT_PERFORMED"),
        (4,"Expose relevant placental receiver cells to hUC-MSC CM and EV","Primary placental stromal cells first; macrophage/Hofbauer model second; vehicle and donor-matched preparations","Pregnancy/cell-preparation-level IFN program score and viability","Establish exposure response before mechanism","NOT_PERFORMED"),
        (5,"Perturb ENPP1 specifically","ENPP1 depletion/neutralization in CM/EV and recombinant/add-back controls","Loss and rescue of molecular effect","Separates ENPP1 dependence from correlated cargo","NOT_PERFORMED"),
        (6,"Measure receiver IFN-related state","Frozen Module07 primary and Modules01/04 sensitivity; targeted RNA/protein readouts","Signed change opposite PE-associated direction without toxicity","Use frozen program membership; no outcome-selected targets","NOT_PERFORMED"),
        (7,"Test whether ENPP1 dependence is causal","Factorial donor×compartment×ENPP1 perturbation with blinded analysis","Specific loss-of-effect and add-back restoration across donors","Required before mediator language","NOT_PERFORMED"),
        (8,"Assess PE-relevant functional phenotypes only after molecular causality","Trophoblast/stromal/immune functional assays selected prospectively","Predefined functional endpoint with dose and safety window","Functional benefit cannot substitute for receiver-direction mechanism","NOT_PERFORMED"),
    ]
    df=pd.DataFrame(rows,columns=["step","future_experiment","design_unit_and_controls","decision_readout","interpretive_rule","completion_status"])
    df["source_url"]="docs/FINAL_SYNTHESIS_REPORT.md|results/final_synthesis/final_17_candidate_evidence.csv"
    return df


def build_risks() -> pd.DataFrame:
    rows=[
        ("FS_RISK01","CRITICAL","RECEIVER_PROVENANCE","Admati true raw UMI counts were not publicly recoverable; corrected analysis used donor-mean normalized/ceiled expression","Use only corrected Phase2A2 statistics; label Phase2A legacy","OPEN_FROZEN","results/02_phase2a2/provenance/admati_expression_layer_audit.csv"),
        ("FS_RISK02","CRITICAL","CLINICAL_CONFOUNDING","EOPE delivery mode/induction and IUGR are structurally non-estimable; early controls n=3","No adjusted causal disease interpretation","OPEN_FROZEN","results/02_phase2a1/clinical_confounding_audit.csv"),
        ("FS_RISK03","HIGH","LOCALIZATION","Most receiver cell-type localizations depend on Admati alone","Separate localization evidence from bulk program support","OPEN","results/final_synthesis/final_receiver_evidence.csv"),
        ("FS_RISK04","HIGH","BULK_VALIDATION","Zero receiver modules met robust bulk criteria","Describe only directional tissue-level support","OPEN_FROZEN","results/02_phase2b/meta/program_module_validation.csv"),
        ("FS_RISK05","HIGH","SENDER_REPLICATION","Licensing arm contains only two independent donors and passage-dependent strata","Treat licensing as context, not benefit","OPEN_FROZEN","results/03_phase3/licensing/licensing_ligand_classification.csv"),
        ("FS_RISK06","HIGH","PROTEIN_SCOPE","Only the frozen 17 underwent systematic topology/proteomics audit; 197 senders remain transcriptomic-only for protein inference","Do not generalize protein evidence to all 214","OPEN","results/final_synthesis/final_sender_evidence.csv"),
        ("FS_RISK07","CRITICAL","COMPARTMENT","Phase4B1 evidence is EV-based and no complete soluble-CM proteome was evaluable","Never equate EV detection with soluble secretion","OPEN","results/04_phase4b1/corrected_hucmsc_protein_source_evidence.csv"),
        ("FS_RISK08","CRITICAL","DIRECTION","NicheNet compatibility is unsigned and the Phase4A signed layer is generic","Use reversal language only for explicitly signed evidence and retain insufficiency","OPEN_FROZEN","results/04_phase4a/signed/signed_reversal_analysis.csv"),
        ("FS_RISK09","CRITICAL","ENPP1_CONTEXT","ENPP1 reversal evidence is non-placental and candidate has mixed receiver-module directions","Keep context-dependent hypothesis wording","OPEN","results/final_synthesis/final_17_candidate_evidence.csv"),
        ("FS_RISK10","CRITICAL","CAUSALITY","No wet-lab PE or placental validation was performed","No therapeutic mediator, efficacy or target-established claims","OPEN","results/final_synthesis/experimental_validation_roadmap.csv"),
        ("FS_RISK11","MEDIUM","HISTORICAL_DOCUMENTATION","Phase2A and Phase4B historical files contain superseded method/protein statements","Preserve history; current synthesis explicitly points to Phase2A2 and Phase4B1 corrections","RESOLVED_BY_VERSIONING","docs/PHASE2A2_RECEIVER_FRAMEWORK_CORRECTION_REPORT.md|docs/PHASE4B1_PROTEIN_SOURCE_COMPLETENESS_REPORT.md"),
        ("FS_RISK12","HIGH","NEGATIVE_RESULT","PROGRAM_MODULE_10 has zero frozen target-compatible sender axes","Do not create a retrospective connection from literature","OPEN_FROZEN","results/04_phase4a/targets/nichenet_target_compatibility.csv"),
    ]
    return pd.DataFrame(rows,columns=["risk_id","severity","domain","risk","mitigation_or_claim_ceiling","status","source_url"])


def validate(outputs: dict[str, pd.DataFrame]) -> None:
    assert len(outputs["phase_registry.csv"]) == 13
    receiver=outputs["final_receiver_evidence.csv"]
    assert len(receiver)==11
    assert receiver["final_receiver_level"].value_counts().to_dict()=={"R2A":4,"NOT_PROMOTED":4,"R1":1,"R2B":1,"HOLD":1}
    sender=outputs["final_sender_evidence.csv"]
    assert len(sender)==214 and sender["gene"].nunique()==214
    assert sender["sender_evidence_level"].value_counts().to_dict()=={"S2":176,"S1":38}
    candidates=outputs["final_17_candidate_evidence.csv"]
    assert candidates["candidate"].tolist()==PHASE4B17
    enpp=candidates.set_index("candidate").loc["ENPP1"]
    assert enpp["final_frozen_status"]=="MECHANICALLY_TRIANGULATED_LEAD_HYPOTHESIS"
    assert enpp["context_dependent"]=="YES" and enpp["ev_only_protein_source"]=="YES"
    assert enpp["non_placental_empirical_reversal"]=="YES" and enpp["direct_pe_experimental_validation"]=="NO"
    assert (candidates["final_deterministic_category"]=="TRIANGULATED_HIGH_PRIORITY").sum()==1
    assert len(outputs["negative_results_register.csv"])>=8
    for name,df in outputs.items():
        assert len(df)>0, name
        assert "source_url" in df.columns, name
        assert df["source_url"].astype(str).str.len().gt(0).all(), name


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args=parser.parse_args()
    root=args.root.resolve()
    outdir=root/"results"/"final_synthesis"
    outputs={
        "phase_registry.csv":build_phase_registry(),
        "final_receiver_evidence.csv":build_receiver(root),
        "final_sender_evidence.csv":build_sender(root),
        "final_17_candidate_evidence.csv":build_candidates(root),
        "negative_results_register.csv":build_negative_results(),
        "claim_ceiling_matrix.csv":build_claims(),
        "figure_source_registry.csv":build_figures(),
        "supplementary_table_registry.csv":build_supplementary(root),
        "experimental_validation_roadmap.csv":build_roadmap(),
        "final_synthesis_risk_flags.csv":build_risks(),
    }
    validate(outputs)
    for name,df in outputs.items(): write_csv(df,outdir/name)
    manifest={
        "built_from_existing_outputs_only":True,
        "input_sha256":{
            rel:sha256(root/rel) for rel in [
                "results/04_phase4a/freeze/phase4_receiver_hierarchy.csv",
                "results/03_phase3/sender/frozen_phase4_sender_candidates.csv",
                "results/04_phase4a/integration/phase4a_candidate_hierarchy.csv",
                "results/04_phase4b1/corrected_phase4b_candidate_evidence_matrix.csv",
                "results/04_phase4b1/corrected_phase4b_candidate_classification.csv",
            ]
        },
        "output_rows":{name:len(df) for name,df in outputs.items()},
    }
    (outdir/"build_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    print(json.dumps(manifest["output_rows"],indent=2))


if __name__=="__main__":
    main()
