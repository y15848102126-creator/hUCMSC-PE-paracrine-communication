#!/usr/bin/env python3
"""Freeze Phase 4A sender scopes, receiver hierarchy, and signed resources."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "04_phase4a"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_gmt(path: Path, collection: str) -> dict[str, list[str]]:
    answer: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        answer[f"{collection}::{fields[0]}"] = sorted(set(fields[2:]))
    return answer


def main() -> int:
    cfg = json.loads((ROOT / "config" / "phase4a_analysis.json").read_text(encoding="utf-8"))
    assert cfg["outcome_inspection_status_at_freeze"] == "NO_PHASE4_COMPATIBILITY_OR_SIGNED_RESULTS_INSPECTED"
    for group in ("upstream_sha256", "resource_sha256", "msigdb_sha256"):
        for filename, expected in cfg[group].items():
            assert sha256(ROOT / filename).lower() == expected.lower(), filename

    for directory in ("freeze", "signed"):
        (OUT / directory).mkdir(parents=True, exist_ok=True)

    sender = pd.read_csv(ROOT / "results/03_phase3/sender/frozen_phase4_sender_candidates.csv")
    assert len(sender) == 214 and set(sender.sender_evidence_level) == {"S1", "S2"}
    scopes = sender[[
        "gene", "candidate_rank", "sender_evidence_level", "baseline_classification",
        "licensing_classification", "omnipath_secreted_annotation",
        "omnipath_extracellular_annotation", "protein_secretome_annotation",
    ]].copy()
    scopes["P1_PARACRINE_CORE"] = scopes.omnipath_secreted_annotation.eq("YES").map({True: "YES", False: "NO"})
    scopes["P2_EXTRACELLULAR_EXTENDED"] = (
        scopes.omnipath_secreted_annotation.eq("YES") | scopes.omnipath_extracellular_annotation.eq("YES")
    ).map({True: "YES", False: "NO"})
    scopes["P3_FULL_LR_SENSITIVITY"] = "YES"
    scopes["primary_scope_membership"] = "P3_FULL_LR_SENSITIVITY"
    scopes.loc[scopes.P2_EXTRACELLULAR_EXTENDED.eq("YES"), "primary_scope_membership"] = "P2_EXTRACELLULAR_EXTENDED"
    scopes.loc[scopes.P1_PARACRINE_CORE.eq("YES"), "primary_scope_membership"] = "P1_PARACRINE_CORE"
    scopes["scope_freeze_status"] = "FROZEN_BEFORE_PHASE4_COMPATIBILITY_OUTCOMES"
    scopes["interpretation"] = "SENDER_SCOPE_NOT_THERAPEUTIC_PRIORITY"
    scopes["source_url"] = "results/03_phase3/sender/frozen_phase4_sender_candidates.csv|config/phase4a_analysis.json"
    scopes = scopes.sort_values("candidate_rank")
    assert (scopes.P1_PARACRINE_CORE == "YES").sum() == 148
    assert (scopes.P2_EXTRACELLULAR_EXTENDED == "YES").sum() == 190
    scopes.to_csv(OUT / "freeze/phase4_sender_scopes.csv", index=False)

    gmt = {}
    gmt.update(read_gmt(ROOT / "data/raw/phase2a_resources/h.all.v2026.1.Hs.symbols.gmt", "HALLMARK"))
    gmt.update(read_gmt(ROOT / "data/raw/phase2a_resources/c2.cp.reactome.v2026.1.Hs.symbols.gmt", "REACTOME"))
    gmt.update(read_gmt(ROOT / "data/raw/phase2a_resources/c5.go.bp.v2026.1.Hs.symbols.gmt", "GOBP"))
    modules = pd.read_csv(ROOT / "results/02_phase2a2/corrected_analysis/corrected_program_modules.csv")
    modules = modules.loc[modules.record_type.eq("ORIGINAL_GENE_SET")].copy()
    evidence = pd.read_csv(ROOT / "results/02_phase2b/evidence/updated_receiver_evidence_hierarchy.csv")
    rows = []
    for module, group in modules.groupby("program_module", sort=True):
        memberships = []
        hashes = []
        for pathway in group.pathway:
            genes = gmt[pathway]
            memberships.extend(genes)
            hashes.append(hashlib.sha256(";".join(genes).encode()).hexdigest())
        union = sorted(set(memberships))
        first = group.iloc[0]
        rows.append({
            "program_module": module, "module_label": first.module_label, "celltype": first.celltype,
            "frozen_direction": first.frozen_direction, "constituent_gene_set_n": len(group),
            "constituent_pathways": ";".join(group.pathway), "constituent_membership_sha256": ";".join(hashes),
            "module_union_gene_n": len(union), "module_union_gene_membership": ";".join(union),
            "module_union_sha256": hashlib.sha256(";".join(union).encode()).hexdigest(),
        })
    receiver = pd.DataFrame(rows).merge(evidence[[
        "program_module", "CORRECTED_ADMATI_SUPPORT", "EXTERNAL_SCRNA_SUPPORT",
        "INDEPENDENT_BULK_PROGRAM_SUPPORT", "celltype_localization_claim",
    ]], on="program_module", how="left")
    receiver["receiver_level"] = "NOT_PROMOTED"
    receiver.loc[receiver.program_module.eq("PROGRAM_MODULE_07"), "receiver_level"] = "R1"
    receiver.loc[receiver.program_module.isin(["PROGRAM_MODULE_01", "PROGRAM_MODULE_04", "PROGRAM_MODULE_05", "PROGRAM_MODULE_10"]), "receiver_level"] = "R2A"
    receiver.loc[receiver.program_module.eq("PROGRAM_MODULE_06"), "receiver_level"] = "R2B"
    receiver.loc[receiver.program_module.eq("PROGRAM_MODULE_08"), "receiver_level"] = "HOLD"
    receiver["phase4a_analysis_scope"] = receiver.receiver_level.map({
        "R1": "PRIMARY", "R2A": "PRIMARY", "R2B": "SECONDARY_SENSITIVITY",
        "HOLD": "HOLD_NOT_TESTED", "NOT_PROMOTED": "NOT_TESTED_NOT_PROMOTED",
    })
    receiver["molecular_evidence_counting_rule"] = "IDENTICAL_MEMBERSHIP_SHA256_NOT_INDEPENDENT;CELL_LOCALIZATION_SEPARATE"
    receiver["freeze_status"] = "FROZEN_BEFORE_PHASE4_COMPATIBILITY_OUTCOMES"
    receiver["source_url"] = "results/02_phase2a2/corrected_analysis/corrected_program_modules.csv|results/02_phase2b/evidence/updated_receiver_evidence_hierarchy.csv|MSigDB:2026.1.Hs|config/phase4a_analysis.json"
    assert len(receiver) == 11 and (receiver.receiver_level == "R1").sum() == 1 and (receiver.receiver_level == "R2A").sum() == 4
    receiver.sort_values("program_module").to_csv(OUT / "freeze/phase4_receiver_hierarchy.csv", index=False)

    manifest = pd.read_csv(ROOT / "data/raw/phase4a/download_manifest.csv")
    registry = pd.DataFrame({
        "resource_id": ["NICHENET_LIGAND_TARGET", "OMNIPATH_LR_CROSSCHECK", "OMNIPATH_SIGNED_ACTIVITY_FLOW", "COLLECTRI_SIGNED_TF_TARGET", "CYTOSIG_CORE_PERTURBATION"],
        "filename": manifest.filename,
        "version": ["NicheNet-v2 2021-12-21 / Zenodo 7074291", "OmniPath API snapshot 2026-08-15", "OmniPath activity-flow snapshot 2026-08-15", "CollecTRI via OmniPath snapshot 2026-08-15", "CytoSig core signature GitHub snapshot 2026-08-15"],
        "sha256": manifest.sha256, "bytes": manifest.bytes, "url": manifest.url,
        "phase4a_role": ["UNSIGNED_TARGET_COMPATIBILITY", "LR_EDGE_CROSSCHECK", "SIGNED_RECEPTOR_TO_TF_PROPAGATION", "SIGNED_TF_TO_TARGET_REGULON", "DIRECT_SIGNED_PERTURBATION_SENSITIVITY"],
        "directionality": ["UNSIGNED", "MIXED_EDGE_ANNOTATION", "DIRECTED_UNAMBIGUOUS_SIGN_ONLY", "DIRECTED_UNAMBIGUOUS_SIGN_ONLY", "SIGNED_EXPRESSION_RESPONSE"],
        "limitation": ["regulatory potential has no activation/inhibition sign", "database support does not prove signaling in placenta", "generic prior; paths limited to three edges", "generic TF-target prior; exact signs required", "exact ligand-column symbol matches only; context-generic"],
        "freeze_status": "FROZEN_BEFORE_PHASE4_COMPATIBILITY_OUTCOMES",
        "source_url": "data/raw/phase4a/download_manifest.csv|config/phase4a_analysis.json",
    })
    registry.to_csv(OUT / "signed/signed_network_resource_registry.csv", index=False)
    print("PHASE4A_FREEZE_OK: P1=148 P2=190 P3=214; primary receivers=5; sensitivity receivers=1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
