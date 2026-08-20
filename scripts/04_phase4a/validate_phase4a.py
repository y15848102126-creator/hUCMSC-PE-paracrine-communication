#!/usr/bin/env python3
"""Mechanical validation for frozen Phase 4A outputs."""

from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/04_phase4a"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> int:
    cfg = json.loads((ROOT / "config/phase4a_analysis.json").read_text(encoding="utf-8"))
    expected = [
        "freeze/phase4_sender_scopes.csv", "freeze/phase4_receiver_hierarchy.csv",
        "receptors/receiver_receptor_competence.csv", "lr/sender_receiver_lr_compatibility.csv",
        "targets/nichenet_target_compatibility.csv", "targets/ligand_receiver_target_edges.csv",
        "signed/signed_network_resource_registry.csv", "signed/signed_reversal_analysis.csv",
        "integration/sender_receiver_evidence_matrix.csv", "integration/phase4a_candidate_hierarchy.csv",
        "integration/disease_concordant_candidates.csv", "qc/phase4a_risk_flags.csv",
        "qc/phase4a_method_diagnostics.csv",
    ]
    for relative in expected:
        path = OUT / relative
        require(path.is_file() and path.stat().st_size > 20, f"missing/empty: {relative}")
        frame = pd.read_csv(path)
        require(len(frame.columns) >= 2, f"too few columns: {relative}")
        require(any(c in frame.columns for c in ("source_url", "url")), f"source field missing: {relative}")

    sender = pd.read_csv(OUT / "freeze/phase4_sender_scopes.csv")
    receiver = pd.read_csv(OUT / "freeze/phase4_receiver_hierarchy.csv")
    receptor = pd.read_csv(OUT / "receptors/receiver_receptor_competence.csv")
    lr = pd.read_csv(OUT / "lr/sender_receiver_lr_compatibility.csv")
    target = pd.read_csv(OUT / "targets/nichenet_target_compatibility.csv")
    signed = pd.read_csv(OUT / "signed/signed_reversal_analysis.csv")
    evidence = pd.read_csv(OUT / "integration/sender_receiver_evidence_matrix.csv")
    hierarchy = pd.read_csv(OUT / "integration/phase4a_candidate_hierarchy.csv")
    disease = pd.read_csv(OUT / "integration/disease_concordant_candidates.csv")

    require(len(sender) == 214, "sender freeze must have 214 rows")
    require((sender.P1_PARACRINE_CORE == "YES").sum() == 148, "P1 count")
    require((sender.P2_EXTRACELLULAR_EXTENDED == "YES").sum() == 190, "P2 count")
    require((sender.P3_FULL_LR_SENSITIVITY == "YES").sum() == 214, "P3 count")
    tested = receiver[receiver.phase4a_analysis_scope.isin(["PRIMARY", "SECONDARY_SENSITIVITY"])]
    require(set(tested.program_module) == {"PROGRAM_MODULE_01", "PROGRAM_MODULE_04", "PROGRAM_MODULE_05", "PROGRAM_MODULE_06", "PROGRAM_MODULE_07", "PROGRAM_MODULE_10"}, "receiver test set")
    require(receiver.loc[receiver.program_module.eq("PROGRAM_MODULE_08"), "phase4a_analysis_scope"].iloc[0] == "HOLD_NOT_TESTED", "module 08 must stay hold")
    require(len(target) == 214 * 6 and len(signed) == 214 * 6 and len(evidence) == 214 * 6, "complete ligand-module universe")
    require(target.nichenet_directionality.eq("UNSIGNED_COMPATIBILITY_ONLY").all(), "NicheNet must remain unsigned")
    require(set(signed.signed_reversal_class).issubset({"REVERSAL_SUPPORTED", "DISEASE_CONCORDANT_POTENTIAL", "SIGNED_EVIDENCE_INSUFFICIENT"}), "signed classes")
    require(len(hierarchy) == 214 and hierarchy.ligand.nunique() == 214, "one hierarchy row per sender")
    require("mixed_signed_direction_across_modules" in hierarchy.columns, "mixed signed direction flag missing")
    require(len(disease) == signed.signed_reversal_class.eq("DISEASE_CONCORDANT_POTENTIAL").sum(), "disease-concordant axes lost")
    require(lr[lr.receptor_competence.eq("RECEPTOR_COMPETENT")].ligand.nunique() == 192, "competent ligand count")
    require(receptor.receptor_competence.isin(["RECEPTOR_COMPETENT", "RECEPTOR_WEAK", "RECEPTOR_NOT_SUPPORTED"]).all(), "receptor class")
    require((target.target_compatibility_class.eq("TARGET_COMPATIBLE")).sum() == 791, "frozen target-compatible result changed")
    require((signed.signed_reversal_class.eq("REVERSAL_SUPPORTED")).sum() == 87, "frozen signed reversal result changed")
    require((signed.signed_reversal_class.eq("DISEASE_CONCORDANT_POTENTIAL")).sum() == 285, "frozen disease concordance changed")
    require((hierarchy.best_phase4a_tier.eq("TIER_A_DIRECTIONAL_RESCUE_CANDIDATE")).sum() == 17, "Tier A count changed")
    require(not (ROOT / "results/04_phase4b").exists(), "Phase 4B output must not exist")
    report = (ROOT / "docs/PHASE4A_SENDER_RECEIVER_INTEGRATION_REPORT.md").read_text(encoding="utf-8")
    require("GO_TO_PHASE4B_WITH_RESTRICTIONS" in report, "gate missing")
    require("not therapeutic proof" in report.lower(), "interpretive warning missing")

    delivery = ROOT / "Phase4A_Delivery.zip"
    if delivery.exists():
        allowed_suffix = {".md", ".csv", ".png", ".txt"}
        with zipfile.ZipFile(delivery) as archive:
            names = archive.namelist()
            require(all(Path(name).suffix.lower() in allowed_suffix for name in names), "delivery contains disallowed format")
            forbidden = ("data/raw", "data/interim", ".git", ".rds", ".bin", ".mtx", ".fastq")
            require(not any(any(token in name.lower() for token in forbidden) for name in names), "delivery contains raw/interim/git data")
            require("PHASE4A_SENDER_RECEIVER_INTEGRATION_REPORT.md" in names, "report missing from delivery")
            require("git_log_oneline.txt" in names, "git log missing from delivery")
    print("PHASE4A_VALIDATION_OK: 13 CSV outputs; P1/P2/P3=148/190/214; 1,284 axes complete; Phase4B locked")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
