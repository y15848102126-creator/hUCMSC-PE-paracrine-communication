# DEPRECATED / HISTORICAL ONLY: Phase 2A count-likelihood receiver inference is superseded by Phase 2A.2 pregnancy-level continuous-expression analysis. Excluded from the default execution workflow.
#!/usr/bin/env python3
"""Validate Phase 2A identity, model, program and delivery invariants."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/02_phase2a"


def require_columns(frame: pd.DataFrame, columns: set[str], name: str) -> None:
    missing = columns - set(frame.columns)
    assert not missing, f"{name}: missing {sorted(missing)}"


def main() -> int:
    required = [
        "metadata/patient_registry.csv", "metadata/celltype_annotation_registry.csv", "metadata/pseudobulk_eligibility.csv",
        "pseudobulk/pseudobulk_registry.csv", "DE/celltype_DE_summary.csv", "programs/pe_cellstate_programs.csv",
        "programs/shared_pe_programs.csv", "programs/eope_specific_programs.csv", "programs/lope_specific_programs.csv",
        "regulons/cellstate_regulon_activity.csv", "qc/phase2a_qc_summary.csv", "qc/phase2a_risk_flags.csv",
        "qc/pseudobulk_mds_coordinates.csv", "qc/phase2a_session_info.txt"
    ]
    for name in required:
        assert (OUT / name).exists(), name
    patients = pd.read_csv(OUT / "metadata/patient_registry.csv")
    assert len(patients) == patients.patient_id.nunique() == 26
    assert patients.cell_count.sum() == 86752
    assert patients.library_count.sum() == 31
    assert patients.pe_subtype_or_control_group.value_counts().to_dict() == {"EOPE":10,"LOPE":7,"LATE_CONTROL":6,"EARLY_CONTROL":3}
    assert patients.donor_identity_validated.eq("YES").all()
    ann = pd.read_csv(OUT / "metadata/celltype_annotation_registry.csv")
    assert len(ann) == 46 and ann.harmonized_annotation.nunique() == 15
    assert ann.marker_validation_status.eq("PASS").all()
    assert not any(".x" in c or ".y" in c for c in ann.columns)
    pb = pd.read_csv(OUT / "pseudobulk/pseudobulk_registry.csv")
    assert len(pb) == 26 * 15 and pb.patient_id.nunique() == 26
    assert pb.pseudobulk_library_umi.sum() == 962152952
    eligibility = pd.read_csv(OUT / "metadata/pseudobulk_eligibility.csv")
    assert len(eligibility) == 3 * 26 * 15
    assert set(eligibility.contrast) == {"EOPE","LOPE","COMBINED_PE_SECONDARY"}
    de_summary = pd.read_csv(OUT / "DE/celltype_DE_summary.csv")
    assert de_summary.full_rank.eq("YES").all()
    assert len(de_summary[de_summary.contrast.eq("EOPE")]) == 11
    assert len(de_summary[de_summary.contrast.eq("LOPE")]) == 12
    assert len(de_summary[de_summary.contrast.eq("COMBINED_PE_SECONDARY")]) == 13
    de_fields = {"gene","logFC","SE","statistic","P","BH_FDR","n_PE","n_control"}
    for path in sorted((OUT / "DE").glob("*_DE.csv")):
        frame = pd.read_csv(path, low_memory=False)
        require_columns(frame, de_fields | {"source_url","source_accession"}, path.name)
        assert frame.gene.notna().all() and frame.gene.is_unique, path.name
        assert frame.P.between(0,1).all() and frame.BH_FDR.between(0,1).all(), path.name
        assert frame.logFC.notna().all() and frame.statistic.notna().all(), path.name
        assert frame.n_PE.nunique() == frame.n_control.nunique() == 1, path.name
    programs = pd.read_csv(OUT / "programs/pe_cellstate_programs.csv")
    require_columns(programs,{"celltype","contrast","gene_set","NES","P","BH_FDR","direction","classification","all_loo_directions_same"},"programs")
    allowed = {"SHARED_PE","EOPE_ENRICHED","LOPE_ENRICHED","UNSTABLE","NOT_SIGNIFICANT","SECONDARY_COMBINED_SUPPORT_ONLY"}
    assert set(programs.classification) <= allowed
    assert not programs.duplicated(["celltype","contrast","pathway"]).any()
    assert programs.P.between(0,1).all() and programs.BH_FDR.between(0,1).all()
    assert programs.source_url.notna().all()
    shared = pd.read_csv(OUT / "programs/shared_pe_programs.csv")
    assert len(shared) == 20
    assert (shared.BH_FDR_EOPE < 0.05).all() and (shared.BH_FDR_LOPE < 0.05).all()
    assert ((shared.NES_EOPE > 0) == (shared.NES_LOPE > 0)).all()
    assert shared.all_loo_directions_same_EOPE.eq("YES").all() and shared.all_loo_directions_same_LOPE.eq("YES").all()
    regs = pd.read_csv(OUT / "regulons/cellstate_regulon_activity.csv")
    require_columns(regs,{"celltype","contrast","tf","target_n","activity_effect","SE","statistic","P","BH_FDR","network","method"},"regulons")
    assert regs.target_n.ge(10).all()
    assert regs.P.between(0,1).all() and regs.BH_FDR.between(0,1).all()
    mds = pd.read_csv(OUT / "qc/pseudobulk_mds_coordinates.csv")
    assert set(mds.outlier_flag) <= {"RETAIN","FLAG_REVIEW"}
    assert not patients.include_phase2a.eq("NO").any()
    qc = pd.read_csv(OUT / "qc/phase2a_qc_summary.csv")
    gate = qc.loc[qc.metric.eq("phase2a_gate"),"value"].iloc[0]
    assert gate == "GO_TO_PHASE2B_WITH_RESTRICTIONS"
    report = (ROOT / "docs/PHASE2A_PE_CELLSTATE_DISEASE_PROGRAM_REPORT.md").read_text(encoding="utf-8")
    for phrase in [gate,"26","SHARED_PE","Phase 1B","Phase 2B was not started"]:
        assert phrase in report
    config = json.loads((ROOT / "config/phase2a_analysis.json").read_text(encoding="utf-8"))
    assert config["eligibility"]["rules_frozen_before_outcome_analysis"] is True
    assert config["phase1b_integration"].startswith("Phase 1B is frozen negative")
    print("Phase 2A validation passed")
    return 0


if __name__ == "__main__": raise SystemExit(main())
