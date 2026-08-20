#!/usr/bin/env python3
"""Validate Phase 2A.2 required outputs, frozen history, and analysis gates."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/02_phase2a2"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(rel: str) -> list[dict[str, str]]:
    with (OUT / rel).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    required = [
        "provenance/admati_expression_layer_audit.csv", "provenance/admati_raw_count_search.csv",
        "corrected_analysis/corrected_gene_statistics.csv", "corrected_analysis/frozen20_corrected_retest.csv",
        "corrected_analysis/corrected_program_rediscovery.csv", "corrected_analysis/corrected_program_modules.csv",
        "external_validation/zheng_eope_dataset_audit.csv", "external_validation/zheng_eope_targeted_validation.csv",
        "external_validation/yang_lope_updated_evidence.csv", "evidence/receiver_module_evidence_hierarchy.csv",
        "evidence/phase2a2_risk_flags.csv", "figures/A_frozen20_corrected_direction.png",
        "figures/B_tier1_classification.png", "figures/C_corrected_module_status.png",
    ]
    for rel in required:
        p = OUT / rel
        assert p.exists() and p.stat().st_size > 0, rel
    assert (ROOT / "docs/PHASE2A2_CORRECTED_RECEIVER_ANALYSIS_PLAN.md").exists()
    report = (ROOT / "docs/PHASE2A2_RECEIVER_FRAMEWORK_CORRECTION_REPORT.md").read_text(encoding="utf-8")
    assert "GO_TO_PHASE2B_WITH_RESTRICTIONS" in report and "Phase 2B was not started" in report

    config = json.loads((ROOT / "config/phase2a2_analysis.json").read_text(encoding="utf-8"))
    for rel, expected in config["history_policy"]["frozen_sha256"].items():
        assert sha(ROOT / rel) == expected, f"history changed: {rel}"
    assert config["admati_provenance_decision"] == "PUBLIC_MATRIX_NORMALIZED_RAW_NOT_PUBLIC"
    search = rows("provenance/admati_raw_count_search.csv")
    assert any(r["result"] == "PUBLIC_MATRIX_NORMALIZED_RAW_NOT_PUBLIC" for r in search)

    genes = rows("corrected_analysis/corrected_gene_statistics.csv")
    assert len(genes) == 368370
    assert {r["contrast"] for r in genes} == {"EOPE", "LOPE"}
    assert all(r["source_url"] for r in genes)
    t1 = rows("corrected_analysis/frozen20_corrected_retest.csv")
    assert len(t1) == 20
    assert all(r["direction_agrees_EOPE"] == "YES" and r["direction_agrees_LOPE"] == "YES" for r in t1)
    assert all(float(r["BH_FDR_frozen20_EOPE"]) < .05 and float(r["BH_FDR_frozen20_LOPE"]) < .05 for r in t1)
    assert all(r["classification"] == "CORRECTED_SHARED_SUPPORT" for r in t1)
    red = rows("corrected_analysis/corrected_program_rediscovery.csv")
    assert len(red) == 55560 and all(r["interpretation"] == "EXPLORATORY_CORRECTED_REDISCOVERY" for r in red)
    modules = [r for r in rows("corrected_analysis/corrected_program_modules.csv") if r["record_type"] == "PROGRAM_MODULE"]
    assert len(modules) == 11 and all(r["corrected_module_status"] == "CORRECTED_ADMATI_SUPPORT" for r in modules)

    zheng = rows("external_validation/zheng_eope_dataset_audit.csv")
    assert len(zheng) == 5 and len({r["gsm"] for r in zheng}) == 5
    assert {r["gse"] for r in zheng} == {"GSE282038", "GSE267340", "GSE298119"}
    ztests = rows("external_validation/zheng_eope_targeted_validation.csv")
    assert len(ztests) == 20 and all(r["validation_status"] == "NOT_RUN_EOPE_ESTIMAND_CONTRADICTION" for r in ztests)
    assert all(not r["patient_level_effect"] and not r["exact_or_permutation_P"] for r in ztests)
    yang = rows("external_validation/yang_lope_updated_evidence.csv")
    assert len(yang) == 20 and all(r["rerun_status"] == "NOT_RERUN_VALUES_PRESERVED" for r in yang)
    hierarchy = rows("evidence/receiver_module_evidence_hierarchy.csv")
    assert len(hierarchy) == 11
    assert sum(r["evidence_level"] == "LEVEL_A" for r in hierarchy) == 2
    assert sum(r["evidence_level"] == "LEVEL_B" for r in hierarchy) == 8
    assert sum(r["phase2b_program_validation_candidate"].startswith("YES") for r in hierarchy) == 10
    print("Phase 2A.2 validation PASSED: history frozen; 368370 genes; 20 retests; 11 modules; Zheng NOT_RUN; 10 Phase2B candidates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
