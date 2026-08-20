#!/usr/bin/env python3
"""Validate Phase 4B.1 outputs, correction boundary, and deterministic gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/04_phase4b1"
CFG = json.loads((ROOT / "config/phase4b1_analysis.json").read_text(encoding="utf-8"))
CANDIDATES = set(CFG["frozen_candidates"])


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    expected = {
        "hucmsc_proteomics_dataset_registry.csv",
        "candidate_protein_detection_matrix.csv",
        "corrected_hucmsc_protein_source_evidence.csv",
        "protein_identifier_mapping_audit.csv",
        "corrected_phase4b_candidate_evidence_matrix.csv",
        "corrected_phase4b_candidate_classification.csv",
        "phase4b1_risk_flags.csv",
    }
    assert {p.name for p in OUT.glob("*.csv")} == expected
    tables = {name: pd.read_csv(OUT / name) for name in expected}
    for name, table in tables.items():
        assert len(table) > 0, name
        assert "source_url" in table.columns, name
        assert table.source_url.fillna("").str.len().gt(0).all(), name

    for rel, expected_hash in CFG["immutable_upstream_sha256"].items():
        assert digest(ROOT / rel) == expected_hash, rel

    registry = tables["hucmsc_proteomics_dataset_registry.csv"]
    detection = tables["candidate_protein_detection_matrix.csv"]
    source = tables["corrected_hucmsc_protein_source_evidence.csv"]
    mapping = tables["protein_identifier_mapping_audit.csv"]
    evidence = tables["corrected_phase4b_candidate_evidence_matrix.csv"]
    classes = tables["corrected_phase4b_candidate_classification.csv"]

    assert registry.technical_evaluability.eq("TECHNICALLY_EVALUABLE").sum() == 4
    assert len(detection) == 17 * 7
    assert set(source.candidate) == set(mapping.candidate) == set(evidence.candidate) == set(classes.candidate) == CANDIDATES
    assert source.EV_direct.eq("YES").sum() == 15
    assert source.soluble_CM_direct.eq("YES").sum() == 0
    assert set(source.loc[source.EV_direct.eq("NO"), "candidate"]) == {"ADAM17", "FURIN"}
    for gene in ("NID1", "DCN", "ENPP1"):
        row = source[source.candidate.eq(gene)].iloc[0]
        assert row.corrected_protein_source_classification == "HUCMSC_EV_PROTEIN_DIRECT"
    assert source.phase4b_false_negative.eq("YES").sum() == 11
    enpp1 = classes[classes.candidate.eq("ENPP1")].iloc[0]
    assert enpp1.corrected_primary_classification == "TRIANGULATED_HIGH_PRIORITY"
    assert classes.TRIANGULATED_HIGH_PRIORITY.eq("YES").sum() == 1

    # All non-protein evidence columns must equal frozen Phase 4B exactly.
    old = pd.read_csv(ROOT / "results/04_phase4b/integration/phase4b_candidate_evidence_matrix.csv").sort_values("candidate").reset_index(drop=True)
    new = evidence.sort_values("candidate").reset_index(drop=True)
    mutable = {"HUCMSC_PROTEIN_SOURCE", "classification_flags", "source_url"}
    for column in old.columns:
        if column not in mutable:
            assert old[column].fillna("").astype(str).equals(new[column].fillna("").astype(str)), column

    report = (ROOT / "docs/PHASE4B1_PROTEIN_SOURCE_COMPLETENESS_REPORT.md").read_text(encoding="utf-8")
    assert "GO_TO_FINAL_SYNTHESIS_WITH_RESTRICTIONS" in report
    assert "Phase 4A was not rerun" in report
    print("PHASE4B1_VALIDATION_OK: files=7; candidates=17; direct_EV=15; high_priority=1; immutable_dimensions_verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
