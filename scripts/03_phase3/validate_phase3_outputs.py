#!/usr/bin/env python3
"""Independent structural and mechanical validation of Phase 3 outputs."""

from __future__ import annotations

import csv
import math
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "03_phase3"


def rows(relative: str) -> list[dict[str, str]]:
    path = OUT / relative
    assert path.is_file() and path.stat().st_size > 0, f"missing/empty: {path}"
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def num(value: str) -> float:
    return float(value) if value not in ("", None) else math.nan


def main() -> int:
    required = [
        "metadata/hucmsc_donor_registry.csv", "metadata/hucmsc_dataset_role_registry.csv",
        "metadata/licensing_strata_registry.csv", "ligand_universe/frozen_ligand_universe.csv",
        "baseline/baseline_ligand_expression_by_donor.csv", "baseline/baseline_sender_robustness.csv",
        "baseline/cross_dataset_sender_concordance.csv", "licensing/licensing_effect_by_donor_passage.csv",
        "licensing/licensing_effect_by_donor.csv", "licensing/licensing_ligand_classification.csv",
        "licensing/licensing_programs.csv", "sender/sender_evidence_hierarchy.csv",
        "sender/frozen_phase4_sender_candidates.csv", "qc/phase3_qc_summary.csv", "qc/phase3_risk_flags.csv",
    ]
    for name in required:
        rows(name)
    assert (ROOT / "docs" / "PHASE3_HUCMSC_SENDER_PROGRAM_REPORT.md").is_file()
    assert (ROOT / "docs" / "PHASE3_INDEPENDENT_SENDER_ANALYSIS_PLAN.md").is_file()

    donor = rows("metadata/hucmsc_donor_registry.csv")
    sender_donors = [r for r in donor if r["is_hucmsc_sender_donor"] == "YES"]
    assert len(sender_donors) == 8
    assert Counter(r["dataset"] for r in sender_donors) == {"GSE182158": 2, "GSE199071": 4, "GSE117837": 2}
    assert all(r["biological_unit"] == "DONOR" for r in sender_donors)

    universe = rows("ligand_universe/frozen_ligand_universe.csv")
    genes = [r["gene"] for r in universe]
    assert len(genes) == len(set(genes)) == 1226
    assert all(r["freeze_status"] == "FROZEN_BEFORE_EXPRESSION_OUTCOMES" for r in universe)

    baseline = rows("baseline/baseline_ligand_expression_by_donor.csv")
    assert len(baseline) == 1226 * 15
    assert set(r["biological_unit"] for r in baseline) == {"DONOR"}
    core = [r for r in baseline if r["analysis_role"] == "CORE_BASELINE"]
    assert len(core) == 1226 * 6
    assert len(set((r["dataset"], r["donor_id"]) for r in core)) == 6
    assert set(r["expression_layer"] for r in baseline) == {"TRUE_CELLRANGER_RAW_COUNTS"}

    robust = rows("baseline/baseline_sender_robustness.csv")
    assert len(robust) == 1226
    for r in robust:
        s182, s199 = int(r["gse182_supported_donor_n"]), int(r["gse199_supported_donor_n"])
        l182, l199 = int(r["gse182_low_detected_donor_n"]), int(r["gse199_low_detected_donor_n"])
        if s182 == 2 and s199 >= 3:
            expected = "ROBUST_BASELINE_SENDER"
        elif l182 >= 1 and l199 >= 1:
            expected = "MULTIDATASET_LOW_EXPRESSION"
        elif (s182 >= 2 and l199 == 0) or (s199 >= 2 and l182 == 0):
            expected = "DATASET_SPECIFIC"
        elif s182 + s199 >= 2:
            expected = "DONOR_VARIABLE"
        else:
            expected = "NOT_RELIABLY_EXPRESSED"
        assert r["baseline_classification"] == expected, (r["gene"], expected)
    assert Counter(r["baseline_classification"] for r in robust)["ROBUST_BASELINE_SENDER"] == 214

    strata = rows("metadata/licensing_strata_registry.csv")
    valid = [r for r in strata if r["valid_within_stratum_contrast"] == "YES"]
    assert [r["stratum_id"] for r in valid] == ["Donor1_P5", "Donor2_P2", "Donor2_P5"]
    assert next(r for r in strata if r["stratum_id"] == "Donor2_P0")["valid_within_stratum_contrast"] == "NO_UNPAIRED"
    passage = rows("licensing/licensing_effect_by_donor_passage.csv")
    assert len(passage) == 1226 * 3
    assert set(r["biological_unit"] for r in passage) == {"DONOR_PASSAGE_STRATUM_NOT_CELL"}
    donor_effect = rows("licensing/licensing_effect_by_donor.csv")
    assert len(donor_effect) == 1226 * 2
    assert set(r["biological_unit"] for r in donor_effect) == {"DONOR"}

    licensing = rows("licensing/licensing_ligand_classification.csv")
    assert len(licensing) == 1226
    assert Counter(r["licensing_classification"] for r in licensing) == {
        "NO_CLEAR_LICENSING_EFFECT": 778, "PASSAGE_DEPENDENT": 281,
        "DONOR_DEPENDENT": 65, "LICENSING_CONSISTENT_UP": 63,
        "LICENSING_CONSISTENT_DOWN": 39,
    }

    hierarchy = rows("sender/sender_evidence_hierarchy.csv")
    assert len(hierarchy) == 1226
    for r in hierarchy:
        if (r["baseline_classification"] == "ROBUST_BASELINE_SENDER"
                and r["licensing_classification"] in {"LICENSING_CONSISTENT_UP", "LICENSING_CONSISTENT_DOWN"}
                and r["retained_after_stimulation_all_strata"] == "YES"):
            expected = "S1"
        elif r["baseline_classification"] == "ROBUST_BASELINE_SENDER":
            expected = "S2"
        elif r["baseline_classification"] in {"MULTIDATASET_LOW_EXPRESSION", "DATASET_SPECIFIC", "DONOR_VARIABLE"}:
            expected = "S3"
        else:
            expected = "S4"
        assert r["sender_evidence_level"] == expected, (r["gene"], expected)
        assert r["ligand_claim"] == "PUTATIVE_TRANSCRIPTOMIC_LIGAND"
        assert "THERAPEUTIC" not in r["sender_competence"]
    assert Counter(r["sender_evidence_level"] for r in hierarchy) == {"S4": 793, "S3": 219, "S2": 176, "S1": 38}
    candidates = rows("sender/frozen_phase4_sender_candidates.csv")
    assert len(candidates) == 214
    assert set(r["sender_evidence_level"] for r in candidates) == {"S1", "S2"}
    assert all(r["candidate_label"] == "HIGH_CONFIDENCE_SENDER_CANDIDATE_NOT_THERAPEUTIC" for r in candidates)

    programs = rows("licensing/licensing_programs.csv")
    assert len(programs) == 5316
    assert set(r["MSigDB_version"] for r in programs) == {"2026.1.Hs"}
    assert all(r["interpretation"] == "HUCMSC_INFLAMMATORY_LICENSING_ONLY_NOT_PE_THERAPEUTIC_PROGRAM" for r in programs)

    for name in required:
        table = rows(name)
        source_cols = [c for c in table[0] if c.startswith("source_url")]
        assert source_cols, f"no provenance column: {name}"
        assert all(any(r.get(c, "") for c in source_cols) for r in table), f"blank provenance: {name}"

    for script in (ROOT / "scripts" / "03_phase3").glob("*"):
        if script.suffix.lower() not in {".r", ".py", ".ps1"}:
            continue
        text = script.read_text(encoding="utf-8", errors="replace").lower()
        if script.name != Path(__file__).name:
            assert "results/02_phase2" not in text, f"receiver path in sender script: {script.name}"
            assert "nichenet ligand activity" not in text
            assert "cellchat" not in text

    figures = list((OUT / "figures").glob("*.png"))
    assert len(figures) >= 5 and all(p.stat().st_size > 10_000 for p in figures)
    print("PHASE3_VALIDATION_OK: 15 required tables; 1226 ligands; 8 sender donors; 214 S1/S2 candidates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
