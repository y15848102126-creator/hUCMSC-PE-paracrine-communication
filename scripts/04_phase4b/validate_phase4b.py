#!/usr/bin/env python3
"""Validate Phase 4B freeze, tables, conclusions and delivery constraints."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/04_phase4b"
CFG = json.loads((ROOT / "config/phase4b_analysis.json").read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load(rel: str) -> pd.DataFrame:
    path = OUT / rel
    assert path.exists() and path.stat().st_size > 10, f"missing/empty: {path}"
    frame = pd.read_csv(path, keep_default_na=False)
    assert len(frame) > 0, f"no rows: {path}"
    assert "source_url" in frame.columns, f"source_url missing: {path}"
    assert frame.source_url.astype(str).str.strip().ne("").all(), f"blank source_url: {path}"
    return frame


def main() -> int:
    for rel, expected in CFG["phase4a_upstream_sha256"].items():
        assert sha256(ROOT / rel) == expected, f"Phase4A hash changed: {rel}"

    frozen = load("freeze/phase4b_frozen_candidates.csv")
    assert frozen.ligand.tolist() == CFG["frozen_candidates"]
    assert frozen.phase4b_freeze_status.eq("FROZEN_BEFORE_EXTERNAL_EVIDENCE_REVIEW").all()
    assert frozen.manuscript_facing_label.eq("GENERIC_SIGNED_PRIOR_REVERSAL_CANDIDATE").all()
    assert frozen.phase4a_internal_label.eq("TIER_A_DIRECTIONAL_RESCUE_CANDIDATE").all()

    topology = load("topology/candidate_secretion_topology.csv")
    protein = load("protein/hucmsc_protein_source_evidence.csv")
    proteomics = load("protein/public_secretome_proteomics_registry.csv")
    perturb = load("perturbation/empirical_signed_perturbation_evidence.csv")
    disease = load("disease/pe_candidate_context_evidence.csv")
    novelty = load("novelty/direct_msc_pe_overlap.csv")
    searches = load("novelty/literature_search_log.csv")
    matrix = load("integration/phase4b_candidate_evidence_matrix.csv")
    classes = load("integration/phase4b_candidate_classification.csv")
    mixed = load("integration/mixed_direction_stress_test.csv")
    risks = load("qc/phase4b_risk_flags.csv")

    expected = set(CFG["frozen_candidates"])
    for frame, col in [(topology, "candidate"), (protein, "candidate"), (proteomics, "candidate"), (perturb, "candidate"), (disease, "candidate"), (novelty, "candidate"), (matrix, "candidate"), (classes, "candidate"), (mixed, "candidate")]:
        assert set(frame[col]) == expected, f"candidate coverage mismatch: {col}"
    assert len(topology) == len(matrix) == len(classes) == len(mixed) == 17
    assert len(perturb) == 38
    assert not perturb.program_module.eq("PROGRAM_MODULE_10").any(), "Module10 axis was retroactively created"
    assert (perturb.empirical_signed_classification == "EMPIRICAL_REVERSAL_SUPPORTED").sum() == 2
    assert (perturb.empirical_signed_classification == "EMPIRICAL_DISEASE_CONCORDANT").sum() == 2
    assert set(perturb.loc[perturb.empirical_signed_classification == "EMPIRICAL_REVERSAL_SUPPORTED", "candidate"]) == {"ENPP1"}
    assert set(perturb.loc[perturb.empirical_signed_classification == "EMPIRICAL_DISEASE_CONCORDANT", "candidate"]) == {"DCN"}

    direct = protein.loc[protein.protein_source_classification == "HUCMSC_SECRETOME_DIRECT"]
    assert set(direct.candidate) == {"TIMP1", "WNT5A"}
    assert direct.extracellular_compartment.eq("EV_ONLY").all()
    assert classes.TRIANGULATED_HIGH_PRIORITY.eq("NO").all()
    assert classes.TRIANGULATED_CONTEXT_DEPENDENT.eq("NO").all()
    assert (classes.COMPUTATIONAL_ONLY == "YES").sum() == 13
    assert (classes.BIOPHYSICALLY_WEAK_PARACRINE == "YES").sum() == 3
    assert novelty.novelty_classification.isin(CFG["novelty_categories"]).all()
    assert not novelty.novelty_classification.isin(["DIRECT_HUCMSC_PE_MECHANISM_ALREADY_SHOWN", "DIRECT_MSC_PE_MECHANISM_ALREADY_SHOWN"]).any()
    assert searches.search_date.eq(CFG["freeze_date"]).all()
    assert len(searches) >= 109
    assert len(risks) >= 9

    report = (ROOT / "docs/PHASE4B_EXTERNAL_TRIANGULATION_REPORT.md").read_text(encoding="utf-8")
    required_text = ["REVISE_CANDIDATE_FRAMEWORK", "TRIANGULATED_HIGH_PRIORITY", "PROGRAM_MODULE_10", "GENERIC_SIGNED_PRIOR_REVERSAL_CANDIDATE", "Zero"]
    for token in required_text:
        assert token in report, f"report token missing: {token}"
    assert "therapeutically validated" not in report.lower()

    figures = sorted((OUT / "figures").glob("*.png"))
    assert len(figures) >= 3 and all(p.stat().st_size > 1000 for p in figures)
    print("PHASE4B_VALIDATION_OK: freeze=17 axes=38 reversal=2 concordant=2 high=0 context=0 figures=3")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
