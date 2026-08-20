#!/usr/bin/env python3
"""Hash-check and freeze the exact 17 Phase 4A Tier A candidates."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    cfg = json.loads((ROOT / "config/phase4b_analysis.json").read_text(encoding="utf-8"))
    assert cfg["external_evidence_review_status_at_freeze"] == "NOT_STARTED"
    for filename, expected in cfg["phase4a_upstream_sha256"].items():
        assert sha256(ROOT / filename).lower() == expected.lower(), filename
    hierarchy = pd.read_csv(ROOT / "results/04_phase4a/integration/phase4a_candidate_hierarchy.csv")
    frozen = hierarchy[hierarchy.best_phase4a_tier.eq("TIER_A_DIRECTIONAL_RESCUE_CANDIDATE")].copy()
    assert sorted(frozen.ligand) == sorted(cfg["frozen_candidates"])
    assert len(frozen) == 17
    assert "paracrine_scope" in frozen.columns
    assert frozen.paracrine_scope.eq("P1_PARACRINE_CORE").all()
    frozen["phase4a_internal_label"] = frozen.best_phase4a_tier
    frozen["manuscript_facing_label"] = cfg["manuscript_facing_phase4a_label"]
    frozen["phase4b_freeze_status"] = "FROZEN_BEFORE_EXTERNAL_EVIDENCE_REVIEW"
    frozen["phase4a_hierarchy_sha256"] = cfg["phase4a_upstream_sha256"]["results/04_phase4a/integration/phase4a_candidate_hierarchy.csv"]
    frozen["source_url"] = "results/04_phase4a/integration/phase4a_candidate_hierarchy.csv|config/phase4b_analysis.json"
    out = ROOT / "results/04_phase4b/freeze"
    out.mkdir(parents=True, exist_ok=True)
    frozen.sort_values("ligand").to_csv(out / "phase4b_frozen_candidates.csv", index=False)
    print("PHASE4B_FREEZE_OK: 17 candidates; Phase4A hashes verified; external review not started")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
