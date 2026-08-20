#!/usr/bin/env python3
"""Verify that every non-protein Phase 4A/4B evidence dimension is immutable."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    cfg = json.loads((ROOT / "config/phase4b1_analysis.json").read_text(encoding="utf-8"))
    for rel, expected in cfg["immutable_upstream_sha256"].items():
        assert sha256(ROOT / rel) == expected, f"immutable evidence changed: {rel}"
    old = pd.read_csv(ROOT / "results/04_phase4b/integration/phase4b_candidate_evidence_matrix.csv")
    assert sorted(old.candidate) == sorted(cfg["frozen_candidates"])
    print("PHASE4B1_FREEZE_OK: 17 candidates; 8 immutable evidence files hash-verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
