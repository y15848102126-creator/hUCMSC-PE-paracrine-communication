#!/usr/bin/env python3
"""Structural and scientific invariants for the Phase 2A.1 stress test."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/02_phase2a1"


def rows(name: str) -> list[dict[str, str]]:
    path = OUT / name
    assert path.exists() and path.stat().st_size > 20, f"missing/empty {path}"
    with path.open(encoding="utf-8-sig", newline="") as handle:
        data = list(csv.DictReader(handle))
    assert data and data[0], f"no data rows in {path}"
    assert "source_url" in data[0], f"source_url missing in {path}"
    assert all(row.get("source_url", "").strip() for row in data), f"blank source_url in {path}"
    return data


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    config = json.loads((ROOT / "config/phase2a1_analysis.json").read_text(encoding="utf-8"))
    for relative, expected in config["frozen_phase2a"]["sha256"].items():
        assert sha256(ROOT / relative) == expected, f"frozen Phase 2A file changed: {relative}"

    clinical = rows("clinical_confounding_audit.csv")
    iugr = rows("shared_program_iugr_sensitivity.csv")
    delivery = rows("shared_program_delivery_sensitivity.csv")
    scores = rows("patient_program_score_sensitivity.csv")
    yang = rows("yang_lope_replication.csv")
    redundancy = rows("shared_program_redundancy_matrix.csv")
    modules = rows("shared_program_modules.csv")
    camera = rows("cameraPR_statistic_audit.csv")
    risks = rows("phase2a1_risk_flags.csv")

    assert len(clinical) == 12
    assert len(iugr) == 40 and all(r["estimability"] == "ESTIMABLE" for r in iugr)
    assert len(delivery) == 80
    assert len(scores) == 40
    assert len(yang) == 20
    assert len(redundancy) == 400
    assert sum(r["record_type"] == "ORIGINAL_GENE_SET" for r in modules) == 20
    assert sum(r["record_type"] == "PROGRAM_MODULE" for r in modules) == 11
    assert len(camera) == 40
    assert len(risks) >= 8

    assert sum(r["direction_agrees"] == "YES" for r in scores if r["contrast"] == "EOPE") == 19
    assert sum(r["direction_agrees"] == "YES" for r in scores if r["contrast"] == "LOPE") == 20
    assert all(r["direction_agrees"] == "YES" for r in iugr)
    assert all(r["estimability"].startswith("NON_ESTIMABLE") for r in delivery if r["contrast"] == "EOPE")
    lope_csection = [r for r in delivery if r["contrast"] == "LOPE" and r["context"] == "C_SECTION_ONLY"]
    assert len(lope_csection) == 20 and all(r["direction_agrees"] == "YES" for r in lope_csection)

    assert sum(r["mapping_status"] == "APPROXIMATE_MATCH" and r["direction_agrees"] == "YES" for r in yang) == 2
    assert sum(r["mapping_status"] == "DIRECT_MATCH" for r in yang) == 7
    assert all(r["direction_agrees"].startswith("NOT_EVALUABLE") for r in yang if r["mapping_status"] == "DIRECT_MATCH")
    assert all(r["replacement_display_name"] == "SIGNED_MEAN_STATISTIC_SD_SCALED" for r in camera)
    assert all("preserved exactly" in r["significance_status"] for r in camera)

    summary = rows_from_interim(ROOT / "data/interim/phase2a1/count_layer_provenance_summary.csv")
    assert summary[0]["resolution"] == "RESOLVED_EXPRESSION_IS_CEILED_LIBRARY_SIZE_NORMALIZED_NOT_RAW_UMI"
    assert summary[0]["ceiling_10000_consistent_n"] == summary[0]["cell_n"] == "86752"

    report = ROOT / "docs/PHASE2A1_RECEIVER_PROGRAM_STRESS_TEST_REPORT.md"
    assert report.exists() and "**REVISE_RECEIVER_FRAMEWORK**" in report.read_text(encoding="utf-8")
    print("Phase 2A.1 validation passed: frozen Phase 2A hashes, 9 CSVs, and scientific invariants verified")
    return 0


def rows_from_interim(path: Path) -> list[dict[str, str]]:
    assert path.exists(), f"missing {path}"
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


if __name__ == "__main__":
    raise SystemExit(main())
