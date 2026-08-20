#!/usr/bin/env python3
"""Validate Phase 0B registries and high-risk design invariants."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "00_dataset_audit_phase0b"
REQUIRED = {
    "pe_scrna_extended_registry.csv": {"study_id", "pe_donor_n", "control_donor_n", "recommended_role", "source_url", "source_accession"},
    "pe_scrna_data_access.csv": {"study_id", "component", "access_class", "source_url", "source_accession"},
    "phs001886_version_overlap.csv": {"earlier_version", "later_version", "shared_subject_n", "shared_sample_n", "independence_conclusion"},
    "hucmsc_sender_redundancy_registry.csv": {"dataset", "independent_donor_n", "passage", "donor_identity", "source_url"},
    "phase0b_risk_flags.csv": {"risk_id", "scope", "severity", "phase1_blocking", "source_url"},
    "revised_dataset_roles.csv": {"dataset", "phase0a_role", "phase0b_role", "change", "source_url"},
}
ALLOWED_PE_ROLES = {"PRIMARY_PE_SCRNA", "SCRNA_REPLICATION", "PENDING", "EXCLUDE"}


def read_csv(name: str) -> tuple[list[str], list[dict[str, str]]]:
    path = OUT / name
    assert path.exists() and path.stat().st_size > 0, f"missing/empty {name}"
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        rows = list(reader)
    assert rows, f"no rows in {name}"
    assert len(fields) == len(set(fields)), f"duplicate columns in {name}"
    assert REQUIRED[name] <= set(fields), f"missing columns in {name}: {REQUIRED[name] - set(fields)}"
    return fields, rows


def main() -> int:
    tables = {name: read_csv(name)[1] for name in REQUIRED}
    pe = {row["study_id"]: row for row in tables["pe_scrna_extended_registry.csv"]}
    assert pe["Admati_2023_FIGSHARE"]["recommended_role"] == "PRIMARY_PE_SCRNA"
    assert pe["PMID41472684_phs001886v6"]["recommended_role"] == "PENDING"
    assert pe["GSE290578"]["recommended_role"] == "SCRNA_REPLICATION"
    assert pe["Jiao2023"]["exclusion_reason"] == "EXCLUDE_DATA_UNAVAILABLE"
    assert {row["recommended_role"] for row in pe.values()} <= ALLOWED_PE_ROLES
    for name, rows in tables.items():
        for index, row in enumerate(rows, 2):
            source_fields = [value for key, value in row.items() if key.startswith("source_")]
            assert any(value and value not in {"NOT_FOUND", "UNRESOLVED"} for value in source_fields), f"missing source in {name}:{index}"

    overlaps = tables["phs001886_version_overlap.csv"]
    assert len(overlaps) == 15, "expected all 15 version pairs"
    for row in overlaps:
        assert row["independence_conclusion"] == "NOT_INDEPENDENT_CUMULATIVE"
        assert int(row["retired_subject_n"]) == 0 and int(row["retired_sample_n"]) == 0
        assert int(row["shared_subject_n"]) == int(row["earlier_subject_n"])
        assert int(row["shared_sample_n"]) == int(row["earlier_sample_n"])

    huc = {row["dataset"]: row for row in tables["hucmsc_sender_redundancy_registry.csv"]}
    assert "2 UC donors" in huc["GSE182158"]["independent_donor_n"]
    assert "4 HUCMSC donors" in huc["GSE199071"]["independent_donor_n"]
    assert "D2 P0 has no stimulated match" in huc["GSE117837_LICENSED"]["passage"]

    with (ROOT / "results" / "00_dataset_audit" / "proposed_dataset_roles.csv").open(encoding="utf-8-sig", newline="") as handle:
        phase0a_roles = {row["geo_accession"]: row for row in csv.DictReader(handle)}
    assert phase0a_roles["GSE272342"]["proposed_role"] == "SUBTYPE_VALIDATION"
    assert "DESIGN_SPECIFIC_SENSITIVITY" in phase0a_roles["GSE272342"]["mandatory_restrictions"]

    report = (ROOT / "docs" / "DATASET_AUDIT_PHASE0B_REPORT.md").read_text(encoding="utf-8")
    for phrase in ["GO_TO_PHASE1_WITH_RESTRICTIONS", "PRIMARY_PE_SCRNA", "PREPRINT", "GSE298602", "NORMAL_GESTATIONAL_REFERENCE", "DESIGN_SPECIFIC_SENSITIVITY"]:
        assert phrase in report, f"report missing {phrase}"

    matrix_zip = ROOT / "data" / "raw" / "phase0b" / "sc_PE_allcells_with_metadata_29-May-2023.txt.zip"
    if matrix_zip.exists():
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "00_dataset_registry_phase0b" / "inspect_admati_processed.py"), "--zip", str(matrix_zip)],
            check=True, capture_output=True, text=True, encoding="utf-8",
        )
        inspected = json.loads(result.stdout)
        assert inspected["cell_columns"] == 86752
        assert len(inspected["metadata"]["donorID"]) == 26

    for forbidden in [ROOT / "results" / "01_deg", ROOT / "results" / "01_cellchat", ROOT / "results" / "01_nichenet"]:
        assert not forbidden.exists(), f"Phase 1 output unexpectedly exists: {forbidden}"
    print("Phase 0B validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
