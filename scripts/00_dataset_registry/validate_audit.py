#!/usr/bin/env python3
"""Fail fast on Phase 0 structural, count, provenance, and leakage invariants."""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "00_dataset_audit"
CFG = json.loads((ROOT / "config" / "audit_accessions.json").read_text(encoding="utf-8"))
ALLOWED_ROLES = {
    "PRIMARY_BULK", "BULK_REPLICATION", "SUBTYPE_VALIDATION", "PRIMARY_PE_SCRNA",
    "SCRNA_REPLICATION", "PRIMARY_HUCMSC_ATLAS", "HUCMSC_LICENSING",
    "SUPPLEMENTARY", "EXCLUDE", "PENDING",
}


def read_csv(name: str) -> list[dict[str, str]]:
    path = OUT / name
    assert path.exists() and path.stat().st_size > 0, f"missing/empty {path}"
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def require_sources(rows: list[dict[str, str]], name: str) -> None:
    for index, row in enumerate(rows, 2):
        assert row.get("source_url", "").strip(), f"{name}:{index} missing source_url"
        assert row.get("source_accession", "").strip(), f"{name}:{index} missing source_accession"


def main() -> int:
    datasets = read_csv("dataset_registry.csv")
    samples = read_csv("sample_registry.csv")
    overlaps = read_csv("dataset_overlap_matrix.csv")
    risks = read_csv("dataset_risk_flags.csv")
    availability = read_csv("data_availability.csv")
    roles = read_csv("proposed_dataset_roles.csv")

    expected = set(CFG["main_accessions"])
    assert len(datasets) == 13, f"expected 13 dataset rows, got {len(datasets)}"
    assert {row["geo_accession"] for row in datasets} == expected
    assert len(roles) == 13 and {row["geo_accession"] for row in roles} == expected
    assert len(availability) == 13 and {row["geo_accession"] for row in availability} == expected
    assert all(row["proposed_role"] in ALLOWED_ROLES for row in datasets)
    assert all(row["proposed_role"] in ALLOWED_ROLES for row in roles)

    for rows, name in [(datasets, "dataset_registry"), (samples, "sample_registry"),
                       (overlaps, "dataset_overlap_matrix"), (risks, "dataset_risk_flags"),
                       (availability, "data_availability"), (roles, "proposed_dataset_roles")]:
        require_sources(rows, name)

    direct = Counter(row["dataset_accession"] for row in samples if row["analytical_membership"] == "direct_gsm")
    expected_direct = {row["geo_accession"]: int(row["geo_gsm_count"]) for row in datasets}
    assert direct == Counter(expected_direct), f"direct GSM counts differ: {direct - Counter(expected_direct)}"

    gse75010 = [row for row in overlaps if row["dataset_a"] == "GSE75010" and row["overlap_type"] == "exact GSM reanalysis"]
    assert len(gse75010) == 7
    assert sum(int(row["overlap_n"]) for row in gse75010) == 173
    expected_overlap_counts = {"GSE30186": 12, "GSE10588": 43, "GSE24129": 16, "GSE25906": 60, "GSE43942": 12, "GSE4707": 14, "GSE44711": 16}
    assert {row["dataset_b"]: int(row["overlap_n"]) for row in gse75010} == expected_overlap_counts

    g234 = [row for row in samples if row["dataset_accession"] == "GSE234729"]
    assert sum(row["analytical_membership"] == "direct_gsm" for row in g234) == 111
    assert sum(row["analytical_membership"] == "reused_external_gsm" for row in g234) == 12
    assert len({row["sample_accession"] for row in g234}) == 123
    assert sum("severe" in row["disease_group"].lower() for row in g234) == 50
    assert sum("control" in row["disease_group"].lower() for row in g234) == 73

    g173 = [row for row in samples if row["dataset_accession"] == "GSE173193" and row["analytical_membership"] == "direct_gsm"]
    groups = Counter(row["disease_group"] for row in g173)
    assert groups == Counter({"control group": 2, "GDM group": 2, "PE group": 2, "elderly group": 2})

    g282 = [row for row in samples if row["dataset_accession"] == "GSE282038"]
    assert sum(row["analytical_membership"] == "direct_gsm" for row in g282) == 1
    assert sum(row["analytical_membership"] == "external_series_reference" for row in g282) == 4
    assert {row["origin_gse"] for row in g282} == {"GSE282038", "GSE267340", "GSE298119"}

    g329 = [row for row in samples if row["dataset_accession"] == "GSE329173"]
    assert len(g329) == 3 and all("eclamps" in row["disease_group"].lower() for row in g329)

    g182 = [row for row in samples if row["dataset_accession"] == "GSE182158"]
    expected_donors = {
        "A01": ("38", "Female"), "A02": ("46", "Female"), "A03": ("32", "Female"),
        "B01": ("33", "Male"), "B02": ("43", "Male"), "B03": ("27", "Female"),
        "D01": ("33", "Male"), "D02": ("28", "Male"), "D03": ("22", "Male"),
        "U01": ("28", "Female"), "U02": ("37", "Female"),
    }
    assert len(g182) == 11
    assert {row["donor_id"]: (row["donor_age"], row["donor_sex"]) for row in g182} == expected_donors
    assert all(row["passage"] == "P1 or P2; exact donor mapping UNRESOLVED" for row in g182)

    g117 = [row for row in samples if row["dataset_accession"] == "GSE117837"]
    observed = Counter(row["treatment"] for row in g117)
    assert observed == Counter({"naïve": 203, "IFNγ+TNFα": 158}), observed
    assert len(g117) == 361
    assert len({row["donor_id"] for row in g117}) == 2

    assert any(row["flag_id"] == "OVERLAP_7_SOURCE_GSE" and row["severity"] == "CRITICAL" for row in risks)
    assert any(row["flag_id"] == "NO_INTERNAL_CONTROL" and row["geo_accession"] == "GSE282038" for row in risks)
    assert any(row["flag_id"] == "CASE_ONLY" and row["geo_accession"] == "GSE329173" for row in risks)
    assert any(row["flag_id"] == "PASSAGE_MAPPING_GAP" and row["geo_accession"] == "GSE182158" for row in risks)

    report = ROOT / "docs" / "DATASET_AUDIT_REPORT.md"
    text = report.read_text(encoding="utf-8")
    assert "**GO_WITH_MODIFICATIONS**" in text
    assert "PRIMARY_PE_SCRNA" not in {row["proposed_role"] for row in roles}
    for accession in expected:
        assert accession in text, f"report missing {accession}"

    print(f"PASS: 13 datasets, {len(samples)} sample/dependency rows, {len(overlaps)} overlap rows, {len(risks)} risk flags")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except AssertionError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
