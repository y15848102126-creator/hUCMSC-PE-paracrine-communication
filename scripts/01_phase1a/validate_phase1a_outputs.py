#!/usr/bin/env python3
"""Validate Phase 1A freeze tables, design invariants, and phase-boundary safety."""

from __future__ import annotations

import csv
import hashlib
import subprocess
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "01_phase1a"
REQUIRED = {
    "bulk_sample_freeze.csv": {"dataset", "GSM/sample ID", "patient/pregnancy ID if known", "PE/control", "PE subtype", "GA", "delivery GA", "FGR", "singleton/twin", "fetal sex", "ancestry", "platform", "batch", "include_phase1b", "exclusion_reason", "independence_group", "source"},
    "bulk_cohort_registry.csv": {"dataset", "role_phase1b", "frozen_independent_unit_n", "phase1b_main_pe_n", "phase1b_main_control_n", "source"},
    "bulk_overlap_freeze.csv": {"dataset_a", "dataset_b", "overlap_n", "independence_decision", "source"},
    "gene_mapping_registry.csv": {"original_probe_id", "original_gene_id", "original_symbol", "mapped_symbol", "mapped_ensembl_gene_id", "mapping_source", "mapping_version", "source"},
    "bulk_processing_registry.csv": {"dataset", "expression_unit", "prespecified_qc_transform", "gene_matrix_path", "gene_matrix_sha256", "source"},
    "bulk_exclusion_log.csv": {"dataset", "sample_id", "exclusion_scope", "exclusion_reason", "source"},
    "phase1a_risk_flags.csv": {"risk_id", "severity", "dataset", "risk", "required_action", "status", "source"},
    "bulk_sample_qc.csv": {"dataset", "sample_id", "qc_flag", "action", "method", "source"},
}


def read_csv(name: str) -> list[dict[str, str]]:
    path = OUT / name
    assert path.exists() and path.stat().st_size > 0, f"missing/empty {path}"
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        rows = list(reader)
    assert len(fields) == len(set(fields)), f"duplicate columns in {name}"
    assert REQUIRED[name] <= set(fields), f"missing columns in {name}: {REQUIRED[name] - set(fields)}"
    assert rows, f"no rows in {name}"
    for line, row in enumerate(rows, 2):
        assert row.get("source") not in {None, "", "UNRESOLVED"}, f"missing source in {name}:{line}"
    return rows


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    tables = {name: read_csv(name) for name in REQUIRED}
    freeze = tables["bulk_sample_freeze.csv"]
    assert len(freeze) == 538
    assert len({(row["dataset"], row["GSM/sample ID"]) for row in freeze}) == len(freeze)
    assert {row["include_phase1b"] for row in freeze} <= {"YES", "NO"}

    expected = {
        "GSE234729_ANALYTIC_123": (123, 50, 73, 0, 0),
        "GSE75010_BIOBANK": (157, 80, 77, 63, 53),
        "GSE30186": (12, 6, 6, 6, 6),
        "GSE10588": (43, 17, 26, 17, 26),
        "GSE24129": (24, 8, 8, 8, 8),
        "GSE25906": (60, 23, 37, 23, 37),
        "GSE43942": (12, 5, 7, 5, 7),
        "GSE4707": (14, 10, 4, 0, 0),
        "GSE44711": (16, 8, 8, 0, 0),
        "GSE190639": (32, 19, 13, 0, 0),
        "GSE272342": (45, 20, 9, 0, 0),
    }
    grouped = defaultdict(list)
    for row in freeze:
        grouped[row["dataset"]].append(row)
    for dataset, (total, pe, control, main_pe, main_control) in expected.items():
        rows = grouped[dataset]
        assert len(rows) == total, (dataset, len(rows), total)
        units = defaultdict(set)
        main = defaultdict(set)
        for row in rows:
            units[row["PE/control"]].add(row["independence_group"])
            if row["include_phase1b"] == "YES":
                main[row["PE/control"]].add(row["independence_group"])
        assert (len(units["PE"]), len(units["CONTROL"])) == (pe, control), dataset
        assert (len(main["PE"]), len(main["CONTROL"])) == (main_pe, main_control), dataset

    g234 = grouped["GSE234729_ANALYTIC_123"]
    assert Counter(row["source_provenance"] for row in g234) == Counter({"GSE234729 direct GSM": 111, "GSE186257 reanalysis included once in GSE234729 analytical matrix": 12})
    assert len({row["GSM/sample ID"] for row in g234}) == 123
    assert all(row["batch"] == "UNRESOLVED" and row["GA"] == "UNRESOLVED" for row in g234)

    assert all(row["dataset"] != "GSE75010_COMPOSITE_330" for row in freeze)
    biobank_ch = [row for row in grouped["GSE75010_BIOBANK"] if "-CH" in row["sample_title"]]
    assert len(biobank_ch) == 41 and all(row["include_phase1b"] == "NO" for row in biobank_ch)
    fgr_only = [row for row in grouped["GSE24129"] if row["PE/control"] == "OTHER"]
    assert len(fgr_only) == 8 and all(row["include_phase1b"] == "NO" for row in fgr_only)
    g272 = grouped["GSE272342"]
    assert len({row["independence_group"] for row in g272}) == 29
    assert all(row["include_phase1b"] == "NO" for row in g272)

    excluded = {(row["dataset"], row["sample_id"]): row for row in tables["bulk_exclusion_log.csv"]}
    no_rows = [row for row in freeze if row["include_phase1b"] == "NO"]
    assert len(excluded) == len(no_rows)
    assert all((row["dataset"], row["GSM/sample ID"]) in excluded for row in no_rows)
    assert all(row["exclusion_reason"] == "NONE" for row in freeze if row["include_phase1b"] == "YES")

    overlaps = tables["bulk_overlap_freeze.csv"]
    assert sum(int(row["overlap_n"]) for row in overlaps if row["dataset_a"] == "GSE75010_COMPOSITE_330") == 330
    assert any(row["dataset_a"] == "GSE234729_ANALYTIC_123" and row["dataset_b"] == "GSE186257" and row["overlap_n"] == "12" for row in overlaps)
    assert any(row["dataset_a"] == "GSE272342" and row["overlap_n"] == "3" for row in overlaps)

    processing = tables["bulk_processing_registry.csv"]
    assert len(processing) == 11
    for row in processing:
        if row["gene_matrix_path"].startswith("data/"):
            path = ROOT / row["gene_matrix_path"]
            assert path.exists() and sha256(path) == row["gene_matrix_sha256"]
    mappings = tables["gene_mapping_registry.csv"]
    assert len(mappings) > 250000
    assert all(row["collapse_rule"] in {"MEDIAN_AFTER_TRANSFORM_NO_OUTCOME", "GENE_INPUT_MEDIAN_IF_DUPLICATE_NO_OUTCOME"} for row in mappings)

    qc = tables["bulk_sample_qc.csv"]
    assert len(qc) == 493
    assert all(row["action"] == "FLAG_ONLY_NO_AUTOMATIC_EXCLUSION" for row in qc if row["qc_flag"] != "PASS")
    qc_flags = {(row["dataset"], row["sample_id"]) for row in qc if row["qc_flag"] != "PASS"}
    excluded_ids = {(row["dataset"], row["sample_id"]) for row in tables["bulk_exclusion_log.csv"]}
    assert not (qc_flags - excluded_ids) or all(
        next(x for x in freeze if x["dataset"] == key[0] and x["GSM/sample ID"] == key[1])["include_phase1b"] == "YES"
        for key in qc_flags - excluded_ids
    ), "QC flags may remain included; they must not be silently deleted"

    for doc in ["PHASE1B_STATISTICAL_ANALYSIS_PLAN.md", "PHASE1A_BULK_DATA_FREEZE_REPORT.md", "PHASE1A_CHANGELOG.md"]:
        path = ROOT / "docs" / doc
        assert path.exists() and path.stat().st_size > 0, f"missing {doc}"
    report = (ROOT / "docs" / "PHASE1A_BULK_DATA_FREEZE_REPORT.md").read_text(encoding="utf-8")
    plan = (ROOT / "docs" / "PHASE1B_STATISTICAL_ANALYSIS_PLAN.md").read_text(encoding="utf-8")
    for phrase in ["GO_TO_PHASE1B_WITH_RESTRICTIONS", "GSE75010_BIOBANK", "GSE234729_ANALYTIC_123", "GSE272342", "No Phase 1B analysis was run"]:
        assert phrase in report
    for phrase in ["random-effects", "I²", "leave-one-cohort-out", "at least four independent core cohorts", "ComBat"]:
        assert phrase in plan

    for phase0b in ["pe_scrna_extended_registry.csv", "pe_scrna_data_access.csv", "phs001886_version_overlap.csv", "hucmsc_sender_redundancy_registry.csv", "phase0b_risk_flags.csv", "revised_dataset_roles.csv"]:
        assert (ROOT / "results" / "00_dataset_audit_phase0b" / phase0b).exists()

    tracked = subprocess.run(
        ["git", "-c", f"safe.directory={ROOT.as_posix()}", "ls-files", "data/raw", "data/interim"],
        cwd=ROOT, check=True, capture_output=True, text=True,
    ).stdout.splitlines()
    assert set(tracked) <= {"data/raw/.gitkeep", "data/interim/.gitkeep"}, tracked
    for forbidden in [ROOT / "results" / "01_deg", ROOT / "results" / "01_cellchat", ROOT / "results" / "01_nichenet", ROOT / "results" / "01_wgcna", ROOT / "results" / "01_ml"]:
        assert not forbidden.exists(), f"Phase 1B/later output unexpectedly exists: {forbidden}"
    print("Phase 1A validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
