#!/usr/bin/env python3
"""Validate Phase 1A.1 registries, matrices, preregistration and phase boundary."""

from __future__ import annotations

import csv
import hashlib
import subprocess
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "01_phase1a1"
REQUIRED = {
    "preprocessing_comparison_registry.csv": {"dataset", "common_feature_n", "common_sample_n", "median_sample_spearman", "formal_selection", "source"},
    "formal_phase1b_matrix_registry.csv": {"dataset", "formal_matrix_choice", "feature_matrix_path", "feature_matrix_sha256", "gene_matrix_path", "gene_matrix_sha256", "arbitrary_shift_log_in_formal_matrix", "source"},
    "qc_flag_reinterpretation.csv": {"dataset", "sample_id", "batch", "retain_phase1b", "phase1a1_interpretation", "source"},
    "cohort_model_formula_registry.csv": {"dataset", "primary_formula", "mandatory_covariates", "fallback_formula", "estimability_checks", "source"},
    "phase1a1_risk_flags.csv": {"risk_id", "severity", "dataset", "risk", "required_action", "status", "source"},
}


def rows(name: str) -> list[dict[str, str]]:
    path = OUT / name
    assert path.exists() and path.stat().st_size > 0, f"missing/empty {path}"
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        data = list(reader)
    assert len(fields) == len(set(fields)), f"duplicate columns: {name}"
    assert REQUIRED[name] <= set(fields), f"missing columns: {name}"
    assert data and all(row.get("source") for row in data), f"empty rows/source: {name}"
    return data


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    tables = {name: rows(name) for name in REQUIRED}
    comparison = {row["dataset"]: row for row in tables["preprocessing_comparison_registry.csv"]}
    assert set(comparison) == {"GSE30186", "GSE10588", "GSE43942"}
    assert {d: int(r["common_sample_n"]) for d, r in comparison.items()} == {"GSE30186": 12, "GSE10588": 43, "GSE43942": 12}
    assert all(float(row["median_sample_spearman"]) > 0.97 for row in comparison.values())
    assert comparison["GSE30186"]["formal_selection"] == "RECONSTRUCTED_FORMAL"
    assert all(row["comparison_is_outcome_blind"] == "YES_NO_DISEASE_LABEL_USED" for row in comparison.values())

    formal = tables["formal_phase1b_matrix_registry.csv"]
    assert len(formal) == 6
    expected = {"GSE75010_BIOBANK", "GSE30186", "GSE10588", "GSE24129", "GSE25906", "GSE43942"}
    assert {row["dataset"] for row in formal} == expected
    assert all(row["arbitrary_shift_log_in_formal_matrix"] == "NO" for row in formal)
    assert all(row["cross_cohort_batch_correction"] == "PROHIBITED" for row in formal)
    for row in formal:
        for path_col, hash_col in [("feature_matrix_path", "feature_matrix_sha256"), ("gene_matrix_path", "gene_matrix_sha256")]:
            path = ROOT / row[path_col]
            assert path.exists(), path
            assert sha256(path) == row[hash_col], (path, sha256(path), row[hash_col])

    qc = tables["qc_flag_reinterpretation.csv"]
    assert len(qc) == 18
    assert Counter(row["dataset"] for row in qc) == Counter({"GSE25906": 15, "GSE30186": 2, "GSE75010_BIOBANK": 1})
    assert Counter(row["PE_control"] for row in qc if row["dataset"] == "GSE25906") == Counter({"CONTROL": 8, "PE": 7})
    assert all(row["batch"] == "A" for row in qc if row["dataset"] == "GSE25906")
    assert all(row["retain_phase1b"] == "YES" and row["exclusion_reason"].startswith("NONE") for row in qc)

    models = {row["dataset"]: row for row in tables["cohort_model_formula_registry.csv"]}
    assert set(models) == expected
    assert models["GSE25906"]["primary_formula"] == "~ disease + batch + GA_c + fetal_sex"
    assert models["GSE75010_BIOBANK"]["primary_formula"] == "~ disease + GA_c + fetal_sex"
    assert all(row["primary_design_rank"] == row["primary_design_column_n"] for row in models.values())
    assert all(float(row["primary_max_vif"]) < 5 for row in models.values())
    assert all(row["primary_estimability_decision"] == "FULL_RANK_NO_ZERO_CELL_VIF_LT5" for row in models.values())
    assert all(row["outcome_blind_formula_lock"] == "YES_BEFORE_DEG" for row in models.values())

    plan = (ROOT / "docs" / "PHASE1B_STATISTICAL_ANALYSIS_PLAN.md").read_text(encoding="utf-8")
    amendment = (ROOT / "docs" / "PHASE1A1_PREPROCESSING_AMENDMENT.md").read_text(encoding="utf-8")
    report = (ROOT / "docs" / "PHASE1A_BULK_DATA_FREEZE_REPORT.md").read_text(encoding="utf-8")
    for phrase in ["at least four independent core cohorts", "75%", "I² ≤ 60%", "Every valid leave-one-cohort-out", "80%", "no universal absolute 0.25", "Hedges' g"]:
        assert phrase in plan, phrase
    assert "The absolute random-effects summary log2 fold change is at least 0.25" not in plan
    for phrase in ["GO_TO_PHASE1B_WITH_RESTRICTIONS", "All 18", "GSE30186", "GSE43942", "Phase 1B has not been run"]:
        assert phrase in amendment, phrase
    assert "Phase 1A.1 mandatory amendment" in report

    tracked = subprocess.run(
        ["git", "-c", f"safe.directory={ROOT.as_posix()}", "ls-files", "data/raw", "data/interim"],
        cwd=ROOT, check=True, capture_output=True, text=True, encoding="utf-8", errors="replace",
    ).stdout.splitlines()
    assert set(tracked) <= {"data/raw/.gitkeep", "data/interim/.gitkeep"}, tracked
    forbidden = ["01_deg", "01_meta_analysis", "01_cellchat", "01_nichenet", "01_wgcna", "01_ml"]
    for name in forbidden:
        assert not (ROOT / "results" / name).exists(), f"Phase 1B/later output exists: {name}"
    print("Phase 1A.1 validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
