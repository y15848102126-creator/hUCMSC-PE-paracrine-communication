#!/usr/bin/env python3
"""Validate Phase 1B statistical outputs, frozen rules, and phase boundaries."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "01_phase1b"
COHORTS = ["GSE75010_BIOBANK", "GSE30186", "GSE10588", "GSE24129", "GSE25906", "GSE43942"]


def bh(p: np.ndarray) -> np.ndarray:
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1)
    return result


def assert_csv(path: Path, required: set[str], allow_empty: bool = False) -> pd.DataFrame:
    assert path.exists() and path.stat().st_size > 0, path
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
    assert len(header) == len(set(header)), f"duplicate columns: {path}"
    assert required <= set(header), f"missing columns {required - set(header)}: {path}"
    frame = pd.read_csv(path)
    assert allow_empty or len(frame), f"unexpected empty table: {path}"
    if "source" in frame.columns and len(frame):
        assert frame["source"].notna().all() and (frame["source"].astype(str).str.len() > 0).all(), path
    return frame


def main() -> int:
    config = json.loads((ROOT / "config" / "phase1b_analysis.json").read_text(encoding="utf-8"))
    assert config["meta_analysis"]["method"] == "REML" and config["meta_analysis"]["inference"] == "MODIFIED_KNHA_T"
    assert config["stable_rule"]["absolute_meta_log2fc_cutoff"] is None
    registry = pd.read_csv(ROOT / "results" / "01_phase1a1" / "formal_phase1b_matrix_registry.csv")
    expected_rows = dict(zip(registry["dataset"], registry["mapped_gene_n"]))
    expected_counts = {
        "GSE75010_BIOBANK": (63, 53), "GSE30186": (6, 6), "GSE10588": (17, 26),
        "GSE24129": (8, 8), "GSE25906": (23, 37), "GSE43942": (5, 7),
    }
    de_tables = {}
    for cohort in COHORTS:
        de = assert_csv(
            OUT / "cohort_DE" / f"{cohort}_DE.csv",
            {"gene", "log2FC", "SE", "t_statistic", "raw_P", "BH_FDR", "n_PE", "n_control", "model_formula", "source"},
        )
        assert len(de) == int(expected_rows[cohort]), (cohort, len(de), expected_rows[cohort])
        assert de["gene"].is_unique
        assert np.isfinite(de[["log2FC", "SE", "t_statistic", "raw_P", "BH_FDR"]].to_numpy()).all()
        assert (de["SE"] > 0).all() and de["raw_P"].between(0, 1).all() and de["BH_FDR"].between(0, 1).all()
        assert np.max(np.abs(bh(de["raw_P"].to_numpy()) - de["BH_FDR"].to_numpy())) < 1e-10
        assert (int(de["n_PE"].iloc[0]), int(de["n_control"].iloc[0])) == expected_counts[cohort]
        expected_formula = "~ disease + batch" if cohort == "GSE25906" else "~ disease"
        assert set(de["model_formula"]) == {expected_formula}
        de_tables[cohort] = de

    availability = assert_csv(
        OUT / "qc" / "phase1b_gene_availability.csv",
        {"gene", "cohort", "platform_available", "mapping_status", "estimable_primary", "absence_interpretation", "mapping_version", "source"},
    )
    assert len(availability) == 24403 * 6
    assert not availability.duplicated(["gene", "cohort"]).any()
    assert set(availability["cohort"]) == set(COHORTS)
    assert set(availability["mapping_version"]) == {"HGNC_2026-08-09_faaeb6ae1e2a596b"}

    meta = assert_csv(
        OUT / "meta" / "pe_gene_meta_analysis.csv",
        {"gene", "k_cohorts", "pooled_log2FC", "pooled_SE", "CI_lower", "CI_upper", "test_statistic", "test_df", "raw_meta_P", "BH_FDR", "tau2", "I2", "Cochran_Q", "Q_P", "direction_consistency", "all_LOCO_same_direction", "LOCO_FDR_lt_0_10_proportion", "category", "meta_method", "source"},
    )
    assert len(meta) == 17731 and meta["gene"].is_unique
    assert (meta["k_cohorts"] >= 4).all()
    assert set(meta["meta_method"]) == {"REML_MODIFIED_KNHA_T"}
    assert np.max(np.abs(bh(meta["raw_meta_P"].to_numpy()) - meta["BH_FDR"].to_numpy())) < 1e-10
    assert meta["direction_consistency"].between(0, 1).all() and meta["I2"].between(0, 100).all()
    expected_stable = (
        (meta["k_cohorts"] >= 4) & (meta["BH_FDR"] < 0.05) & (meta["direction_consistency"] >= 0.75)
        & (meta["I2"] <= 60) & meta["all_LOCO_same_direction"].astype(bool)
        & (meta["LOCO_FDR_lt_0_10_proportion"] >= 0.80)
    )
    assert (meta["category"].eq("STABLE") == expected_stable).all()
    allowed = {"STABLE", "ROBUST_BUT_HETEROGENEOUS", "DIRECTION_CONSISTENT_NON_SIGNIFICANT", "COHORT_SPECIFIC", "UNSTABLE"}
    assert set(meta["category"]) <= allowed
    assert not expected_stable.any() and not (meta["BH_FDR"] < 0.05).any()
    assert "absolute_meta_log2fc_cutoff" not in meta.columns
    for cohort in COHORTS:
        assert (meta[f"{cohort}_log2FC"].notna().sum() >= 0)
    availability_k = availability[availability["estimable_primary"] == "YES"].groupby("gene").size()
    assert set(meta["gene"]) == set(availability_k[availability_k >= 4].index)

    stable = assert_csv(OUT / "meta" / "stable_pe_genes.csv", set(meta.columns), allow_empty=True)
    assert len(stable) == int(expected_stable.sum()) == 0
    heterogeneous = assert_csv(OUT / "meta" / "heterogeneous_pe_genes.csv", {"gene", "I2", "category", "source"}, allow_empty=True)
    assert (heterogeneous["I2"] > 60).all()
    direction = assert_csv(OUT / "meta" / "direction_consistency.csv", {"gene", "direction_consistency", "category", "source"})
    assert len(direction) == len(meta)

    loco = assert_csv(
        OUT / "robustness" / "leave_one_cohort_out.csv",
        {"gene", "omitted_cohort", "remaining_cohort_n", "pooled_log2FC", "BH_FDR", "I2", "same_direction_as_full", "BH_family", "source"},
        allow_empty=True,
    )
    assert len(loco) == 0
    standardized = assert_csv(
        OUT / "robustness" / "standardized_effect_sensitivity.csv",
        {"gene", "pooled_Hedges_g", "BH_FDR", "standardized_direction_consistency", "direction_matches_primary", "is_primary_STABLE", "interpretation", "source"},
    )
    assert len(standardized) == len(meta)
    assert standardized["standardized_direction_consistency"].between(0, 1).all()
    assert set(standardized["interpretation"]) == {"PLATFORM_DYNAMIC_RANGE_SENSITIVITY_NOT_COVARIATE_ADJUSTED"}
    assert not standardized["is_primary_STABLE"].astype(bool).any()
    cov = assert_csv(
        OUT / "robustness" / "covariate_sensitivity.csv",
        {"gene", "cohort", "GA_SENSITIVE", "SEX_SENSITIVE", "LABOR_SENSITIVE", "interpretation", "source"},
        allow_empty=True,
    )
    assert len(cov) == 0

    diagnostics = assert_csv(
        OUT / "qc" / "phase1b_model_diagnostics.csv",
        {"cohort", "model_id", "model_formula", "analysis_role", "complete_case_n", "design_rank", "design_column_n", "max_VIF", "estimability", "source"},
    )
    models = diagnostics[diagnostics["analysis_role"] != "TECHNICAL_VALIDATION"]
    assert len(models) == 11
    assert (models["design_rank"] == models["design_column_n"]).all() and (models["max_VIF"] < 5).all()
    assert (models["sample_exclusions"] == "NONE").all()
    engine = diagnostics[diagnostics["model_id"] == "REML_ENGINE_VALIDATION"]
    assert len(engine) == 1 and "MAX_TAU2_DIFF_METAFOR" in engine["estimability"].iloc[0]
    risks = assert_csv(OUT / "qc" / "phase1b_risk_flags.csv", {"risk_id", "severity", "scope", "risk", "mitigation", "status", "source"})
    assert len(risks) >= 6

    flags = pd.read_csv(ROOT / "results" / "01_phase1a1" / "qc_flag_reinterpretation.csv")
    freeze = pd.read_csv(ROOT / "results" / "01_phase1a" / "bulk_sample_freeze.csv")
    assert len(flags) == 18 and flags["retain_phase1b"].eq("YES").all()
    frozen = set(zip(freeze.loc[freeze["include_phase1b"] == "YES", "dataset"], freeze.loc[freeze["include_phase1b"] == "YES", "GSM/sample ID"]))
    assert all((row.dataset, row.sample_id) in frozen for row in flags.itertuples())

    session = (OUT / "qc" / "phase1b_session_info.txt").read_text(encoding="utf-8")
    for phrase in ["Random seed: 20260809", "metafor_5.0-1", "limma_3.66.0", "HGNC_2026-08-09_faaeb6ae1e2a596b"]:
        assert phrase in session
    figures = sorted((OUT / "figures").glob("*.png"))
    assert len(figures) >= 6
    for figure in figures:
        assert figure.stat().st_size > 5000 and figure.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", figure

    report = (ROOT / "docs" / "PHASE1B_PE_DISEASE_SIGNATURE_REPORT.md").read_text(encoding="utf-8")
    for phrase in ["Final gate: **NO_GO**", "17,731", "STABLE = 0", "Phase 1C was not started", "no stable signature exists"]:
        assert phrase in report, phrase
    plan = (ROOT / "docs" / "PHASE1B_STATISTICAL_ANALYSIS_PLAN.md").read_text(encoding="utf-8")
    for phrase in ["~ disease + batch", "Sensitivity-adjusted coefficients never replace", "Descriptive category precedence", "no stable gene is found"]:
        assert phrase in plan, phrase

    tracked = subprocess.run(
        ["git", "-c", f"safe.directory={ROOT.as_posix()}", "ls-files", "data/raw", "data/interim"],
        cwd=ROOT, check=True, capture_output=True, text=True, encoding="utf-8", errors="replace",
    ).stdout.splitlines()
    assert set(tracked) <= {"data/raw/.gitkeep", "data/interim/.gitkeep"}, tracked
    forbidden = ["02_phase1c", "02_pathway", "02_wgcna", "02_ml", "02_cellchat", "02_nichenet", "02_sender"]
    for name in forbidden:
        assert not (ROOT / "results" / name).exists(), name
    print("Phase 1B validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
