#!/usr/bin/env python3
"""Structural and frozen-rule validation for Phase 2B outputs."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "02_phase2b"


def rows(rel: str) -> list[dict[str, str]]:
    with (ROOT / rel).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def sha256(rel: str) -> str:
    h = hashlib.sha256()
    with (ROOT / rel).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def num(value: str) -> float:
    return float(value)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    required = [
        "config/phase2b_analysis.json",
        "docs/PHASE2B_STATISTICAL_ANALYSIS_PLAN.md",
        "docs/PHASE2B_BULK_PROGRAM_VALIDATION_REPORT.md",
        "results/02_phase2b/hypotheses/frozen_phase2b_gene_sets.csv",
        "results/02_phase2b/hypotheses/frozen_phase2b_modules.csv",
        "results/02_phase2b/scores/cohort_program_scores_summary.csv",
        "results/02_phase2b/scores/program_gene_coverage.csv",
        "results/02_phase2b/meta/program_gene_set_meta_analysis.csv",
        "results/02_phase2b/meta/program_module_validation.csv",
        "results/02_phase2b/robustness/program_leave_one_cohort_out.csv",
        "results/02_phase2b/robustness/cameraPR_method_concordance.csv",
        "results/02_phase2b/robustness/program_covariate_sensitivity.csv",
        "results/02_phase2b/evidence/updated_receiver_evidence_hierarchy.csv",
        "results/02_phase2b/evidence/phase2b_risk_flags.csv",
        "results/02_phase2b/phase2b_session_info.txt",
    ]
    for rel in required:
        require((ROOT / rel).is_file(), f"missing required output: {rel}")

    cfg = json.loads((ROOT / "config/phase2b_analysis.json").read_text(encoding="utf-8"))
    for rel, expected in cfg["receiver_history"]["upstream_sha256"].items():
        require(sha256(rel) == expected, f"upstream freeze hash changed: {rel}")
    require(cfg["outcome_inspection_status_at_freeze"] == "NOT_INSPECTED", "freeze provenance altered")

    hyp = rows("results/02_phase2b/hypotheses/frozen_phase2b_gene_sets.csv")
    modules = rows("results/02_phase2b/hypotheses/frozen_phase2b_modules.csv")
    require(len(hyp) == 19, "frozen family is not exactly 19 rows")
    require(len(modules) == 10, "frozen module family is not exactly 10 rows")
    require(len({r["gene_set"] for r in hyp}) == 16, "expected 16 unique memberships across 19 hypotheses")
    require(all(r["program_module"] != "PROGRAM_MODULE_08" for r in hyp), "held Module 08 entered hypothesis family")
    require({r["program_module"] for r in hyp} == set(cfg["candidate_modules"]), "module family differs from config")
    for r in hyp:
        genes = r["original_gene_membership"].split(";")
        require(len(genes) == int(r["original_gene_n"]), f"membership count mismatch: {r['hypothesis_id']}")
        require(hashlib.sha256(";".join(genes).encode()).hexdigest() == r["membership_sha256"], f"membership hash mismatch: {r['hypothesis_id']}")

    cohorts = cfg["core_cohorts"]
    coverage = rows("results/02_phase2b/scores/program_gene_coverage.csv")
    scores = rows("results/02_phase2b/scores/cohort_program_scores_summary.csv")
    require(len(coverage) == 19 * 6, "coverage table must be 19 x 6")
    require(len(scores) == 19 * 6, "model table must be 19 x 6")
    require({r["cohort"] for r in coverage} == set(cohorts), "wrong core cohorts in coverage")
    require(all(r["estimable"] == "YES" for r in coverage), "unexpected non-estimable frozen hypothesis")
    for r in coverage:
        require(int(r["measured_gene_n"]) >= 15, "measured-gene rule violated")
        require(num(r["membership_coverage"]) >= 0.60, "coverage fraction rule violated")
    for r in scores:
        require(r["model_formula"] == cfg["primary_formulas"][r["cohort"]], f"formula mismatch: {r['cohort']}")
        require(math.isfinite(num(r["beta_disease"])) and num(r["SE"]) > 0, "invalid cohort effect")

    meta = rows("results/02_phase2b/meta/program_gene_set_meta_analysis.csv")
    loco = rows("results/02_phase2b/robustness/program_leave_one_cohort_out.csv")
    mv = rows("results/02_phase2b/meta/program_module_validation.csv")
    require(len(meta) == 19 and len(loco) == 19 * 6 and len(mv) == 10, "meta/module/LOCO row count mismatch")
    expected_class = Counter({"BULK_ROBUST_SUPPORT": 0, "BULK_DIRECTIONAL_SUPPORT": 9, "BULK_HETEROGENEOUS_SUPPORT": 0, "NOT_BULK_SUPPORTED": 10})
    actual_class = Counter(r["classification"] for r in meta)
    for key, val in expected_class.items():
        require(actual_class[key] == val, f"unexpected constituent classification count: {key}")
    require(min(num(r["BH_FDR"]) for r in meta) >= 0.05, "report expects no multiplicity-significant constituent")
    directional = [r for r in meta if r["classification"] == "BULK_DIRECTIONAL_SUPPORT"]
    require(all(r["all_valid_LOCO_expected_direction"] == "YES" for r in directional), "directional support lost LOCO direction")
    require(all(r["GSE75010_driver"] == "NO" for r in directional), "supported set is BioBank-driven")

    expected_mod = Counter({"BULK_MODULE_SUPPORTED": 0, "BULK_MODULE_DIRECTIONAL": 5, "BULK_MODULE_DISCORDANT": 1, "NOT_BULK_SUPPORTED": 4})
    actual_mod = Counter(r["module_classification"] for r in mv)
    for key, val in expected_mod.items():
        require(actual_mod[key] == val, f"unexpected module classification count: {key}")
    directional_modules = {r["program_module"] for r in mv if r["module_classification"] == "BULK_MODULE_DIRECTIONAL"}
    require(directional_modules == {"PROGRAM_MODULE_01", "PROGRAM_MODULE_04", "PROGRAM_MODULE_05", "PROGRAM_MODULE_07", "PROGRAM_MODULE_10"}, "directional module set changed")

    camera = rows("results/02_phase2b/robustness/cameraPR_method_concordance.csv")
    cov = rows("results/02_phase2b/robustness/program_covariate_sensitivity.csv")
    require(len(camera) == 19 * 6, "cameraPR sensitivity must be 19 x 6")
    require(len(cov) == 19 * 3, "covariate sensitivity must be 19 x 3")
    require(sum(r["primary_camera_direction_agreement"] == "YES" for r in camera) == 95, "camera primary agreement changed")
    require(sum(r["expected_direction_agreement"] == "YES" for r in camera) == 73, "camera expected agreement changed")

    by_model: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in cov:
        by_model[r["model_id"]].append(r)
    expected_material = {"GSE75010_BIOBANK_GA_SEX": 1, "GSE25906_GA_SEX": 14, "GSE25906_GA_SEX_LABOR": 16}
    expected_flips = {"GSE75010_BIOBANK_GA_SEX": 0, "GSE25906_GA_SEX": 1, "GSE25906_GA_SEX_LABOR": 8}
    for key in expected_material:
        require(sum(r["material_change_flag"] == "YES" for r in by_model[key]) == expected_material[key], f"material sensitivity count changed: {key}")
        require(sum(r["sign_concordant"] == "NO" for r in by_model[key]) == expected_flips[key], f"sign sensitivity count changed: {key}")

    evidence = rows("results/02_phase2b/evidence/updated_receiver_evidence_hierarchy.csv")
    require(len(evidence) == 11, "updated hierarchy must preserve all 11 historical modules")
    held = next(r for r in evidence if r["program_module"] == "PROGRAM_MODULE_08")
    require(held["INDEPENDENT_BULK_PROGRAM_SUPPORT"] == "NOT_TESTED_HOLD_EXTERNAL_DISCORDANCE", "Module 08 hold was lost")
    require(all(r["celltype_localization_claim"] == "DEPENDENT_ON_SINGLE_CELL_EVIDENCE_NOT_VALIDATED_BY_BULK" for r in evidence), "bulk overclaims cell localization")

    report = (ROOT / "docs/PHASE2B_BULK_PROGRAM_VALIDATION_REPORT.md").read_text(encoding="utf-8")
    for token in ["STABLE = 0", "GO_TO_PHASE3_WITH_RESTRICTIONS", "BULK_ROBUST_SUPPORT = 0", "Bulk placenta cannot establish"]:
        require(token in report, f"report missing mandatory statement: {token}")

    print("Phase 2B structural/frozen-rule validation: PASS")
    print("19 hypotheses; 10 modules; Module 08 held; 0 robust; 9 directional; 5 directional modules")


if __name__ == "__main__":
    main()
