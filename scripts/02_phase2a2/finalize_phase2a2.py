#!/usr/bin/env python3
"""Assemble Phase 2A.2 evidence hierarchy, risks, report, and README summary."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/02_phase2a2"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def f(x: str) -> float:
    return float(x) if x else float("nan")


def main() -> int:
    t1 = read_csv(OUT / "corrected_analysis/frozen20_corrected_retest.csv")
    redisc = read_csv(OUT / "corrected_analysis/corrected_program_rediscovery.csv")
    module_rows = read_csv(OUT / "corrected_analysis/corrected_program_modules.csv")
    modules = [r for r in module_rows if r["record_type"] == "PROGRAM_MODULE"]
    yang = read_csv(OUT / "external_validation/yang_lope_updated_evidence.csv")
    yang_by_module: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in yang:
        yang_by_module[row["program_module"]].append(row)

    hierarchy = []
    for row in modules:
        y = yang_by_module[row["program_module"]]
        yes = sum(r["direction_agrees"] == "YES" for r in y)
        no = sum(r["direction_agrees"] == "NO" for r in y)
        evaluable = yes + no
        corrected = row["corrected_module_status"] == "CORRECTED_ADMATI_SUPPORT"
        if corrected and yes > 0 and no == 0:
            level = "LEVEL_A"
            rationale = "Corrected Admati support plus independent Yang LOPE directional support."
            candidate = "YES"
        elif corrected and evaluable == 0:
            level = "LEVEL_B"
            rationale = "Corrected Admati support; Yang is not evaluable and Zheng targeted EOPE validation did not pass its audit gate."
            candidate = "YES_WITH_RESTRICTIONS"
        elif corrected and no > 0:
            level = "UNCLASSIFIED_FROZEN_HIERARCHY_GAP"
            rationale = "Corrected Admati support but evaluable Yang evidence is directionally discordant; frozen A-D definitions contain no honest category for this combination."
            candidate = "NO_HOLD_EXTERNAL_DISCORDANCE"
        elif row["corrected_module_status"] == "DIRECTIONAL_CORRECTED_SUPPORT_ONLY":
            level = "LEVEL_C"; rationale = "Directional corrected support only."; candidate = "NO"
        else:
            level = "LEVEL_D"; rationale = "Not supported after corrected analysis."; candidate = "NO"
        hierarchy.append({
            "program_module": row["program_module"], "module_label": row["module_label"], "celltype": row["celltype"], "frozen_direction": row["frozen_direction"],
            "constituent_gene_set_n": row["constituent_gene_set_n"], "corrected_module_status": row["corrected_module_status"],
            "corrected_direction_agree_both_proportion": row["direction_agree_both_proportion"], "corrected_any_constituent_fdr05": row["any_constituent_fdr05"],
            "yang_evaluable_constituent_n": evaluable, "yang_direction_agree_n": yes, "yang_direction_disagree_n": no,
            "zheng_eope_status": "NOT_RUN_EOPE_ESTIMAND_CONTRADICTION", "evidence_level": level, "phase2b_program_validation_candidate": candidate,
            "rationale": rationale, "bulk_evidence_used": "NO_PHASE2B_LOCKED",
            "source_url": "results/02_phase2a2/corrected_analysis/corrected_program_modules.csv|results/02_phase2a2/external_validation/yang_lope_updated_evidence.csv|results/02_phase2a2/external_validation/zheng_eope_dataset_audit.csv",
        })
    write_csv(OUT / "evidence/receiver_module_evidence_hierarchy.csv", hierarchy)

    risks = [
        ("P2A2-001","CRITICAL","RESOLVED_METHOD_ROUTE","Admati true raw UMI counts are not publicly recoverable","The only public scRNA table is normalized to ~10,000 per cell and ceiled.","Use donor means plus limma Route B; never use edgeR/DESeq2 NB on this layer."),
        ("P2A2-002","HIGH","OPEN_PROVENANCE_LABEL_CONTRADICTION","Figshare calls the matrix UMI counts","Numerical values and author code contradict the repository description; the reason for the label is undocumented.","Preserve evidence and request true Cell Ranger outputs from authors if future access is possible."),
        ("P2A2-003","CRITICAL","NON_ESTIMABLE","EOPE delivery-mode confounding","Early controls are all vaginal while EOPE is predominantly C-section.","Do not fit a delivery adjustment without positivity; require external matched validation."),
        ("P2A2-004","HIGH","NON_ESTIMABLE","Induction confounding","Disease and induction lack adequate cross-class positivity in the frozen cohort.","No forced regression adjustment."),
        ("P2A2-005","HIGH","OPEN","IUGR/FGR confounding","IUGR occurs only in PE; restriction sensitivity does not identify a disease effect among FGR cases.","Carry Phase 2A.1 restriction results and interpret receiver programs as PE-associated, not PE-specific."),
        ("P2A2-006","CRITICAL","UNRESOLVED","Zheng EOPE identity contradiction","The paper labels all three cases EOPE, while GSE298119 explicitly describes PE002/PE003 as term PE; sample-level GA is absent.","Do not run targeted EOPE statistics until an authoritative sample-level crosswalk resolves the estimand."),
        ("P2A2-007","HIGH","OPEN","Zheng public annotation gap","10x matrices are public, but the 13-type cell annotation-to-barcode map was not identified.","Require public/publisher annotation object or perform a separately preregistered outcome-blind reannotation after clinical identity is resolved."),
        ("P2A2-008","HIGH","OPEN","Zheng series and batch structure","Five subjects are split across three BioProjects/GSEs; one control comes from a term GDM/macrosomia study.","Treat series as technical structure, not independent replication; n=3+2 cannot support complex adjustment."),
        ("P2A2-009","MEDIUM","OPEN_FRAMEWORK_GAP","External discordance is not represented by frozen Level A-D labels","Module 08 has corrected Admati support but an evaluable opposite Yang direction.","Keep an explicit unclassified status and withhold it from Phase 2B rather than relabel post hoc."),
        ("P2A2-010","HIGH","OPEN","Clinical limitations survive expression-model correction","EOPE GA/sex/age are limited by early-control n=3; LOPE induction is non-estimable and IUGR positivity incomplete.","Retain restrictions in every downstream interpretation."),
        ("P2A2-011","MEDIUM","OPEN","Full corrected rediscovery is exploratory and LOPE-heavy","Prior results were known; corrected full-universe results cannot be called independent confirmation.","Keep Tier 2 separate from frozen Tier 1 and validate only pre-specified modules in Phase 2B."),
    ]
    risk_rows = [{"risk_id":a,"severity":b,"status":c,"risk":d,"impact":e,"mitigation":g,"source_url":"config/phase2a2_analysis.json|docs/PHASE2A2_CORRECTED_RECEIVER_ANALYSIS_PLAN.md"} for a,b,c,d,e,g in risks]
    write_csv(OUT / "evidence/phase2a2_risk_flags.csv", risk_rows)

    t1_counts = Counter(r["classification"] for r in t1)
    red_counts = Counter(r["classification"] for r in redisc)
    legacy_paths = {(r["celltype"],r["pathway"]) for r in t1}
    rediscovery_shared_legacy_overlap = sum(r["classification"] == "CORRECTED_SHARED_PE" and (r["celltype"],r["pathway"]) in legacy_paths for r in redisc)
    level_counts = Counter(r["evidence_level"] for r in hierarchy)
    candidates = [r for r in hierarchy if r["phase2b_program_validation_candidate"].startswith("YES")]
    candidate_labels = ", ".join(r["program_module"] for r in candidates)
    report = f"""# Phase 2A.2 receiver framework correction report

**Audit/analysis date:** 2026-08-09  
**Historical boundary:** Phase 2A is preserved as `LEGACY_COUNT_MODEL_DISCOVERY`; Phase 2A.1 is unchanged.  
**Final gate:** `GO_TO_PHASE2B_WITH_RESTRICTIONS`

## Executive conclusion

The Admati public scRNA expression layer is not an unmodified raw-UMI matrix. Figshare article 23264102.v1 contains one ZIP whose 86,752 cell columns behave as per-cell library normalization to approximately 10,000 followed by `ceil()`. The author GitHub code independently contains the same transformation, while its Cell Ranger loader points to internal paths that are not publicly deposited. Searches of Figshare, GitHub, GEO/GDS, SRA, BioProject and BioSample found no donor-linked true raw Admati UMI matrix or FASTQ. The formal provenance conclusion is `PUBLIC_MATRIX_NORMALIZED_RAW_NOT_PUBLIC`.

The corrected analysis therefore uses Route B: frozen eligible cells are averaged within pregnancy × cell type from the public normalized layer, transformed as `log2(mean + 1)`, and analyzed by pregnancy-level limma empirical-Bayes models. No negative-binomial likelihood is applied. All 20 legacy hypotheses retain their frozen direction in both EOPE and LOPE; all {t1_counts['CORRECTED_SHARED_SUPPORT']} meet the predeclared `CORRECTED_SHARED_SUPPORT` definition and both subtype FDRs are <0.05 within the frozen 20-test family.

This is method correction in the same discovery cohort, not independent validation. The full frozen MSigDB rerun is explicitly exploratory: {red_counts['CORRECTED_SHARED_PE']:,} sets are corrected-shared, {red_counts['CORRECTED_EOPE_ENRICHED']:,} EOPE-enriched, {red_counts['CORRECTED_LOPE_ENRICHED']:,} LOPE-enriched, {red_counts['CORRECTED_UNSTABLE']:,} unstable and {red_counts['CORRECTED_NOT_SIGNIFICANT']:,} not significant. {rediscovery_shared_legacy_overlap} of the 20 legacy hypotheses also pass the much larger full-universe shared correction.

## Provenance resolution

The exact Phase 2A file was `sc_PE_allcells_with_metadata_29-May-2023.txt.zip` (256,501,211 bytes; SHA-256 `1b317645daf5331bb6ac4d7d858fea3f78deb712dcfbce2df5f386c161d88446`), containing only `sc_PE_allcells_with_metadata_29-May-2023.txt` (5,871,738,353 bytes). The file-level inventory, ZIP member, checksum, numerical signature, loader evidence and linked nuclei article are recorded in `admati_expression_layer_audit.csv`.

The repository description's phrase “UMI counts” cannot be reconciled with the deposited values from public evidence. The most specific evidence-supported statement is that the deposited integer table is the result of the public author transformation `ceil(data / column_sum × 10,000)`. Why Figshare used the raw-UMI label is **UNRESOLVED**; no motive or undocumented processing step is inferred.

True Cell Ranger `matrix.mtx`, `barcodes.tsv`, `features.tsv` and FASTQ are not publicly recoverable as of the audit date. The internal loader proves that such files existed in the authors' workflow, not that they were deposited. The linked Figshare article 23264165.v1 is a separate trophoblast single-nucleus dataset.

## Corrected method and gene-level receiver statistics

The corrected input is the same public Figshare table, interpreted as normalized/ceiled continuous expression. For each frozen eligible donor-cell-type stratum, gene values are summed across cells, divided by the frozen eligible cell count, and transformed by `log2(mean + 1)`. EOPE is compared only with early controls; LOPE only with late controls. Each cell type uses `~ disease` with limma `eBayes(trend=TRUE, robust=TRUE)`. The gene filter is donor mean >0 in at least three eligible donor summaries, independent of disease direction.

There are 368,370 corrected gene-statistic rows across 23 eligible subtype × cell-type models. These statistics—not the legacy edgeR outputs—are the frozen future receiver gene layer. They are not yet NicheNet targets.

## Tier 1: frozen 20-hypothesis retest

All 20 hypotheses agree with their historical direction in both subtypes and survive BH correction within the frozen 20-hypothesis family in both EOPE and LOPE. This supports the direction of the historical themes under an appropriate continuous-expression model, while leaving the original Phase 2A P/FDR values untouched and historically labeled.

The 20 rows reduce to 11 legacy modules: four Hofbauer interferon/antiviral modules, one macrophage interferon module, three placental-stromal interferon modules, and three SCT mitochondrial/replication modules. All 11 have corrected Admati support under the frozen constituent rule.

## Tier 2: exploratory corrected rediscovery

The full-universe rerun does not reproduce the legacy result as a simple 20-row confirmation. It yields only {red_counts['CORRECTED_SHARED_PE']} full-universe shared sets and remains strongly LOPE-heavy ({red_counts['CORRECTED_LOPE_ENRICHED']} LOPE-enriched versus {red_counts['CORRECTED_EOPE_ENRICHED']} EOPE-enriched). This tier is not preregistered confirmation because the historical outcome was already known. The difference between Tier 1 and Tier 2 is expected from their different multiplicity families: 20 frozen tests versus thousands of measured sets per cell type.

## Zheng elective-C-section EOPE rescue audit

The reported five subjects are not three independent GEO cohorts. The reconstructable crosswalk is PE001 (`GSM8634701`, GSE282038), PE002/PE003 (`GSM9008678/79`, GSE298119), CONTROL (`GSM8264272`, GSE267340), and CTL2 (`GSM9008680`, GSE298119). Cell Ranger-style matrices and SRA reads are public for all five.

Targeted validation was not run. GSE298119 describes PE002/PE003 as **term PE**, whereas the Frontiers paper calls all three cases EOPE; sample-level GA/FGR data for PE002, PE003 and CTL2 are absent. The exact public cell annotation-to-barcode map was also not found. Moreover, the five samples span three GSEs/BioProjects, and CONTROL is embedded in a term GDM/macrosomia series. These unresolved issues fail the frozen EOPE-estimand gate despite the paper-level statement that all five underwent elective cesarean section. Consequently no EOPE receiver theme can honestly be claimed to replicate in Zheng.

## Yang LOPE evidence and receiver hierarchy

Phase 2A.1 Yang results were not rerun. They retain partial approximate-match directional support for placental-stromal IFN-alpha response (Module 06) and IFN-alpha/beta signaling (Module 07). The stromal IFN-stimulated host-response Module 08 is directionally discordant. SCT is not evaluable under the frozen minimum-cell rule, and Yang's macrophage label cannot be forced into Hofbauer equivalence.

The hierarchy contains {level_counts['LEVEL_A']} Level A modules and {level_counts['LEVEL_B']} Level B modules. Module 08 is explicitly `UNCLASSIFIED_FROZEN_HIERARCHY_GAP`: it has corrected Admati support but evaluable, directionally discordant Yang evidence, a combination not covered by the frozen A-D definitions. It is withheld rather than relabeled post hoc. The ten modules eligible for independent bulk program-level validation are: {candidate_labels}.

## Clinical confounding that remains non-identifiable

- EOPE delivery mode remains `NON_ESTIMABLE`: early controls are all vaginal and no adequate C-section early-control support exists.
- Induction remains `NON_ESTIMABLE`; disease and induction lack positivity.
- IUGR/FGR remains structurally disease-linked. Non-IUGR sensitivity cannot identify effects among FGR pregnancies.
- EOPE GA, fetal sex and maternal age remain weakly estimable at best because early-control n=3.
- LOPE C-section restriction remains estimable from Phase 2A.1, but induction is non-estimable and IUGR positivity is incomplete.
- Corrected expression modeling does not repair any of these clinical design limitations.

## Explicit answers

1. **Are true raw Admati UMI counts publicly recoverable?** No. Status: `PUBLIC_MATRIX_NORMALIZED_RAW_NOT_PUBLIC`.
2. **What exact layer is used?** Figshare 23264102.v1's per-cell ~10,000-normalized, ceiled expression, aggregated to donor means and transformed `log2(mean+1)`.
3. **Do the 20 retain direction?** Yes, 20/20 in both EOPE and LOPE.
4. **Which survive frozen-family correction?** All 20 survive BH FDR<0.05 in both subtypes within the frozen 20-test family.
5. **What changes in full rediscovery?** Only {red_counts['CORRECTED_SHARED_PE']} sets are full-universe corrected-shared; subtype-enriched results remain LOPE-heavy. It is exploratory and not a replacement for Tier 1.
6. **Which modules remain credible?** All 11 have corrected Admati support; Modules 06/07 reach Level A, eight reach Level B, and Module 08 is held because Yang points oppositely. Ten enter the Phase 2B candidate list.
7. **Is Zheng independently reusable?** Not as an EOPE targeted validation cohort until its term-versus-EOPE identity and annotation crosswalk are resolved. Its three GSEs are parts of one reported cohort.
8. **Do EOPE themes replicate there?** Not assessable; no targeted statistics were run.
9. **Which retain Yang LOPE support?** Approximate placental-stromal IFN-alpha and IFN-alpha/beta signaling; not SCT or Hofbauer/macrophage, and the stromal host-response module is discordant.
10. **Which modules can enter bulk program validation?** Modules 01-07 and 09-11, with Level B restrictions where external scRNA evidence is not evaluable. Module 08 is held.
11. **Which confounders remain fundamentally non-identifiable?** EOPE delivery mode and induction; disease-linked IUGR/FGR; severely limited EOPE GA/sex/age; LOPE induction and incomplete IUGR positivity.

## Final gate

`GO_TO_PHASE2B_WITH_RESTRICTIONS`

Phase 2B may test only the ten frozen candidate modules at program level across independent bulk cohorts. It must not reinterpret legacy edgeR gene statistics, treat Zheng as completed EOPE validation, or unlock sender/communication analyses. Phase 2B was not started here.

## Key artifacts

- [Corrected plan](PHASE2A2_CORRECTED_RECEIVER_ANALYSIS_PLAN.md)
- [Admati layer audit](../results/02_phase2a2/provenance/admati_expression_layer_audit.csv)
- [Frozen 20 retest](../results/02_phase2a2/corrected_analysis/frozen20_corrected_retest.csv)
- [Corrected modules](../results/02_phase2a2/corrected_analysis/corrected_program_modules.csv)
- [Zheng audit](../results/02_phase2a2/external_validation/zheng_eope_dataset_audit.csv)
- [Receiver hierarchy](../results/02_phase2a2/evidence/receiver_module_evidence_hierarchy.csv)

![Frozen 20 corrected direction](../results/02_phase2a2/figures/A_frozen20_corrected_direction.png)

![Tier-1 classification](../results/02_phase2a2/figures/B_tier1_classification.png)

![Module support](../results/02_phase2a2/figures/C_corrected_module_status.png)
"""
    (ROOT / "docs/PHASE2A2_RECEIVER_FRAMEWORK_CORRECTION_REPORT.md").write_text(report,encoding="utf-8")
    print({"hierarchy": dict(level_counts), "phase2b_candidates": len(candidates), "risks": len(risks)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
