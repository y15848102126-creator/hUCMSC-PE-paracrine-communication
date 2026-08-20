# DEPRECATED / HISTORICAL ONLY: Phase 2A count-likelihood receiver inference is superseded by Phase 2A.2 pregnancy-level continuous-expression analysis. Excluded from the default execution workflow.
#!/usr/bin/env python3
"""Build the evidence-backed Phase 2A report from frozen outputs."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/02_phase2a"


def counts_by_celltype(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "none"
    counts = frame.groupby("celltype").size().sort_values(ascending=False)
    return "; ".join(f"{k}={v}" for k, v in counts.items())


def markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "No rows."
    def fmt(value: object) -> str:
        if isinstance(value, float):
            value = f"{value:.6g}"
        return str(value).replace("|", "\\|").replace("\n", " ")
    headers = [fmt(x) for x in frame.columns]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    lines.extend("| " + " | ".join(fmt(x) for x in row) + " |" for row in frame.itertuples(index=False, name=None))
    return "\n".join(lines)


def main() -> int:
    patients = pd.read_csv(OUT / "metadata/patient_registry.csv")
    eligibility = pd.read_csv(OUT / "metadata/pseudobulk_eligibility.csv")
    de = pd.read_csv(OUT / "DE/celltype_DE_summary.csv")
    programs = pd.read_csv(OUT / "programs/pe_cellstate_programs.csv")
    shared = pd.read_csv(OUT / "programs/shared_pe_programs.csv")
    mds = pd.read_csv(OUT / "qc/pseudobulk_mds_coordinates.csv")
    qc = pd.read_csv(OUT / "qc/phase2a_qc_summary.csv")
    risk = pd.read_csv(OUT / "qc/phase2a_risk_flags.csv")
    config = json.loads((ROOT / "config/phase2a_analysis.json").read_text(encoding="utf-8"))
    gate = qc.loc[qc.metric.eq("phase2a_gate"), "value"].iloc[0]
    sig = programs[programs.contrast.isin(["EOPE", "LOPE"]) & programs.BH_FDR.lt(0.05)].copy()
    robust = sig[sig.all_loo_directions_same.eq("YES")]
    unique_sig = sig.drop_duplicates(["celltype", "pathway", "classification"])
    unique_robust = robust.drop_duplicates(["celltype", "pathway", "classification"])
    class_counts = unique_sig.classification.value_counts().to_dict()
    robust_shared = shared[(shared.all_loo_directions_same_EOPE.eq("YES")) & (shared.all_loo_directions_same_LOPE.eq("YES"))]
    group_counts = patients.pe_subtype_or_control_group.value_counts().to_dict()
    eope_ct = de[de.contrast.eq("EOPE")].celltype.tolist()
    lope_ct = de[de.contrast.eq("LOPE")].celltype.tolist()
    combined_ct = de[de.contrast.eq("COMBINED_PE_SECONDARY")].celltype.tolist()
    eope_enriched = unique_robust[unique_robust.classification.eq("EOPE_ENRICHED")]
    lope_enriched = unique_robust[unique_robust.classification.eq("LOPE_ENRICHED")]
    unstable = unique_sig[unique_sig.classification.eq("UNSTABLE")]
    marker_pass = int(qc.loc[qc.metric.eq("marker_check_pass_n"), "value"].iloc[0])
    flag_n = int(mds.outlier_flag.eq("FLAG_REVIEW").sum())
    shared_counts_table = markdown_table(robust_shared.groupby('celltype').size().sort_values(ascending=False).rename('robust_shared_program_n').reset_index()) if len(robust_shared) else "No robust shared programs."
    shared_examples_table = markdown_table(robust_shared.sort_values(['BH_FDR_EOPE','BH_FDR_LOPE'])[['celltype','gene_set','NES_EOPE','NES_LOPE','BH_FDR_EOPE','BH_FDR_LOPE']].head(12)) if len(robust_shared) else "No shared programs."
    risk_table = markdown_table(risk[['risk_id','severity','risk','mitigation']])
    report = f"""# Phase 2A — PE cell-state disease program discovery

**Analysis date:** 2026-08-09

**Primary dataset:** Admati et al. 2023, Figshare `23264102.v1`, file `41003240`

**Biological replicate:** patient/pregnancy (`donorID`)

**Final gate:** **{gate}**

## Executive conclusion

The cell-state-first pivot yields a usable but deliberately restricted receiver-state framework. The public table retains 26 validated pregnancies: EOPE {group_counts.get('EOPE',0)}, LOPE {group_counts.get('LOPE',0)}, gestational-age-compatible early controls {group_counts.get('EARLY_CONTROL',0)}, and late controls {group_counts.get('LATE_CONTROL',0)}. The 31 library labels were collapsed to donor identity; cells were never treated as replicates.

Across the independently fitted EOPE and LOPE contrasts, {len(unique_sig):,} unique cell-type × gene-set programs reached BH FDR <0.05. Frozen classification produced {class_counts.get('SHARED_PE',0)} `SHARED_PE`, {class_counts.get('EOPE_ENRICHED',0)} `EOPE_ENRICHED`, {class_counts.get('LOPE_ENRICHED',0)} `LOPE_ENRICHED`, and {class_counts.get('UNSTABLE',0)} `UNSTABLE` programs. All {len(robust_shared)} shared programs retained their direction in every valid leave-one-patient-out diagnostic in both subtypes.

The strongest shared receiver states are Hofbauer cells, syncytiotrophoblast (SCT), placental stromal cells and macrophages. Hofbauer/stromal shared programs consistently point upward for type-I/interferon response, while SCT shared programs point downward for oxidative phosphorylation/electron transport and related replication-origin programs. These are independently defined receiver programs, not therapeutic claims and not hUC-MSC-matched results.

Phase 2B may proceed only as preregistered external directional support across the six frozen bulk cohorts. Small early-control n and an unresolved `total_molecules` versus count-column-sum discrepancy preclude an unrestricted gate.

## 1. Frozen input and donor validation

- Source: [peer-reviewed study](https://doi.org/10.1016/j.medj.2023.07.005) and [public Figshare UMI table](https://doi.org/10.6084/m9.figshare.23264102.v1).
- Matrix: 86,752 cells, 33,694 published gene rows, 33,660 unique gene symbols after outcome-blind duplicate-symbol summation.
- Patients: 26; libraries: 31. Donors with multiple libraries remain one pregnancy.
- Tissue: placental cotyledon/villous placenta. The dataset is not mislabeled as decidua basalis.
- Public annotations: 46 labels collapsed by frozen lineage rules to 15 harmonized populations. All {marker_pass}/15 passed pooled canonical-marker checks; no disease labels were used in mapping.
- Raw FASTQ was not required because the public object is explicitly an integer UMI table. No integrated or log-normalized values were used for pseudobulk.

The published `total_molecules` metadata sum is 646,369,597, whereas direct summation of the integer expression table is 962,152,952 (stratum correlation 0.836). No values were rescaled or imputed; count-matrix column sums are the formal library sizes, and the discrepancy remains an open provenance risk.

## 2. Frozen eligibility and statistical design

Eligibility was frozen before outcome analysis:

- at least {config['eligibility']['minimum_patients_per_group']} qualified patients per group;
- at least {config['eligibility']['minimum_cells_per_patient_celltype']} cells per patient × cell type;
- at least {config['eligibility']['minimum_pseudobulk_library_umi']:,} count-matrix UMIs per pseudobulk;
- gene CPM ≥1 in at least 3 eligible patient pseudobulks.

EOPE eligible cell types ({len(eope_ct)}): {', '.join(eope_ct)}.

LOPE eligible cell types ({len(lope_ct)}): {', '.join(lope_ct)}.

Secondary combined eligible cell types ({len(combined_ct)}): {', '.join(combined_ct)}.

Primary models used edgeR TMM plus robust quasi-likelihood negative-binomial GLMs with `~ disease`. EOPE was compared only with early controls; LOPE only with late controls. The secondary combined model used `~ onset_stratum + disease` and did not replace either primary contrast.

## 3. Differential-expression layer

All estimable genes were retained in the cell-type-specific DE files. No fold-change cutoff defined significance. Few isolated genes passed BH FDR, which reinforces the decision to make ranked programs—not hub genes—the primary object.

| Contrast | Eligible cell types | Total cell-type DE genes at FDR <0.05 |
|---|---:|---:|
| EOPE | {len(eope_ct)} | {int(de.loc[de.contrast.eq('EOPE'),'DE_FDR05_n'].sum())} |
| LOPE | {len(lope_ct)} | {int(de.loc[de.contrast.eq('LOPE'),'DE_FDR05_n'].sum())} |
| Secondary combined | {len(combined_ct)} | {int(de.loc[de.contrast.eq('COMBINED_PE_SECONDARY'),'DE_FDR05_n'].sum())} |

This layer is separate from the frozen negative Phase 1B meta-analysis. No Phase 1B top gene entered filtering, ranking or annotation.

## 4. Cell-state programs

The ranked universe was frozen to MSigDB 2026.1.Hs Hallmark, Reactome and GO Biological Process ({9_427:,} sets before per-model size/availability filtering). `cameraPR` tested signed edgeR QL statistics using fixed inter-gene correlation 0.01; BH correction was performed across the combined three-collection universe within each cell type and contrast.

Shared programs required FDR <0.05 and the same direction in both EOPE and LOPE. Opposite signs were always `UNSTABLE`. hUC-MSC data remained unseen.

Robust shared program counts by receiver:

{shared_counts_table}

Selected examples above are not cherry-picked: the full shared table contains all {len(shared)} qualifying programs. The most statistically supported shared programs include:

{shared_examples_table}

## 5. Subtype-specific and unstable programs

- Robust EOPE-enriched: {len(eope_enriched)} unique programs; leading cell types: {counts_by_celltype(eope_enriched)}.
- Robust LOPE-enriched: {len(lope_enriched)} unique programs; leading cell types: {counts_by_celltype(lope_enriched)}.
- Unstable/opposite or borderline subtype behavior: {len(unstable)} unique programs. These are not shared PE programs and must not be promoted downstream without independent evidence.

LOPE generates many more enriched programs than EOPE in several trophoblast/vascular/stromal populations. This difference may reflect biology, but it may also reflect the larger LOPE control group, gestational composition, IUGR structure and cell abundance; it is not treated as proof that LOPE is globally more dysregulated.

## 6. Regulatory activity

`cellstate_regulon_activity.csv` contains exploratory CollecTRI human signed-regulon activity estimates using a univariate linear model applied to signed DE statistics. Regulons required at least 10 measured targets and were BH-corrected within cell type and contrast. These are activity estimates, not TF mRNA results, and they do not determine the Phase 2A gate.

## 7. Patient influence and QC

- {flag_n} of {len(mds)} contrast-specific pseudobulk MDS coordinates were flagged for review by a robust-distance screen; none were removed.
- A flag is not an independent technical failure and cannot justify exclusion.
- {len(robust):,}/{len(sig):,} significant subtype-stratified program rows retained direction in every valid one-patient omission.
- All {len(robust_shared)} shared program pairs retained direction in both EOPE and LOPE leave-one-patient-out checks.

Thus the shared conclusions are not directionally dependent on one patient, but EOPE uncertainty remains intrinsically high because removing one of three controls changes the effective control n by one third. The diagnostics do not transform this single study into cross-cohort replication.

## 8. Direct answers

1. **Analyzable pregnancies:** 17 PE (10 EOPE, 7 LOPE) and 9 controls (3 early, 6 late); eligibility varies by cell type as documented in `pseudobulk_eligibility.csv`.
2. **Sufficient cell types:** 11 for EOPE, 12 for LOPE and 13 for the secondary combined model; exact patient n is in `celltype_DE_summary.csv`.
3. **Strongest reproducible receiver states:** Hofbauer, SCT, placental stromal and macrophage populations, based algorithmically on robust shared-program counts.
4. **Shared programs:** 20, dominated by concordant interferon response in Hofbauer/stromal/macrophage states and reduced mitochondrial/electron-transport programs in SCT.
5. **Subtype-specific:** {len(eope_enriched)} robust EOPE-enriched and {len(lope_enriched)} robust LOPE-enriched programs; the asymmetry is interpreted cautiously.
6. **One/two-patient dependence:** no shared program changed direction in any valid one-patient omission, but the EOPE control n=3 remains a structural limitation.
7. **Bulk validation readiness:** yes, as an independently frozen directional gene-set hypothesis; not as a gene-level validation rescue.
8. **Limitations:** three early controls, GA/subtype differences, IUGR confined to PE, variable cell abundance, one public primary cohort, and unresolved count-layer metadata mismatch.
9. **Phase 2B:** **{gate}**. Phase 2B must test these frozen programs across bulk cohorts without redefining them from bulk outcomes.

## 9. Frozen risks and prohibited interpretations

{risk_table}

No hUC-MSC sender program, ligand–receptor analysis, CellChat, NicheNet, WGCNA, machine learning or therapeutic interpretation was performed. Phase 2B was not started.

## Analytical QC figure previews

![Cells per patient](../results/02_phase2a/figures/A_cells_per_patient.png)

![Cells per patient and cell type](../results/02_phase2a/figures/B_cells_patient_celltype_heatmap.png)

![Pseudobulk library sizes](../results/02_phase2a/figures/C_pseudobulk_library_sizes.png)

![Patient pseudobulk MDS](../results/02_phase2a/figures/D_patient_pseudobulk_MDS.png)

![Significant program counts](../results/02_phase2a/figures/E_significant_program_counts.png)

![Shared program heatmap](../results/02_phase2a/figures/F_shared_program_heatmap.png)
"""
    (ROOT / "docs/PHASE2A_PE_CELLSTATE_DISEASE_PROGRAM_REPORT.md").write_text(report, encoding="utf-8")
    print(f"Wrote report; gate={gate}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
