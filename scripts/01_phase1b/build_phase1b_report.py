#!/usr/bin/env python3
"""Build the Phase 1B report from frozen result tables without pathway inference."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "01_phase1b"
DOC = ROOT / "docs" / "PHASE1B_PE_DISEASE_SIGNATURE_REPORT.md"


def fnum(value: float, digits: int = 3) -> str:
    if pd.isna(value):
        return "NOT_APPLICABLE"
    return f"{value:.{digits}f}"


def pct(value: float) -> str:
    if pd.isna(value):
        return "NOT_APPLICABLE"
    return f"{100 * value:.1f}%"


def main() -> int:
    summary_df = pd.read_csv(ROOT / "data" / "interim" / "phase1b" / "phase1b_summary.csv")
    summary = dict(zip(summary_df["metric"], summary_df["value"]))
    meta = pd.read_csv(OUT / "meta" / "pe_gene_meta_analysis.csv")
    standardized = pd.read_csv(OUT / "robustness" / "standardized_effect_sensitivity.csv")
    diagnostics = pd.read_csv(OUT / "qc" / "phase1b_model_diagnostics.csv")
    risk = pd.read_csv(OUT / "qc" / "phase1b_risk_flags.csv")

    cohorts = ["GSE75010_BIOBANK", "GSE30186", "GSE10588", "GSE24129", "GSE25906", "GSE43942"]
    cohort_lines = []
    for cohort in cohorts:
        de = pd.read_csv(OUT / "cohort_DE" / f"{cohort}_DE.csv")
        cohort_lines.append(
            f"| {cohort} | {len(de):,} | {int(de['n_PE'].iloc[0])} | {int(de['n_control'].iloc[0])} | "
            f"{int((de['BH_FDR'] < 0.05).sum()):,} | `{de['model_formula'].iloc[0]}` |"
        )

    model_lines = []
    for _, row in diagnostics[diagnostics["analysis_role"] != "TECHNICAL_VALIDATION"].iterrows():
        model_lines.append(
            f"| {row['model_id']} | `{row['model_formula']}` | {int(row['complete_case_n'])} | "
            f"{int(row['design_rank'])}/{int(row['design_column_n'])} | {float(row['max_VIF']):.3f} | {row['estimability']} |"
        )

    category_counts = meta["category"].value_counts().to_dict()
    raw_p_n = int((meta["raw_meta_P"] < 0.05).sum())
    meta_fdr_n = int(summary["meta_fdr_lt_0_05_n"])
    stable_n = int(summary["stable_n"])
    robust_hetero_n = int(summary["robust_but_heterogeneous_n"])
    estimable_n = int(summary["estimable_ge4_n"])
    union_n = int(summary["gene_union_n"])
    i2_high_n = int(summary["i2_gt_60_n"])
    minimum_fdr = float(summary["minimum_meta_BH_FDR"])

    gate = "NO_GO" if stable_n == 0 else ("REMAIN_IN_PHASE1B" if stable_n < 10 else "GO_TO_PHASE1C_WITH_RESTRICTIONS")
    assert gate == "NO_GO", "Frozen gate logic changed unexpectedly"

    text = f"""# Phase 1B — multi-cohort PE disease-signature report

Analysis date: 2026-08-09

Evidence/mapping freeze: 2026-08-08 / HGNC 2026-08-09

Final gate: **{gate}**

Phase 1B is the first formal outcome-analysis phase. It fitted cohort-wise placental PE contrasts and combined evidence without pooling expression matrices. It did not run pathway enrichment, WGCNA, machine learning, single-cell disease-state analysis, CellChat, NicheNet, hUC-MSC sender analysis or therapeutic inference. Phase 1C was not started.

## Executive conclusion

The frozen multi-cohort framework was technically executable, but it did **not** identify a gene meeting the preregistered stable PE criteria. Of {union_n:,} genes in the union, {estimable_n:,} were estimable in at least four independent cohorts. Although {raw_p_n:,} had an unadjusted meta-analysis P below 0.05, none survived the single BH family; the smallest meta BH FDR was {minimum_fdr:.3f}. Consequently, `STABLE = 0` and `ROBUST_BUT_HETEROGENEOUS = 0`.

The negative gate is not caused by adding an absolute log2FC threshold—none exists. It follows directly from the frozen BH, direction, heterogeneity and LOCO requirements. Hedges' g and covariate sensitivities were not used to rescue genes.

## Estimand amendment frozen before outcomes

The primary synthesis used minimally adjusted, comparable coefficients:

- `~ disease` for GSE75010 BioBank, GSE30186, GSE10588, GSE24129 and GSE43942.
- `~ disease + batch` for GSE25906, because batch is a technical design variable.

GA, fetal sex and labor were evaluated only in nested sensitivity models. Their adjusted coefficients were never substituted into only part of the primary meta-analysis. This amendment is recorded in `docs/PHASE1B_STATISTICAL_ANALYSIS_PLAN.md` and `config/phase1b_analysis.json` before outcome fitting.

## Cohort-wise differential expression

No cohort-level DEG filter was applied before synthesis. Every estimable mapped gene was passed forward.

| Cohort | Genes modeled | PE n | Control n | Cohort BH FDR<0.05 (descriptive) | Primary model |
|---|---:|---:|---:|---:|---|
{chr(10).join(cohort_lines)}

Cohort-specific significance is descriptive and is not a substitute for the cross-cohort stable rule.

## Model diagnostics

| Model | Formula | Complete n | Rank/columns | Maximum VIF | Status |
|---|---|---:|---:|---:|---|
{chr(10).join(model_lines)}

All primary and sensitivity designs were full rank, had residual degrees of freedom, maximum VIF below 5 and no disease-by-included-categorical-covariate zero cell. All frozen analytical samples—including the 18 Phase 1A review flags—were retained.

The vectorized REML engine was checked against `metafor 5.0-1` on the first 40 alphabetically ordered eligible genes. Maximum absolute differences were {float(summary['meta_engine_max_tau2_diff']):.2e} for tau² and {float(summary['meta_engine_max_effect_diff']):.2e} for the pooled effect.

## Gene universe and availability

- Union of mapped genes: **{union_n:,}**.
- Estimable in at least four cohorts: **{estimable_n:,}**.
- Platform/mapping absence is encoded as unavailable, never as zero expression or a biological negative result.
- Mapping version: `HGNC_2026-08-09_faaeb6ae1e2a596b`.

The complete long-form gene × cohort availability matrix is `results/01_phase1b/qc/phase1b_gene_availability.csv`.

## Primary REML meta-analysis

The primary analysis used intercept-only random-effects REML with conservative modified Knapp–Hartung t inference. BH correction covered all {estimable_n:,} genes with at least four effects.

| Outcome category | Gene count |
|---|---:|
| Meta BH FDR <0.05 | {meta_fdr_n:,} |
| STABLE | {stable_n:,} |
| ROBUST_BUT_HETEROGENEOUS | {robust_hetero_n:,} |
| DIRECTION_CONSISTENT_NON_SIGNIFICANT | {category_counts.get('DIRECTION_CONSISTENT_NON_SIGNIFICANT', 0):,} |
| COHORT_SPECIFIC | {category_counts.get('COHORT_SPECIFIC', 0):,} |
| UNSTABLE | {category_counts.get('UNSTABLE', 0):,} |
| I² >60% (regardless of significance) | {i2_high_n:,} |

The {category_counts.get('DIRECTION_CONSISTENT_NON_SIGNIFICANT', 0):,} direction-consistent genes are descriptive signals only. Their meta FDR values fail the stable rule and they must not be treated as a frozen disease signature.

## Leave-one-cohort-out robustness

No gene passed the preliminary meta-FDR and direction requirements, so there was no candidate-stable set for formal LOCO membership testing. `leave_one_cohort_out.csv` therefore intentionally contains its schema and zero data rows.

As a non-membership diagnostic across all eligible genes, omitting GSE75010 BioBank retained the full-analysis direction for {pct(float(summary['biobank_omission_all_sign_agreement']))} and produced a pooled-effect correlation of {fnum(float(summary['biobank_omission_all_effect_correlation']))}; still, zero genes reached BH FDR <0.05 in that omission. Thus BioBank is influential for some individual estimates but is not “driving” a stable signature—no stable signature exists with or without it.

## Standardized-effect sensitivity

Unadjusted Hedges' g was estimable for all {estimable_n:,} primary-eligible genes. Relative to the primary pooled log2 effects:

- Direction agreement: **{pct(float(summary['standardized_all_direction_support']))}**.
- Effect correlation: **{fnum(float(summary['standardized_effect_correlation']))}**.
- Rank Spearman correlation: **{fnum(float(summary['standardized_rank_spearman']))}**.
- Standardized meta BH FDR <0.05: **{int((standardized['BH_FDR'] < 0.05).sum()):,}**.

This supports broad directional concordance across platform scales, but it cannot support or rescue stable membership. Unadjusted Hedges' g is not a substitute for covariate-adjusted regression coefficients.

## Biological-covariate sensitivity

Because `STABLE = 0`, the preregistered stable-gene `covariate_sensitivity.csv` intentionally contains zero data rows. There are no “core genes” for GA-, sex- or labor-sensitivity classification.

Across all modeled genes, used only as model diagnostics:

- GSE75010 primary versus GA/sex-adjusted effects: correlation {fnum(float(summary['GSE75010_GA_sex_effect_correlation']))}, sign agreement {pct(float(summary['GSE75010_GA_sex_sign_agreement']))}, median attenuation ratio {fnum(float(summary['GSE75010_GA_sex_median_attenuation']))}.
- GSE25906 batch-adjusted primary versus GA/sex-adjusted effects: correlation {fnum(float(summary['GSE25906_GA_sex_effect_correlation']))}, sign agreement {pct(float(summary['GSE25906_GA_sex_sign_agreement']))}, median attenuation ratio {fnum(float(summary['GSE25906_GA_sex_median_attenuation']))}.
- GSE25906 primary versus the full GA/sex/labor model: correlation {fnum(float(summary['GSE25906_labor_effect_correlation']))}, sign agreement {pct(float(summary['GSE25906_labor_sign_agreement']))}, median attenuation ratio {fnum(float(summary['GSE25906_labor_median_attenuation']))}.

Direction flips among near-zero effects are not interpreted as biological findings and do not trigger outcome-driven deletion.

## Heterogeneity interpretation

{i2_high_n:,} eligible genes had I² above 60%, but none met meta FDR <0.05. Heterogeneity may reflect severe-PE composition, unmeasured GA/FGR differences, platform dynamic range or cohort size. These explanations remain descriptive hypotheses; no uncontrolled post-hoc subtype search was performed.

## Analytical figures

Figures use algorithmic selection only. Because there are no stable genes, forest, heatmap and LOCO panels display an explicit no-stable-gene message rather than substituting familiar biological candidates. The volcano, standardized-effect concordance and heterogeneity-distribution previews remain populated.

## Answers to the ten Phase 1B questions

1. **Genes estimable in ≥4 cohorts:** {estimable_n:,}.
2. **Meta BH FDR <0.05:** {meta_fdr_n:,}.
3. **All STABLE criteria:** {stable_n:,}.
4. **ROBUST_BUT_HETEROGENEOUS:** {robust_hetero_n:,}.
5. **Is BioBank driving the signature?** No stable signature exists. Its omission preserved {pct(float(summary['biobank_omission_all_sign_agreement']))} of eligible-gene directions, but still produced zero BH-significant genes.
6. **LOCO robustness:** Not applicable to stable membership because no preliminary stable candidate exists; the formal LOCO table is correctly empty.
7. **Standardized-effect support:** Broad direction agreement is {pct(float(summary['standardized_all_direction_support']))}, but standardized FDR also yields zero significant genes and cannot rescue membership.
8. **Effect of GA/sex adjustment on core genes:** There are no core stable genes to test. All-gene diagnostic correlations are high in GSE75010 and lower in GSE25906, consistent with material cohort-specific confounding sensitivity.
9. **Enough stable genes for pathway/cell-state mapping?** No. Zero stable genes cannot support a frozen gene-level substrate for Phase 1C.
10. **Unresolved risks:** missing GA in three cohorts, unobserved BioBank processing batch, GSE25906 technical/bio-covariate structure, heterogeneous PE subtype/FGR composition, platform dynamic range and limited small-cohort precision.

## Final gate

**{gate}**

The preregistered gate defines zero stable genes as `NO_GO`. This does not prove that PE lacks placental transcriptomic abnormalities; it means the current six-cohort, genome-wide, multiple-testing-controlled framework did not establish a reproducible gene-level signature strong enough for Phase 1C. Thresholds were not relaxed after viewing outcomes, and Phase 1C was not started.

## Reproducibility

- Random seed: `20260809`.
- Primary models and thresholds: `config/phase1b_analysis.json`.
- Formal matrices and hashes: `results/01_phase1a1/formal_phase1b_matrix_registry.csv`.
- Software/session record: `results/01_phase1b/qc/phase1b_session_info.txt`.
- Result risks: `results/01_phase1b/qc/phase1b_risk_flags.csv`.
- Primary evidence sources: GEO records for [GSE75010](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75010), [GSE30186](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30186), [GSE10588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE10588), [GSE24129](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24129), [GSE25906](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25906) and [GSE43942](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43942).
"""
    DOC.write_text(text, encoding="utf-8", newline="\n")
    print(f"Phase 1B report written: {DOC}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
