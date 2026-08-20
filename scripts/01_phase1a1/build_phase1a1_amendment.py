#!/usr/bin/env python3
"""Build the Phase 1A.1 preprocessing freeze without disease-outcome analysis."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
P1A = ROOT / "results" / "01_phase1a"
OUT = ROOT / "results" / "01_phase1a1"
INTERIM = ROOT / "data" / "interim" / "phase1a1"
FORMAL = INTERIM / "formal_matrices"
SUBMITTED = INTERIM / "submitted"
RECON = INTERIM / "reconstructed"

GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
SOURCES = {
    "GSE75010_BIOBANK": GEO + "GSE75010|https://pubmed.ncbi.nlm.nih.gov/27160201/",
    "GSE30186": GEO + "GSE30186|https://pubmed.ncbi.nlm.nih.gov/22702245/",
    "GSE10588": GEO + "GSE10588|https://pubmed.ncbi.nlm.nih.gov/19249095/|https://www.bioconductor.org/packages/devel/bioc/vignettes/ABarray/inst/doc/ABarray.pdf",
    "GSE24129": GEO + "GSE24129|https://pmc.ncbi.nlm.nih.gov/articles/PMC3199758/",
    "GSE25906": GEO + "GSE25906|https://pmc.ncbi.nlm.nih.gov/articles/PMC3039036/",
    "GSE43942": GEO + "GSE43942|https://pubmed.ncbi.nlm.nih.gov/23544093/|https://www.bioconductor.org/packages/release/bioc/manuals/oligo/man/oligo.pdf",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gzip_content_sha256(path: Path) -> str:
    """Hash decompressed matrix bytes so gzip timestamps do not affect the registry."""
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_matrix(path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(path, sep="\t", index_col=0)
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    return matrix


def write_matrix(path: Path, matrix: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # Fix gzip mtime so repeated builds produce byte-identical matrices.
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                matrix.to_csv(text, sep="\t", index=True, index_label="feature_id", lineterminator="\n", float_format="%.10g")


def paired_corr(a: np.ndarray, b: np.ndarray, axis: int) -> np.ndarray:
    a = a - np.nanmean(a, axis=axis, keepdims=True)
    b = b - np.nanmean(b, axis=axis, keepdims=True)
    num = np.nansum(a * b, axis=axis)
    den = np.sqrt(np.nansum(a * a, axis=axis) * np.nansum(b * b, axis=axis))
    return np.divide(num, den, out=np.full_like(num, np.nan, dtype=float), where=den > 0)


def rank_array(x: np.ndarray, axis: int) -> np.ndarray:
    return pd.DataFrame(x).rank(axis=axis, method="average", na_option="keep").to_numpy(dtype=float)


def vector_spearman(a: np.ndarray, b: np.ndarray) -> float:
    ar = rank_array(np.asarray(a).reshape(-1, 1), axis=0).ravel()
    br = rank_array(np.asarray(b).reshape(-1, 1), axis=0).ravel()
    return float(paired_corr(ar.reshape(1, -1), br.reshape(1, -1), axis=1)[0])


def sample_distances(x: np.ndarray) -> np.ndarray:
    # Condensed Euclidean distance vector, equivalent to scipy.pdist.
    gram = x.T @ x
    squared = np.maximum(np.diag(gram)[:, None] + np.diag(gram)[None, :] - 2 * gram, 0)
    upper = np.triu_indices(squared.shape[0], 1)
    return np.sqrt(squared[upper])


def matrix_comparison(submitted: pd.DataFrame, reconstructed: pd.DataFrame) -> dict[str, object]:
    common_features = submitted.index.intersection(reconstructed.index, sort=False)
    common_samples = submitted.columns.intersection(reconstructed.columns, sort=False)
    a = submitted.loc[common_features, common_samples].to_numpy(dtype=float)
    b = reconstructed.loc[common_features, common_samples].to_numpy(dtype=float)
    complete = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
    a, b = a[complete], b[complete]

    sample_pearson = paired_corr(a, b, axis=0)
    sample_spearman = paired_corr(rank_array(a, axis=0), rank_array(b, axis=0), axis=0)
    gene_pearson = paired_corr(a, b, axis=1)
    gene_spearman = paired_corr(rank_array(a, axis=1), rank_array(b, axis=1), axis=1)

    # Outcome-blind PCA-structure surrogate: correlation between all pairwise
    # sample distances after feature-wise centering/scaling on the 2,000 most
    # variable common features. No PE/control label enters feature selection.
    av = np.nanvar(a, axis=1)
    bv = np.nanvar(b, axis=1)
    top = np.argsort(np.nan_to_num(av + bv, nan=-np.inf))[-min(2000, len(av)):]
    def scaled_sample_distance(x: np.ndarray) -> np.ndarray:
        y = x[top]
        sd = np.nanstd(y, axis=1, ddof=1)
        keep = sd > 0
        y = (y[keep] - np.nanmean(y[keep], axis=1, keepdims=True)) / sd[keep, None]
        return sample_distances(y)
    pca_structure = vector_spearman(scaled_sample_distance(a), scaled_sample_distance(b))

    quantiles = np.linspace(0.01, 0.99, 99)
    qa = np.nanquantile(a, quantiles, axis=0).ravel()
    qb = np.nanquantile(b, quantiles, axis=0).ravel()
    distribution = vector_spearman(qa, qb)
    return {
        "common_feature_n": int(a.shape[0]),
        "common_sample_n": int(a.shape[1]),
        "median_sample_pearson": float(np.nanmedian(sample_pearson)),
        "median_sample_spearman": float(np.nanmedian(sample_spearman)),
        "pca_structure_distance_spearman": float(pca_structure),
        "distribution_quantile_spearman": float(distribution),
        "median_gene_pearson": float(np.nanmedian(gene_pearson)),
        "median_gene_spearman": float(np.nanmedian(gene_spearman)),
    }


def collapse_to_gene(dataset: str, feature: pd.DataFrame, mapping_rows: list[dict[str, str]]) -> pd.DataFrame:
    mapping = {
        row["original_probe_id"]: row["mapped_symbol"]
        for row in mapping_rows
        if dataset in row["dataset"].split("|") and row["mapping_status"] == "MAPPED_UNAMBIGUOUS"
    }
    keep = [x for x in feature.index if x in mapping]
    temp = feature.loc[keep].copy()
    temp.insert(0, "mapped_symbol", [mapping[x] for x in keep])
    return temp.groupby("mapped_symbol", sort=True).median(numeric_only=True)


def make_formal_matrices(mapping_rows: list[dict[str, str]]) -> tuple[dict[str, dict[str, object]], dict[str, dict[str, object]]]:
    comparisons: dict[str, dict[str, object]] = {}
    formal: dict[str, dict[str, object]] = {}
    decisions = {
        "GSE30186": {
            "formal": "reconstructed",
            "transform": "non-normalized BeadStudio AVG_Signal -> limma normexp(saddle, offset=16) -> quantile -> log2",
            "unit": "reconstructed log2 Illumina intensity",
            "reason": "The submitted matrix contains negative background-adjusted values; the Phase 1A cohort-minimum shift is arbitrary. The public non-normalized summary is rebuildable, although negative-control intensities are absent.",
        },
        "GSE10588": {
            "formal": "submitted",
            "transform": "log2(submitted quantile-normalized VALUE); no +1 pseudocount",
            "unit": "submitted ABarray quantile-normalized log2 intensity",
            "reason": "GEO defines VALUE as the submitted quantile-normalized signal. Raw Signal and quality fields are present, but historical ABarray flag filtering/imputation is not fully recoverable; the submitted matrix is the more faithful formal input.",
        },
        "GSE43942": {
            "formal": "submitted",
            "transform": "log2(submitted positive RMA/quantile-normalized VALUE); no +1 pseudocount",
            "unit": "submitted NimbleScan/RMA log2 intensity",
            "reason": "PAIR PM values and design files permit a transparent modern RMA-like reconstruction, but archived PAIR is not a turnkey input to current oligo and the exact historical NimbleScan 2.5 implementation is not reproducibly certified. Use reconstruction as preprocessing sensitivity.",
        },
    }
    recon_names = {
        "GSE30186": "GSE30186_normexp16_quantile_log2_feature_level.tsv.gz",
        "GSE10588": "GSE10588_raw_signal_quantile_log2_feature_level.tsv.gz",
        "GSE43942": "GSE43942_pair_rma_quantile_medianpolish_log2_transcript_level.tsv.gz",
    }
    for dataset in decisions:
        submitted_linear = read_matrix(SUBMITTED / f"{dataset}_submitted_linear.tsv.gz")
        reconstructed = read_matrix(RECON / recon_names[dataset])
        if dataset == "GSE30186":
            submitted_candidate = submitted_linear
        else:
            if (submitted_linear <= 0).any().any():
                raise AssertionError(f"{dataset} submitted matrix is not strictly positive")
            submitted_candidate = np.log2(submitted_linear)
        comparisons[dataset] = matrix_comparison(submitted_candidate, reconstructed)
        feature = reconstructed if decisions[dataset]["formal"] == "reconstructed" else submitted_candidate
        feature_path = FORMAL / f"{dataset}_formal_feature_level.tsv.gz"
        gene_path = FORMAL / f"{dataset}_formal_gene_level.tsv.gz"
        write_matrix(feature_path, feature)
        gene = collapse_to_gene(dataset, feature, mapping_rows)
        write_matrix(gene_path, gene)
        formal[dataset] = {
            **decisions[dataset],
            "feature_path": feature_path.relative_to(ROOT).as_posix(),
            "feature_sha256": sha256(feature_path),
            "gene_path": gene_path.relative_to(ROOT).as_posix(),
            "gene_sha256": sha256(gene_path),
            "feature_n": len(feature),
            "gene_n": len(gene),
            "sample_n": feature.shape[1],
            "reconstructed_path": (RECON / recon_names[dataset]).relative_to(ROOT).as_posix(),
            "reconstructed_content_sha256": gzip_content_sha256(RECON / recon_names[dataset]),
        }
    return comparisons, formal


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    FORMAL.mkdir(parents=True, exist_ok=True)
    mapping_rows = read_csv(P1A / "gene_mapping_registry.csv")
    processing = {row["dataset"]: row for row in read_csv(P1A / "bulk_processing_registry.csv")}
    comparisons, formal = make_formal_matrices(mapping_rows)

    comparison_rows = []
    audit = {
        "GSE30186": ("YES_FROM_NON_NORMALIZED_SUMMARY", "NO_NEGATIVE_CONTROL_DEFINITIONS_ONLY", "limma normexp saddle offset 16; quantile normalization; log2", "RECONSTRUCTED_FORMAL"),
        "GSE10588": ("YES_WITH_HISTORICAL_LIMITATIONS", "YES_RAW_SIGNAL_FLAGS_AND_CONTROLS_PRESENT", "ABI regular-probe Signal; natural-scale quantile; log2; flag variants audited", "SUBMITTED_FORMAL_RAW_SENSITIVITY"),
        "GSE43942": ("YES_RMA_LIKE_NOT_EXACT_VENDOR_CERTIFIED", "YES_PAIR_PM_NDF_NGD_AND_CONTROLS_PRESENT", "PAIR PM; RMA background; natural-scale quantile; log2; transcript median polish", "SUBMITTED_FORMAL_RAW_SENSITIVITY"),
    }
    for dataset, metrics in comparisons.items():
        feasible, controls, method, selection = audit[dataset]
        reason = formal[dataset]["reason"]
        comparison_rows.append({
            "dataset": dataset,
            "raw_reprocessing_feasibility": feasible,
            "raw_control_probe_sufficiency": controls,
            "scientifically_appropriate_raw_method": method,
            "submitted_candidate_transform": "untransformed submitted background-adjusted VALUE (comparison only)" if dataset == "GSE30186" else "log2(submitted VALUE), without +1",
            "reconstructed_matrix_path": formal[dataset]["reconstructed_path"],
            "reconstructed_matrix_content_sha256": formal[dataset]["reconstructed_content_sha256"],
            **{key: f"{value:.6f}" if isinstance(value, float) else value for key, value in metrics.items()},
            "comparison_is_outcome_blind": "YES_NO_DISEASE_LABEL_USED",
            "formal_selection": selection,
            "selection_reason": reason,
            "source": SOURCES[dataset],
        })
    comparison_fields = list(comparison_rows[0])
    write_csv(OUT / "preprocessing_comparison_registry.csv", comparison_rows, comparison_fields)

    formal_rows = []
    for dataset in ["GSE75010_BIOBANK", "GSE30186", "GSE10588", "GSE24129", "GSE25906", "GSE43942"]:
        if dataset in formal:
            x = formal[dataset]
            row = {
                "dataset": dataset, "formal_matrix_choice": x["formal"].upper(), "expression_unit": x["unit"],
                "formal_transform": x["transform"], "feature_matrix_path": x["feature_path"], "feature_matrix_sha256": x["feature_sha256"],
                "gene_matrix_path": x["gene_path"], "gene_matrix_sha256": x["gene_sha256"], "feature_n": x["feature_n"],
                "mapped_gene_n": x["gene_n"], "sample_n": x["sample_n"], "raw_reconstruction_role": "FORMAL" if dataset == "GSE30186" else "PREPROCESSING_SENSITIVITY",
                "arbitrary_shift_log_in_formal_matrix": "NO", "cross_cohort_batch_correction": "PROHIBITED", "decision_reason": x["reason"], "source": SOURCES[dataset],
            }
        else:
            p = processing[dataset]
            row = {
                "dataset": dataset, "formal_matrix_choice": "PHASE1A_SUBMITTED_FROZEN", "expression_unit": p["expression_unit"],
                "formal_transform": "none; already defensible submitted log scale", "feature_matrix_path": p["feature_matrix_path"], "feature_matrix_sha256": p["feature_matrix_sha256"],
                "gene_matrix_path": p["gene_matrix_path"], "gene_matrix_sha256": p["gene_matrix_sha256"], "feature_n": p["submitted_feature_n"],
                "mapped_gene_n": p["mapped_gene_n"], "sample_n": {"GSE75010_BIOBANK": 157, "GSE24129": 24, "GSE25906": 60}[dataset],
                "raw_reconstruction_role": "NOT_NEEDED", "arbitrary_shift_log_in_formal_matrix": "NO", "cross_cohort_batch_correction": "PROHIBITED",
                "decision_reason": "Submitted processed values are explicitly normalized log-scale measurements; unnecessary reprocessing would not improve provenance.", "source": SOURCES[dataset],
            }
        formal_rows.append(row)
    write_csv(OUT / "formal_phase1b_matrix_registry.csv", formal_rows, list(formal_rows[0]))

    freeze_rows = read_csv(P1A / "bulk_sample_freeze.csv")
    freeze_lookup = {(row["dataset"], row["GSM/sample ID"]): row for row in freeze_rows}
    flagged = [row for row in read_csv(P1A / "bulk_sample_qc.csv") if row["qc_flag"] == "FLAG_REVIEW"]
    qc_rows = []
    for old in flagged:
        dataset, sample = old["dataset"], old["sample_id"]
        frozen = freeze_lookup[(dataset, sample)]
        if dataset == "GSE25906":
            interpretation = "BATCH_A_ASSOCIATED_QC_STRUCTURE_NOT_15_INDEPENDENT_OUTLIERS"
            group = "GSE25906_BATCH_A_15_FLAGGED_BOTH_PE_AND_CONTROL"
            model = "Retain; disease model considers batch + GA + fetal_sex, with labor sensitivity if estimable"
        elif dataset == "GSE30186":
            interpretation = "TWO_PE_SAMPLES_REQUIRE_PROVENANCE_REVIEW_NO_FAILURE_EVIDENCE"
            group = "GSE30186_TWO_FLAGGED_PE"
            model = "Retain; influence diagnostics only; no disease-separation-based exclusion"
        else:
            interpretation = "SINGLE_LOW_CORRELATION_CONTROL_REVIEW_TRIGGER_NO_FAILURE_EVIDENCE"
            group = "GSE75010_SINGLE_FLAGGED_CONTROL"
            model = "Retain; influence diagnostics only; no disease-separation-based exclusion"
        qc_rows.append({
            "dataset": dataset, "sample_id": sample, "PE_control": frozen["PE/control"], "batch": frozen["batch"],
            "phase1a_qc_flag": old["qc_flag"], "phase1a_median_sample_correlation": old["median_sample_correlation"],
            "reinterpretation_group": group, "phase1a1_interpretation": interpretation,
            "retain_phase1b": "YES", "exclusion_reason": "NONE_NO_INDEPENDENT_TECHNICAL_OR_PROVENANCE_FAILURE",
            "model_or_review_action": model, "source": f"results/01_phase1a/bulk_sample_qc.csv|{frozen['source']}",
        })
    write_csv(OUT / "qc_flag_reinterpretation.csv", qc_rows, list(qc_rows[0]))

    model_specs = [
        ("GSE75010_BIOBANK", "~ disease + GA_c + fetal_sex", "GA_c (cohort-mean-centered GA);fetal_sex", "maternal_age;ancestry;BMI;parity in separately reported complete-case sensitivity", "FGR/IUGR;labor;delivery_mode;birth outcomes excluded from primary adjustment as potential downstream/mediator variables; processing batch unavailable", "~ disease + GA_c"),
        ("GSE30186", "~ disease", "disease", "maternal_age as a single-covariate sensitivity only if full rank", "ancestry constant; GA unavailable; additional covariates excluded to preserve degrees of freedom", "~ disease"),
        ("GSE10588", "~ disease", "disease", "NONE", "sample-level GA and clinical covariates unavailable", "~ disease"),
        ("GSE24129", "~ disease", "disease", "NONE", "isolated-FGR samples excluded from primary contrast; caesarean status constant; sample-level GA unavailable", "~ disease"),
        ("GSE25906", "~ disease + batch + GA_c + fetal_sex", "batch;GA_c (cohort-mean-centered GA);fetal_sex", "labor induction added only as a prespecified sensitivity if design remains full rank", "No additional weak covariates; 15 QC flags are represented as one batch-A structure, not deletion indicators", "~ disease + batch + GA_c"),
        ("GSE43942", "~ disease", "disease", "NONE", "caesarean status constant; sample-level GA unavailable; n=12 prohibits weak multi-covariate adjustment", "~ disease"),
    ]

    def ga_number(value: str) -> float:
        if "w+" in value and value.endswith("d"):
            weeks, days = value[:-1].split("w+")
            return float(weeks) + float(days) / 7.0
        return float(value)

    def design_diagnostics(dataset: str) -> dict[str, object]:
        cohort = [row for row in freeze_rows if row["dataset"] == dataset and row["include_phase1b"] == "YES"]
        columns = [np.ones(len(cohort)), np.array([1.0 if row["PE/control"] == "PE" else 0.0 for row in cohort])]
        labels = ["intercept", "disease"]
        categoricals: list[tuple[str, list[str]]] = []
        if dataset in {"GSE75010_BIOBANK", "GSE25906"}:
            ga = np.array([ga_number(row["GA"]) for row in cohort])
            columns.append(ga - ga.mean())
            labels.append("GA_c")
        if dataset == "GSE25906":
            columns.append(np.array([1.0 if row["batch"] == "B" else 0.0 for row in cohort]))
            labels.append("batch_B")
            categoricals.append(("batch", [row["batch"] for row in cohort]))
        if dataset in {"GSE75010_BIOBANK", "GSE25906"}:
            columns.append(np.array([1.0 if row["fetal sex"] == "M" else 0.0 for row in cohort]))
            labels.append("fetal_sex_M")
            categoricals.append(("fetal_sex", [row["fetal sex"] for row in cohort]))
        matrix = np.column_stack(columns)
        zero_cells = []
        diseases = [row["PE/control"] for row in cohort]
        for name, values in categoricals:
            for disease in sorted(set(diseases)):
                for level in sorted(set(values)):
                    if not any(d == disease and v == level for d, v in zip(diseases, values)):
                        zero_cells.append(f"{name}:{disease}x{level}")
        rank = int(np.linalg.matrix_rank(matrix))
        vifs = []
        for j in range(1, matrix.shape[1]):
            y = matrix[:, j]
            others = np.delete(matrix, j, axis=1)
            fitted = others @ np.linalg.lstsq(others, y, rcond=None)[0]
            total = np.sum((y - y.mean()) ** 2)
            residual = np.sum((y - fitted) ** 2)
            r_squared = 1.0 - residual / total if total > 0 else 1.0
            vifs.append(1.0 / max(1.0 - r_squared, np.finfo(float).eps))
        return {
            "primary_complete_case_n": len(cohort), "primary_design_columns": ";".join(labels),
            "primary_design_rank": rank, "primary_design_column_n": matrix.shape[1],
            "primary_condition_number_centered": f"{np.linalg.cond(matrix):.4f}",
            "primary_max_vif": f"{max(vifs):.4f}",
            "disease_categorical_zero_cells": ";".join(zero_cells) if zero_cells else "NONE",
            "primary_estimability_decision": "FULL_RANK_NO_ZERO_CELL_VIF_LT5" if rank == matrix.shape[1] and not zero_cells and max(vifs) < 5 else "USE_FALLBACK_OR_REVIEW",
        }

    model_rows = []
    for d, f, m, o, e, fb in model_specs:
        model_rows.append({
            "dataset": d, "analytical_n": {"GSE75010_BIOBANK": 116, "GSE30186": 12, "GSE10588": 43, "GSE24129": 16, "GSE25906": 60, "GSE43942": 12}[d],
            "primary_formula": f, "mandatory_covariates": m, "optional_covariates": o, "excluded_covariates_and_reason": e, "fallback_formula": fb,
            **design_diagnostics(d),
            "estimability_checks": "model.matrix full rank; disease-by-covariate cross-tabs; zero-variance levels; pairwise correlation; condition number; VIF target <5; retain simpler fallback if unstable",
            "missing_covariate_rule": "NO_IMPUTATION; optional complete-case loss must be enumerated", "outcome_blind_formula_lock": "YES_BEFORE_DEG",
            "source": SOURCES[d],
        })
    write_csv(OUT / "cohort_model_formula_registry.csv", model_rows, list(model_rows[0]))

    risks = [
        ("P1A1-R01", "HIGH", "GSE30186", "Non-normalized AVG_Signal is public but negative-control intensities are absent; reconstruction is control-free.", "Use the frozen normexp/quantile/log2 route; report submitted-matrix sensitivity; do not restore cohort-minimum shift.", "RESTRICTED_RESOLVED_FOR_GATE"),
        ("P1A1-R02", "MODERATE", "GSE10588", "Historical ABarray flag filtering/imputation details are not completely recoverable from the public workflow.", "Use log2 submitted quantile-normalized VALUE formally; raw reconstruction only as sensitivity.", "RESOLVED_BY_FORMAL_CHOICE"),
        ("P1A1-R03", "MODERATE", "GSE43942", "Modern RMA-like PAIR reconstruction is transparent but not certified identical to NimbleScan 2.5.", "Use log2 submitted RMA VALUE formally; use rebuilt PAIR matrix as preprocessing sensitivity.", "RESOLVED_BY_FORMAL_CHOICE"),
        ("P1A1-R04", "HIGH", "GSE25906", "Fifteen review flags cluster in batch A and include both PE and control samples.", "Retain all; adjust for batch with GA and fetal sex when estimable; perform influence review without label-driven deletion.", "OPEN_PHASE1B_RESTRICTION"),
        ("P1A1-R05", "HIGH", "GSE30186;GSE10588;GSE43942", "Small or metadata-poor cohorts lack sample-level GA and cannot support complex adjustment.", "Use parsimonious frozen formulas and interpret cross-cohort heterogeneity/LOO results.", "OPEN_PHASE1B_RESTRICTION"),
        ("P1A1-R06", "MODERATE", "CROSS_PLATFORM", "Microarray effect scales have heterogeneous dynamic ranges and preprocessing histories.", "No universal absolute meta-log2FC threshold; run standardized-effect-size sensitivity without replacing the primary gate.", "METHOD_LOCKED"),
        ("P1A1-R07", "MODERATE", "GSE75010_BIOBANK", "Sample-level processing groups are not public.", "Use GA/fetal-sex model and report latent-batch limitation; do not infer batch labels.", "OPEN_PHASE1B_RESTRICTION"),
        ("P1A1-R08", "CRITICAL_GATE", "META_ANALYSIS", "A stable label is impossible if fewer than four independent core cohorts yield estimable effects or LOO/heterogeneity gates fail.", "Return no stable label for that gene; do not relax the preregistered rule after outcome inspection.", "FROZEN_STOP_RULE"),
    ]
    risk_rows = [{"risk_id": i, "severity": s, "dataset": d, "risk": r, "required_action": a, "status": st, "source": "docs/PHASE1A1_PREPROCESSING_AMENDMENT.md|" + (SOURCES.get(d, "docs/PHASE1B_STATISTICAL_ANALYSIS_PLAN.md"))} for i, s, d, r, a, st in risks]
    write_csv(OUT / "phase1a1_risk_flags.csv", risk_rows, list(risk_rows[0]))

    print("Phase 1A.1 registries and formal matrices built")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
