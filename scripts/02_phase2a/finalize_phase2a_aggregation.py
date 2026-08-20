#!/usr/bin/env python3
"""Reconcile streamed counts and freeze pseudobulk eligibility registries."""

from __future__ import annotations

import csv
import gzip
import json
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONFIG = json.loads((ROOT / "config/phase2a_analysis.json").read_text(encoding="utf-8"))
OUT = ROOT / "results/02_phase2a"
INTERIM = ROOT / "data/interim/phase2a"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main() -> int:
    strata = pd.read_csv(INTERIM / "pseudobulk_strata.tsv", sep="\t")
    totals = pd.read_csv(INTERIM / "pseudobulk_observed_totals.tsv", sep="\t")
    check = strata.merge(totals, on=["group_id", "matrix_column"], validate="one_to_one")
    exact_reconciliation = bool((check.expected_total_umi == check.observed_total_umi).all())
    total_molecule_correlation = float(check.expected_total_umi.corr(check.observed_total_umi))
    if total_molecule_correlation < 0.75:
        raise AssertionError("Published total_molecules is insufficiently concordant with streamed gene counts")
    patients = pd.read_csv(OUT / "metadata/patient_registry.csv", dtype=str).set_index("patient_id")
    min_cells = int(CONFIG["eligibility"]["minimum_cells_per_patient_celltype"])
    min_umi = int(CONFIG["eligibility"]["minimum_pseudobulk_library_umi"])
    min_patients = int(CONFIG["eligibility"]["minimum_patients_per_group"])
    registry = []
    for row in check.itertuples(index=False):
        p = patients.loc[row.patient_id]
        qualified = row.expected_cell_count >= min_cells and row.observed_total_umi >= min_umi
        reasons = []
        if row.expected_cell_count < min_cells: reasons.append(f"CELL_COUNT_LT_{min_cells}")
        if row.observed_total_umi < min_umi: reasons.append(f"LIBRARY_UMI_LT_{min_umi}")
        registry.append({"dataset": "Admati_2023_FIGSHARE", "patient_id": row.patient_id, "harmonized_annotation": row.harmonized_annotation, "disease_status": p.disease_status, "pe_subtype_or_control_group": p.pe_subtype_or_control_group, "cell_count": int(row.expected_cell_count), "published_total_molecules_sum": int(row.expected_total_umi), "pseudobulk_library_umi": int(row.observed_total_umi), "count_to_published_molecule_ratio": round(row.observed_total_umi / row.expected_total_umi, 6) if row.expected_total_umi else "", "stratum_passes_numeric_thresholds": "YES" if qualified else "NO", "exclusion_reason": "" if qualified else ";".join(reasons), "matrix_column": row.matrix_column, "source_url": "https://doi.org/10.6084/m9.figshare.23264102.v1"})
    write_csv(OUT / "pseudobulk/pseudobulk_registry.csv", registry)
    contrasts = {"EOPE": ("EOPE", "EARLY_CONTROL"), "LOPE": ("LOPE", "LATE_CONTROL"), "COMBINED_PE_SECONDARY": ("PE", "CONTROL")}
    eligibility = []
    for contrast, groups in contrasts.items():
        for celltype in sorted(check.harmonized_annotation.unique()):
            subset = [r for r in registry if r["harmonized_annotation"] == celltype]
            group_field = "disease_status" if contrast == "COMBINED_PE_SECONDARY" else "pe_subtype_or_control_group"
            count_by_group = {g: sum(r["stratum_passes_numeric_thresholds"] == "YES" and r[group_field] == g for r in subset) for g in groups}
            contrast_pass = all(n >= min_patients for n in count_by_group.values())
            for r in subset:
                belongs = r[group_field] in groups
                include = belongs and r["stratum_passes_numeric_thresholds"] == "YES" and contrast_pass
                if not belongs: reason = "NOT_IN_CONTRAST"
                elif r["stratum_passes_numeric_thresholds"] != "YES": reason = r["exclusion_reason"]
                elif not contrast_pass: reason = "CELLTYPE_FAILS_MIN_PATIENTS_PER_GROUP"
                else: reason = ""
                eligibility.append({"contrast": contrast, "patient_id": r["patient_id"], "harmonized_annotation": celltype, "group": r[group_field], "cell_count": r["cell_count"], "pseudobulk_library_umi": r["pseudobulk_library_umi"], "group_qualified_patient_n": count_by_group.get(r[group_field], 0), "celltype_contrast_eligible": "YES" if contrast_pass else "NO", "include_in_contrast": "YES" if include else "NO", "exclusion_reason": reason, "thresholds": f"patients/group>={min_patients};cells/stratum>={min_cells};UMI>={min_umi}", "source_url": "config/phase2a_analysis.json"})
    write_csv(OUT / "metadata/pseudobulk_eligibility.csv", eligibility)
    matrix = INTERIM / "admati_harmonized_pseudobulk_counts.csv.gz"
    with gzip.open(matrix, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split(",")
        gene_n = sum(1 for _ in handle)
    uncompressed_matrix = INTERIM / "admati_harmonized_pseudobulk_counts.csv"
    with gzip.open(matrix, "rb") as source, uncompressed_matrix.open("wb") as target:
        shutil.copyfileobj(source, target, length=1024 * 1024)
    summary = {"gene_rows": gene_n, "matrix_columns": len(header)-1, "harmonized_celltypes": check.harmonized_annotation.nunique(), "patients": check.patient_id.nunique(), "source_cell_count": int(check.expected_cell_count.sum()), "published_total_molecules_sum": int(check.expected_total_umi.sum()), "pseudobulk_count_matrix_sum": int(check.observed_total_umi.sum()), "total_molecule_stratum_correlation": total_molecule_correlation, "umi_reconciliation": "EXACT" if exact_reconciliation else "NOT_EXACT_FLAGGED; COUNT_MATRIX_COLUMN_SUM_USED_AS_LIBRARY_SIZE", "matrix_gzip": str(matrix), "matrix_uncompressed_for_R": str(uncompressed_matrix)}
    (INTERIM / "pseudobulk_build_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__": raise SystemExit(main())
