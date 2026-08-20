#!/usr/bin/env python3
"""Stream the published UMI table into donor-by-harmonized-cell-type counts."""

from __future__ import annotations

import csv
import gzip
import json
import time
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from freeze_admati_metadata import ZIP, harmonize, read_metadata

ROOT = Path(__file__).resolve().parents[2]
CONFIG = json.loads((ROOT / "config/phase2a_analysis.json").read_text(encoding="utf-8"))
OUT = ROOT / "results/02_phase2a"
INTERIM = ROOT / "data/interim/phase2a"
MATRIX = INTERIM / "admati_harmonized_pseudobulk_counts.csv.gz"


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    start = time.time()
    md, _, _, _ = read_metadata()
    patients = pd.read_csv(OUT / "metadata/patient_registry.csv", dtype=str).set_index("patient_id")
    donors = sorted(patients.index)
    celltypes = sorted({harmonize(x) for x in md["celltype"]})
    strata = [(d, c) for c in celltypes for d in donors]
    stratum_index = {x: i for i, x in enumerate(strata)}
    cell_group = np.fromiter((stratum_index[(d, harmonize(c))] for d, c in zip(md["donorID"], md["celltype"])), dtype=np.int32, count=len(md["donorID"]))
    cell_counts = np.bincount(cell_group, minlength=len(strata)).astype(np.int64)
    expected_umi = np.bincount(cell_group, weights=np.asarray(md["total_molecules"], dtype=np.int64), minlength=len(strata)).round().astype(np.int64)
    INTERIM.mkdir(parents=True, exist_ok=True)
    gene_n = 0
    observed_umi = np.zeros(len(strata), dtype=np.int64)
    with zipfile.ZipFile(ZIP) as archive, archive.open(archive.infolist()[0]) as handle, gzip.open(MATRIX, "wt", encoding="utf-8", newline="") as out:
        for _ in range(23):
            handle.readline()
        writer = csv.writer(out)
        writer.writerow(["gene"] + [f"{d}__{c}" for d, c in strata])
        for raw in handle:
            gene_raw, values_raw = raw.rstrip(b"\r\n").split(b"\t", 1)
            values = np.fromstring(values_raw, sep="\t", dtype=np.int64)
            if values.size != cell_group.size:
                raise AssertionError(f"{gene_raw!r}: {values.size} values, expected {cell_group.size}")
            sums = np.bincount(cell_group, weights=values, minlength=len(strata)).round().astype(np.int64)
            observed_umi += sums
            writer.writerow([gene_raw.decode("utf-8")] + sums.tolist())
            gene_n += 1
            if gene_n % 5000 == 0:
                print(f"aggregated {gene_n} genes in {time.time()-start:.1f}s", flush=True)
    if not np.array_equal(observed_umi, expected_umi):
        bad = np.flatnonzero(observed_umi != expected_umi)
        raise AssertionError(f"UMI reconciliation failed for {len(bad)} strata")
    min_cells = int(CONFIG["eligibility"]["minimum_cells_per_patient_celltype"])
    min_umi = int(CONFIG["eligibility"]["minimum_pseudobulk_library_umi"])
    min_patients = int(CONFIG["eligibility"]["minimum_patients_per_group"])
    registry = []
    for i, (donor, celltype) in enumerate(strata):
        p = patients.loc[donor]
        qualified = cell_counts[i] >= min_cells and observed_umi[i] >= min_umi
        reasons = []
        if cell_counts[i] < min_cells:
            reasons.append(f"CELL_COUNT_LT_{min_cells}")
        if observed_umi[i] < min_umi:
            reasons.append(f"LIBRARY_UMI_LT_{min_umi}")
        registry.append({"dataset": "Admati_2023_FIGSHARE", "patient_id": donor, "harmonized_annotation": celltype, "disease_status": p.disease_status, "pe_subtype_or_control_group": p.pe_subtype_or_control_group, "cell_count": int(cell_counts[i]), "pseudobulk_library_umi": int(observed_umi[i]), "stratum_passes_numeric_thresholds": "YES" if qualified else "NO", "exclusion_reason": "" if qualified else ";".join(reasons), "matrix_column": f"{donor}__{celltype}", "source_url": "https://doi.org/10.6084/m9.figshare.23264102.v1"})
    write_csv(OUT / "pseudobulk/pseudobulk_registry.csv", registry, list(registry[0]))
    contrasts = {"EOPE": ("EOPE", "EARLY_CONTROL"), "LOPE": ("LOPE", "LATE_CONTROL"), "COMBINED_PE_SECONDARY": ("PE", "CONTROL")}
    eligibility = []
    for contrast, groups in contrasts.items():
        for celltype in celltypes:
            subset = [r for r in registry if r["harmonized_annotation"] == celltype]
            if contrast == "COMBINED_PE_SECONDARY":
                count_by_group = {g: sum(r["stratum_passes_numeric_thresholds"] == "YES" and r["disease_status"] == g for r in subset) for g in groups}
            else:
                count_by_group = {g: sum(r["stratum_passes_numeric_thresholds"] == "YES" and r["pe_subtype_or_control_group"] == g for r in subset) for g in groups}
            contrast_pass = all(n >= min_patients for n in count_by_group.values())
            for r in subset:
                belongs = r["disease_status"] in groups if contrast == "COMBINED_PE_SECONDARY" else r["pe_subtype_or_control_group"] in groups
                include = belongs and r["stratum_passes_numeric_thresholds"] == "YES" and contrast_pass
                reason = ""
                if not belongs:
                    reason = "NOT_IN_CONTRAST"
                elif r["stratum_passes_numeric_thresholds"] != "YES":
                    reason = r["exclusion_reason"]
                elif not contrast_pass:
                    reason = "CELLTYPE_FAILS_MIN_PATIENTS_PER_GROUP"
                eligibility.append({"contrast": contrast, "patient_id": r["patient_id"], "harmonized_annotation": celltype, "group": r["disease_status"] if contrast == "COMBINED_PE_SECONDARY" else r["pe_subtype_or_control_group"], "cell_count": r["cell_count"], "pseudobulk_library_umi": r["pseudobulk_library_umi"], "group_qualified_patient_n": count_by_group.get(r["disease_status"] if contrast == "COMBINED_PE_SECONDARY" else r["pe_subtype_or_control_group"], 0), "celltype_contrast_eligible": "YES" if contrast_pass else "NO", "include_in_contrast": "YES" if include else "NO", "exclusion_reason": reason, "thresholds": f"patients/group>={min_patients};cells/stratum>={min_cells};UMI>={min_umi}", "source_url": "config/phase2a_analysis.json"})
    write_csv(OUT / "metadata/pseudobulk_eligibility.csv", eligibility, list(eligibility[0]))
    summary = {"gene_rows": gene_n, "matrix_columns": len(strata), "harmonized_celltypes": len(celltypes), "patients": len(donors), "source_cell_count": len(cell_group), "source_total_umi": int(expected_umi.sum()), "pseudobulk_total_umi": int(observed_umi.sum()), "elapsed_seconds": round(time.time()-start, 2), "matrix": str(MATRIX)}
    (INTERIM / "pseudobulk_build_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
