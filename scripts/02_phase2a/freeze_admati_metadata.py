#!/usr/bin/env python3
"""Freeze Admati donor identities and published annotation mapping before outcomes."""

from __future__ import annotations

import csv
import hashlib
import json
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ZIP = ROOT / "data/raw/phase0b/sc_PE_allcells_with_metadata_29-May-2023.txt.zip"
OUT = ROOT / "results/02_phase2a/metadata"
INTERIM = ROOT / "data/interim/phase2a"

MARKERS = {
    "EVT": "HLA-G;MMP2;ITGA5",
    "VCT": "EGFR;TP63;KRT7",
    "SCT": "CGA;CGB3;ERVW-1",
    "ENDOTHELIAL": "PECAM1;VWF;KDR",
    "PERICYTE": "RGS5;CSPG4;PDGFRB",
    "VASCULAR_SMOOTH_MUSCLE": "ACTA2;MYH11;TAGLN",
    "PLACENTAL_STROMAL": "COL1A1;COL3A1;DCN",
    "HOFBAUER": "C1QC;CD163;FOLR2",
    "MACROPHAGE": "CD68;LST1;FCGR3A",
    "MONOCYTE": "S100A8;S100A9;FCN1",
    "NK": "NKG7;GNLY;KLRD1",
    "T_CELL": "CD3D;CD3E;TRAC",
    "B_CELL": "CD79A;MS4A1;CD37",
    "NEUTROPHIL": "FCGR3B;CSF3R;S100A8",
    "IMMUNE_PROGENITOR_OR_PROLIFERATING": "CD34;MKI67;STMN1",
}


def harmonize(label: str) -> str:
    if label.startswith("TB_EVT"):
        return "EVT"
    if label.startswith("TB_VCT"):
        return "VCT"
    if label.startswith("TB_SCT"):
        return "SCT"
    if label.startswith("VASCULAR_CAPILLARY") or label.startswith("VASCULAR_ARTERIAL") or label.startswith("VASCULAR_VENOUS") or label.startswith("VASCULAR_EC"):
        return "ENDOTHELIAL"
    if label.startswith("VASCULAR_PERICYTE"):
        return "PERICYTE"
    if label == "VASCULAR_SM":
        return "VASCULAR_SMOOTH_MUSCLE"
    if label.startswith("STROMAL_"):
        return "PLACENTAL_STROMAL"
    if label.startswith("IMMUNE_HOF"):
        return "HOFBAUER"
    if label.startswith("IMMUNE_MACROPHAGE"):
        return "MACROPHAGE"
    if label.startswith("IMMUNE_MONOCYTE"):
        return "MONOCYTE"
    if label.startswith("IMMUNE_NK"):
        return "NK"
    if label.startswith("IMMUNE_T-"):
        return "T_CELL"
    if label == "IMMUNE_B":
        return "B_CELL"
    if label == "IMMUNE_NEURTOPHIL":
        return "NEUTROPHIL"
    if label in {"IMMUNE_HSC", "IMMUNE_proliferating"}:
        return "IMMUNE_PROGENITOR_OR_PROLIFERATING"
    raise ValueError(f"Unmapped published annotation: {label}")


def read_metadata() -> tuple[dict[str, list[str]], str, int, int]:
    selected = {"cellID", "celltype", "sample", "donorID", "total_molecules", "early_control", "late_control", "early_PE", "late_PE", "female_fetus", "IUGR", "C-section_birth", "vaginal_birth", "induction", "non-induction", "magnesium", "spinal_anaesthesia", "epidural_anaesthesia", "general_anaesthesia", "delivery_week", "weight", "wieght_percentile-Dolberg", "donor_age"}
    rows: dict[str, list[str]] = {}
    with zipfile.ZipFile(ZIP) as archive:
        member = archive.infolist()[0]
        with archive.open(member) as handle:
            while True:
                raw = handle.readline()
                if not raw:
                    break
                fields = raw.decode("utf-8").rstrip("\r\n").split("\t")
                label = fields[0]
                if label in selected:
                    rows[label] = fields[1:]
                elif len(rows) >= len(selected):
                    break
        return rows, member.filename, member.file_size, member.compress_size


def singleton(values: list[str], indices: list[int], field: str, donor: str) -> str:
    found = sorted({values[i] for i in indices})
    if len(found) != 1:
        raise AssertionError(f"{donor} has inconsistent {field}: {found}")
    return found[0]


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    INTERIM.mkdir(parents=True, exist_ok=True)
    md, member, uncompressed, compressed = read_metadata()
    n = len(md["cellID"])
    assert n == 86752 and all(len(v) == n for v in md.values())
    flags = ["early_control", "late_control", "early_PE", "late_PE"]
    for i in range(n):
        assert sum(int(md[x][i]) for x in flags) == 1
    donor_indices: dict[str, list[int]] = defaultdict(list)
    for i, donor in enumerate(md["donorID"]):
        donor_indices[donor].append(i)
    assert len(donor_indices) == 26
    patients: list[dict[str, object]] = []
    flag_to_group = {"early_control": ("CONTROL", "EARLY_CONTROL", "EARLY"), "late_control": ("CONTROL", "LATE_CONTROL", "LATE"), "early_PE": ("PE", "EOPE", "EARLY"), "late_PE": ("PE", "LOPE", "LATE")}
    for donor in sorted(donor_indices):
        ix = donor_indices[donor]
        active = [f for f in flags if singleton(md[f], ix, f, donor) == "1"]
        assert len(active) == 1
        disease, subtype, ga_group = flag_to_group[active[0]]
        samples = sorted({md["sample"][i] for i in ix})
        patients.append({
            "dataset": "Admati_2023_FIGSHARE", "patient_id": donor, "pregnancy_id": donor,
            "disease_status": disease, "pe_subtype_or_control_group": subtype,
            "gestational_age_group": ga_group, "delivery_gestational_age_weeks": float(singleton(md["delivery_week"], ix, "delivery_week", donor)),
            "tissue_compartment": "PLACENTAL_COTYLEDON_VILLI", "cell_count": len(ix), "library_count": len(samples),
            "library_ids": ";".join(samples), "female_fetus": singleton(md["female_fetus"], ix, "female_fetus", donor),
            "iugr": singleton(md["IUGR"], ix, "IUGR", donor), "maternal_age_years": singleton(md["donor_age"], ix, "donor_age", donor),
            "delivery_mode": "C_SECTION" if singleton(md["C-section_birth"], ix, "C-section_birth", donor) == "1" else "VAGINAL",
            "induction": singleton(md["induction"], ix, "induction", donor),
            "raw_counts_available": "YES_PUBLISHED_UMI_TABLE", "published_annotations_available": "YES",
            "donor_identity_validated": "YES", "include_phase2a": "YES",
            "source_url": "https://doi.org/10.6084/m9.figshare.23264102.v1|https://doi.org/10.1016/j.medj.2023.07.005",
            "source_accession": "Figshare:23264102.v1:file41003240", "audit_date": "2026-08-09"
        })
    write_csv(OUT / "patient_registry.csv", patients, list(patients[0]))
    counts = Counter(md["celltype"])
    donors_by_type: dict[str, set[str]] = defaultdict(set)
    for donor, ct in zip(md["donorID"], md["celltype"]):
        donors_by_type[ct].add(donor)
    annotations = []
    for original in sorted(counts):
        h = harmonize(original)
        annotations.append({
            "original_annotation": original, "harmonized_annotation": h,
            "published_cell_count": counts[original], "donors_represented": len(donors_by_type[original]),
            "marker_evidence": MARKERS[h], "marker_validation_status": "PENDING_COUNT_LEVEL_CHECK",
            "annotation_confidence": "PUBLISHED_HIGH; HARMONIZATION_RULE_FROZEN",
            "mapping_rationale": "Prefix/taxonomy collapse within the published lineage; no disease labels used",
            "source_url": "https://doi.org/10.1016/j.medj.2023.07.005|https://doi.org/10.6084/m9.figshare.23264102.v1",
            "audit_date": "2026-08-09"
        })
    write_csv(OUT / "celltype_annotation_registry.csv", annotations, list(annotations[0]))
    checksum = hashlib.md5(ZIP.read_bytes()).hexdigest()
    freeze = {"zip_path": str(ZIP), "zip_md5": checksum, "member": member, "uncompressed_bytes": uncompressed, "compressed_bytes": compressed, "cell_columns": n, "donor_n": len(donor_indices), "library_n": len(set(md["sample"])), "original_annotation_n": len(counts), "harmonized_annotation_n": len({harmonize(x) for x in counts}), "metadata_rows": sorted(md), "counts_begin_after_row": 23}
    (INTERIM / "admati_metadata_freeze.json").write_text(json.dumps(freeze, indent=2), encoding="utf-8")
    print(json.dumps(freeze, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
