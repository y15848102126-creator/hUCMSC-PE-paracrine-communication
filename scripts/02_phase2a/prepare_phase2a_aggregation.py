#!/usr/bin/env python3
"""Prepare an outcome-blind cell-to-pseudobulk grouping map."""

from __future__ import annotations

import csv
import struct
from pathlib import Path

import pandas as pd

from freeze_admati_metadata import harmonize, read_metadata

ROOT = Path(__file__).resolve().parents[2]
INTERIM = ROOT / "data/interim/phase2a"


def main() -> int:
    md, _, _, _ = read_metadata()
    patients = pd.read_csv(ROOT / "results/02_phase2a/metadata/patient_registry.csv", dtype=str).set_index("patient_id")
    donors = sorted(patients.index)
    celltypes = sorted({harmonize(x) for x in md["celltype"]})
    strata = [(d, c) for c in celltypes for d in donors]
    index = {x: i for i, x in enumerate(strata)}
    group_ids = [index[(d, harmonize(c))] for d, c in zip(md["donorID"], md["celltype"])]
    INTERIM.mkdir(parents=True, exist_ok=True)
    with (INTERIM / "cell_group_ids.int32").open("wb") as handle:
        handle.write(struct.pack(f"<{len(group_ids)}i", *group_ids))
    with (INTERIM / "pseudobulk_strata.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["group_id", "patient_id", "harmonized_annotation", "matrix_column", "expected_cell_count", "expected_total_umi"])
        for i, (donor, celltype) in enumerate(strata):
            indices = [j for j, gid in enumerate(group_ids) if gid == i]
            writer.writerow([i, donor, celltype, f"{donor}__{celltype}", len(indices), sum(int(md["total_molecules"][j]) for j in indices)])
    print(f"Prepared {len(group_ids)} cells -> {len(strata)} frozen strata")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
