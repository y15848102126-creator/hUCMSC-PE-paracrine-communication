#!/usr/bin/env python3
"""Export untransformed submitted matrices and platform/sample crosswalks for the R audit."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INTERIM = ROOT / "data" / "interim" / "phase1a1"
SUBMITTED = INTERIM / "submitted"


def load_phase1a_builder():
    path = ROOT / "scripts" / "01_phase1a" / "build_phase1a_freeze.py"
    spec = importlib.util.spec_from_file_location("phase1a_builder", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> int:
    builder = load_phase1a_builder()
    SUBMITTED.mkdir(parents=True, exist_ok=True)
    for accession in ["GSE30186", "GSE10588", "GSE43942"]:
        matrix = builder.read_soft_expression(accession)
        matrix.to_csv(SUBMITTED / f"{accession}_submitted_linear.tsv.gz", sep="\t", compression="gzip", lineterminator="\n")
    annotation = builder.read_platform_annotation("GSE10588")
    annotation.to_csv(INTERIM / "GPL2986_annotation.tsv.gz", sep="\t", compression="gzip", index=False, lineterminator="\n")

    with (ROOT / "results" / "01_phase1a" / "bulk_sample_freeze.csv").open(encoding="utf-8-sig", newline="") as handle:
        rows = [row for row in csv.DictReader(handle) if row["dataset"] in {"GSE30186", "GSE10588", "GSE43942"}]
    fields = ["dataset", "GSM/sample ID", "sample_title", "PE/control", "batch", "GA"]
    with (INTERIM / "priority_sample_crosswalk.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("Prepared submitted Phase 1A.1 inputs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
