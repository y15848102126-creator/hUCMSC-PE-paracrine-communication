#!/usr/bin/env python3
"""Freeze Phase 2B module hypotheses and exact MSigDB membership before outcomes."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/02_phase2b/hypotheses"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader(); writer.writerows(rows)


def gmt(path: Path, collection: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        out[f"{collection}::{fields[0]}"] = list(dict.fromkeys(fields[2:]))
    return out


def main() -> int:
    cfg = json.loads((ROOT / "config/phase2b_analysis.json").read_text(encoding="utf-8"))
    for rel, expected in cfg["receiver_history"]["upstream_sha256"].items():
        assert sha(ROOT / rel) == expected, f"Frozen upstream changed: {rel}"
    resource = ROOT / "data/raw/phase2a_resources"
    mapping = {}
    mapping.update(gmt(resource / "h.all.v2026.1.Hs.symbols.gmt", "HALLMARK"))
    mapping.update(gmt(resource / "c2.cp.reactome.v2026.1.Hs.symbols.gmt", "REACTOME"))
    mapping.update(gmt(resource / "c5.go.bp.v2026.1.Hs.symbols.gmt", "GOBP"))
    modules_all = read_csv(ROOT / "results/02_phase2a2/corrected_analysis/corrected_program_modules.csv")
    members = [r for r in modules_all if r["record_type"] == "ORIGINAL_GENE_SET" and r["program_module"] in cfg["candidate_modules"]]
    assert len(members) == cfg["expected_constituent_hypothesis_n"] == 19
    hypothesis_rows = []
    for i, row in enumerate(sorted(members, key=lambda r: (r["program_module"], r["gene_set"])), 1):
        genes = mapping[row["pathway"]]
        membership = ";".join(genes)
        hypothesis_rows.append({
            "hypothesis_id": f"P2B_H{i:02d}", "program_module": row["program_module"], "module_label": row["module_label"],
            "celltype_origin": row["celltype"], "scRNA_direction": row["frozen_direction"], "collection": row["collection"],
            "gene_set": row["gene_set"], "pathway": row["pathway"], "original_gene_n": len(genes),
            "original_gene_membership": membership, "membership_sha256": hashlib.sha256(membership.encode()).hexdigest(),
            "MSigDB_version": cfg["gene_sets"]["version"], "validation_status": "FROZEN_BEFORE_BULK_OUTCOME",
            "source_url": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/|results/02_phase2a2/corrected_analysis/corrected_program_modules.csv"
        })
    write_csv(OUT / "frozen_phase2b_gene_sets.csv", hypothesis_rows)

    hierarchy = {r["program_module"]: r for r in read_csv(ROOT / "results/02_phase2a2/evidence/receiver_module_evidence_hierarchy.csv")}
    module_rows = []
    for module in cfg["candidate_modules"]:
        h = hierarchy[module]
        subset = [r for r in hypothesis_rows if r["program_module"] == module]
        n = len(subset)
        module_rows.append({
            "program_module": module, "module_label": h["module_label"], "celltype_origin": h["celltype"], "scRNA_direction": h["frozen_direction"],
            "constituent_gene_set_n": n, "required_robust_or_directional_n": int(n * 0.60 + 0.999999),
            "hypothesis_ids": ";".join(r["hypothesis_id"] for r in subset), "constituent_gene_sets": ";".join(r["gene_set"] for r in subset),
            "corrected_admati_support": h["corrected_module_status"], "external_scRNA_evidence_level": h["evidence_level"],
            "phase2a2_candidate_status": h["phase2b_program_validation_candidate"], "validation_status": "FROZEN_BEFORE_BULK_OUTCOME",
            "bulk_celltype_claim_allowed": "NO_TISSUE_LEVEL_PROGRAM_ONLY", "source_url": "results/02_phase2a2/evidence/receiver_module_evidence_hierarchy.csv|config/phase2b_analysis.json"
        })
    write_csv(OUT / "frozen_phase2b_modules.csv", module_rows)
    print(json.dumps({"modules": len(module_rows), "constituent_hypotheses": len(hypothesis_rows), "unique_pathways": len({r['pathway'] for r in hypothesis_rows})}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
