#!/usr/bin/env python3
"""Create the lightweight Phase2B_Delivery.zip after the Phase 2B commit."""

from __future__ import annotations

import subprocess
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "Phase2B_Delivery.zip"
FILES = [
    "docs/PHASE2B_BULK_PROGRAM_VALIDATION_REPORT.md",
    "results/02_phase2b/hypotheses/frozen_phase2b_modules.csv",
    "results/02_phase2b/meta/program_gene_set_meta_analysis.csv",
    "results/02_phase2b/meta/program_module_validation.csv",
    "results/02_phase2b/robustness/program_leave_one_cohort_out.csv",
    "results/02_phase2b/robustness/cameraPR_method_concordance.csv",
    "results/02_phase2b/evidence/updated_receiver_evidence_hierarchy.csv",
    "results/02_phase2b/evidence/phase2b_risk_flags.csv",
    "results/02_phase2b/figures/phase2b_module_by_cohort_direction_heatmap.png",
    "results/02_phase2b/figures/phase2b_representative_constituent_forest.png",
    "results/02_phase2b/figures/phase2b_leave_one_cohort_out_robustness.png",
    "results/02_phase2b/figures/phase2b_scrna_bulk_evidence_matrix.png",
    "results/02_phase2a2/external_validation/yang_lope_updated_evidence.csv",
]


def main() -> int:
    for rel in FILES:
        if not (ROOT / rel).is_file():
            raise FileNotFoundError(rel)
    log = subprocess.check_output(
        ["git", "-c", "safe.directory=*", "log", "--oneline"],
        cwd=ROOT,
        text=True,
        encoding="utf-8",
    )
    with zipfile.ZipFile(DEST, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for rel in FILES:
            archive.write(ROOT / rel, rel)
        archive.writestr("git_log_oneline.txt", log)
    with zipfile.ZipFile(DEST) as archive:
        names = archive.namelist()
        assert len(names) == len(FILES) + 1
        assert not any(n.startswith((".git/", "data/", "environment/", "cache/")) for n in names)
        assert not any(n.lower().endswith((".fastq", ".fastq.gz", ".mtx", ".rds", ".h5", ".h5ad")) for n in names)
        assert "results/02_phase2a2/external_validation/yang_lope_updated_evidence.csv" in names
        assert "git_log_oneline.txt" in names
    print(f"Created {DEST.name}: {DEST.stat().st_size} bytes, {len(FILES) + 1} entries")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
