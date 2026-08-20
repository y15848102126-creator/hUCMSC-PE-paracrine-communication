#!/usr/bin/env python3
"""Create the lightweight Phase2A2_Delivery.zip after the Phase 2A.2 commit."""

from __future__ import annotations

import subprocess
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "Phase2A2_Delivery.zip"
FILES = [
    "docs/PHASE2A2_RECEIVER_FRAMEWORK_CORRECTION_REPORT.md",
    "docs/PHASE2A2_CORRECTED_RECEIVER_ANALYSIS_PLAN.md",
    "results/02_phase2a2/provenance/admati_expression_layer_audit.csv",
    "results/02_phase2a2/provenance/admati_raw_count_search.csv",
    "results/02_phase2a2/corrected_analysis/frozen20_corrected_retest.csv",
    "results/02_phase2a2/corrected_analysis/corrected_program_modules.csv",
    "results/02_phase2a2/external_validation/zheng_eope_dataset_audit.csv",
    "results/02_phase2a2/external_validation/zheng_eope_targeted_validation.csv",
    "results/02_phase2a2/evidence/receiver_module_evidence_hierarchy.csv",
    "results/02_phase2a2/evidence/phase2a2_risk_flags.csv",
    "results/02_phase2a2/figures/A_frozen20_corrected_direction.png",
    "results/02_phase2a2/figures/B_tier1_classification.png",
    "results/02_phase2a2/figures/C_corrected_module_status.png",
]


def main() -> int:
    for rel in FILES:
        path = ROOT / rel
        if not path.exists():
            raise FileNotFoundError(rel)
    log = subprocess.check_output(["git", "log", "--oneline"], cwd=ROOT, text=True, encoding="utf-8")
    with zipfile.ZipFile(DEST, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for rel in FILES:
            archive.write(ROOT / rel, rel)
        archive.writestr("git_log_oneline.txt", log)
    with zipfile.ZipFile(DEST) as archive:
        names = archive.namelist()
        assert ".git" not in "|".join(names).lower()
        assert not any("matrix" in n.lower() or "fastq" in n.lower() for n in names)
    print(f"Created {DEST.name}: {DEST.stat().st_size} bytes, {len(FILES)+1} entries")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
