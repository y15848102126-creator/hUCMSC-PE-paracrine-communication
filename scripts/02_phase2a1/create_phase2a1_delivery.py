#!/usr/bin/env python3
"""Build the lightweight, raw-matrix-free Phase 2A.1 delivery archive."""

from __future__ import annotations

import subprocess
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "Phase2A1_Delivery.zip"
FILES = [
    "docs/PHASE2A1_RECEIVER_PROGRAM_STRESS_TEST_REPORT.md",
    "results/02_phase2a1/clinical_confounding_audit.csv",
    "results/02_phase2a1/shared_program_iugr_sensitivity.csv",
    "results/02_phase2a1/shared_program_delivery_sensitivity.csv",
    "results/02_phase2a1/patient_program_score_sensitivity.csv",
    "results/02_phase2a1/yang_lope_replication.csv",
    "results/02_phase2a1/shared_program_redundancy_matrix.csv",
    "results/02_phase2a1/shared_program_modules.csv",
    "results/02_phase2a1/cameraPR_statistic_audit.csv",
    "results/02_phase2a1/phase2a1_risk_flags.csv",
    "results/02_phase2a1/phase2a1_session_info.txt",
    "results/02_phase2a1/figures/A_clinical_overlap.png",
    "results/02_phase2a1/figures/B_shared_program_stress_test.png",
    "results/02_phase2a1/figures/C_yang_replication.png",
    "results/02_phase2a1/figures/D_program_module_overlap.png",
]


def main() -> int:
    for relative in FILES:
        path = ROOT / relative
        if not path.exists():
            raise FileNotFoundError(path)
    log = subprocess.check_output(["git", "log", "--oneline"], cwd=ROOT, text=True, encoding="utf-8")
    with zipfile.ZipFile(ARCHIVE, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for relative in FILES:
            archive.write(ROOT / relative, relative)
        archive.writestr("git_log_oneline.txt", log)
    with zipfile.ZipFile(ARCHIVE) as archive:
        names = set(archive.namelist())
        expected = set(FILES) | {"git_log_oneline.txt"}
        assert names == expected
        forbidden = ("data/raw/", ".git/", ".fastq", ".gz", ".rds", ".ipynb")
        assert not any(any(token in name.lower() for token in forbidden) for name in names)
        assert archive.testzip() is None
    print(f"Created {ARCHIVE.name}: {ARCHIVE.stat().st_size} bytes, {len(FILES) + 1} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
