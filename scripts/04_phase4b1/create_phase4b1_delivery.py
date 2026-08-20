#!/usr/bin/env python3
"""Create the lightweight Phase 4B.1 delivery archive after the phase commit."""

from __future__ import annotations

import subprocess
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ZIP = ROOT / "Phase4B1_Delivery.zip"


def main() -> int:
    requested = [
        "docs/PHASE4B1_PROTEIN_SOURCE_COMPLETENESS_REPORT.md",
        "results/04_phase4b1/hucmsc_proteomics_dataset_registry.csv",
        "results/04_phase4b1/candidate_protein_detection_matrix.csv",
        "results/04_phase4b1/corrected_hucmsc_protein_source_evidence.csv",
        "results/04_phase4b1/protein_identifier_mapping_audit.csv",
        "results/04_phase4b1/corrected_phase4b_candidate_evidence_matrix.csv",
        "results/04_phase4b1/corrected_phase4b_candidate_classification.csv",
        "results/04_phase4b1/phase4b1_risk_flags.csv",
    ]
    files = [(ROOT / rel, Path(rel).name) for rel in requested]
    log = subprocess.run(
        ["git", "-c", f"safe.directory={ROOT}", "log", "--oneline"],
        cwd=ROOT, check=True, capture_output=True, text=True,
    ).stdout
    log_path = ROOT / "data/interim/phase4b1_delivery/git_log_oneline.txt"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(log, encoding="utf-8")
    files.append((log_path, "git_log_oneline.txt"))
    with zipfile.ZipFile(ZIP, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path, name in files:
            if not path.is_file():
                raise FileNotFoundError(path)
            archive.write(path, name)
    with zipfile.ZipFile(ZIP) as archive:
        names = archive.namelist()
        assert len(names) == len(set(names)) == len(files)
        forbidden = (".mtx", ".h5", ".fastq", ".fastq.gz", ".rds", ".raw", ".mztab", ".xlsx", ".xls")
        assert all(not name.lower().endswith(forbidden) for name in names)
        assert all("data/raw" not in name.lower() and ".git" not in name.lower() for name in names)
    print(f"PHASE4B1_DELIVERY_OK: {ZIP.name}; files={len(files)}; bytes={ZIP.stat().st_size}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
