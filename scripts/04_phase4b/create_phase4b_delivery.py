#!/usr/bin/env python3
"""Create the lightweight Phase 4B delivery archive after the phase commit."""

from __future__ import annotations

import subprocess
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ZIP = ROOT / "Phase4B_Delivery.zip"


def main() -> int:
    requested = [
        "docs/PHASE4B_EXTERNAL_TRIANGULATION_REPORT.md",
        "results/04_phase4b/topology/candidate_secretion_topology.csv",
        "results/04_phase4b/protein/hucmsc_protein_source_evidence.csv",
        "results/04_phase4b/perturbation/empirical_signed_perturbation_evidence.csv",
        "results/04_phase4b/disease/pe_candidate_context_evidence.csv",
        "results/04_phase4b/novelty/direct_msc_pe_overlap.csv",
        "results/04_phase4b/integration/phase4b_candidate_evidence_matrix.csv",
        "results/04_phase4b/integration/phase4b_candidate_classification.csv",
        "results/04_phase4b/integration/mixed_direction_stress_test.csv",
        "results/04_phase4b/qc/phase4b_risk_flags.csv",
        "results/04_phase4b/novelty/literature_search_log.csv",
    ]
    files = [(ROOT / rel, Path(rel).name) for rel in requested]
    for figure in sorted((ROOT / "results/04_phase4b/figures").glob("*.png")):
        files.append((figure, f"figure_previews/{figure.name}"))

    git_exe = "git"
    log = subprocess.run([git_exe, "-c", f"safe.directory={ROOT}", "log", "--oneline"], cwd=ROOT, check=True, capture_output=True, text=True).stdout
    log_path = ROOT / "data/interim/phase4b_delivery/git_log_oneline.txt"
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
        assert all(not name.lower().endswith((".mtx", ".h5", ".fastq", ".fastq.gz", ".rds")) for name in names)
        assert all("raw" not in name.lower() for name in names)
    print(f"PHASE4B_DELIVERY_OK: {ZIP.name}; files={len(files)}; bytes={ZIP.stat().st_size}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
