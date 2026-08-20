#!/usr/bin/env python3
"""Create and validate the lightweight Phase 3 review archive."""

from __future__ import annotations

import shutil
import subprocess
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STAGE = ROOT / "data" / "interim" / "phase3_delivery"
ZIP = ROOT / "Phase3_Delivery.zip"

FILES = [
    ("docs/PHASE3_HUCMSC_SENDER_PROGRAM_REPORT.md", "PHASE3_HUCMSC_SENDER_PROGRAM_REPORT.md"),
    ("results/03_phase3/metadata/hucmsc_donor_registry.csv", "hucmsc_donor_registry.csv"),
    ("results/03_phase3/ligand_universe/frozen_ligand_universe.csv", "frozen_ligand_universe.csv"),
    ("results/03_phase3/baseline/baseline_sender_robustness.csv", "baseline_sender_robustness.csv"),
    ("results/03_phase3/licensing/licensing_ligand_classification.csv", "licensing_ligand_classification.csv"),
    ("results/03_phase3/licensing/licensing_programs.csv", "licensing_programs.csv"),
    ("results/03_phase3/sender/sender_evidence_hierarchy.csv", "sender_evidence_hierarchy.csv"),
    ("results/03_phase3/sender/frozen_phase4_sender_candidates.csv", "frozen_phase4_sender_candidates.csv"),
    ("results/03_phase3/qc/phase3_risk_flags.csv", "phase3_risk_flags.csv"),
]


def main() -> int:
    if STAGE.exists():
        shutil.rmtree(STAGE)
    STAGE.mkdir(parents=True)
    for source, name in FILES:
        src = ROOT / source
        assert src.is_file() and src.stat().st_size > 0, source
        shutil.copy2(src, STAGE / name)
    figure_dir = STAGE / "figure_previews"
    figure_dir.mkdir()
    figures = sorted((ROOT / "results" / "03_phase3" / "figures").glob("*.png"))
    assert len(figures) >= 5
    for src in figures:
        shutil.copy2(src, figure_dir / src.name)
    log = subprocess.run(["git", "log", "--oneline"], cwd=ROOT, check=True, capture_output=True, text=True).stdout
    (STAGE / "git_log_oneline.txt").write_text(log, encoding="utf-8")
    if ZIP.exists():
        ZIP.unlink()
    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(STAGE.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(STAGE).as_posix())
    with zipfile.ZipFile(ZIP) as archive:
        names = archive.namelist()
        assert len(names) == len(FILES) + len(figures) + 1
        forbidden = (".git/", "fastq", "matrix.mtx", "data/raw", "environment/", "cache")
        assert not any(any(token in name.lower() for token in forbidden) for name in names)
        assert all(info.file_size > 0 for info in archive.infolist())
    print(f"PHASE3_DELIVERY_OK: {ZIP.name}; {ZIP.stat().st_size} bytes; {len(names)} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
