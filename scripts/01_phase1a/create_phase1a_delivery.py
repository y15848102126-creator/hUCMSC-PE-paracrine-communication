#!/usr/bin/env python3
"""Create a lightweight review ZIP containing reports and CSV registries only."""

from __future__ import annotations

import subprocess
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT / "Phase1A_Bulk_Data_Freeze_Delivery.zip"
DOCS = [
    "docs/DATASET_AUDIT_REPORT.md",
    "docs/DATASET_AUDIT_PHASE0B_REPORT.md",
    "docs/PHASE0B_CHANGELOG.md",
    "docs/PHASE1A_BULK_DATA_FREEZE_REPORT.md",
    "docs/PHASE1B_STATISTICAL_ANALYSIS_PLAN.md",
    "docs/PHASE1A_CHANGELOG.md",
]
PHASE0B_REQUIRED = {
    "pe_scrna_extended_registry.csv", "pe_scrna_data_access.csv", "phs001886_version_overlap.csv",
    "hucmsc_sender_redundancy_registry.csv", "phase0b_risk_flags.csv", "revised_dataset_roles.csv",
}


def main() -> int:
    files = [ROOT / "README.md", *(ROOT / name for name in DOCS)]
    for directory in [
        ROOT / "results" / "00_dataset_audit",
        ROOT / "results" / "00_dataset_audit_phase0b",
        ROOT / "results" / "01_phase1a",
    ]:
        files.extend(sorted(directory.glob("*.csv")))
    missing = [path for path in files if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing delivery inputs: {missing}")
    phase0b_names = {path.name for path in files if path.parent.name == "00_dataset_audit_phase0b"}
    if not PHASE0B_REQUIRED <= phase0b_names:
        raise AssertionError(f"Missing Phase 0B delivery fix: {PHASE0B_REQUIRED - phase0b_names}")

    git_log = subprocess.run(
        ["git", "-c", f"safe.directory={ROOT.as_posix()}", "log", "--oneline"],
        cwd=ROOT, check=True, capture_output=True, text=True, encoding="utf-8",
    ).stdout
    with zipfile.ZipFile(TARGET, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in files:
            archive.write(path, path.relative_to(ROOT).as_posix())
        archive.writestr("git log --oneline.txt", git_log)

    with zipfile.ZipFile(TARGET) as archive:
        names = archive.namelist()
        forbidden = [name for name in names if any(token in name.lower() for token in [".git/", "data/raw/", "data/interim/", "fastq", ".ipynb", "environment/", "cache/"])]
        if forbidden:
            raise AssertionError(f"Forbidden delivery members: {forbidden}")
        if len(names) != len(set(names)):
            raise AssertionError("Duplicate delivery members")
        for required in PHASE0B_REQUIRED:
            expected = f"results/00_dataset_audit_phase0b/{required}"
            if expected not in names:
                raise AssertionError(f"Missing {expected}")
    print(f"Created {TARGET} ({TARGET.stat().st_size} bytes, {len(names)} files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
