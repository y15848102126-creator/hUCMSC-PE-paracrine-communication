#!/usr/bin/env python3
"""Create the lightweight final-synthesis delivery archive."""

from __future__ import annotations

import argparse
import subprocess
import zipfile
from pathlib import Path


CSV_FILES = [
    "phase_registry.csv", "final_receiver_evidence.csv", "final_sender_evidence.csv",
    "final_17_candidate_evidence.csv", "negative_results_register.csv",
    "claim_ceiling_matrix.csv", "figure_source_registry.csv",
    "supplementary_table_registry.csv", "experimental_validation_roadmap.csv",
    "final_synthesis_risk_flags.csv",
]
DOC_FILES = [
    "FINAL_SYNTHESIS_REPORT.md", "MANUSCRIPT_RESULTS_BLUEPRINT.md",
    "MANUSCRIPT_DISCUSSION_BLUEPRINT.md", "FIGURE_PLAN_FINAL.md",
    "SUPPLEMENTARY_PLAN_FINAL.md",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.root.resolve()
    stage = root / "data" / "interim" / "final_synthesis_delivery"
    stage.mkdir(parents=True, exist_ok=True)
    log = subprocess.run(
        ["git", "log", "--oneline"], cwd=root, check=True, capture_output=True, text=True, encoding="utf-8"
    ).stdout
    log_path = stage / "git_log_oneline.txt"
    log_path.write_text(log, encoding="utf-8")

    members = [(root / "results" / "final_synthesis" / name, f"results/final_synthesis/{name}") for name in CSV_FILES]
    members += [(root / "docs" / name, f"docs/{name}") for name in DOC_FILES]
    members += [(log_path, "git_log_oneline.txt")]
    for source, _ in members:
        if not source.exists():
            raise FileNotFoundError(source)

    archive = root / "Final_Synthesis_Delivery.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for source, arcname in members:
            zf.write(source, arcname)

    with zipfile.ZipFile(archive) as zf:
        names = zf.namelist()
        expected = [arcname for _, arcname in members]
        if names != expected or len(names) != 16:
            raise RuntimeError((names, expected))
        forbidden = (".git/", "data/raw/", ".mtx", ".h5", ".h5ad", ".fastq", ".fq", "__pycache__", "environment/")
        if any(any(token in name.lower() for token in forbidden) for name in names):
            raise RuntimeError("Forbidden delivery member detected")
    print(f"Created {archive} with {len(members)} files")


if __name__ == "__main__":
    main()
