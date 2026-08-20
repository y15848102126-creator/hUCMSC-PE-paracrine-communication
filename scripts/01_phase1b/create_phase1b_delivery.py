#!/usr/bin/env python3
"""Create the lightweight Phase 1B delivery archive after the Phase 1B commit."""

from __future__ import annotations

import shutil
import subprocess
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STAGE = ROOT / "data" / "interim" / "phase1b_delivery"
ZIP = ROOT / "Phase1B_Delivery.zip"


def main() -> int:
    if STAGE.exists():
        shutil.rmtree(STAGE)
    STAGE.mkdir(parents=True)
    files = {
        ROOT / "docs" / "PHASE1B_PE_DISEASE_SIGNATURE_REPORT.md": STAGE / "PHASE1B_PE_DISEASE_SIGNATURE_REPORT.md",
        ROOT / "results" / "01_phase1b" / "meta" / "stable_pe_genes.csv": STAGE / "stable_pe_genes.csv",
        ROOT / "results" / "01_phase1b" / "meta" / "pe_gene_meta_analysis.csv": STAGE / "pe_gene_meta_analysis.csv",
        ROOT / "results" / "01_phase1b" / "robustness" / "leave_one_cohort_out.csv": STAGE / "leave_one_cohort_out.csv",
        ROOT / "results" / "01_phase1b" / "robustness" / "standardized_effect_sensitivity.csv": STAGE / "standardized_effect_sensitivity.csv",
        ROOT / "results" / "01_phase1b" / "robustness" / "covariate_sensitivity.csv": STAGE / "covariate_sensitivity.csv",
        ROOT / "results" / "01_phase1b" / "qc" / "phase1b_model_diagnostics.csv": STAGE / "phase1b_model_diagnostics.csv",
        ROOT / "results" / "01_phase1b" / "qc" / "phase1b_risk_flags.csv": STAGE / "phase1b_risk_flags.csv",
    }
    for source, target in files.items():
        assert source.exists(), source
        shutil.copy2(source, target)
    figure_stage = STAGE / "figure_previews"
    figure_stage.mkdir()
    for source in sorted((ROOT / "results" / "01_phase1b" / "figures").glob("*.png")):
        shutil.copy2(source, figure_stage / source.name)
    log = subprocess.run(
        ["git", "-c", f"safe.directory={ROOT.as_posix()}", "log", "--oneline"],
        cwd=ROOT, check=True, capture_output=True, text=True, encoding="utf-8", errors="replace",
    ).stdout
    (STAGE / "git_log_oneline.txt").write_text(log, encoding="utf-8", newline="\n")
    if ZIP.exists():
        ZIP.unlink()
    with zipfile.ZipFile(ZIP, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(STAGE.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(STAGE).as_posix())
    with zipfile.ZipFile(ZIP) as archive:
        names = archive.namelist()
        assert not any("raw" in x.lower() or "matrix" in x.lower() or "fastq" in x.lower() for x in names)
        assert "git_log_oneline.txt" in names and len([x for x in names if x.endswith(".png")]) >= 6
    print(f"Created {ZIP} ({ZIP.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
