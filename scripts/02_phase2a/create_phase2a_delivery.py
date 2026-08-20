# DEPRECATED / HISTORICAL ONLY: Phase 2A count-likelihood receiver inference is superseded by Phase 2A.2 pregnancy-level continuous-expression analysis. Excluded from the default execution workflow.
#!/usr/bin/env python3
"""Create the lightweight Phase 2A delivery ZIP without full matrices."""

from __future__ import annotations

import shutil
import subprocess
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STAGE = ROOT / "data/interim/phase2a_delivery"
ZIP = ROOT / "Phase2A_Delivery.zip"

FILES = {
    ROOT / "docs/PHASE2A_PE_CELLSTATE_DISEASE_PROGRAM_REPORT.md": "PHASE2A_PE_CELLSTATE_DISEASE_PROGRAM_REPORT.md",
    ROOT / "results/02_phase2a/metadata/patient_registry.csv": "metadata/patient_registry.csv",
    ROOT / "results/02_phase2a/metadata/celltype_annotation_registry.csv": "metadata/celltype_annotation_registry.csv",
    ROOT / "results/02_phase2a/metadata/pseudobulk_eligibility.csv": "metadata/pseudobulk_eligibility.csv",
    ROOT / "results/02_phase2a/pseudobulk/pseudobulk_registry.csv": "pseudobulk/pseudobulk_registry.csv",
    ROOT / "results/02_phase2a/DE/celltype_DE_summary.csv": "DE/celltype_DE_summary.csv",
    ROOT / "results/02_phase2a/programs/pe_cellstate_programs.csv": "programs/pe_cellstate_programs.csv",
    ROOT / "results/02_phase2a/programs/shared_pe_programs.csv": "programs/shared_pe_programs.csv",
    ROOT / "results/02_phase2a/programs/eope_specific_programs.csv": "programs/eope_specific_programs.csv",
    ROOT / "results/02_phase2a/programs/lope_specific_programs.csv": "programs/lope_specific_programs.csv",
    ROOT / "results/02_phase2a/regulons/cellstate_regulon_activity.csv": "regulons/cellstate_regulon_activity.csv",
    ROOT / "results/02_phase2a/qc/phase2a_qc_summary.csv": "qc/phase2a_qc_summary.csv",
    ROOT / "results/02_phase2a/qc/phase2a_risk_flags.csv": "qc/phase2a_risk_flags.csv",
}


def main() -> int:
    if STAGE.exists(): shutil.rmtree(STAGE)
    STAGE.mkdir(parents=True)
    for src, rel in FILES.items():
        dst = STAGE / rel; dst.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(src,dst)
    for fig in sorted((ROOT / "results/02_phase2a/figures").glob("*.png")):
        dst=STAGE/"figure_previews"/fig.name; dst.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(fig,dst)
    log = subprocess.check_output(["git","log","--oneline"],cwd=ROOT,text=True,encoding="utf-8")
    (STAGE/"git_log_oneline.txt").write_text(log,encoding="utf-8")
    with zipfile.ZipFile(ZIP,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for path in sorted(STAGE.rglob("*")):
            if path.is_file(): archive.write(path,path.relative_to(STAGE).as_posix())
    print(f"Created {ZIP} ({ZIP.stat().st_size} bytes)")
    return 0


if __name__ == "__main__": raise SystemExit(main())
