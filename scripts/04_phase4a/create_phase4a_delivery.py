#!/usr/bin/env python3
"""Create the lightweight Phase 4A delivery archive after the phase commit."""

from __future__ import annotations

import subprocess
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ZIP = ROOT / "Phase4A_Delivery.zip"


def main() -> int:
    files = [
        (ROOT / "docs/PHASE4A_SENDER_RECEIVER_INTEGRATION_REPORT.md", "PHASE4A_SENDER_RECEIVER_INTEGRATION_REPORT.md"),
        (ROOT / "results/04_phase4a/freeze/phase4_sender_scopes.csv", "phase4_sender_scopes.csv"),
        (ROOT / "results/04_phase4a/freeze/phase4_receiver_hierarchy.csv", "phase4_receiver_hierarchy.csv"),
        (ROOT / "results/04_phase4a/receptors/receiver_receptor_competence.csv", "receiver_receptor_competence.csv"),
        (ROOT / "results/04_phase4a/lr/sender_receiver_lr_compatibility.csv", "sender_receiver_lr_compatibility.csv"),
        (ROOT / "results/04_phase4a/targets/nichenet_target_compatibility.csv", "nichenet_target_compatibility.csv"),
        (ROOT / "results/04_phase4a/signed/signed_reversal_analysis.csv", "signed_reversal_analysis.csv"),
        (ROOT / "results/04_phase4a/integration/phase4a_candidate_hierarchy.csv", "phase4a_candidate_hierarchy.csv"),
        (ROOT / "results/04_phase4a/integration/disease_concordant_candidates.csv", "disease_concordant_candidates.csv"),
        (ROOT / "results/04_phase4a/qc/phase4a_risk_flags.csv", "phase4a_risk_flags.csv"),
    ]
    for figure in sorted((ROOT / "results/04_phase4a/figures").glob("*.png")):
        files.append((figure, f"figure_previews/{figure.name}"))
    log = subprocess.run(["git", "-c", f"safe.directory={ROOT}", "log", "--oneline"], cwd=ROOT, check=True, capture_output=True, text=True).stdout
    log_path = ROOT / "data/interim/phase4a_delivery/git_log_oneline.txt"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(log, encoding="utf-8")
    files.append((log_path, "git_log_oneline.txt"))
    with zipfile.ZipFile(ZIP, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path, name in files:
            if not path.is_file():
                raise FileNotFoundError(path)
            archive.write(path, name)
    print(f"PHASE4A_DELIVERY_OK: {ZIP.name}; files={len(files)}; bytes={ZIP.stat().st_size}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
