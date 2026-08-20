#!/usr/bin/env python3
"""Download the six frozen Yang placenta matrices and auditable public scripts."""

from __future__ import annotations

import hashlib
import gzip
import json
import shutil
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/raw/phase2a1"
YANG = OUT / "yang_2023_placenta"
ADM = OUT / "admati_2023"
UNCOMPRESSED = ROOT / "data/interim/phase2a1/yang_uncompressed"
CONTENTS = ROOT / "data/raw/phase0b/yang2023_github_contents.json"

YANG_FILES = {
    "pla1_placenta0423-2_exon_tagged.dge.txt.gz",
    "pla2_placenta-2018_combined_exon_tagged.dge.txt.gz",
    "pla3_Placenta20190402_L3_1000902_exon_tagged.dge.txt.gz",
    "pla4_placenta-2019_combined_exon_tagged.dge.txt.gz",
    "pla5_placenta508-1_exon_tagged.dge.txt.gz",
    "pla6_placenta508-2_exon_tagged.dge.txt.gz",
}
SCRIPT_URLS = {
    "yang_hfzV16_step1.R": "https://raw.githubusercontent.com/JustMoveOnnn/preeclampsia/main/single_cell_matrix/code/hfzV16_step1.R",
    "yang_step4_cluster.R": "https://raw.githubusercontent.com/JustMoveOnnn/preeclampsia/main/single_cell_matrix/code/step4_cluster.R",
    "yang_step5_anno.R": "https://raw.githubusercontent.com/JustMoveOnnn/preeclampsia/main/single_cell_matrix/code/step5_anno.R",
    "admati_load_sc_PE_data_and_save_v1.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/load_sc_PE_data_and_save_v1.m",
    "admati_sc_step1_split_to_classes_v1.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_step1_split_to_classes_v1.m",
    "admati_sc_step2_immune_v1.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_step2_immune_v1.m",
    "admati_sc_step2_stromal_v1.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_step2_stromal_v1.m",
    "admati_sc_step2_trophoblasts_v1.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_step2_trophoblasts_v1.m",
    "admati_sc_step2_vascular_v1.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_step2_vascular_v1.m",
    "admati_sc_step3_immune_v1.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_step3_immune_v1.m",
    "admati_sc_step3_stromal_v1.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_step3_stromal_v1.m",
    "admati_sc_step3_trophoblasts_v1.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_step3_trophoblasts_v1.m",
    "admati_sc_step3_vascular_v1.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_step3_vascular_v1.m",
    "admati_sc_step4_immune_v3.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_step4_immune_v3.m",
    "admati_sc_step4_stromal_v3.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_step4_stromal_v3.m",
    "admati_sc_step4_trophoblasts_v3.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_step4_trophoblasts_v3.m",
    "admati_sc_step4_vascular_v2.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_step4_vascular_v2.m",
    "admati_sc_plot_scatter_pece_plct_alltypes_v3.m": "https://raw.githubusercontent.com/zeiselamit/PE_2023/main/sc_plot_scatter_pece_plct_alltypes_v3.m",
}


def fetch(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size:
        return
    req = urllib.request.Request(url, headers={"User-Agent": "hUCMSC-PE-reproducibility-audit"})
    with urllib.request.urlopen(req, timeout=120) as response:
        path.write_bytes(response.read())


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    YANG.mkdir(parents=True, exist_ok=True)
    ADM.mkdir(parents=True, exist_ok=True)
    records = json.loads(CONTENTS.read_text(encoding="utf-8"))
    urls = {record["name"]: record["download_url"] for record in records}
    manifest = []
    for name in sorted(YANG_FILES):
        if name not in urls:
            raise RuntimeError(f"Frozen Yang file absent from audited repository listing: {name}")
        path = YANG / name
        fetch(urls[name], path)
        UNCOMPRESSED.mkdir(parents=True, exist_ok=True)
        unpacked = UNCOMPRESSED / name.removesuffix(".gz")
        if not unpacked.exists() or not unpacked.stat().st_size:
            with gzip.open(path, "rb") as source, unpacked.open("wb") as destination:
                shutil.copyfileobj(source, destination, length=1 << 20)
        manifest.append({"source": "Yang2023_GitHub", "file": name, "bytes": path.stat().st_size, "sha256": digest(path), "url": urls[name]})
    for name, url in SCRIPT_URLS.items():
        path = (ADM if name.startswith("admati_") else YANG) / name
        fetch(url, path)
        manifest.append({"source": "Admati2023_GitHub" if name.startswith("admati_") else "Yang2023_GitHub", "file": name, "bytes": path.stat().st_size, "sha256": digest(path), "url": url})
    (OUT / "source_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"downloaded_or_verified": len(manifest), "bytes": sum(x["bytes"] for x in manifest)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
