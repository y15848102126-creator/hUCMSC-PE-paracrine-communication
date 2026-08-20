#!/usr/bin/env python3
"""Download and structurally extract small OA hUC-MSC CM peptidomics tables."""

from __future__ import annotations

import csv
import hashlib
import urllib.request
from pathlib import Path

from docx import Document


ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "data/raw/phase4b/pmc7510303_supplements"
DEST.mkdir(parents=True, exist_ok=True)
FILES = {
    "13287_2020_1931_MOESM3_ESM.docx": "https://static-content.springer.com/esm/art%3A10.1186%2Fs13287-020-01931-0/MediaObjects/13287_2020_1931_MOESM3_ESM.docx",
    "13287_2020_1931_MOESM4_ESM.docx": "https://static-content.springer.com/esm/art%3A10.1186%2Fs13287-020-01931-0/MediaObjects/13287_2020_1931_MOESM4_ESM.docx",
    "13287_2020_1931_MOESM5_ESM.docx": "https://static-content.springer.com/esm/art%3A10.1186%2Fs13287-020-01931-0/MediaObjects/13287_2020_1931_MOESM5_ESM.docx",
    "13287_2020_1931_MOESM6_ESM.docx": "https://static-content.springer.com/esm/art%3A10.1186%2Fs13287-020-01931-0/MediaObjects/13287_2020_1931_MOESM6_ESM.docx",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    return h


def main() -> int:
    manifest, extracted = [], []
    for filename, url in FILES.items():
        path = DEST / filename
        if not path.exists():
            request = urllib.request.Request(url, headers={"User-Agent": "Phase4B-evidence-audit/1.0"})
            with urllib.request.urlopen(request, timeout=90) as response:
                path.write_bytes(response.read())
        document = Document(path)
        for table_i, table in enumerate(document.tables, 1):
            for row_i, row in enumerate(table.rows, 1):
                values = [cell.text.strip().replace("\n", " | ") for cell in row.cells]
                extracted.append({"source_file": filename, "table_index": table_i, "row_index": row_i, "cell_values": "\t".join(values), "source_url": url})
        manifest.append({"filename": filename, "bytes": path.stat().st_size, "sha256": sha256(path), "url": url, "PMID": "32967723", "DOI": "10.1186/s13287-020-01931-0"})
    with (DEST / "manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0])); writer.writeheader(); writer.writerows(manifest)
    with (DEST / "extracted_table_rows.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(extracted[0])); writer.writeheader(); writer.writerows(extracted)
    print(f"PHASE4B_SUPPLEMENT_OK files={len(manifest)} table_rows={len(extracted)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
