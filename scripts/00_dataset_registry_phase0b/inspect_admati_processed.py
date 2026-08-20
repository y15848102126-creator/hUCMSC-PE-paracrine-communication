#!/usr/bin/env python3
"""Inspect only the metadata rows of the public Admati matrix ZIP."""

from __future__ import annotations

import argparse
import json
import zipfile
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ZIP = ROOT / "data" / "raw" / "phase0b" / "sc_PE_allcells_with_metadata_29-May-2023.txt.zip"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    args = parser.parse_args()
    selected = {
        "sample", "donorID", "early_control", "late_control", "early_PE", "late_PE",
        "female_fetus", "IUGR", "C-section_birth", "vaginal_birth", "induction",
        "non-induction", "magnesium", "spinal_anaesthesia", "epidural_anaesthesia",
        "general_anaesthesia", "delivery_week", "weight", "wieght_percentile-Dolberg",
        "donor_age",
    }
    result: dict[str, object] = {"zip": str(args.zip), "metadata": {}}
    with zipfile.ZipFile(args.zip) as archive:
        member = archive.infolist()[0]
        result["member"] = member.filename
        result["uncompressed_bytes"] = member.file_size
        result["compressed_bytes"] = member.compress_size
        with archive.open(member) as handle:
            for line_no in range(80):
                raw = handle.readline()
                if not raw:
                    break
                fields = raw.rstrip(b"\r\n").split(b"\t")
                label = fields[0].decode("utf-8", "replace")
                if line_no == 0:
                    result["cell_columns"] = len(fields) - 1
                    continue
                if label in selected:
                    values = [value.decode("utf-8", "replace") for value in fields[1:]]
                    result["metadata"][label] = dict(sorted(Counter(values).items()))
                if label.startswith("ENSG") or (line_no > 25 and label not in selected):
                    result["first_expression_row"] = label
                    break
    print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
