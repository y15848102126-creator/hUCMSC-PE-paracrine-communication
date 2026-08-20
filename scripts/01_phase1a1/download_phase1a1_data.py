#!/usr/bin/env python3
"""Download the small public raw/non-normalized archives used in Phase 1A.1."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import time
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "phase1a1"
CFG = json.loads((ROOT / "config" / "phase1a1_sources.json").read_text(encoding="utf-8"))
USER_AGENT = "hUCMSC-PE-Phase1A1-preprocessing-audit/1.0"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(url: str, path: Path, force: bool) -> None:
    if path.exists() and path.stat().st_size and not force:
        return
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    error = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as handle:
                shutil.copyfileobj(response, handle)
            return
        except Exception as exc:
            error = exc
            if path.exists():
                path.unlink()
            time.sleep(2**attempt)
    raise RuntimeError(f"Download failed: {url}: {error}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    manifest = []
    for item in CFG["files"]:
        path = RAW / item["name"]
        fetch(item["url"], path, args.force)
        if path.stat().st_size != item["expected_bytes"]:
            raise AssertionError(f"Size mismatch for {path}: {path.stat().st_size}")
        observed_hash = sha256(path)
        if observed_hash != item["expected_sha256"]:
            raise AssertionError(f"SHA-256 mismatch for {path}: {observed_hash}")
        manifest.append({
            "kind": item["kind"], "path": path.relative_to(ROOT).as_posix(),
            "bytes": path.stat().st_size, "sha256": observed_hash, "source_url": item["url"],
        })
    (RAW / "download_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Phase 1A.1 raw cache ready: {RAW}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
