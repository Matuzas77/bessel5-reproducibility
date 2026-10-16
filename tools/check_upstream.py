#!/usr/bin/env python3
"""Verify every upstream member against the exact, published source ZIP."""
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def check():
    provenance = json.loads((ROOT / "provenance/source.json").read_text())
    archive = ROOT / provenance["archive_path"]
    if hashlib.sha256(archive.read_bytes()).hexdigest() != provenance["archive_sha256"]:
        raise RuntimeError("Published source archive SHA-256 mismatch")
    if hashlib.md5(archive.read_bytes()).hexdigest() != provenance["zenodo_archive_md5"]:
        raise RuntimeError("Zenodo source archive MD5 mismatch")
    with ZipFile(archive) as source:
        names = {name for name in source.namelist() if not name.endswith("/")}
        actual = {str(p.relative_to(ROOT / "upstream"))
                  for p in (ROOT / "upstream").rglob("*") if p.is_file()}
        if names != actual:
            raise RuntimeError(f"Upstream file inventory mismatch: missing={names - actual}, extra={actual - names}")
        for name in sorted(names):
            if name.startswith("/") or ".." in Path(name).parts:
                raise RuntimeError("Unsafe archive path")
            if (ROOT / "upstream" / name).read_bytes() != source.read(name):
                raise RuntimeError(f"Upstream bytes differ from published archive: {name}")
    print(f"PASS: {len(names)} upstream files match the published v4 archive byte-for-byte.", flush=True)
    return len(names)


if __name__ == "__main__":
    check()

