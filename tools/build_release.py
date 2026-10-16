#!/usr/bin/env python3
"""Build a deterministic ZIP of the public release with content checksums."""
from pathlib import Path
from zipfile import ZipFile, ZipInfo, ZIP_STORED
import hashlib
import json

from check_upstream import check

ROOT = Path(__file__).resolve().parents[1]
NAME = "bessel5-reproducibility-v1.0.0"
TOP_FILES = ["README.md", "CITATION.cff", "LICENSE", "NOTICE.md", "requirements.lock", ".python-version", ".gitignore"]
TREES = ["upstream", "tools", "docs", "provenance", ".github", "results/2026-10-16-macos-arm64"]


def main():
    check()
    paths = [ROOT / name for name in TOP_FILES]
    for tree in TREES:
        paths.extend(p for p in (ROOT / tree).rglob("*")
                     if p.is_file() and "__pycache__" not in p.parts and p.suffix != ".pyc")
    paths = sorted(paths, key=lambda p: str(p.relative_to(ROOT)))
    payload = {str(p.relative_to(ROOT)): p.read_bytes() for p in paths}
    manifest = "".join(f"{hashlib.sha256(data).hexdigest()}  {name}\n" for name, data in payload.items())
    payload["RELEASE_SHA256SUMS"] = manifest.encode()
    output = ROOT / "dist"
    output.mkdir(exist_ok=True)
    archive = output / f"{NAME}.zip"
    with ZipFile(archive, "w", compression=ZIP_STORED) as target:
        for name, data in sorted(payload.items()):
            info = ZipInfo(f"{NAME}/{name}", date_time=(2026, 10, 16, 0, 0, 0))
            info.compress_type = ZIP_STORED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            target.writestr(info, data)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    (output / "SHA256SUMS").write_text(f"{digest}  {archive.name}\n")
    print(json.dumps({"archive": archive.name, "sha256": digest, "files": len(payload), "bytes": archive.stat().st_size}, indent=2))


if __name__ == "__main__":
    main()

