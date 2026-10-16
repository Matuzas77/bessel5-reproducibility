#!/usr/bin/env python3
"""Check embedded-code parity and optionally execute in a clean directory.

All the work was done by AI.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tempfile

import nbformat
from nbclient import NotebookClient
from build_notebook import module_source, ROOT, DEST


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.execute and (args.output is None or args.report is None):
        parser.error("Execution requires --output notebook.ipynb and --report run.json")
    nb = nbformat.read(DEST, as_version=4)
    nbformat.validate(nb)
    sources = [cell.source for cell in nb.cells if cell.cell_type == "code"]
    for name in ("bessel.py", "lvalue.py", "compare.py"):
        if sources.count(module_source(name)) != 1:
            raise RuntimeError(f"Notebook code differs from independent/{name}")
        digest = hashlib.sha256((ROOT/"independent"/name).read_bytes()).hexdigest()
        if nb.metadata["source_module_sha256"][name] != digest:
            raise RuntimeError(f"Notebook provenance hash differs: {name}")
    print("Notebook schema and embedded-code parity passed.", flush=True)
    if args.execute:
        for cell in nb.cells:
            if cell.cell_type == "code":
                cell.outputs = []
                cell.execution_count = None
        with tempfile.TemporaryDirectory(prefix="bessel5-notebook-") as directory:
            NotebookClient(nb, timeout=600, kernel_name="python3",
                           resources={"metadata": {"path": directory}}).execute()
            reports = list(Path(directory).glob("bessel5-fresh-*.json"))
            if len(reports) != 1:
                raise RuntimeError("Expected exactly one freshly generated report")
            report = json.loads(reports[0].read_text())
            if not report["assessment"]["passed"]:
                raise RuntimeError("Notebook numerical assessment failed")
            args.report.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(reports[0], args.report)
        nbformat.validate(nb)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        nbformat.write(nb, args.output)
        print(f"Clean-kernel execution passed: {args.output}")
        print(f"Actual numerical evidence: {args.report}")


if __name__ == "__main__":
    main()
