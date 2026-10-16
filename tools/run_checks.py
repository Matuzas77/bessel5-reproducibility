#!/usr/bin/env python3
"""Run unchanged scientific verifiers and save actual logs and environment evidence."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import time

from check_upstream import check

ROOT = Path(__file__).resolve().parents[1]
SUITES = {
    "stored": "verify_package.py",
    "replay": "anc/verify.py",
    "primary": "anc/verify_primary.py",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=[*SUITES, "all"], default="all")
    parser.add_argument("--output", type=Path, default=ROOT / "results/current")
    args = parser.parse_args()
    if not __debug__ or os.environ.get("PYTHONOPTIMIZE"):
        raise RuntimeError("Assertions are required: do not use -O or PYTHONOPTIMIZE.")
    if sys.version_info < (3, 11):
        raise RuntimeError("The upstream verifier requires Python >= 3.11.")
    check()
    output = args.output.resolve()
    if ROOT / "upstream" == output or ROOT / "upstream" in output.parents:
        raise RuntimeError("Logs must not be placed inside the preserved upstream tree.")
    if output.exists() and any(output.iterdir()):
        raise RuntimeError("Output directory is not empty; choose a new --output path.")
    output.mkdir(parents=True, exist_ok=True)
    packages = {}
    for package in ["sympy", "mpmath", "pip"]:
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            pass
    record = {
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "packages": packages,
        "runtime_matches_pin": platform.python_version() == (ROOT / ".python-version").read_text().strip(),
        "scientific_files_unchanged_before": True,
        "bessel_conjectures_proved_by_execution": False,
        "jobs": [],
    }
    env = os.environ.copy()
    env.update(PYTHONHASHSEED="0", PYTHONDONTWRITEBYTECODE="1", PYTHONUNBUFFERED="1")
    for name in SUITES if args.suite == "all" else [args.suite]:
        command = [sys.executable, "-B", str(ROOT / "upstream" / SUITES[name])]
        started = time.monotonic()
        print(f"Running {name}: python -B upstream/{SUITES[name]}", flush=True)
        try:
            completed = subprocess.run(command, cwd=ROOT / "upstream", env=env,
                                       text=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, timeout=1800)
            text, returncode = completed.stdout, completed.returncode
        except subprocess.TimeoutExpired as error:
            text = error.stdout or b""
            if isinstance(text, bytes):
                text = text.decode("utf-8", errors="replace")
            text += "\nWRAPPER TIMEOUT after 1800 seconds; this suite did not pass.\n"
            returncode = 124
        text = text.replace(str(ROOT), "<repository>")
        (output / f"{name}.log").write_text(text)
        print(text, end="", flush=True)
        record["jobs"].append({
            "suite": name, "command": f"python -B upstream/{SUITES[name]}",
            "returncode": returncode, "elapsed_seconds": round(time.monotonic() - started, 3),
            "log": f"{name}.log", "log_sha256": hashlib.sha256(text.encode()).hexdigest(),
            "exact_record_matches": text.count("EXACT RECORD MATCH:"),
        })
        if returncode:
            break
    try:
        check()
        record["scientific_files_unchanged_after"] = True
    except Exception as error:
        record["scientific_files_unchanged_after"] = False
        record["integrity_error"] = str(error)
    record["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    record["status"] = "PASS" if all(job["returncode"] == 0 for job in record["jobs"]) and record["scientific_files_unchanged_after"] else "FAIL"
    (output / "run.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"{record['status']}: evidence saved to {output}", flush=True)
    return 0 if record["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())

