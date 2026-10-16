#!/usr/bin/env python3
"""Fresh direct quadrature and eta-Mellin check. All the work was done by AI."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from independent.bessel import compute_bessel
from independent.lvalue import compute_lvalue
from independent.compare import assess_checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--digits", default="40,60")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    digits = [int(x) for x in args.digits.split(",")]
    if len(digits) < 2 or digits != sorted(set(digits)) or min(digits) < 30:
        parser.error("Use at least two increasing, distinct precisions >=30")
    if args.output.exists():
        parser.error("Output already exists; use a new path to preserve evidence")
    report = {"started_utc": datetime.now(timezone.utc).isoformat(),
              "python": sys.version, "platform": platform.platform(),
              "base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "module_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in sorted((ROOT/"independent").glob("*.py"))},
              "input_policy": "Only mathematical definitions; no upstream code, archived results or certificates are read.",
              "disclosure": "All the work was done by AI.", "runs": []}
    for i, dps in enumerate(digits):
        print(f"Computing fresh direct Bessel integral at {dps} digits...", flush=True)
        bessel = compute_bessel(dps, "primary" if i % 2 == 0 else "alternate")
        modular = compute_lvalue(dps, "primary" if i % 2 == 0 else "shifted")
        report["runs"].append({"dps": dps, "bessel": bessel, "modular": modular})
        print(f"  Integral = {bessel['value']}\n  RHS      = {modular['rhs']}", flush=True)
    report["assessment"] = assess_checks(report["runs"])
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report["assessment"], indent=2))
    return 0 if report["assessment"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
