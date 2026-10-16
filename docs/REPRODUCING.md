# Reproduction guide

All the work was done by AI. Jonas Matuzas is responsible for the claims and release. Running these programs provides reproducibility evidence, not independent peer review.

## Environment

The release uses CPython 3.12.12 in an isolated virtual environment, SymPy 1.14.0 and mpmath 1.3.0. The two dependency versions and their pure-Python wheel hashes are pinned in `requirements.lock`; `.python-version` records the interpreter version. CI selects the same interpreter on Ubuntu 24.04. The exact local OS, interpreter, dependency versions and timings are in the release [run.json](../results/2026-10-16-macos-arm64/run.json).

The scientific programs use CPU and RAM. No direct network access is made by the inspected verifier code. The dependency install requires internet access unless wheels are cached. Runtime is hardware dependent; use the measured durations in `run.json` rather than a promised universal completion time. Each upstream subprocess has a 240-second timeout; slow systems may need a separately documented timeout-only patch, which this release does not apply.

The standard-library replay can also be used without installing SymPy or mpmath:

```sh
python3.12 -m venv .venv
.venv/bin/python -B tools/run_checks.py --suite stored --output results/stored-run
.venv/bin/python -B tools/run_checks.py --suite replay --output results/replay-run
```

For the primary reconstruction:

```sh
.venv/bin/python -m pip install --require-hashes --only-binary=:all: -r requirements.lock
.venv/bin/python -B tools/run_checks.py --suite primary --output results/primary-run
```

Use an empty output directory each time. The wrapper fails if a previous log would be overwritten. It also rejects `-O` and `PYTHONOPTIMIZE`, because upstream assertions are required.

## Direct upstream commands

These invoke the original programs without the release wrapper:

```sh
.venv/bin/python -B upstream/verify_package.py
.venv/bin/python -B upstream/anc/verify.py
.venv/bin/python -B upstream/anc/verify_primary.py
```

Do **not** add `--write`: that option deliberately replaces upstream reference outputs and is unnecessary for reproduction. The release wrapper never passes it.

`verify_package.py` reads stored files. `anc/verify.py` freshly recomputes the directed intervals and finite witnesses and compares the whole nine JSON records. `anc/verify_primary.py` reconstructs the primary exact algebra, comparing the critical-norm records as well. Its word “independent” describes an alternative implementation/reconstruction within the supplied scientific package, not an independent research group.

## Understanding a result

`run.json` records return codes, elapsed times, package versions and SHA-256 hashes of the actual logs. A successful three-suite run has:

- `status: "PASS"` and three jobs with zero return codes;
- nine `EXACT RECORD MATCH` lines in `replay.log`;
- two `EXACT RECORD MATCH` lines and fifteen successful jobs in `primary.log`;
- `scientific_files_unchanged_before` and `scientific_files_unchanged_after` both true;
- `bessel_conjectures_proved_by_execution: false`.

A mismatch or failed assertion is a failed reproduction and must not be replaced by an old passing log. A timeout is incomplete execution, not evidence that the mathematics is wrong. Preserve the actual log and environment before investigating. The wrapper stops at the first failed suite.

The logs shipped under `upstream/anc/logs/` belong to the original source release. The separate `results/2026-10-16-macos-arm64/` directory records this release's fresh execution. GitHub Actions artifacts record subsequent CI runs and are subject to GitHub's retention period.

## Integrity and release archive

```sh
.venv/bin/python -B tools/check_upstream.py
.venv/bin/python -B tools/build_release.py
```

The first command compares every upstream file byte-for-byte with the included source ZIP, whose SHA-256 and Zenodo MD5 are checked. The second creates `dist/bessel5-reproducibility-v1.0.0.zip` and `dist/SHA256SUMS`. The ZIP has fixed timestamps, file permissions and ordering; it uses uncompressed entries for deterministic bytes. Its inner `RELEASE_SHA256SUMS` covers all payload files except the manifest itself.

The archive excludes Git internals, virtual environments, fresh `results/current` runs and unrelated local files. Its included dated logs document one specific execution, not every computer or future dependency version.

