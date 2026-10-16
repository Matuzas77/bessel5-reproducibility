# Bessel 5: reproduce the numerical evidence and exact witnesses

[![Reproduce published v4 checks](https://github.com/Matuzas77/bessel5-reproducibility/actions/workflows/verify.yml/badge.svg)](https://github.com/Matuzas77/bessel5-reproducibility/actions/workflows/verify.yml)
[![Paper DOI](https://img.shields.io/badge/paper-10.5281%2Fzenodo.22793737-blue)](https://doi.org/10.5281/zenodo.22793737)

Reproducibility companion to Jonas Matuzas, **[A quadratic critical-value conjecture for the fifth Bessel moment](https://doi.org/10.5281/zenodo.22793737)**, published source version 4 (16 September 2026).

**All the work was done by AI.** Jonas Matuzas is the human author and is responsible for the published claims and this release. The original manuscript's AI disclosure is preserved. This release is an AI-executed reproducibility check; it is not independent human peer review.

The paper proposes a formula for the unweighted fifth Bessel moment

$$A_0=\int_0^\infty K_0(t)^5\,dt.$$

For the paper's specified weight-three, level-60 newform, with $r=\sqrt5$ and $L=L(f,2)$, the proposed evaluation is

$$A_0\;\stackrel{?}{=}\;\frac\pi4\left[(125-13r)\operatorname{Re}(L^2)+(41r-125)\operatorname{Im}(L^2)\right].$$

The measure is **`dt`**, and the square is taken **before** the real or imaginary part. The exact normalization and form are specified in the [manuscript source](upstream/main.tex).

**The Bessel equality remains conjectural.** The manuscript separately presents proofs of modular norm, twist and twisted symmetric-square identities. Computation verifies finite witnesses and numerical enclosures; it does not formalize their analytic proofs or turn the Bessel conjecture into a theorem.

## Run it

Use Python **3.12.12**, the runtime used for the recorded release run and pinned in CI. The upstream programs require Python 3.11 or later; other runtimes have not been validated by this release.

```sh
git clone https://github.com/Matuzas77/bessel5-reproducibility.git
cd bessel5-reproducibility
python3.12 -m venv .venv
.venv/bin/python -m pip install --require-hashes --only-binary=:all: -r requirements.lock
.venv/bin/python -B tools/run_checks.py --suite all --output results/current
```

On Windows, use `.venv\Scripts\python.exe` in place of `.venv/bin/python`. The output folder must be empty; use a new `--output` path for another run. Do not use Python `-O`: upstream assertions are part of the checks.

No credentials, GPU, commercial software, or network access are needed for the verifiers after dependencies are installed. [Full reproduction instructions](docs/REPRODUCING.md) explain each stage and the recorded environment.

## What is actually checked?

| Suite | Computation | What a pass establishes |
|---|---|---|
| `stored` | Original SHA-256 manifest, stored certificate consistency, manuscript equation-status checks | The published package is intact and its saved finite checks are consistent. This does **not** regenerate numerical results. |
| `replay` | Ten jobs in a temporary directory, including decimal and binary directed arithmetic and exact modular witnesses | Nine freshly generated JSON records match the published records exactly. This **does** recompute the published interval certificates; it is not direct Bessel quadrature. |
| `primary` | Fifteen jobs reconstructing modular/Manin-symbol and critical-norm witnesses using SymPy and mpmath | The second exact reconstruction succeeds and its two critical-norm records match. It performs no Bessel quadrature and proves no Bessel equality. |

Both fresh runners use empty temporary output trees and leave the published files untouched. The wrapper checks all **65 upstream files against the original ZIP before and after execution**. The original internal manifest covers 64 files, excluding itself.

The published certificate bounds the individual-period discrepancy by **`10^-358`**, and each determinant discrepancy by **`10^-352`**. These finite bounds support the conjectures; they are not proofs of equality. [Claim boundaries and analytic inputs](docs/CLAIMS.md).

## Evidence and source preservation

The **16 October 2026 local run passed all three suites** on macOS arm64 with CPython 3.12.12: stored checks **0.033 s**, fresh interval replay **13.436 s**, primary reconstruction **15.132 s**. These are measured durations for that machine, not a runtime guarantee. Every upstream file was unchanged after the run.

- [Actual release-run logs and environment](results/2026-10-16-macos-arm64/) contain the stored check, fresh interval replay and primary reconstruction logs.
- [Unmodified published source](upstream/) includes the manuscript, full mathematical arguments, engines, exact endpoints and original logs.
- [Source provenance](provenance/source.json) records the source DOI, exact archive SHA-256, Zenodo MD5 and per-file hashes. The original ZIP and a public Zenodo metadata snapshot are included.
- [GitHub Actions](https://github.com/Matuzas77/bessel5-reproducibility/actions/workflows/verify.yml) runs the same three suites on Ubuntu and saves its actual logs. A badge reports workflow state, not the truth of a conjecture.
- [Release notes](docs/RELEASE_NOTES.md) separate this packaging release **v1.0.0** from the scientific source **v4**. No scientific source or saved numerical record has been edited.

## Read and cite

- [Published v4 manuscript and source, Zenodo](https://doi.org/10.5281/zenodo.22793737)
- [arXiv manuscript](https://arxiv.org/abs/2609.22318) — arXiv and Zenodo use different version numbering; this repository reproduces the specified Zenodo v4 archive.
- [Author's research website](https://jonas-matuzas-research.bronzeharp16.chatgpt.site/)
- [Bessel-moment topic guide](https://jonas-matuzas-research.bronzeharp16.chatgpt.site/research/bessel-moments/)
- [ResearchGate paper](https://www.researchgate.net/publication/414561930)

Use [CITATION.cff](CITATION.cff) or cite the paper DOI. This GitHub packaging release has **no separately minted DOI**. To identify an exact computational run, also record the Git commit/release tag and attach its `run.json` and logs.

## License and attribution

The Zenodo source record specifies **CC BY 4.0**. That license and the original attribution are preserved; the added documentation and orchestration are also released under CC BY 4.0. See [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md). The license does not imply peer review, endorsement or a mathematical correctness guarantee.

