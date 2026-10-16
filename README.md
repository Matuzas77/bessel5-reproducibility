# Bessel 5: reproduce the numerical evidence and exact witnesses

[![Reproduce published v4 checks](https://github.com/Matuzas77/bessel5-reproducibility/actions/workflows/verify.yml/badge.svg)](https://github.com/Matuzas77/bessel5-reproducibility/actions/workflows/verify.yml)
[![Paper DOI](https://img.shields.io/badge/paper-10.5281%2Fzenodo.22793737-blue)](https://doi.org/10.5281/zenodo.22793737)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Matuzas77/bessel5-reproducibility/blob/main/notebooks/Bessel5_independent_check.ipynb)

Reproducibility companion to Jonas Matuzas, **[A quadratic critical-value conjecture for the fifth Bessel moment](https://doi.org/10.5281/zenodo.22793737)**, published source version 4 (16 September 2026).

**All the work was done by AI.** Jonas Matuzas is the human author and is responsible for the published claims and this release. The original manuscript's AI disclosure is preserved. This release is an AI-executed reproducibility check; it is not independent human peer review.

The paper proposes a formula for the unweighted fifth Bessel moment

$$A_0=\int_0^\infty K_0(t)^5\,dt.$$

For the paper's specified weight-three, level-60 newform, with $r=\sqrt5$ and $L=L(f,2)$, the proposed evaluation is

$$A_0\;\stackrel{?}{=}\;\frac\pi4\left[(125-13r)\operatorname{Re}(L^2)+(41r-125)\operatorname{Im}(L^2)\right].$$

The measure is **`dt`**, and the square is taken **before** the real or imaginary part. The exact normalization and form are specified in the [manuscript source](upstream/main.tex).

**The Bessel equality remains conjectural.** The manuscript separately presents proofs of modular norm, twist and twisted symmetric-square identities. Computation verifies finite witnesses and numerical enclosures; it does not formalize their analytic proofs or turn the Bessel conjecture into a theorem.

## Run it

**In your browser:** [open the self-contained notebook in Google Colab](https://colab.research.google.com/github/Matuzas77/bessel5-reproducibility/blob/main/notebooks/Bessel5_independent_check.ipynb) and select **Runtime → Run all**. Use a CPU runtime. No GPU, Drive mount or access token is needed. [View the executed notebook on GitHub](notebooks/Bessel5_independent_check.ipynb).

The new calculation directly integrates `K0(t)^5` and separately integrates the modular form's defining eta expression. It reads no saved scientific answers or upstream numerical routines. Default 40/60-digit runs include changed coordinates, convergence checks and deliberately incorrect phase/normalization controls. The notebook contains all its calculation code and pins mpmath 1.3.0. Colab supplies the Python runtime; the recorded local/CI runtime is pinned separately below.

**Measured fresh result:** `135.268302580868837594226279646192207420305889359…`. The separate 80-digit direct and modular calculations have a relative discrepancy about **9.73 × 10⁻⁸¹** after decimal serialization, passing a predeclared **10⁻⁷⁰** threshold. This is empirical agreement, **not** an interval-certified error bound or a proof. Both implementations use mpmath; “independent” means new calculation routes, not independent software libraries or human review. [Methods and actual evidence](docs/INDEPENDENT_CHECK.md).

**Locally, for the original certificate suites:**

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
| `independent` notebook / `tools/run_independent.py` | Fresh direct K0 quadrature and eta-Mellin quadrature with varied precision and coordinates | Empirical agreement of separately computed values, convergence and rejection of wrong phase/normalization controls. No archived answers are inputs; global errors are not rigorously enclosed. |

Both fresh runners use empty temporary output trees and leave the published files untouched. The wrapper checks all **65 upstream files against the original ZIP before and after execution**. The original internal manifest covers 64 files, excluding itself.

The published certificate bounds the individual-period discrepancy by **`10^-358`**, and each determinant discrepancy by **`10^-352`**. These finite bounds support the conjectures; they are not proofs of equality. [Claim boundaries and analytic inputs](docs/CLAIMS.md).

## Evidence and source preservation

The **16 October 2026 local run passed all three suites** on macOS arm64 with CPython 3.12.12: stored checks **0.033 s**, fresh interval replay **13.436 s**, primary reconstruction **15.132 s**. These are measured durations for that machine, not a runtime guarantee. Every upstream file was unchanged after the run.

- [Actual release-run logs and environment](results/2026-10-16-macos-arm64/) contain the stored check, fresh interval replay and primary reconstruction logs.
- [Unmodified published source](upstream/) includes the manuscript, full mathematical arguments, engines, exact endpoints and original logs.
- [Source provenance](provenance/source.json) records the source DOI, exact archive SHA-256, Zenodo MD5 and per-file hashes. The original ZIP and a public Zenodo metadata snapshot are included.
- [GitHub Actions](https://github.com/Matuzas77/bessel5-reproducibility/actions/workflows/verify.yml) runs the same three suites on Ubuntu and saves its actual logs. A badge reports workflow state, not the truth of a conjecture.
- [New direct-check results and notebook run](results/2026-10-16-independent/) include 40/60/80-digit values, numerical settings, tail expressions, controls, source hashes and an AI review record.
- [Release notes](docs/RELEASE_NOTES.md) separate software releases **v1.0.0 / v1.1.0** from the scientific source **v4**. No scientific source or original saved numerical record has been edited.

## Read and cite

- [Published v4 manuscript and source, Zenodo](https://doi.org/10.5281/zenodo.22793737)
- [arXiv manuscript](https://arxiv.org/abs/2609.22318) — arXiv and Zenodo use different version numbering; this repository reproduces the specified Zenodo v4 archive.
- [Author's research website](https://jonas-matuzas-research.bronzeharp16.chatgpt.site/)
- [Bessel-moment topic guide](https://jonas-matuzas-research.bronzeharp16.chatgpt.site/research/bessel-moments/)
- [ResearchGate paper](https://www.researchgate.net/publication/414561930)

Use [CITATION.cff](CITATION.cff) or cite the paper DOI. This GitHub packaging release has **no separately minted DOI**. To identify an exact computational run, also record the Git commit/release tag and attach its `run.json` and logs.

## License and attribution

The Zenodo source record specifies **CC BY 4.0**. That license and the original attribution are preserved; the added documentation and orchestration are also released under CC BY 4.0. See [LICENSE](LICENSE) and [NOTICE.md](NOTICE.md). The license does not imply peer review, endorsement or a mathematical correctness guarantee.

