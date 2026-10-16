#!/usr/bin/env python3
"""Generate a self-contained notebook from the separate check implementations.

All the work was done by AI. Requires nbformat==5.10.4 to rebuild.
"""
import ast
import hashlib
from pathlib import Path
import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "notebooks/Bessel5_independent_check.ipynb"


def module_source(name):
    source = (ROOT / "independent" / name).read_text()
    tree = ast.parse(source)
    # CLI entry points are not part of the interactive definitions.
    for node in tree.body:
        if isinstance(node, ast.If) and ast.unparse(node.test) == "__name__ == '__main__'":
            source = "\n".join(source.splitlines()[:node.lineno-1]) + "\n"
    return source.strip()


def build():
    md, code = nbf.v4.new_markdown_cell, nbf.v4.new_code_cell
    cells = [md(r"""# A fresh numerical check of the fifth Bessel moment

**All the work was done by AI.** Jonas Matuzas is responsible for the claims.

This notebook independently recomputes both sides of the conjecture in
[the published v4 paper](https://doi.org/10.5281/zenodo.22793737).
It uses new numerical routes rather than the paper's certificate programs or stored answers.
It is an AI-executed calculation, **not independent human peer review**.

**The equality remains conjectural.** Agreement at many digits is numerical evidence, not a proof.

Open in Google Colab, select **Runtime → Run all**, and use the ordinary CPU runtime.
No GPU, mounted Google Drive, private data or access tokens are needed.
The notebook contains its calculation code; only the pinned `mpmath` package may need downloading.
The two default runs use 40 and 60 significant digits, with 20 guard digits internally.
Runtime depends on the hosted CPU. All outputs below must be regenerated to count as your own run.

[GitHub source and method notes](https://github.com/Matuzas77/bessel5-reproducibility/blob/main/docs/INDEPENDENT_CHECK.md)
"""), md(r"""## 1. Define exactly what is being tested

The measure is **\(dt\)**, not \(t\,dt\):

\[
A_0=\int_0^\infty K_0(t)^5\,dt
\stackrel{?}{=}\frac{\pi}{4}\left[(125-13\sqrt5)\Re(L^2)
+(41\sqrt5-125)\Im(L^2)\right].
\]

Here \(L=L(f,2)=4\pi^2\int_0^\infty y f(iy)\,dy\), the classical unshifted value.
The **complex square** is taken before the real and imaginary parts.
With \(\eta_d=\eta(d\tau)\), define

\[
P_1=\eta_2\eta_6^2\eta_{10}^4/\eta_{30},\quad
P_2=\eta_2^4\eta_{10}\eta_{30}^2/\eta_6,
\]
\[
P_3=\eta_6^4\eta_{10}^2\eta_{30}/\eta_2,\quad
P_4=\eta_2^2\eta_6\eta_{30}^4/\eta_{10},
\]
\[
\omega=(1-2i)/\sqrt5,\quad
\epsilon=(5+4\sqrt5+i(2\sqrt5-10))/15,
\quad f=P_1+\omega P_2+3\epsilon P_3+3\epsilon\bar\omega P_4.
\]

Freshly generated Fourier coefficients check the embedding
\(a_1=1, a_3=\sqrt5-2i, a_5=\sqrt5(1+2i), a_7=-8i\).
No saved integral, L-value, coefficient table or certificate is a numerical input.
"""), md("## 2. Pin the numerical library and record the runtime"), code('''import sys, platform, subprocess, tempfile
from importlib.metadata import version, PackageNotFoundError
try:
    needs_install = version("mpmath") != "1.3.0"
except PackageNotFoundError:
    needs_install = True
if needs_install:
    requirement = "mpmath==1.3.0 --hash=sha256:a0b2b9fe80bbcd81a6647ff13108738cfb482d481d826cc0e02f5b35e5c88d2c\\n"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt") as handle:
        handle.write(requirement)
        handle.flush()
        subprocess.check_call([sys.executable, "-m", "pip", "install", "--quiet",
                               "--require-hashes", "--only-binary=:all:", "-r", handle.name])
import mpmath
if mpmath.__version__ != "1.3.0":
    raise RuntimeError("Restart the runtime after installing the pinned mpmath version")
print("Python:", sys.version)
print("Platform:", platform.platform())
print("mpmath:", mpmath.__version__)
print("All the work was done by AI.")'''), md(r"""## 3. Compute the Bessel side directly

The code below calls `mpmath.besselk(0, t)` and quadratures it directly.
The logarithmic singularity near zero is handled by a coordinate change.
An alternate run changes that coordinate and the interval partition.
Analytic endpoint-tail bounds are recorded, but adaptive quadrature and rounding
are **not certified interval arithmetic**.
""")]
    cells.append(code(module_source("bessel.py")))
    cells.extend([md(r"""## 4. Independently compute the modular side

Dedekind eta is evaluated from its Euler pentagonal expansion.
The transformation \(\eta(ix)=\eta(i/x)/\sqrt{x}\) keeps the series rapidly convergent.
We integrate the defining Mellin transform, changing to \(u=\log y\) in the alternate run.
This computation does not use the direct Bessel calculation or its answer.
Both routes share the `mpmath` arithmetic library, so they are not fully independent software implementations.
"""), code(module_source("lvalue.py")), md(r"""## 5. Run both calculations from scratch

Default: 40-digit primary and 60-digit alternate runs. For a longer check,
change `PRECISIONS` to `[40, 60, 80]` and rerun this cell and all cells below it.
Run each side completely before comparing; never use the proposed equality to obtain either input.
"""), code('''from datetime import datetime, timezone
import json, time
PRECISIONS = [40, 60]
if len(PRECISIONS) < 2 or PRECISIONS != sorted(set(PRECISIONS)) or min(PRECISIONS) < 30:
    raise ValueError("Choose at least two increasing precisions >=30")
runs = []
started_utc = datetime.now(timezone.utc).isoformat()
for index, digits in enumerate(PRECISIONS):
    print(f"Computing {digits} digits...", flush=True)
    bessel = compute_bessel(digits, "primary" if index % 2 == 0 else "alternate")
    modular = compute_lvalue(digits, "primary" if index % 2 == 0 else "shifted")
    runs.append({"dps": digits, "bessel": bessel, "modular": modular})
    print("Direct integral:", bessel["value"])
    print("Modular RHS:    ", modular["rhs"])
    print("L real:         ", modular["real"])
    print("L imaginary:    ", modular["imag"])
    print("Seconds:", round(bessel["elapsed_seconds"] + modular["seconds"], 3), "\\n")'''), md(r"""## 6. Convergence and deliberate failure controls

Before seeing the results, set the relative comparison threshold at
\(10^{-(p-10)}\) for a run requesting \(p\) significant digits.
For adjacent runs, test convergence at the threshold of the lower precision.
Decimal strings are compared with extra working precision, not ordinary Python floats.

Three controls deliberately use the wrong complex embedding, omit the Mellin factor,
or replace \(L^2\) by \(|L|^2\). Each must disagree by more than 0.1%.
A synthetic relative perturbation of \(10^{-12}\) must also be detected.
These test sensitivity; they do not prove the original formula.
"""), code(module_source("compare.py")), code('''assessment = assess_checks(runs)
print(json.dumps(assessment, indent=2))
if not assessment["passed"]:
    raise ArithmeticError("A convergence, agreement or negative-control check failed")
print("PASS: empirical agreement and convergence; incorrect controls rejected.")
print("The Bessel equality remains conjectural. All the work was done by AI.")'''), md("## 7. Preserve the complete fresh-run evidence"), code('''from pathlib import Path
report = {
    "started_utc": started_utc,
    "finished_utc": datetime.now(timezone.utc).isoformat(),
    "python": sys.version, "platform": platform.platform(),
    "mpmath": mpmath.__version__, "runs": runs, "assessment": assessment,
    "input_policy": "Definitions only; no upstream numerical data or saved answers read.",
    "disclosure": "All the work was done by AI."
}
output_path = Path("bessel5-fresh-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") + ".json")
output_path.write_text(json.dumps(report, indent=2) + "\\n")
print("Saved", output_path)
print("Download this JSON and your executed notebook to retain the settings and results.")'''), md(r"""## What this establishes—and what remains open

A passing run supplies fresh numerical evidence using direct Bessel quadrature and
the defining eta-Mellin integral, with changed precision and coordinates.
It does not provide a rigorous global error enclosure or a proof of equality.
The original certificate pipeline is a separate method, preserved unmodified in `upstream/`.
No inference of external peer review or human verification is intended.

**All the work was done by AI.** Jonas Matuzas retains responsibility for the claims.

References: [paper and source v4](https://doi.org/10.5281/zenodo.22793737),
[DLMF Bessel integral representations](https://dlmf.nist.gov/10.32),
[DLMF eta transformations](https://dlmf.nist.gov/23.18),
[mpmath integration documentation](https://mpmath.org/doc/current/calculus/integration.html).
""")])
    for index, cell in enumerate(cells):
        if cell.cell_type == "markdown":
            # Colab's Markdown parser requires dollar delimiters for math.
            cell.source = (cell.source.replace(r"\(", "$" ).replace(r"\)", "$")
                           .replace(r"\[", "$$").replace(r"\]", "$$"))
        cell.id = hashlib.sha256(f"{index}:{cell.source}".encode()).hexdigest()[:16]
    nb = nbf.v4.new_notebook(cells=cells, metadata={
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12.12"},
        "colab": {"name": DEST.name},
        "source_module_sha256": {name: hashlib.sha256((ROOT/"independent"/name).read_bytes()).hexdigest()
                                 for name in ("bessel.py", "lvalue.py", "compare.py")}})
    nbf.validate(nb)
    DEST.parent.mkdir(exist_ok=True)
    nbf.write(nb, DEST)
    print(DEST)


if __name__ == "__main__":
    build()
