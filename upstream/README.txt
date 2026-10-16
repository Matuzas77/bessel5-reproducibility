A quadratic critical-value conjecture for the fifth Bessel moment
Jonas Matuzas -- version 4, 16 September 2026

Typeset master: main.tex (pdflatex; compile twice).
The PDF remains four pages. No arXiv or Zenodo upload is performed by this archive.

Mathematical update: Proposition 3.2 proves equation (3.5), the purely modular
twisted symmetric-square critical-value/Petersson-norm identity. Equations
(1.2), (3.3) and the Bessel equality in (3.4) remain conjectural.

All full proofs are under anc/proofs/. Exact rational endpoints and finite
modular certificates are under anc/outputs/.

Python 3.11+, standard library only:
  python -B anc/verify.py
This runs 10 jobs and reproduces 9 JSON records without modifying them.

Optional second exact reconstruction (SymPy and mpmath required):
  python -B anc/verify_primary.py
This runs 15 exact jobs and independently reconstructs the critical-norm
coefficients and witnesses, including boundary-case checks.

  python -B verify_package.py
checks every SHA-256, the equation status and preserved data.

Analytic continuation, Stokes and other all-index arguments are written
mathematics, not proof-assistant formalizations. Verification success is not
an equality proof for the Bessel conjectures. No new precision record is
claimed. The main AI disclosure is retained from version 3.
