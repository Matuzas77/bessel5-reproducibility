# Version 4 audit — the twisted symmetric-square norm bridge

**Manuscript:** A quadratic critical-value conjecture for the fifth Bessel moment  
**Author:** Jonas Matuzas  
**Revision date:** 16 September 2026

## Verdict and exact scope

The uploaded bridge supplies the missing independent modular argument. Equation (3.5) is now an ordinary equality, stated as Proposition 3.2:

\[
\mathscr L(2)=\sqrt{15}\pi^2\langle f,f\rangle_{60}
=\frac{3\sqrt{15}(5-\sqrt5)}{2\pi^2}|L(f,2)|^2,
\qquad \mathscr L(s)=L(\chi_{-4}\operatorname{Sym}^2 f,s).
\]

The new intermediate evaluation is

\[
\mathscr L(3)=\frac{\pi^4}{8}\langle f,f\rangle_{60}.
\]

The Petersson integral has no division by the group index 144. The form, character, coefficient embedding and classical unshifted L-function convention are unchanged. This is the Bessel5 level-60 problem, not the separate seven-Bessel level-525 calculation.

The main pure-moment formula (1.2) and the companion Bessel comparisons (3.3) and (3.4) remain conjectural. Proposition 3.2 proves that the two companion descriptions (3.3) and (3.4) are equivalent; it does not prove either equality with a Bessel integral. It also does not imply the individual-entry formula (1.2).

## What was checked in the analytic proof

A formal Poincare sum at a critical exponent is not sufficient. The audit checked the actual limiting argument in the uploaded report, before interpreting any finite Fourier computation as a proof.

The weight-three Poincare-series normalization gives
\(\langle f,\mathcal P_j\rangle=a_j/(16\pi^2j^2)\).
Consequently, for the square-index kernel
\(\mathcal K_u=\sum_{m\ge1}\chi_{-4}(m)m^{1-2u}\mathcal P_{m^2}\),
its pairing with f equals \(D(3+2u)/(16\pi^2)\) for real u in the initial convergence region. The complex-parameter statement is used with the sesquilinear pairing convention kept explicit.

The full proof now spells out an initial overlap strip \(1/4<\Re u<1/2\) for Poisson rearrangement. The quadratic Gauss bound gives an absolute c-sum there, and the transform supplies the required summability in the Fourier index. The Weber--Schafheitlin constants were independently reduced using gamma reflection and duplication, starting from DLMF 10.22.56.

The important exceptional cases are the principal and zero discriminants. The local factor at 2 has a zero at u=0 when \(\Delta=60v^2\), cancelling the potentially principal zeta pole before the transform is truncated. For \(\Delta=0\), the revised full proof gives an explicit formula rather than relying on the nonzero-discriminant factorization: if \(n=s^2\) and \(\ell=2s\), the local factor is zero identically for even s and, for odd s, equals

\[
2\chi_{-4}(s)\frac{t(1-2t)}{1-t},\qquad t=2^{-1-4u}.
\]

It vanishes at u=0 and has derivative \(8\chi_{-4}(s)\log2\). The remaining Fourier sum has a locally uniform summable bound. Thus the critical limit is a genuine cusp form, not an unproved termwise deletion of infinitely many hyperbolic terms. These limiting arguments are written mathematics; the finite tests below do not replace them.

The critical coefficients then reduce to a finite rational class-number formula. The exact eta witness and Sturm's bound 36 determine the holomorphic kernel's projection onto the four-form space. Hecke adjointness, with \(T_7^*=-T_7\), removes the complementary witness. The exact spectral coefficient is

\[
\delta=\frac{(7-i\sqrt5)(7-i)}{576}.
\]

Finally, the full Euler-factor identity was checked, including the missing primes in the zeta factor:

\[
\zeta^{(30)}(2s-4)D(s)=R_0(s)\mathscr L(s),\quad
R_0(s)=\bigl(1+(1+4i\sqrt5)3^{-s}\bigr)
       \bigl(1+(15+20i)5^{-s}\bigr).
\]

At s=3, the spectral and Euler factors give exactly
\(128\overline\delta/(75R_0(3))=1/8\).
The conductor-960 completion then gives \(\mathscr L(2)=8\sqrt{15}\mathscr L(3)/\pi^2\). The ordinary functional-equation relation, the trivial local factor at 2 and the gamma/conductor normalization were checked against Lim--Tu--Yu, Appendix A.3 and equation (44). Their question-mark Bessel assertion was not used as a theorem.

## Exact reconstructions and additional tests

The uploaded release's original verifier passed before integration. Its primary implementation reconstructs the class-number coefficients through q^200; its separately written standard-library implementation reconstructs them through q^240. All 201 common coefficients, including the constant term, agree exactly. Both construct the eta coefficients through q^1764 and check the rational Hecke witness through the Sturm bound. The same spectral and Euler constants are obtained without numerical periods or recognition of constants.

Version 4 adds 16 direct cyclotomic-ring checks at higher prime powers and at zero, principal and nonprincipal discriminants. Eight have nonzero answers and eight are zero controls. Together with the uploaded primary script's 20 direct cyclotomic checks, this gives 36 direct checks of the finite Gauss/CRT specialization. A separate set of 117 raw character sums tests the zero-discriminant boundary formula. Exact gamma-function and spectral-scalar checks are included. These are checks of specified formulas, not 153 independent analytic theorems.

The standard-library publication replay now runs ten jobs and regenerates nine JSON records. The optional SymPy/mpmath reconstruction runs fifteen jobs and regenerates the additional primary bridge and boundary-audit records. Both routes passed during preparation; their complete logs are included. They were also run from a fresh extraction of the final source ZIP.

All eight original v3 output records are retained byte-for-byte. The source and modular interval calculations were rerun, not merely copied or compared as pre-existing decimals. They reproduce the strict absolute bounds below 10^-358 for the pure-moment discrepancy and below 10^-352 for each displayed determinant discrepancy. Every residual interval contains zero. These are unchanged numerical inequalities, not equality proofs.

## Manuscript and packaging changes

The main article remains four A4 pages in the original 11-point Latin Modern template, with the same author/email block and ChatGPT-6-pro disclosure. Equation labels (1.1)--(3.5) are preserved. Only the last numbered equation changes proof status. Proposition 3.2 includes a proof sketch; the complete argument and rational witness are supplied in `proofs/critical_norm_bridge.md` and `proofs/critical_norm_certificate.md`.

The abstract and explanatory paragraphs now state the exact modular result and the equivalence of the companion conjectures. Earlier statements that the modular bridge was unproved have been superseded. The v2/v3 audits are retained only in `history/` with an explicit warning that they are historical. All eleven prior references are kept, with DLMF added for the transform calculation. Chuang's sole authorship and Yu's appendix collaboration, the Broadhurst--Roberts prior framework and the distinction between the level-15 weighted sector and the present dt integral remain explicit.

The four-page PDF was compiled without undefined citations, cross-references, overfull boxes or LaTeX warnings, and all four rendered pages were inspected. A clean source extraction reproduces the PDF byte-for-byte. The SHA-256 manifest covers all source-package files other than itself. The OpenDocument Text export is an additional editable copy with native MathML objects and resolved references; it reflows to six pages in the tested office renderer. The PDF and LaTeX are the authoritative four-page publication layout.

## Reproduction and qualifications

From the source ZIP's root:

```sh
python -B verify_package.py
python -B anc/verify.py
python -B anc/verify_primary.py
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The first command checks integrity and stored cross-consistency, the second recomputes the standard-library certificates and interval claims, and the third reconstructs the optional exact witnesses. Do not use Python's `-O` option.

No Bessel equality, algebraic correspondence, unrestricted novelty claim, external human review, or proof-assistant formalization is asserted. The public input theorems are cited dependencies, not re-proved here. The new proof is an independently audited written analytic argument together with finite computational witnesses. Nothing was uploaded to Zenodo or arXiv and no email was sent in preparing this version.
