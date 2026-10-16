> HISTORICAL AUDIT — superseded by version 4 for equation (3.5).
> The current mathematical status is in ../revision_audit_v4.md and
> ../proofs/critical_norm_bridge.md. This file records the earlier audit only.

# Audit of the six criticisms

**Article:** *A quadratic critical-value conjecture for the fifth Bessel moment*  
**Author:** Jonas Matuzas  
**Audited revision:** 15 September 2026  
**Result:** revised four-page manuscript and reproducible source package; not submitted to arXiv.

## Executive finding

The criticism correctly identifies two omitted modular proofs, a bibliographic authorship error, and insufficient contextual explanation. Its proposed promotion of the twisted symmetric-square value to a theorem is not justified by the recovered norm and twist certificates. Its statements that the cited author block contains an affiliation and that Lim–Tu–Yu leave the relevant `c_4` blank are contradicted by the public PDFs. The even/odd determinant factor does have a complete exact derivation from Chuang's published period formulas; it is not merely a numerical observation.

The revision therefore adds the genuine modular proposition, gives the determinant derivation in the article, corrects the citation, expands the lineage, and preserves the unproved question marks.

## 1. Proved modular results: recovered and added, with an important limitation

### What was actually recovered

The complete proof `NORM_PROOF.md` and its exact programs were recovered from the author's 5 September archive. The 6 September archive supplies the full critical-value twist argument, all 33 raw rational Manin coefficients and two exact reconstructions. The original norm archive and the exact part of the original twist archive were rerun successfully during this audit. Publication-local versions are now in `proofs/norm_proof.md`, `proofs/twist_proof.md` and `proofs/MODULAR_SYMBOL_CERTIFICATE.md`.

Proposition 2.1 now states, for the fixed embedding and `r=sqrt(5)`,

\[
\langle f,f\rangle_{60}=
\frac{3(5-r)}{2\pi^4}|L(f,2)|^2,
\qquad
L(f^\sigma,2)=\frac{(r-1)(4+3i)}{10}L(f,2),
\]

where `sigma(r)=-r`, `sigma(i)=-i` acts on the Fourier coefficients. The Petersson integral is **not divided by the group index 144**. Coefficient conjugation is not interpreted as an arbitrary automorphism acting on a transcendental integral.

The norm proof uses a weight-one Eisenstein series `g`, a regular weight-two Eisenstein series at infinity, and `H=g E_h-theta(g)/12`. Stokes removes the raising term. The complete unfolded Dirichlet series is

\[
\mathcal D(s)=
\frac{L(f,s)L(\bar f,s)(1-\bar a_3 3^{-s})(1-\bar a_5 5^{-s})}
{\zeta^{(30)}(2s-2)}.
\]

An explicit 14-eta-term plus two-Eisenstein-term witness satisfies

\[
H-H_f=(T_7^2+64)R.
\]

Its coefficient identity is verified beyond the Sturm bound. Since `T_7^*=-T_7` and `T_7 f=-8if`, the witness contributes zero to the Petersson pairing. The exact spectral coefficient and finite Euler factors give the stated norm constant. The written analytic argument and the exact finite witness are both necessary; a count of passing software assertions alone is not the proof.

For the twist, the certificate proves

\[
w_5=w\frac{51-13U_3+3U_3^2-U_3^3}{30}+R_5(T_7^2+64).
\]

The level-300 identity is the **depleted** twist

\[
f\otimes\chi_5=f^\sigma-a_5(f^\sigma)f^\sigma(5\tau),
\qquad a_5(f^\sigma)=-r+2ir.
\]

Keeping this prime-5 term is essential to the multiplier. The current standard-library verifier rebuilds the modular-symbol quotient from 288 generators and 576 relations, using independent projective representatives and ceiling continued fractions. It does not initialize from the primary reduced Hecke matrices.

### Why the old equation (2.5) is not promoted

Write

\[
\mathscr L(s)=L(\chi_{-4}\operatorname{Sym}^2 f,s).
\]

The extra comparison needed is

\[
\mathscr L(2)\stackrel?=\sqrt{15}\pi^2\langle f,f\rangle_{60}.
\]

Only the next step,

\[
\sqrt{15}\pi^2\langle f,f\rangle_{60}
=\frac{3\sqrt{15}(5-r)}{2\pi^2}|L(f,2)|^2,
\]

is supplied by the recovered norm theorem. The degree-three function `mathscr L` is not the Rankin series `mathcal D` displayed above. Neither recovered certificate evaluates `mathscr L(2)` or gives the additional normalized comparison with the norm. Both sides being modular quantities is not an equality proof.

Accordingly the old equation (2.5), now equation (3.5), retains its question mark. This finding is a proof-status boundary, not a refutation of the proposed value identity. The paper is no longer conjecture-only, but its entire set of assertions cannot yet be described as one proved modular identity plus just one open scalar equality.

## 2. Affiliation: the cited public author block has none

The first page of [arXiv:2609.00290](https://arxiv.org/pdf/2609.00290) was inspected as a page image. It prints the name **Jonas Matuzas** and **jonas.matuzas@gmail.com**, with no affiliation. That public paper is titled *Shifted poles and chamber cancellation for classical Witten zeta functions*.

The revised author block matches those two lines. Adding an inferred institution or location would not make it identical to the cited publication. The recovered source supplies the matching typography; this audit also verifies the actual public first page rather than relying on the source filename alone.

## 3. Chuang authorship: corrected

The bibliography now gives:

> P.-H. Chuang, *On the periods of twisted moments of the Kloosterman connection* (with an appendix in collaboration with J.-D. Yu), Ramanujan J. 65 (2024), no. 3, 1227–1277; arXiv:2306.15216v2.

This matches [Lim–Tu–Yu, reference Chu24](https://arxiv.org/pdf/2609.13576v1), printed page 35, and the [Chuang preprint](https://arxiv.org/html/2306.15216v2). The authoritative article and ancillary prose no longer present Yu as coauthor of the complete Ramanujan Journal article. Yu's appendix contribution remains explicitly credited.

## 4. The determinant factor is an exact consequence, not just a decimal test

Let `A_j`, `B_j`, `C_j` be the integrals of `t^(2j) K_0^5`, `t^(2j) I_0 K_0^4` and `t^(2j) I_0^2 K_0^3`. In [Lim–Tu–Yu, equation (42)](https://arxiv.org/pdf/2609.13576v1),

\[
D_{5,\mathrm{odd}}=B_0,
\qquad D_{5,\mathrm{even}}=A_0C_1-A_1C_0.
\]

Take the three-by-three matrix `X` with rows `A_j`, `pi B_j`, `pi^2 C_j`, indexed from zero. Chuang's Propositions 13, 18, 24 and 32, together with Corollary 39, specialize to

\[
XJX^t=64\pi^6 H_B,
\quad H_B=\begin{pmatrix}1/4&0&1/8\\0&-1/10&0\\1/8&0&3/40\end{pmatrix},
\quad J_{2,*}=(120,0,0),
\quad\det X=\frac{2\pi^9}{75\sqrt{15}}.
\]

Thus

\[
(X^{-1})_{2,1}=-\frac{75 B_0}{4\pi^5}.
\]

The same inverse entry, evaluated by its cofactor, is

\[
(X^{-1})_{2,1}=-\frac{\pi^2D_{5,\mathrm{even}}}{\det X}.
\]

Equating these expressions gives

\[
D_{5,\mathrm{even}}=\frac{\pi^2}{2\sqrt{15}}D_{5,\mathrm{odd}}.
\]

All phase changes, row and column scalings and the determinant sign are recorded in `verification_notes.md`, section 4. The raw Chuang Betti matrix, inverse metric and cofactor constants are checked again in `engine/verify_modular_bridge.py`. Corollary 3.1 now gives this proof in the article itself and labels it as a consequence of Chuang. No Bessel-to-modular conjecture is used.

## 5. Motivation: useful, but the proposed Table 1 explanation is false

The revised introduction explains the limited and relevant distinction: a quadratic expression involving `L(f,2)^2` retains phase information, while `|L(f,2)|^2` supplies a norm. The determinant relation alone does not determine the individual pure entry `A_0`.

However [Lim–Tu–Yu, Table 1](https://arxiv.org/pdf/2609.13576v1), printed page 34, does **not** leave `c_4` blank. Substituting `k=5=4r+1`, `r=1`, in its displayed row gives

\[
c_3=\pi^{-2}D_{5,\mathrm{odd}},
\qquad c_4=\pi^{-2}D_{5,\mathrm{even}}.
\]

The table specifies critical-period representatives, not a finite individual-entry evaluation of `A_0`. The revised text says exactly this.

The statement that every published five-Bessel evaluation is a norm or determinant is also too broad. [Zhou's level-15 work](https://arxiv.org/abs/1706.08308) includes weighted five-Bessel evaluations in the `t dt` sector. The article explicitly distinguishes that sector from the present `dt` integral. The motivation is not presented as a proof that a PSLQ origin is impossible, a proof of novelty, or a mathematical explanation for all earlier non-solutions.

## 6. Bibliographic lineage: expanded and dates corrected

The article now has ten references. The four requested additions are:

- D. Broadhurst, *Feynman integrals, L-series and Kloosterman moments*, Commun. Number Theory Phys. 10 (2016), no. 3, 527–569; [arXiv:1604.03057](https://arxiv.org/abs/1604.03057).
- R. Evans, *Hypergeometric 3F2(1/4) evaluations over finite fields and Hecke eigenforms*, Proc. Amer. Math. Soc. 138 (2010), no. 2, 517–531. The relevant modularity conjectures are explicitly attributed to this paper by Lim–Tu–Yu, introduction and Theorem 1.2.2.
- J. Fresán, C. Sabbah and J.-D. Yu, *Quadratic relations between Bessel moments*, Algebra Number Theory 17 (2023), no. 3, 541–602; [arXiv:2006.02702](https://arxiv.org/abs/2006.02702).
- Y. Zhou, *Wick rotations, Eichler integrals, and multi-loop Feynman diagrams*, Commun. Number Theory Phys. 12 (2018), no. 1, 127–192; [arXiv:1706.08308](https://arxiv.org/abs/1706.08308), a 2017 preprint.

The prose credits Broadhurst's research program, rather than claiming without a pinpoint source that he posed this exact newly displayed formula. Zhou's preprint year and journal year are both specified.

## Reproduction and remaining scope

The current default replay is

```sh
python -B anc/verify.py
```

It requires only Python's standard library, regenerates the source and modular interval calculations, rebuilds the exact proof witnesses, and checks all eight distributed JSON records. The optional second reconstruction is

```sh
python -B anc/verify_primary.py
```

and requires SymPy and mpmath. It runs thirteen exact reconstruction jobs. Both replays completed successfully during preparation; their logs are under `logs/`. All numerical residual intervals contain zero. The advertised strict absolute bounds remain `10^-358` for the pure integral and `10^-352` for each of the two determinant comparisons. Neither the high-precision test nor its successful replay proves exact equality.

The article compiles to four pages with resolved citations and cross-references; all four rendered pages were visually inspected. Exact modular proof certificates are not formal proof-assistant verification, and the public input theorems have not been independently re-proved in this package. No external human review, unrestricted novelty conclusion or arXiv submission is claimed.
