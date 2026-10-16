# Exact modular-symbol certificate

This certificate proves the modular-symbol identity used in Theorem 1 of `twist_proof.md`. It proves no equality with a Bessel integral.

The primary scripts named below are in `../primary/twist/` and regenerate
`results/` in a temporary directory when `../verify_primary.py` is run.
The default replay also runs `../engine/verify_twist_symbols.py`, an I/O-adapted
copy of the independent verifier; it loads only the raw witness in
`../engine/twist_witness.json`, not a reduction matrix.

## 1. Definitions and conventions

The level is 60, weight is 3, and the character is
\[
 \chi(a)=\chi_{-3}(a)\chi_5(a)
\]
on units modulo 60. In particular \(\chi(-1)=-1\).

A Manin symbol \([P,(c,d)]\), with \(P=X\) or \(Y\), means
\[
 [P,g]=g\bigl(P\{0,\infty\}\bigr),\qquad
 g=\begin{pmatrix}a&b\\c&d\end{pmatrix}\in SL_2(\mathbb Z).
\]
The polynomial left action is
\[
 gP(X,Y)=P(dX-bY,-cX+aY).
\]
A primitive bottom row modulo 60 can always be lifted to a coprime integer row. The expression depends on its character-normalized coset, not the particular integer lift. Unit scaling is
\[
 [P,(tc,td)]=\chi(t)[P,(c,d)].
\]

For a determinant-one matrix \(h\), the right action is
\[
 [P,g]h=[h^{-1}P,gh].
\]
With
\[
 S=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad
 V=\begin{pmatrix}0&-1\\1&-1\end{pmatrix},
\]
we impose \(s+sS=0\) and \(s+sV+sV^2=0\). The \(-I\) relation is automatic. There are 144 projective rows, two monomials per row, and 576 such relation rows. The rational quotient has dimension 48.

All calculations can instead be read simply as reductions using valid path and character relations. For the displayed identity, exhibiting zero as a rational combination of those relations is sufficient; a numerical rank is not a substitute for that reduction.

For prime \(p\), the Hecke action is the sum of left actions of
\[
 \begin{pmatrix}1&a\\0&p\end{pmatrix},\quad 0\le a<p,
\]
plus \(\chi(p)\) times the left action of \(\left(\begin{smallmatrix}p&0\\0&1\end{smallmatrix}\right)\) when \(p\nmid60\). The latter is omitted when \(p\mid60\). For determinant \(p\), the polynomial uses the adjugate matrix, not its inverse divided by \(p\). This gives the classical coefficient action \(a_{pn}+\chi(p)p^2a_{n/p}\).

## 2. The identity

Set
\[
 w=X\{0,\infty\},\qquad
 w_5=\sum_{a=1}^4\chi_5(a)(X-aY/5)\{a/5,\infty\},
\]
\[
 U=U_3,\qquad p(U)=\frac{51-13U+3U^2-U^3}{30}.
\]
Then
\[
 \boxed{w_5-wp(U)=R(T_7^2+64).}
\]

The symbol \(R\) is \(1/13320000\) times the following integer combination. Every row means the displayed integer times \([P,(c,d)]\); no matrix-basis index is needed to interpret the table.

| P | c | d | Integer coefficient |
|---|---:|---:|---:|
| Y | 1 | 0 | 103500 |
| X | 5 | 3 | -210610 |
| Y | 5 | 4 | -81770 |
| X | 5 | 4 | 248975 |
| Y | 5 | 8 | 58992 |
| X | 5 | 8 | 245940 |
| Y | 5 | 9 | -111324 |
| X | 5 | 9 | -210610 |
| Y | 5 | 22 | -22778 |
| X | 5 | 22 | 3035 |
| Y | 6 | 1 | 177657 |
| X | 6 | 1 | -13758 |
| X | 6 | 5 | 44322 |
| Y | 6 | 13 | -44273 |
| X | 6 | 13 | -190296 |
| Y | 6 | 19 | -113641 |
| X | 6 | 19 | 248376 |
| Y | 10 | 3 | 111324 |
| Y | 10 | 11 | -202596 |
| Y | 12 | 1 | -62807 |
| X | 12 | 1 | -480264 |
| X | 12 | 5 | -2526 |
| Y | 12 | 7 | 102075 |
| X | 12 | 7 | 50190 |
| Y | 12 | 13 | -38059 |
| X | 12 | 13 | 312324 |
| Y | 12 | 19 | -184062 |
| X | 12 | 19 | 115224 |
| X | 15 | 1 | -107820 |
| Y | 15 | 2 | -22778 |
| X | 15 | 2 | -178095 |
| Y | 15 | 4 | 22778 |
| X | 15 | 4 | 285915 |

## 3. Two independent exact reductions

`derive_symbols.py` chooses lexicographically minimal projective rows, uses ordinary/floor continued fractions, and reduces the full relation matrix over the rationals. `build_witness.py` reconstructs the coefficient vector above and verifies the displayed identity exactly. `results/modsym_data.json` retains every representative, the complete reduction matrix, the Hecke matrices, the winding and twist vectors, and this witness.

`verify_symbols_independent.py` uses lexicographically maximal representatives, reverse sparse Gaussian elimination, and negative/ceiling continued fractions. It takes only the displayed raw witness data from the first representation. It does **not** use the primary reduction or Hecke matrices. It verifies all 288 path reconstructions, all 576 relation rows under both U3 and T7, and the complete identity in its different quotient.

The result is therefore an identity of rational modular symbols, not a numerical comparison of their integrals. The Hecke-compatible integration pairing and the exact eigenform identities are separate, explicitly stated inputs to its critical-value consequence.

## 4. Optional eight-dimensional diagnostic

`derive_projection.py` independently recomputes the kernel of T7²+64 from the full matrix and verifies that w5 and wp(U) have the same restriction. Its output `results/modsym_projection.json` includes an exact basis, the restricted 8×8 U3 matrix, and both restricted vectors. This is an additional exact check. The main certificate is the full quotient identity with R above, so no unspecified eigenprojection is required.
