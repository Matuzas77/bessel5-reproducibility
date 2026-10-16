# Fifth Bessel moment: an exact critical-value twist relation and a two-period reduction

Jonas Matuzas. Proof developed 6 September 2026; audited for this note 15 September 2026.

**Target:**
\[
 A=\int_0^\infty K_0(t)^5\,dt.
\]

**Outcome.** This report proves an exact relation between critical values of the fixed modular form and one of its coefficient conjugates. An explicit rational modular-symbol certificate fixes the factor, including the Euler factor at 5. Consequently the four real eta-product critical values used to define the candidate are expressed in two positive real periods. The proposed fifth-Bessel value reduces to a specified real quadratic expression in those two periods. **The equality of that expression with the Bessel integral is not proved.** The new proof is entirely on the modular side; no Bessel value or critical-value decimal enters it. No worldwide novelty or priority claim is made.

## File layout in this release

`MODULAR_SYMBOL_CERTIFICATE.md` contains the raw 33-term rational witness.
The current standard-library verifier `../engine/verify_twist_symbols.py`
reads only that raw witness, not previously reduced matrices or period values.
The second exact reconstruction is in `../primary/twist/` and is replayed
by `../verify_primary.py`. Script names in the proof refer to that directory.

## 1. The fixed form and the unchanged conjecture

Let \(r=\sqrt5>0\), \(q=e^{2\pi i\tau}\), and \(\eta_j=\eta(j\tau)\). Put
\[
 h=\eta_2^2\eta_{10}^2,\qquad
 g=\frac{\eta_3^2\eta_5^2}{\eta_1\eta_{15}},\qquad
 d=\frac{\eta_1^2\eta_{15}^2}{\eta_3\eta_5},
\]
\[
 \omega=\frac{1-2i}{r},\qquad
 \epsilon=\frac{5+4r+i(2r-10)}{15}.
\]
The four rational-coefficient cusp forms used here are
\[
\begin{array}{ll}
 P_1=h(\tau)g(2\tau)=\eta_2\eta_6^2\eta_{10}^4/\eta_{30},&
 P_2=h(\tau)d(2\tau)=\eta_2^4\eta_{10}\eta_{30}^2/\eta_6,\\
 P_3=h(3\tau)g(2\tau)=\eta_6^4\eta_{10}^2\eta_{30}/\eta_2,&
 P_4=h(3\tau)d(2\tau)=\eta_2^2\eta_6\eta_{30}^4/\eta_{10}.
\end{array}                                                     \tag{1}
\]
Define, exactly,
\[
 \boxed{f=P_1+\omega P_2+3\epsilon P_3+3\epsilon\bar\omega P_4.} \tag{2}
\]
Thus \(3\epsilon\bar\omega=r+2i\), and
\[
 f=q+(r-2i)q^3+r(1+2i)q^5-8iq^7+\cdots .
\]
The classical, unshifted normalization is
\[
 L(F,s)=\sum_{n\ge1}a_n(F)n^{-s},
 \qquad L(F,2)=4\pi^2\int_0^\infty yF(iy)\,dy
\]
for a weight-three cusp form, with the Dirichlet series continued when necessary. Every Mellin integral used here converges at both endpoints.

The unchanged candidate is
\[
 \boxed{A_* =\frac\pi4\left[(125-13r)\Re L(f,2)^2
 +(41r-125)\Im L(f,2)^2\right].}                              \tag{3}
\]
Real and imaginary parts are taken **after squaring**. The assertion \(A=A_*\) remains conjectural.

## 2. Main exact result

Let \(\sigma\) be the automorphism of \(\mathbb Q(r,i)\) specified by
\[
 \sigma(r)=-r,\qquad \sigma(i)=-i.
\]
It acts **on the Fourier coefficients**: \(f^\sigma=\sum\sigma(a_n(f))q^n\). It is not an operation being applied to an arbitrary transcendental number.

### Theorem 1
\[
 \boxed{L(f^\sigma,2)=\kappa L(f,2),\qquad
 \kappa=\frac{(r-1)(4+3i)}{10}.}                             \tag{4}
\]
The factor in (4) is proved using finite modular and modular-symbol identities, not selected by numerical recognition. Sections 3–5 provide the complete proof and its executable rational witness.

Set
\[
 l_j=L(P_j,2)=4\pi^2\int_0^\infty yP_j(iy)\,dy.
\]
All four \(l_j\) are positive real numbers: the eta functions on the positive imaginary axis are positive, and all four forms in (1) are cuspidal.

### Corollary 2: two exact real linear relations
\[
 \boxed{
 l_2=\frac{(3+r)l_1-9(1+r)l_3}{4},\qquad
 l_4=\frac{(1-r)l_1+(5+r)l_3}{4}.
 }                                                           \tag{5}
\]
Indeed substitute (2) in (4), use the realness of the \(l_j\), and separate real and imaginary parts. The coefficient determinant for \((l_2,l_4)\) is \(4-4r/5\ne0\). Solving gives (5). Both exact implementations verify this algebra independently.

In particular,
\[
 \boxed{\frac{9(r-1)}2<\frac{l_1}{l_3}<\frac{5+3r}{2}.}       \tag{6}
\]
These strict inequalities follow from \(l_2,l_4>0\). They are bounds for modular periods, not Bessel-period comparisons. No assertion of algebraic linear independence of \(l_1,l_3\) is made.

## 3. Finite modular identities, with the prime-5 factor retained

Let \(\chi=\chi_{-15}\), and let \(\chi_5(n)=(n/5)\), extended by zero on multiples of 5.

### 3.1 Modularity and the required eigenvalues

Newman's eta quotient criterion, including its Nebentypus formula, gives
\[
 P_j\in S_3(\Gamma_0(60),\chi).
\]
For each form the exponent sum is 6, both eta congruences hold modulo 24, and the character squareclass is \(-15\). In the cusp-order formula
\[
 \operatorname{ord}_{a/c}\prod_{\delta\mid60}\eta(\delta\tau)^{r_\delta}
 =\frac{60}{24\gcd(c^2,60)}
 \sum_{\delta\mid60}\frac{r_\delta\gcd(c,\delta)^2}{\delta},
\]
all values are positive. The products determining the character are respectively \(24000,24000,1944000,1944000\), namely \(15\) times the squares of \(40,40,360,360\). The complete table, indexed by \(c\mid60\), is regenerated as `results/modular_coefficients.json` in the temporary twist workspace.

The level-60 index is 144, hence the weight-three characteristic-zero Sturm bound is 36. Exact coefficient calculations through this bound prove
\[
 U_3f=(r-2i)f,\quad U_5f=r(1+2i)f,\quad T_7f=-8if.          \tag{7}
\]
The Hecke normalization is
\[
 [q^n]T_7F=a_{7n}(F)+\chi(7)7^2a_{n/7}(F),
\]
where the second term is zero unless \(7\mid n\). Coefficients through 252 suffice. These are finite proofs only because modularity of every side is first established.

### 3.2 An exact depleted twist, not an undepleted equality

The identity needed at the ramified prime is
\[
 \boxed{f\otimes\chi_5=f^\sigma-a_5(f^\sigma)f^\sigma(5\tau),
 \qquad a_5(f^\sigma)=-r+2ir.}                             \tag{8}
\]
The rightmost term is essential. The coefficient twist on the left is zero at every index divisible by 5; the right side has exactly this depletion.

Here is a direct level proof, so (8) does not depend on a guessed newform classification. The Gauss sum of \(\chi_5\), with the usual positive choice, is \(r\), so
\[
 f\otimes\chi_5=\frac1r\sum_{a=1}^4\chi_5(a)f(\tau+a/5).   \tag{9}
\]
For \(\gamma=\left(\begin{smallmatrix}A&B\\C&D\end{smallmatrix}\right)
\in\Gamma_0(300)\), choose \(b\equiv aD^2\pmod5\). With
\(\delta_a=\left(\begin{smallmatrix}1&a/5\\0&1\end{smallmatrix}\right)\),
\[
 \delta_a\gamma\delta_b^{-1}=
 \begin{pmatrix}
 A+aC/5&B+(aD-bA)/5-abC/25\\
 C&D-bC/5
 \end{pmatrix}\in\Gamma_0(60).
\]
Its lower-right entry has the same character value as \(D\), and
\(\chi_5(b)=\chi_5(a)\). This proves that (9) transforms at level 300 with character \(\chi\). Rational translates of a cusp form vanish at rational cusps in their local coordinates; hence (9) is cuspidal. Both terms on the right of (8) are cusp forms at level 300 as well.

The index at level 300 is 720, giving Sturm bound 180. The two independently reconstructed coefficient lists verify (8) exactly through that bound, proving it globally.

Taking critical values in (8) gives
\[
 \boxed{L(f\otimes\chi_5,2)=
 \left(1-\frac{a_5(f^\sigma)}{25}\right)L(f^\sigma,2).}       \tag{10}
\]
This follows from the Mellin scaling for \(f^\sigma(5\tau)\); no Euler-product hypothesis is necessary.

## 4. The explicit rational modular-symbol certificate

Use weight-three modular symbols with character \(\chi\). A polynomial symbol is paired by
\[
 \langle F,P\{\alpha,\beta\}\rangle
 =\int_\alpha^\beta P(\tau,1)F(\tau)\,d\tau.
\]
For \(g=\left(\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\right)\), the polynomial action is
\[
 gP(X,Y)=P(dX-bY,-cX+aY).
\]
For determinant different from one this is the adjugate action; it is the convention compatible with the classical Hecke coefficients in (7).

Let
\[
 w=X\{0,\infty\},\qquad
 w_5=\sum_{a=1}^4\chi_5(a)(X-aY/5)\{a/5,\infty\}.
\]
Let \(U=U_3\) and put
\[
 \boxed{p(U)=\frac{51-13U+3U^2-U^3}{30}.}                 \tag{11}
\]

### Proposition 3: an equality of modular symbols

There is an explicitly printed rational symbol \(R\) such that
\[
 \boxed{w_5=w\,p(U)+R(T_7^2+64).}                        \tag{12}
\]
The common denominator of \(R\) is 13,320,000. Its 33 nonzero Manin coefficients, including their bottom-row representatives and polynomial labels, are printed in `MODULAR_SYMBOL_CERTIFICATE.md`. That file also specifies the full finite presentation and action conventions. Thus (12) is a checkable finite certificate, not merely a dimension assertion or an eigenvalue fit.

For clarity, the presentation has 144 primitive projective bottom rows modulo60 and two monomials \(X,Y\), giving 288 generators. Scalar bottom-row changes satisfy
\([P,(tc,td)]=\chi(t)[P,(c,d)]\). For
\[
 S=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad
 V=\begin{pmatrix}0&-1\\1&-1\end{pmatrix},
\]
one imposes \(s+sS=0\), \(s+sV+sV^2=0\). The \(-I\) relation is automatic because the polynomial degree and character are both odd. The 576 relation rows reduce to a 48-dimensional quotient. Exact continued-fraction path reduction and rational elimination prove (12).

The independent verifier uses different projective representatives, the reverse elimination order, and negative/ceiling continued fractions instead of ordinary/floor continued fractions. It reads only the raw coefficient list of \(R\), not the primary reduction or Hecke matrices. It verifies (12) in its independently constructed quotient and checks the descent of \(U_3,T_7\) on the complete576-row relation set.

Hecke compatibility of the integration pairing gives
\[
 \langle f,w_5\rangle=p(r-2i)\langle f,w\rangle,
\]
since \((T_7^2+64)f=0\). Direct evaluation of the polynomial is
\[
 \boxed{p(r-2i)=\frac{9-r+i(8-2r)}5.}                     \tag{13}
\]
Neither a modular critical value nor a Bessel integral is involved in obtaining this coefficient.

## 5. Deduce the exact critical-value relation

The normalization of the ordinary winding symbol is
\[
 \langle f,w\rangle=-\frac{L(f,2)}{4\pi^2}.
\]
For \(w_5\), translate each integration path by \(a/5\), and apply the finite Gauss-sum identity
\[
 \sum_{a=1}^4\chi_5(a)e^{2\pi ina/5}=r\chi_5(n).
\]
Cusp decay justifies the translated Mellin integrals and gives
\[
 \langle f,w_5\rangle=-\frac{rL(f\otimes\chi_5,2)}{4\pi^2}.
\]
Therefore (12)–(13) imply
\[
 rL(f\otimes\chi_5,2)=\frac{9-r+i(8-2r)}5 L(f,2).         \tag{14}
\]
Combining (10) and (14) yields
\[
 \frac{p(r-2i)}{r(1-(-r+2ir)/25)}
 =\frac{(r-1)(4+3i)}{10},
\]
which proves Theorem 1.

The \(\sigma\) in the theorem denotes coefficient conjugation of a cusp form. We have not assumed that an extension of \(\sigma\) to all complex numbers commutes with the Mellin integral. In particular, one should not apply \(\sigma\) formally to both sides of a transcendental period identity. All statements here follow from separately defined forms, fixed complex Gauss sums, and the integration pairing.

## 6. Exact simplification of the proposed Bessel value

From (5),
\[
 \boxed{L(f,2)=\frac{2r}{5}(1-2i)l_1+
 \left[\frac{8r}{5}+i\left(5+\frac{9r}{5}\right)\right]l_3.} \tag{15}
\]
Substitution in (3), entirely on the modular side, gives
\[
 \boxed{A_*=25\pi\left[(1-r)l_1^2+36l_1l_3
 -27(1+r)l_3^2\right].}                                   \tag{16}
\]
This eliminates two eta critical values and every complex part from the candidate. It does not prove that this candidate is the Bessel integral.

There is also a rational quadratic identity among the four eta periods:
\[
\begin{split}
0={}&-13l_1^2+150l_1l_2-90l_1l_3+414l_1l_4-25l_2^2\\
&-90l_2l_3-90l_2l_4+675l_3^2+270l_3l_4-513l_4^2.
\end{split}                                                 \tag{17}
\]
It follows by substituting (5), and is not an identity guessed from rounded period values.

The earlier companion *candidates* simplify to
\[
 B_*=30\sqrt3\left[(r-1)l_1^2-8l_1l_3+9(1+r)l_3^2\right],
\]
\[
 C_*=\frac{30}{\pi}\left[(1-r)l_1^2+18l_1l_3-9(1+r)l_3^2\right]. \tag{18}
\]
Their quadratic conic relation holds identically under (18), as checked independently. The corresponding equalities with the mixed Bessel moments remain unproved; consistency on a conic does not identify its physical point.

## 7. What was computed, and what the computation proves

The optional `../verify_primary.py` performs these seven exact twist jobs, in addition to six norm jobs:

1. Construct the exact modular-symbol quotient and Hecke actions from the 576 Manin relations, verifying the actions on every raw generator and exact commutation.
2. Construct the rational witness in (12).
3. Verify its reduction on the eight-dimensional \(T_7^2+64\)-kernel; this is an additional check, not a replacement for the full quotient identity.
4. Reconstruct and verify the witness independently using rational arithmetic and different continued fractions.
5. Compute the eta coefficients via logarithmic derivative recurrence, verify all cusp orders and the exact eigenform and level 300 twist identities.
6. Independently reconstruct the coefficients by pentagonal products and direct binary theta enumeration, and recheck those identities and the two-period algebra.
7. Verify the exact algebra in (5), (15)–(18), including the nonzero linear-system determinant.

The coefficient proof uses q-series only after proving that both sides are modular forms in the same space, with the stated bound 36 or 180. The modular-symbol proof uses exact relations and path operations, not numerical integration. The general modularity and integration-compatibility theorems remain logical inputs; a count of successful software assertions is not substituted for them.

Numerical evidence is supplied separately by the current interval verifier. No historical quadrature is needed by the proof.

## 8. The still-unproved equality

The exact remaining statement is now
\[
 \boxed{
 \int_0^\infty K_0(t)^5\,dt
 \stackrel?=25\pi\left[(1-\sqrt5)l_1^2+36l_1l_3
 -27(1+\sqrt5)l_3^2\right],
 }                                                           \tag{19}
\]
where the two periods have the independent convergent definitions
\[
 l_1=4\pi^2\int_0^\infty y\,
 \frac{\eta(2iy)\eta(6iy)^2\eta(10iy)^4}{\eta(30iy)}\,dy,
\]
\[
 l_3=4\pi^2\int_0^\infty y\,
 \frac{\eta(6iy)^4\eta(10iy)^2\eta(30iy)}{\eta(2iy)}\,dy.
\]
The new theorem proves equivalence of this two-period target with the original specified candidate; it does not prove either equality with the integral.

No direct integral transformation, independent modular solution of the Bessel recurrence, or comparison of the relevant forms and integration cycles establishing (19) was obtained in this attempt. The inspected primary literature supplies the modular-symbol, twist, and finite-coefficient tools used above, not this final specialized identity. The completed result is an exact relation on the modular side and an explicit reduction of its period coordinates. **The finite arithmetic evaluation of the original fifth Bessel integral remains unproved.**

See `standard_sources.md` for standard inputs and the scope of the literature check. No private correspondence, external article PDF, or font is redistributed.
