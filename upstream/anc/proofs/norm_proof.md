# An exact Petersson-norm identity for the level-60 critical-value candidate

Jonas Matuzas. Proof developed 5 September 2026; audited for this note 15 September 2026.

This note proves an identity entirely on the modular side. It does **not** prove the Bessel-to-modular conjecture. There is no worldwide novelty or priority claim. The proof is explicit: modular functions, a finite Hecke witness, the inner-product convention, and all bad Euler factors are specified. The computational identities are independently replayed with exact rational arithmetic.

## File layout in this release

The readable witness table is in this file. `../engine/verify_norm.py` checks
it using only integer and rational arithmetic. The second construction is
`../primary/norm/rankin_projection.py`, with `eta_basis.json` alongside it;
other historical script names below are located under `../primary/norm/`.
`../verify_primary.py` runs those exact jobs in a temporary directory.
The paper's independent four-eta specification is identified with (1)–(2)
by `../engine/verify_modular_bridge.py`. Current numerical certification is
separate, in `../verification_notes.md`; the historical 73-digit numerical
engine is not required or included here.

## 1. Statement and conventions

Write \(r=\sqrt5>0\), \(q=e^{2\pi i\tau}\), and \(\eta_d=\eta(d\tau)\). Put
\[
\begin{aligned}
F={}&-8\eta_{12}^3\eta_{20}^3
+8\eta_4\eta_6^2\eta_{10}^2\eta_{60}
-4\eta_4^2\eta_6\eta_{10}\eta_{60}^2
-6\eta_3\eta_5\eta_{12}^2\eta_{20}^2\\
&+3\eta_3^2\eta_5^2\eta_{12}\eta_{20}
+\eta_3^3\eta_5^3
-3\eta_2\eta_3^2\eta_5^2\eta_{30}
+2\eta_2^2\eta_3\eta_5\eta_{30}^2.
\end{aligned}
\tag{1}
\]
For a Fourier series define \(U(\sum a_nq^n)=\sum a_{3n}q^n\). Define \(f=\sum_{j=0}^3c_jU^jF\), where
\[
\begin{aligned}
c_0&=5/12+r/20+i(1+19r/60),\\
c_1&=-7/36-13r/540+i(29/108+11r/90),\\
c_2&=-1/12-r/20+ir/60,\\
c_3&=-1/36-7r/540-i(1/108+r/90).
\end{aligned}
\tag{2}
\]
This is the form used in the previous candidate, externally labelled `60.3.b.a`; the proof below uses (1)–(2), not the database classification. Its first coefficients include
\[
a_1=1,\quad a_2=0,\quad a_3=r-2i,\quad a_5=r(1+2i),\quad a_7=-8i,\quad a_{11}=4ir.
\]
Use the unshifted classical \(L\)-function \(L(f,s)=\sum a_n n^{-s}\), continued from its half-plane of absolute convergence. Let \(\bar f(\tau)=\sum\bar a_nq^n\) denote coefficient conjugation, not pointwise conjugation.

Our Petersson convention is linear in the first variable:
\[
\langle f,h\rangle_{60}=\int_{\Gamma_0(60)\backslash\mathbb H}
f(\tau)\overline{h(\tau)}y^3\frac{dx\,dy}{y^2}.
\tag{3}
\]
**There is no division by the group index 144.**

### Theorem
\[
\boxed{\langle f,f\rangle_{60}=\frac{3(5-\sqrt5)}{2\pi^4}|L(f,2)|^2.}
\tag{4}
\]

## 2. Exact construction and finite Hecke control

Each eta monomial in (1) has exponent sum 6, both level-60 congruences modulo 24, character squareclass \(-15\), and positive order at every cusp. Thus \(F\in S_3(\Gamma_0(60),\chi_{-15})\). These properties are checked term by term with the eta cusp-order formula.

The index is 144 and the weight-three Sturm bound is 36. Exact Fourier computations prove
\[
(U^4-2U^2+81)F=0,
\qquad \dim\langle F,UF,U^2F,U^3F\rangle=4.
\tag{5}
\]
For every \(1\le n\le36\), the retained verifier computes \(T_nF\) in this cyclic span. Since \(T_n\) commutes with \(U\), the span is stable under these Hecke generators. Sturm duality shows that they span the complete Hecke algebra acting on this cusp character space. The four distinct roots of \(X^4-2X^2+81\) therefore give simultaneous Hecke eigenlines. The root \(r-2i\) with first coefficient 1 gives exactly (2).

The full-space computation used for the Fricke statement is also explicit. In the effective modular curve, \(\Gamma_0(60)\) is torsion-free: neither \(x^2+1\) nor \(x^2+x+1\) has a zero modulo 60. There are 12 cusps and the genus is
\[
1+144/12-12/2=7.
\]
Every cusp is singular for the character. For denominator \(c\mid60\) and width \(w=60/\gcd(c^2,60)\), the character of the parabolic stabilizer is \(\chi_{-15}(1+acw)=1\). Thus the character introduces no fractional cusp orders. The weight-three automorphic line bundle has degree 36. Since \(36>2\cdot7-2\), Riemann–Roch gives
\[
\dim M_3(\Gamma_0(60),\chi_{-15})=36+1-7=30,\qquad \dim S_3=18.
\tag{6}
\]

The 18 eta quotients in `eta_basis.json` are independently verified cuspidal and linearly independent. Add the 12 Eisenstein forms
\[
E_{3,\phi,\psi}(d\tau),\qquad d=1,2,4,
\]
for \((\phi,\psi)=(1,\chi_{-15}),(\chi_{-15},1),(\chi_{-3},\chi_5),(\chi_5,\chi_{-3})\). Their 30 coefficient vectors have rank 30 through the Sturm bound, hence form a basis of the full space by (6).

The exact \(T_7\)-matrix on this basis has characteristic polynomial
\[
X^6(X-48)^6(X+48)^6(X^2+34)^4(X^2+64)^2.
\]
More directly, exact row reduction gives
\[
\operatorname{rank}(T_7^2+64)=26,
\quad \ker(T_7^2+64)=\langle F,UF,U^2F,U^3F\rangle.
\tag{7}
\]
On that four-dimensional space,
\[
T_7=\frac29(7U+U^3),\qquad T_{11}=1-U^2.
\tag{8}
\]
The four roots \(\pm r\pm2i\) have distinct joint \((T_7,T_{11})\)-eigenvalues. In particular the joint eigenspace with \((-8i,4ir)\) is one-dimensional in the **full** modular space. None of these conclusions is inferred from finite Frobenius agreement with a Bessel motive.

## 3. A weight-one Eisenstein input and a regular weight-two Eisenstein value

Here \(\chi=\chi_{-15}\) is the primitive character of conductor 15. In particular \(\chi(2)=1\); this is distinct from the induced level-60 Hecke character, which is zero at 2.

Define
\[
g=E_{1,1,\chi}=1+\sum_{n\ge1}b_nq^n,
\qquad b_n=\sum_{d\mid n}\chi(d).
\tag{9}
\]
The constant is \(L(\chi,0)/2=1\), verified by the Bernoulli sum. General Eisenstein modularity places \(g\) in \(M_1(\Gamma_0(15),\chi)\). An independent useful description is
\[
g(\tau)=\frac12\sum_{m,n\in\mathbb Z}
\left(q^{m^2+mn+4n^2}+q^{2m^2+mn+2n^2}\right).
\tag{10}
\]
Both positive binary theta series have level 15 and character \(\chi\). Equality follows by the weight-one level-15 Sturm bound 2 (the verifier checks through 40).

Let \(E_2=1-24\sum\sigma_1(n)q^n\) and \(E_2^*=E_2-3/(\pi y)\). The weight-two Eisenstein series at infinity for \(\Gamma_0(60)\), with constant term 1, has regular value at its analytic parameter zero
\[
\begin{aligned}
\mathcal E_{2,60}^{\infty}(\tau,0)
&=\frac{25}{16}\sum_{d\mid30}\frac{\mu(d)}{d^2}E_2^*((60/d)\tau)\\
&=E_{\rm hol}(\tau)-\frac1{48\pi y},
\end{aligned}
\tag{11}
\]
where \(E_{\rm hol}\) is the same finite combination of the holomorphic \(E_2\)'s. One derivation sums primitive bottom rows \((c,d)\), \(60\mid c\), by Möbius inversion. Removing the prime factors 2,3,5 gives \(\prod_{p\mid60}(1-p^{-2})^{-1}=25/16\). The nonholomorphic coefficient follows from
\[
\frac{25}{16}\sum_{d\mid30}\frac{\mu(d)}{d^2}=1,
\qquad
\frac{25}{16}\sum_{d\mid30}\frac{\mu(d)}{60d}=\frac1{144}.
\]
The initial series for the unfolding is \(\sum_{\Gamma_\infty\backslash\Gamma_0(60)}\Im(\gamma\tau)^s j(\gamma,\tau)^{-2}\), continued to \(s=0\). Its regular value, not an unregularized conditionally summed lattice series, is used in (11).

## 4. Unfolding and the raising term

Put \(\theta=q\,d/dq\) and
\[
H=gE_{\rm hol}-\frac1{12}\theta g.
\tag{12}
\]
The weight-one Maass raising operator is \(\delta_1g=\theta g-g/(4\pi y)\). Equation (11) gives
\[
g\mathcal E_{2,60}^{\infty}=H+\frac1{12}\delta_1g.
\]
Both terms other than \(H\) are weight three, so \(H\) is holomorphic modular of weight three and character \(\chi\). Its holomorphy at the cusps follows as usual from the cancellation in this expression; equivalently it follows by applying the same raising-operator identity in each cusp coordinate.

The raising term is orthogonal to \(f\). Explicitly,
\[
f\overline{\delta_1g}\,y
=-\frac1{2\pi i}\partial_{\bar\tau}(f\bar g\,y).
\tag{13}
\]
The one-form \(f\bar g\,y\,d\tau\) is invariant. Paired sides of a fundamental domain cancel, and the cusp boundaries vanish because \(f\) is cuspidal and \(g\) is holomorphic at the cusps. Stokes gives \(\langle f,\delta_1g\rangle=0\).

Unfold first where the analytic parameter has large real part, then continue to zero. The Mellin integral \(\int_0^\infty y e^{-4\pi ny}dy=(16\pi^2n^2)^{-1}\) gives
\[
\boxed{\langle f,H\rangle=\frac1{16\pi^2}D(2),\qquad
D(s)=\sum_{n\ge1}a_nb_n n^{-s}.}
\tag{14}
\]
Here \(D(2)\) means analytic continuation of the Rankin series. No assertion of absolute convergence of that raw series at 2 is needed.

## 5. All Euler factors, including 2,3,5

For \(p\nmid30\), write \(\varepsilon=\chi(p)\). Hecke recurrence and the geometric sum defining \(b_{p^j}\) imply
\[
\sum_{j\ge0}a_{p^j}b_{p^j}T^j
=\frac{1-p^2T^2}{(1-a_pT+\varepsilon p^2T^2)
(1-\varepsilon a_pT+\varepsilon p^2T^2)}.
\tag{15}
\]
The identity is verified separately for \(\varepsilon=1\) and \(-1\). Adjointness of the good Hecke operators gives \(\bar a_p=\chi(p)a_p\).

At 2, \(a_{2^j}=0\) for \(j>0\), so the local factor of \(D\) is 1. At 3 and 5, \(b_{p^j}=1\), and the factors are \((1-a_pp^{-s})^{-1}\). Consequently
\[
D(s)=\frac{L(f,s)L(\bar f,s)(1-\bar a_3 3^{-s})(1-\bar a_5 5^{-s})}
{\zeta^{(30)}(2s-2)},
\tag{16}
\]
where \(\zeta^{(30)}\) omits the three primes. In particular
\[
\zeta^{(30)}(2)=\frac{8\pi^2}{75},
\quad
\beta=\left(1-\frac{r+2i}{9}\right)
\left(1-\frac{r(1-2i)}{25}\right),
\]
so (14) becomes
\[
\boxed{\langle f,H\rangle=\frac{75\beta}{128\pi^4}|L(f,2)|^2.}
\tag{17}
\]

## 6. A finite, exact Hecke witness

Define
\[
H_f=\frac5{48}F-\frac{35}{432}UF-\frac5{144}U^2F-\frac5{432}U^3F.
\tag{18}
\]
The explicit witness in the appendix satisfies
\[
\boxed{H-H_f=(T_7^2+64)R.}
\tag{19}
\]
Every function is proved holomorphic modular. Equality of their coefficients through \(q^{36}\) proves (19) by Sturm; the independent verifier checks \(q^0\) through \(q^{60}\).

For clarity, this is not a statement that finite Fourier agreement proves an identity of arbitrary functions. It proves (19) because all functions have first been placed in the same finite-dimensional modular space.

Since \(\chi(7)=-1\), \(T_7^*=-T_7\) and \(T_7^2+64\) is self-adjoint. It annihilates \(f\). Thus \(\langle f,H\rangle=\langle f,H_f\rangle\).

Put \(\lambda=r-2i\). Let \(h_k\) be the coefficient of \(q\) in \(U^kH_f\). Exact series arithmetic gives
\[
(h_0,h_1,h_2,h_3)=(25/72,5/18,5/72,5/9).
\]
The \(\lambda\)-spectral projector for (5) gives the coefficient of normalized \(f\) in \(H_f\):
\[
\delta=
\frac{h_3+\lambda h_2+(\lambda^2-2)h_1+(\lambda^3-2\lambda)h_0}
{4\lambda^3-4\lambda}
=\frac{25+2r+i(5-r)}{288}.
\tag{20}
\]
The other three conjugate eigenforms are orthogonal to \(f\): their joint good-Hecke eigenvalues in (8) differ, and those operators are normal. Therefore
\[
\langle f,H\rangle=\bar\delta\langle f,f\rangle.
\]
Finally, exact arithmetic in \(\mathbb Q(\sqrt5,i)\) proves
\[
\frac{75\beta}{128\bar\delta}=\frac{3(5-r)}2.
\]
Together with (17), this proves (4).

No Bessel integral, numerical L-value, period-fitting coefficient, or conjectural arithmetic correspondence enters this proof.

## 7. The complete Fricke normalization without a database assumption

For \(W_{60}=\left(\begin{smallmatrix}0&-1\\60&0\end{smallmatrix}\right)\), use the usual determinant-normalized weight-three slash operation. Fricke conjugation changes good-Hecke eigenvalues from \(a_p\) to \(\chi(p)a_p=\bar a_p\). The full-space one-dimensional joint eigenspace result (7)–(8) proves
\[
f|_3W_{60}=\lambda_W\bar f.
\]
Thus proportionality does not require an external newform classification. The retained Dedekind-eta transformation calculation finds the first Fourier coefficients of \(U^jF|_3W_{60}\), for \(j=0,1,2,3\), to be
\[
-i,\quad0,\quad3i,\quad-36i.
\]
Substitution of (2) gives
\[
\lambda_W=\frac{10-2r+i(5+4r)}{15}=i\epsilon,
\quad
\boxed{\epsilon=\frac{5+4r+i(2r-10)}{15},\quad |\epsilon|=1.}
\tag{21}
\]
The exact transformation program works in \(\mathbb Q[\xi]/(\xi^{216}-\xi^{108}+1)\), \(\xi=e^{2\pi i/648}\), using rational Dedekind sums. It has been rerun in this package.

On the imaginary axis,
\[
f(i/(60y))=\epsilon60^{3/2}y^3\bar f(iy).
\]
Splitting the Mellin integral at \(1/\sqrt{60}\) consequently gives, with \(h=2\pi/\sqrt{60}\),
\[
\boxed{
L(f,2)=\sum_{n\ge1}\frac{a_ne^{-hn}(1+hn)}{n^2}
+\epsilon h\sum_{n\ge1}\frac{\bar a_ne^{-hn}}n.
}
\tag{22}
\]
Both series converge exponentially. This is the formula enclosed rigorously in the current modular interval engine.

## Appendix: the complete rational witness R

Order the eta divisors as
\[
\mathcal D=(1,2,3,4,5,6,10,12,15,20,30,60),
\quad E_j=\prod_{d\in\mathcal D}\eta(d\tau)^{e_{j,d}}.
\]
The following table supplies every nonzero eta term of \(R\).

| j | coefficient in R | exponent row in the order \(\mathcal D\) |
|---:|---:|---|
| 0 | 47/32256 | 3,0,0,0,0,0,0,0,3,0,0,0 |
| 1 | -9883/48384 | 2,0,0,1,0,0,0,0,2,0,0,1 |
| 2 | 11779/48384 | 2,0,0,0,0,1,1,0,2,0,0,0 |
| 3 | -1367/3024 | 1,0,0,2,0,0,0,0,1,0,0,2 |
| 4 | 18059/24192 | 1,0,0,1,0,1,1,0,1,0,0,1 |
| 5 | -505/1344 | 1,0,0,0,0,2,2,0,1,0,0,0 |
| 6 | -329/1728 | 0,2,0,0,0,0,0,1,0,1,2,0 |
| 7 | -137/288 | 0,1,0,0,0,0,0,2,0,2,1,0 |
| 8 | 47/1792 | 0,0,3,0,3,0,0,0,0,0,0,0 |
| 9 | 83/672 | 0,0,2,0,2,0,0,1,0,1,0,0 |
| 12 | -2/315 | -1,3,3,0,0,0,1,0,0,0,0,0 |
| 13 | -61/1890 | -1,2,3,0,0,1,0,0,0,0,1,0 |
| 14 | -1/1260 | -1,1,4,0,1,1,0,0,0,0,0,0 |
| 15 | 2/21 | -1,1,3,0,0,1,0,1,0,1,0,0 |

To their sum add
\[
\frac1{132608}E_{3,1,\chi}(2\tau)-\frac1{16576}E_{3,1,\chi}(4\tau),
\quad
E_{3,1,\chi}=-8+\sum_{n\ge1}\left(\sum_{d\mid n}\chi(d)d^2\right)q^n.
\]
This completely specifies \(R\), including its constant term \(1/2368\).

For the 18-dimensional cusp basis used in Section 2, add the four eta rows
\[
\begin{aligned}
e_{10}&=(0,0,0,3,0,0,0,0,0,0,0,3),\\
e_{11}&=(0,0,0,0,0,0,0,3,0,3,0,0),\\
e_{16}&=(-1,0,3,1,0,0,0,0,0,3,0,0),\\
e_{17}&=(2,-1,0,0,0,2,0,0,2,0,1,0).
\end{aligned}
\]

## References and proof dependencies

- W. Stein, *Modular Forms: A Computational Approach*, Eisenstein-series chapter, Theorems 5.8–5.10, for the standard Eisenstein modularity and basis statements; dimension chapter for the standard dimension framework. The dimension in (6) is also derived above by Riemann–Roch.
- A. Agashe and W. Stein, *Generating the Hecke algebra*, for the finite Sturm/Hecke method. In the character space, the same characteristic-zero Sturm duality applies: a functional killing all the first Sturm-bounded coefficients corresponds to a zero modular form.
- Classical eta-quotient transformation and cusp-order criteria, Petersson Hecke adjointness, and the analytic continuation of weight-two Eisenstein series. Their normalizations needed here are printed explicitly.
- Retained preceding package: exact eta/Hecke definition and Dedekind-eta Fricke coefficient, freshly replayed. The preceding appeal to a newform classification for Fricke proportionality is replaced by Section 2 of this note.

The calculations of the finite witness are computational proof components, not high-precision numerical tests. The general analytic and geometric arguments are given in the text; passing code does not replace them. The theorem concerns \(f\) alone. The Bessel-to-norm equality remains a separate conjecture.
