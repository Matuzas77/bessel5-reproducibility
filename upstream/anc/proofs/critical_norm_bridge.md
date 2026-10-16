# Proof of the twisted symmetric-square critical-norm identity

**Supporting proof for version 4, 16 September 2026.**

This is the full analytic proof of Proposition 3.2 and equation (3.5) of
*A quadratic critical-value conjecture for the fifth Bessel moment*.
It uses no Bessel integral or numerical L-value. The existing norm theorem
is proved separately in `norm_proof.md`; the new argument evaluates the
twisted symmetric-square L-function in terms of that norm.

The conclusion is
\[
 L(\chi_{-4}\operatorname{Sym}^2f,3)=\frac{\pi^4}{8}\langle f,f\rangle_{60},
 \qquad
 L(\chi_{-4}\operatorname{Sym}^2f,2)=\sqrt{15}\pi^2\langle f,f\rangle_{60}.
\]
The measure in the Petersson norm is not divided by the index 144.
Neither the pure nor the mixed Bessel-to-modular conjecture is proved here.
The finite calculations below are reproducible; the analytic arguments are
written mathematics, not a formal proof-assistant certificate.

## 1. Fixed modular objects and conventions

Write epsilon=chi_{-15}, chi=chi_{-4}, and eta_d=eta(d tau), with q=exp(2 pi i tau). Put
\[
 P_1=\eta_2\eta_6^2\eta_{10}^4/\eta_{30},\quad
 P_2=\eta_2^4\eta_{10}\eta_{30}^2/\eta_6,
\]
\[
 P_3=\eta_6^4\eta_{10}^2\eta_{30}/\eta_2,\quad
 P_4=\eta_2^2\eta_6\eta_{30}^4/\eta_{10}.
\]
For r=sqrt5>0, let
\[
 w=(1-2i)/r,\qquad e=(5+4r+i(2r-10))/15,
 \qquad f=P_1+wP_2+3eP_3+3e\bar wP_4.
\]
Their modularity, the identification of this embedding, and the norm theorem are inherited from the retained sources. The current scripts recheck the eta modularity criteria and the good Hecke identities used here; they do not infer identification with an entire newform orbit from a finite numerical prefix.

Set
\[
 \mathscr L(s)=L(\chi\operatorname{Sym}^2 f,s),\qquad
 D_f(s)=\sum_{m\ge1}\frac{\chi(m)a_f(m^2)}{m^s}
\]
on their initial convergence half-planes. The good local factors imply
\[
 \zeta^{(30)}(2s-4)D_f(s)=R(s)\mathscr L(s),
 \tag{1.1}
\]
where
\[
 R(s)=(1+(1+4ir)3^{-s})(1+(15+20i)5^{-s}).
 \tag{1.2}
\]
Indeed, if A+B=a_f(p), AB=epsilon(p)p^2, then
\[
 \sum_{j\ge0}a_f(p^{2j})X^j=
 \frac{1+ABX}{(1-A^2X)(1-B^2X)}.
\]
Substitute X=chi(p)p^{-s}; relative to the twisted symmetric-square local factor the numerator is 1-p^{4-2s}. At 2 both the depleted square-index series and the adopted symmetric-square factor are 1. At 3 and 5 the square-index series retains just the root chi(p)a_f(p)^2; the other root supplies (1.2). The adopted full inverse factors are
\[
 1+2T+81T^2\quad(p=3),\qquad
 1+30T+625T^2\quad(p=5).
\]
These local factors and the conductor-960 completion are modular-side inputs from the new arithmetic theorem and its normalization, not consequences of the conjectural Bessel formula.

At s=3,
\[
 \zeta^{(30)}(2)=\frac{8\pi^2}{75},\qquad
 R(3)=\frac{16(7+ir)(7+i)}{675}.
 \tag{1.3}
\]

We use a Petersson pairing linear in its first slot:
\[
 \langle g,h\rangle_{60}=\int_{\Gamma_0(60)\backslash\mathbb H}
 g(z)\overline{h(z)}y^3\,\frac{dx\,dy}{y^2}.
\]
For p not dividing 60,
\[
 (T_pg)_n=g_{pn}+\epsilon(p)p^2g_{n/p},
 \quad T_p^*=\overline{\epsilon(p)}T_p.
\]
Coefficients with noninteger indices are zero.

## 2. The critical square-index kernel

Let P_m be the holomorphic weight-three Poincare series with seed exp(2 pi i m z), level 60 and character epsilon, using Gamma_infinity including -I. Then
\[
 \langle f,P_m\rangle_{60}=\frac{a_f(m)}{16\pi^2m^2}.
 \tag{2.1}
\]
This follows by unfolding and the integral of y exp(-4 pi m y). It fixes all index and scalar conventions.

For a complex parameter sigma initially with Re sigma>1/4, define
\[
 \mathcal K_\sigma=\sum_{m\ge1}\chi(m)m^{1-2\sigma}P_{m^2}.
 \tag{2.2}
\]
The sum converges in the finite-dimensional cusp space in this half-plane. One justification, requiring no target period, is as follows. For a fixed positive Fourier index n, the Poincare expansion and the classical fixed-character Kloosterman bound give
\[
 [q^n]P_{m^2}=O_{n,\varepsilon}(m^{-3/2+\varepsilon})
\]
apart from its one possible seed term. Split the c sum at c=m, use
J_2(v)=O(min(v^2,v^{-1/2})), and
\(|S_\epsilon(m^2,n;c)|\ll_{n,\varepsilon}c^{1/2+\varepsilon}\).
The latter is the usual Weil bound at good primes and the elementary quadratic-character/prime-power bound at the fixed bad primes; the fixed n absorbs its gcd factor. The first 36 Fourier coefficients inject the weight-three space by Sturm. Norm equivalence in a finite-dimensional space therefore proves (2.2)'s convergence and holomorphic dependence. For real sigma in that domain, unfolding gives
\[
 \langle f,\mathcal K_\sigma\rangle_{60}
 =\frac{D_f(3+2\sigma)}{16\pi^2}.
 \tag{2.3}
\]
For complex sigma the conjugate parameter is placed in the second slot; we use real sigma for (2.3).

The required critical kernel is the analytic continuation to sigma=0, **not an assertion that the raw m sum converges absolutely there**. Sections 3–5 prove that continuation and compute it. This is the main analytic issue in the proof.

## 3. Fourier analysis and the critical support cutoff

Define
\[
 S_\epsilon(u,n;c)=\sum_{d\pmod c}^{*}\epsilon(d)
 e((u\bar d+nd)/c),\quad e(z)=e^{2\pi iz}.
\]
In our weight-three convention the Poincare Fourier expansion is
\[
 [q^n]P_{m^2}=1_{n=m^2}+
 \frac{2\pi i n}{m^2}\sum_{60\mid c}\frac{S_\epsilon(m^2,n;c)}c
 J_2(4\pi m\sqrt n/c).
 \tag{3.1}
\]
The phase is **+i**. It can be fixed without an odd-weight convention guess: expanding the factor exp(-2 pi i m/(c^2z)) and integrating z^{-k-j}exp(-2 pi i n z) along a horizontal line gives (-2 pi i)^k, hence 2 pi i^{-k} in the Bessel form. For k=3 this is +2 pi i. Formulae stated only for even weights must not be copied with an unchecked odd-weight sign.

For real b, set
\[
 S_\sigma(b)=\int_0^\infty J_2(v)v^{-1-2\sigma}\sin(bv)\,dv.
\]
Near sigma=0 this is an ordinary convergent integral. Mellin evaluation of J_2, or the Weber–Schafheitlin formula, gives for 0<b<1
\[
 S_\sigma(b)=2^{-2\sigma}b
 \frac{\Gamma(3/2-\sigma)}{\Gamma(3/2+\sigma)}
 {}_2F_1(3/2-\sigma,-1/2-\sigma;3/2;b^2),
 \tag{3.2}
\]
and for b>1
\[
 S_\sigma(b)=\frac{\Gamma(2-2\sigma)\sin(\pi\sigma)}8
 b^{2\sigma-2}{}_2F_1(1-\sigma,3/2-\sigma;3;b^{-2}).
 \tag{3.3}
\]
These formulas also follow by expanding the sine in (3.2) and the Bessel function in the Laplace-regularized version of (3.3), then continuing the convergent resulting hypergeometric functions. Gauss summation gives their common value at b=1; it is O(sigma) near zero. The function is odd in b.

Consequently
\[
 S_0(b)=b\sqrt{1-b^2}\quad(0<b<1),\qquad
 S_0(b)=0\quad(b\ge1).
 \tag{3.4}
\]
Equivalently the Fourier transform of J_2(ax)/x is -2i S_0(2 pi t/a). This is the source of the finite discriminant range below.

Poisson summation in m in (3.1), first in the nonempty strip 1/4<Re sigma<1/2 with a small convergence margin, gives
\[
 [q^n]\mathcal K_\sigma
 =1_{n=\square}\chi(\sqrt n)n^{(1-2\sigma)/2}
 +2\pi n(4\pi\sqrt n)^{2\sigma}
 \sum_{l\in\mathbb Z} S_\sigma\left(\frac l{2\sqrt n}\right)
 Z_{n,l}(\sigma),
 \tag{3.5}
\]
where
\[
 Z_{n,l}(\sigma)=\sum_{60\mid c}\frac{A_{n,l}(c)}{c^{2+2\sigma}},
\]
\[
 A_{n,l}(c)=\sum_{r\pmod c}\chi(r)
 \sum_{d\pmod c}^{*}\chi_{60}(d)e(d(r^2+lr+n)/c).
 \tag{3.6}
\]
For the character in (3.6), substitute r -> d r in the finite quadratic Fourier transform: epsilon(d)chi(d)=chi_60(d). The original positive-m sum is half of the full integer sum, and the Fourier transform contributes -2i. These factors yield exactly (3.5).

For completeness, the initial double rearrangement also has a direct
majorant. A quadratic Gauss sum in r, with the fixed modulo-four character,
is O(sqrt(c)); summing the units d gives |A_{n,l}(c)|=O(c^(3/2)),
uniformly in l. Thus the c-sum of its absolute majorant converges for
Re(sigma)>1/4. Formula (3.3) gives a summable l-majorant on compact
substrips with Re(sigma)<1/2. The near-zero singularity of the Fourier
kernel is locally integrable in this strip. Poisson summation may first be
applied with smooth truncations and then passed to the limit. It is not a
formal interchange at the critical boundary.

## 4. Exact local factors, including the exceptional real-square discriminants

Put Q(r)=r^2+lr+n and Delta=l^2-4n. For Delta !=0 define
\[
 t_{2,a}=\sum_{\substack{r\bmod2^{a+2}\\2^a\mid Q(r)}}
 \chi(r)\chi(Q(r)/2^a),
\]
\[
 t_{p,a}=\sum_{\substack{r\bmod p^{a+1}\\p^a\mid Q(r)}}
 \left(\frac{Q(r)/p^a}{p}\right),\quad p=3,5.
\]
Let
\[
 L_2(\sigma)=\sum_{a\ge0}t_{2,a}2^{-a(1+2\sigma)},\quad
 L_3(\sigma)=\sum_{a\ge0}t_{3,a}3^{-a(1+2\sigma)},
\]
\[
 L_5(\sigma)=\sum_{a\ge0}(-1)^at_{5,a}5^{-a(1+2\sigma)}.
 \tag{4.1}
\]
These sums are finite when Delta !=0: t_{p,a}=0 for a>v_p(Delta), p odd, and t_{2,a}=0 for a>=v_2(Delta)+2.

For odd p this follows by completing the square, (2r+l)^2-4Q=Delta. More explicitly, if v=v_p(Delta), then
\[
 t_{p,a}=\begin{cases}
 (p-1)p^{a/2},&a<v,\ a\text{ even},\\
 0,&a<v,\ a\text{ odd},\\
 -p^{v/2},&a=v,\ v\text{ even},\\
 p^{(v+1)/2}\left(\frac{-\Delta/p^v}{p}\right),&a=v,\ v\text{ odd},\\
 0,&a>v.
 \end{cases}
 \tag{4.2}
\]
For p=2, at a>=v_2(Delta)+2 any root has derivative valuation h=v_2(Delta)/2, or there are no roots. Translating r by 2^{a-h} leaves chi(r) unchanged and cycles Q(r)/2^a through all four residues modulo 4. The quadratic error vanishes modulo 4 because a-2h>=2. Summing chi over this orbit gives zero. This proves the finite cutoff, rather than assuming a computed prefix has ended.

For p not dividing 60 let R_a be the number of roots of Q modulo p^a, R_0=1, and t=chi_60(p)p^{-1-2 sigma}. The good local factor is
\[
 E_p(\sigma)=(1-t)\sum_{a\ge0}R_at^a.
 \tag{4.3}
\]
If v=v_p(Delta), then R_a=p^{floor(a/2)} for a<=v. Beyond v there are no roots when v is odd or Delta/p^v is a nonsquare, and R_a=2p^{v/2} otherwise. Hence (4.3) is an explicit rational function. When p does not divide Delta it is 1+chi_{60Delta}(p)p^{-1-2 sigma}.

The primitive Gauss sum of chi_60 is +sqrt60. Chinese remaindering gives
\[
 Z_{n,l}(\sigma)=\frac{\sqrt{60}}{60^{2+2\sigma}}
 L_2(\sigma)L_3(\sigma)L_5(\sigma)
 \prod_{p\nmid60}E_p(\sigma).
 \tag{4.4}
\]
This formula does not divide by a potentially zero first local coefficient. One derivation is to write c=60hk, with h supported on 2,3,5 and (k,60)=1. The part at 60h is sqrt60 times h times the product of the three finite character sums. The CRT factors contributed by h are 1 at 2, 1 at 3 and (-1)^{v_5(h)} at 5. At p^a||k the untwisted double sum is p^a(R_a-R_{a-1}), and CRT supplies chi_60(p)^a. This gives (4.4).

Write 60Delta=delta f^2, delta a fundamental discriminant. For Delta !=0, (4.4) is
\[
 \frac{\sqrt{60}}{60^{2+2\sigma}}L_2L_3L_5
 \frac{L(1+2\sigma,\chi_\delta)}{\zeta(2+4\sigma)}
 \prod_{p\mid60}\frac{1-\chi_\delta(p)p^{-1-2\sigma}}{1-p^{-2-4\sigma}}
 \prod_{\substack{p\mid\Delta\\p\nmid60}}
 \frac{E_p(\sigma)}{(1-p^{-2-4\sigma})/(1-\chi_\delta(p)p^{-1-2\sigma})}.
 \tag{4.5}
\]
Initially this is an absolutely convergent Euler-product calculation. It then gives meromorphic continuation. The only possible pole at sigma=0 is delta=1.

### 4.1 The potentially dangerous Pell contribution vanishes

A positive Delta with 60Delta a square must be Delta=60u^2, because Delta=l^2-4n is a discriminant. In particular l is even and
\[
 Q(r)=(r+l/2)^2-15u^2.
\]
For any nonzero norm v^2-15u^2 in Q_2, its odd part is 1 modulo 4. This can be checked without class field theory: factor the common power of 2 from v,u. If their remaining parities differ, the norm is odd and 1 modulo 4; if both are odd, it is 2 modulo 8.

With Haar measure of Z_2 normalized to one, (4.1) therefore gives
\[
 L_2(0)=4\int_{\mathbb Z_2}\chi(r)
 \chi(\operatorname{oddpart}Q(r))\,dr
 =4\int_{\mathbb Z_2}\chi(r)\,dr=0.
 \tag{4.6}
\]
The character sum decomposition into exact valuations justifies the first equality. The zero set has measure zero. The local series and its derivative converge locally uniformly for Re sigma>-1/4, from the elementary quadratic root-count bound O(2^{a/2}). Thus the zero in (4.6) cancels the simple zeta pole in (4.5).

The boundary case Delta=0 is treated separately, not sent through a fundamental-discriminant formula. Outside 60 its factor is
\[
 \prod_{p\nmid60}\frac{1-p^{-2-4\sigma}}{1-p^{-1-4\sigma}}.
\]
Its only pole near zero is the simple pole of zeta(1+4 sigma). Here Q is a square, its odd part is again 1 modulo 4, and the same L2(0)=0 cancels the pole.

More explicitly, write n=s^2 and l=2s. If s is even then L2(sigma) is
identically zero. If s is odd, decomposition according to
b=v2(r+s)>=1 gives
\[
 L_2(\sigma)=2\chi_{-4}(s)\frac{t(1-2t)}{1-t},
 \qquad t=2^{-1-4\sigma}.
\]
The b=1 shell has character chi4(s), and every b>=2 shell has character
-chi4(s); their Haar measures are 2^(-b-1). In particular
L2(0)=0 and L2'(0)=8 chi4(s) log(2). This checks the square-discriminant
boundary independently of the nonzero-discriminant cutoff formula.

**Therefore Z_{n,l}(sigma) is holomorphic at zero for every integer l.** In particular there is no unaccounted residue from a positive-square discriminant multiplying the vanishing factor in (3.3). This is why an infinite Pell correction is absent in this particular level and character. Such a correction cannot be discarded in a general theta-kernel problem.

## 5. The critical limit is a cusp form and has finite class-number coefficients

The limit in (3.5) may be taken through the l sum. Here are the needed uniform estimates. For fixed n and |sigma| sufficiently small, the explicit local factors in (4.5) are bounded by C_eta(1+|Delta|)^eta times the Dirichlet L-factor, for any fixed eta>0, after slightly increasing the exponent. This follows from the quadratic root counts: the sums have geometric majorants with ratio at most p^{-1/2+2|Re sigma|}. Products of their normalized exceptional factors over primes dividing Delta are bounded by an arbitrarily small positive power of |Delta|. Finite small primes are absorbed in the constant.

For a nonprincipal character of conductor q, periodic partial sums give, uniformly near s=1,
\[
 L(s,chi)=O_\eta(q^{\max(0,1-\Re s)+\eta}).
\]
One can obtain this crude bound by splitting at q and applying partial summation to the tail; no subconvex estimate is needed. The principal-character cases above are bounded after cancellation, including derivatives, by the same kind of bound. Since q<=60|Delta|, choosing a sufficiently narrow neighborhood of sigma=0 yields, for example,
\[
 |Z_{n,l}(\sigma)|\le C_n(1+|l|)^{1/4}.
\]
For |l|>=4sqrt(n), equation (3.3) uniformly gives
\[
 |S_\sigma(l/(2\sqrt n))|\le C_n(1+|l|)^{-2+1/8}
\]
on a smaller neighborhood. Their product has exponent -13/8, less than -1. It is summable uniformly. The finitely many other l cause no problem; (4.6) already removes their possible pole. The same estimates on compact substrips connect this continuation to the initial Poisson strip.

Consequently (3.4) proves
\[
 k_n:=[q^n]\mathcal K_0
 =1_{n=\square}\chi(\sqrt n)\sqrt n
 +\frac\pi2\sum_{l^2<4n}l\sqrt{4n-l^2}\,Z_{n,l}(0).
 \tag{5.1}
\]
All first 36 coefficients have finite analytic continuations. The Sturm coefficient map is injective on the finite-dimensional cusp space; applying a left inverse and continuing the linear relations in its image gives a unique continued cusp form \(\mathcal K_0\in S_3(60,epsilon)\). Identity of every other Fourier coefficient follows by analytic continuation from the initial strip. This establishes modularity before using any finite coefficient comparison below.

For negative Delta, use the classical imaginary-quadratic class-number formula
\[
 L(1,\chi_\delta)=\frac{\pi H(\delta)}{\sqrt{|\delta|}},
 \quad H(\delta)=h(\delta)\ (\delta<-4),\quad
 H(-3)=1/3,\ H(-4)=1/2.
\]
Let \(60\Delta=\delta f^2\), and set \(L_p=L_p(0)\),
\[
 E_{0,p}=\frac{1-p^{-2}}{1-\chi_\delta(p)/p}.
\]
The positive-l version of (5.1) is the entirely rational finite formula
\[
 \boxed{
 k_n=1_{n=\square}\chi(\sqrt n)\sqrt n+
 \frac1{600}\sum_{l=1}^{\lfloor\sqrt{4n-1}\rfloor}
 l f H(\delta)L_2L_3L_5
 \prod_{p=2,3,5}\frac{1-\chi_\delta(p)/p}{1-p^{-2}}
 \prod_{\substack{p\mid\Delta\\p\nmid60}}\frac{E_p(0)}{E_{0,p}}.
 }
 \tag{5.2}
\]
The factors pi and square roots cancel using zeta(2)=pi^2/6. The l and -l terms are equal, giving the denominator 600 rather than 1200. There is no l=0 contribution.

The first nonzero coefficients are
\[
 \mathcal K_0=\frac{17}{16}q+\frac7{48}q^3+\frac1{24}q^5
 -\frac{11}{12}q^7-\frac{55}{16}q^9+\frac{17}{6}q^{11}+\cdots.
\]

## 6. A finite modular certificate and the f component

The index of Gamma_0(60) is 144, so the weight-three Sturm bound is 36. The eighteen eta quotients e_0,...,e_17 used in the retained norm proof are listed afresh, with all exponents, in `critical_norm_certificate.md` and both scripts. They are all cusp forms with character epsilon; the current programs recheck weight, both eta congruences, character squareclass and positivity at every cusp.

Equation (5.2) through q^36 proves the following identity by Sturm:
\[
 \mathcal K_0=\sum_{j=0}^{17}K_j e_j,
\]
where the complete rational vector is
\[
 K=(0,-15/8,269/128,-5/2,127/32,-235/64,-4,35/8,
 29/128,11/16,5/4,1/6,37/64,15/64,33/128,-3/16,13/48,87/128).
 \tag{6.1}
\]
The primary calculation checks this through q^200, and the separate standard-library calculation through q^240. Those extra coefficients are controls; **modularity plus Sturm**, not a long matching prefix alone, proves the identity.

The exact f-orbit component needed for the Petersson pairing is
\[
 \boxed{H_f=\frac{49}{144}P_1-\frac5{48}P_2
                   +\frac5{16}P_3+\frac1{16}P_4.}
 \tag{6.2}
\]
An explicit cusp form R=sum R_j e_j satisfies
\[
 \boxed{\mathcal K_0-H_f=(T_7^2+64)R.} \tag{6.3}
\]
Its coefficients are
\[
 R=(1/192,-43/55296,17/27648,-619/27648,-871/27648,
 -11/384,-269/3456,1/6,-1/6144,-1/3072,0,0,
 49/2880,-13/8640,1/72,-1/160,0,0).
 \tag{6.4}
\]
The certificate verifies (6.3) directly from eta coefficients through q^1764= q^{49*36}, not just by trusting a previously stored Hecke matrix. For a coefficient sequence R_n,
\[
 [(T_7^2+64)R]_n
 =R_{49n}+(15-49\,1_{7\mid n})R_n+2401R_{n/49}.
\]
Since T7*=-T7 and T7f=-8if, equation (6.3) contributes zero to <f,K0>.

In the four-form space, U3 satisfies U3^4-2U3^2+81=0. Its f eigenvalue is lambda=sqrt5-2i. If h_j is the coefficient of q in U3^j H_f, direct computation gives
\[
 (h_0,h_1,h_2,h_3)=(49,-19,-91,-119)/144.
\]
The primitive spectral projector therefore gives the coefficient of normalized f in H_f as
\[
 \delta_f=
 \frac{h_3+\lambda h_2+(\lambda^2-2)h_1+(\lambda^3-2\lambda)h_0}
 {4\lambda^3-4\lambda}
 =\boxed{\frac{(7-i\sqrt5)(7-i)}{576}.}
 \tag{6.5}
\]
An independent calculation expands H_f as the sum of delta_f^sigma f^sigma over the four coefficient embeddings and verifies all four P_j coefficients. The four conjugate forms are Petersson-orthogonal: the pairs of good eigenvalues a7=-8i and a11=4i sqrt5 distinguish all four embeddings, and these Hecke operators are normal. The current scripts verify both good eigenidentities by Sturm. Hence
\[
 \boxed{\langle f,\mathcal K_0\rangle_{60}
 =\overline{\delta_f}\,N_f.} \tag{6.6}
\]

## 7. Exact evaluation of the modular critical value

Continue (2.3) to sigma=0 and use (1.1)–(1.3):
\[
 \overline{\delta_f}N_f
 =\frac{D_f(3)}{16\pi^2}
 =\frac{75R(3)}{128\pi^4}\mathscr L(3).
\]
Now
\[
 \overline{\delta_f}=\frac{(7+i\sqrt5)(7+i)}{576},\qquad
 R(3)=\frac{16(7+i\sqrt5)(7+i)}{675},
\]
so exactly
\[
 \frac{128\overline{\delta_f}}{75R(3)}=\frac18.
\]
This proves the stated value at 3. No numerical critical value selects the factor 1/8.

The adopted completion is, up to an s-independent factor,
\[
 q_L^s\Gamma(s)\Gamma((s-1)/2)\mathscr L(s),
 \quad q_L=\frac{\sqrt{960}}{2\pi^{3/2}},
\]
with sign + and s -> 5-s. Therefore
\[
 \mathscr L(2)=\frac{2q_L}{\sqrt\pi}\mathscr L(3)
 =\frac{8\sqrt{15}}{\pi^2}\mathscr L(3),
\]
which proves the value at 2. Finally `norm_proof.md` gives the critical degree-two expression in equation (3.5) of the article.


## 8. Exact files and mathematical scope

`../engine/verify_critical_norm.py` independently constructs 241 coefficients
(including the constant term) through q^240, the Euler-product eta expansions
through q^1764, the Hecke witness through the Sturm bound 36, the Galois
spectral coefficient, and the exact constants. It uses only Python's
standard library. `../primary/critical_norm/proofs/verify_exact.py` uses a
different class-number enumeration and eta algorithm and reconstructs
coefficients through q^200. Their common 201 coefficients agree exactly.
The additional primary audit checks sixteen higher-prime-power cyclotomic
Gauss sums and 117 raw Delta=0 character sums, and verifies the
Weber--Schafheitlin coefficients by gamma duplication and reflection.
These finite tests supplement the all-index proofs; they do not replace them.

## 9. Sources and normalization dictionary

1. Z. L. Lim, F.-T. Tu and J.-D. Yu, *Galois representations and modularity
   for twisted moments of Kloosterman sums*, arXiv:2609.13576v1.
   Theorem 1.2.2(iii) identifies the arithmetic form. Appendix A.3,
   printed p.31, supplies the trivial Euler factor at 2, the conductor
   2^6*15=960, and the gamma factors. The motive's L5(s+1) equals the
   classical twisted symmetric-square function used here. Equation (44),
   printed p.34, contains the ordinary functional-equation relation
   pi^2*L(2)=8*sqrt(15)*L(3). Its Bessel equality retains a question mark.
2. NIST DLMF 10.22.56, Weber--Schafheitlin integrals, used with
   sin(bv)=sqrt(pi*b*v/2)*J_(1/2)(bv); DLMF 15.4.20, Gauss summation;
   DLMF 5.5.3 and 5.5.5, reflection and duplication.
3. Standard holomorphic Poincare series and their Fourier expansion,
   fixed-character Kloosterman bounds, primitive quadratic Gauss sums,
   the imaginary-quadratic class-number formula, Hecke adjointness and
   Sturm's theorem. The exact conventions and specializations are given above.
4. The complete norm and coefficient-twist proofs are supplied as
   `norm_proof.md` and `twist_proof.md`. They do not depend on this bridge.

This supplement is adapted from the user-supplied 16 September bridge
bundle. The new square-discriminant formula and explicit initial
rearrangement majorant above were added in the v4 audit. No external article,
private correspondence or font is redistributed. No novelty, independent
human-review or proof-assistant-verification claim is made.
