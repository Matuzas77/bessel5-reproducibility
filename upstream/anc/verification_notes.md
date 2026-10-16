# Verification notes: a quadratic critical-value conjecture for the fifth Bessel moment

Jonas Matuzas. Numerical notes retained from 15 September 2026; v4 scope updated 16 September 2026.

Version 4 proves the purely modular identity (3.5) by the separate argument in `proofs/critical_norm_bridge.md`. None of the numerical comparisons below proves a Bessel equality. The period-enclosure algorithms and all eight original numerical/proof JSON records are unchanged.

These notes document the analytic inputs, normalizations, and complete error bounds used by `verify.py`. The exact Bessel-to-modular evaluations in the paper are conjectures. A containing-zero interval with radius below 10^(-358) is not a proof that its enclosed number vanishes. The computations below prove smallness bounds, subject to the stated classical analytic inputs and correctness of the supplied integer-arithmetic programs; they are not proof-assistant formalizations.

## 1. Quantities and precise claim boundary

All square roots of positive real numbers are positive. The integration measure is **dt**, not t dt. For j >= 0, set

\[
 A_j=\int_0^\infty t^{2j}K_0(t)^5\,dt,\quad
 B_j=\int_0^\infty t^{2j}I_0(t)K_0(t)^4\,dt,\quad
 C_j=\int_0^\infty t^{2j}I_0(t)^2K_0(t)^3\,dt.
\]

All these integrals converge. At zero their singularities are integrable powers of log(t); at infinity their exponential decay rates are respectively 5, 3, and 1. Let r=sqrt(5), and let f and L=L(f,2) be defined by the embedding, Mellin definition, and eta formula in the paper. The square in Re(L^2) or Im(L^2) is taken before the real or imaginary part. Define the modular targets

\[
 A_*={\pi\over4}[(125-13r)\Re(L^2)+(41r-125)\Im(L^2)],
\]
\[
 B_*={3\sqrt{15}(5-r)\over2}|L|^2,\qquad
 E_*={3\pi^2(5-r)\over4}|L|^2,
\]
\[
 C_*={1\over4\pi}[(42+6r)\Re(L^2)+(33r-69)\Im(L^2)].
\]

The C_* formula is an additional companion check retained in the ancillary outputs; it is not needed to state the short note. In Lim–Tu–Yu's notation [LTY, (42)],

\[
 D_{5,\mathrm{odd}}=B_0,\qquad
 D_{5,\mathrm{even}}=A_0C_1-A_1C_0.
\]

The paper reports conservative bounds for A_0-A_*, B_0-B_*, and D_{5,even}-E_*. The program also reports C_0-C_*. The even and odd conjectures are equivalent by the exact source identity in Section 4 below. The pure-moment conjecture is a separate assertion.

[LTY, (44)] already conjectures

\[
 B_0=\pi^2L(\chi_{-4}\operatorname{Sym}^2 f,2).
\]

Its Bessel-to-L equality carries a question mark in v1. We credit this conjecture to Lim–Tu–Yu and do not use it to initialize any numerical computation. The present 300-plus-digit checks do **not** independently evaluate the symmetric-square L-function. They compare the Bessel source with the quadratic expressions in the degree-two L-value L(f,2).

## 2. Bessel source: positive series, not stored Bessel decimals

### 2.1 Classical triple moments

Let T_n=int_0^infinity t^(2n) K_0(t)^3 dt. The classical evaluations and recurrence in [BBBG, Table 1 and (8)] give

\[
 T_0={\sqrt3\pi^3\over8}H_0,\qquad
 T_1={\sqrt3\pi^3\over288}H_1,
\]
\[
 H_0=\sum_{n\ge0}{\binom{2n}{n}^3\over256^n},\qquad
 H_1=\sum_{n\ge0}{\binom{2n}{n}^3\over256^n(n+1)^2},
\]
\[
 9T_{n+2}=2(20n^2+40n+21)T_{n+1}-(2n+1)^4T_n.
\]

Both H-series have positive rational terms and consecutive-term ratio below 1/4. After the term with index M, the complete remainder is smaller than the last retained term divided by 3. No gamma value or floating special function is used.

Squaring the positive I_0 series and applying Tonelli gives the exact identity

\[
 C_j=\sum_{n\ge0}{\binom{2n}{n}\over4^n(n!)^2}T_{n+j}.
\]

The decimal implementation evolves h_n=T_n/[4^n(n!)^2]; the binary implementation evolves r_n=T_n/(2n)!. Their recurrences are obtained by substituting these definitions into the displayed T-recurrence. The initial values and every subtraction are outward enclosed. Forward recurrence is ill-conditioned here; working precision alone is not an accuracy claim. The programs check that every triple-moment interval used has positive lower endpoint.

### 2.2 The differentiated Struve companion and its contact terms

Use the Hilbert-transform convention H g(x)=(1/pi)PV int g(y)/(x-y) dy. With even I_0(|x|), K_0(|x|), [Zhou, (3.1)] gives

\[
 H(\pi^2 I_0^2K_0^2-K_0^4)=-2\pi I_0K_0^3\operatorname{sgn}(x).
\]

The Fourier transform of K_0(c|x|), or its elementary cosine representation, gives

\[
 H[K_0(c|x|)]={\pi\over2}\operatorname{sgn}(x)
 [I_0(c|x|)-\mathbf L_0(c|x|)].
\]

The latter difference also has the convergent Laplace representation in [DLMF, 11.5]. Both even functions paired here belong to L^2(R). Skew-adjointness of H therefore yields, for c>0,

\[
 \mathcal A(c)-\pi^2\mathcal C(c)
 =\pi^2\int_0^\infty I_0(t)K_0(t)^3
 [I_0(ct)-\mathbf L_0(ct)]\,dt,
\]

where mathcal A(c)=int K_0(t)^4 K_0(ct) dt and mathcal C(c)=int I_0(t)^2 K_0(t)^2 K_0(ct) dt. Separate the right-hand terms only for 0<c<2.

Set D_c=partial_c^2+c^(-1)partial_c and V(c)=int I_0 K_0^3 boldL_0(ct) dt. The modified Struve equation [DLMF, 11.2] implies

\[
 D_c\mathbf L_0(ct)=t^2\mathbf L_0(ct)+{2t\over\pi c}.
\]

A finite number of derivatives can be passed under the integrals on any compact neighborhood of c=1 within (0,2): a polynomial times exp(-(2-c_max)t) dominates at infinity, and integrable logarithmic powers dominate at zero. Since D_c acts on I_0(ct) and K_0(ct) as multiplication by t^2, evaluation at c=1 gives

\[
 A_j=\pi^2(2C_j-W_j),\qquad W_j=(D_c^jV)(1).
\]

Define Domb integers

\[
 D_n=\sum_{k=0}^n\binom nk^2\binom{2k}{k}\binom{2n-2k}{n-k}.
\]

The classical quartic moments in [BBBG] give int t^(2n+1) I_0 K_0^3 dt=pi^2 D_n(n!)^2/16^(n+1). Substituting the modified-Struve series and differentiating gives

\[
 W_j={\pi\over8}\sum_{n\ge0}
 {D_n\prod_{h=0}^{j-1}(2n+1-2h)^2
 \over4^n(2n+1)^2\binom{2n}{n}^2}.
\]

The empty product is 1. Terms with n<j must not be discarded. W_j is **not** simply the t^(2j)-weighted Struve integral. For example,

\[
 W_1-\int_0^\infty t^2 I_0K_0^3\mathbf L_0(t)\,dt={\pi\over8},
\]

and the analogous correction for j=2 is 5pi/32. These contact terms are automatically included by the differentiated series. A_j and C_j require no Bessel period pairing and no modular input.

### 2.3 Complete source tails

The positive hyperbolic representation of K_0 and cosh(u)>=1+u^2/2 imply

\[
 K_0(t)\le\sqrt{\pi/(2t)}e^{-t}\quad(t>0).
\]

After retaining n=0,...,N in the C_j series, put n=N+1 and m=n+j. For j=0,1,2 and N>=20, the complete omitted sum is less than

\[
 {25\over2}{(2n)!(4m-2)!
 \over4^{n+2m-1}(n!)^4(2m-1)!9^m}.
\]

Indeed, integrating the K-bound makes the first omitted term at most pi^2 sqrt(6)/4 times this rational expression without its prefactor 25/2. The consecutive majorant ratio is

\[
 {(2n+1)(16(n+j)^2-1)\over72(n+1)^3}
 <{4\over9}\left({n+j\over n+1}\right)^2
 \le {4\over9}(22/21)^2<{1\over2}.
\]

Sum that geometric bound and use pi^2<10 and sqrt(6)<5/2.

For the Struve series, D_n<=16^n and binom(2n,n)>=4^n/(2n+1). Its tail, including the prefactor pi/8, is bounded by

\[
 \sum_{n>N}(\text{term of }W_j)<{(2N+3)^{2j}\over4^{N+1}}.
\]

To check this, the unscaled majorant is 4^(-n)(2n+1)^(2j). Its ratio is at most (1/4)(43/41)^4<1/2 for the indices used; pi/4<1 completes the bound. Thus every source truncation, not only a difference between two truncations, is bounded.

## 3. Modular value: an independent four-eta Mellin calculation

The eta quotients P_1,...,P_4 and the coefficients omega,epsilon are exactly those in the paper's eta formula. Define eta(tau)=q^(1/24) product_(n>=1)(1-q^n), q=exp(2pi i tau). The implementation independently constructs each integer coefficient array by logarithmic differentiation and Euler multiplication.

For each exponent vector e_d, the standard eta-quotient criterion uses sum e_d=6, sum d e_d=0 mod 24, sum (60/d)e_d=0 mod 24, and character squareclass (-1)^3 product d^(e_d)=-15 modulo rational squares. The cusp order at denominator c|60 is

\[
 {60\over24\gcd(c^2,60)}\sum_{d\mid60}{e_d\gcd(c,d)^2\over d}.
\]

Every such order is positive for all four products. The exact audit records these rational orders. The level index is 144, giving the characteristic-zero Sturm bound 3*144/12=36. The supplied exact Q(sqrt(5),i) calculation checks

\[
 a_{pn}+\chi_{-15}(p)p^2a_{n/p}=a_p a_n\qquad(1\le n\le36)
\]

for every prime p<=31, with the second term omitted for p|60 or p not dividing n. These 396 coefficient witnesses certify the indicated Hecke relations by the standard Sturm theorem; they are not 396 independent theorems. The orbit label is used as in [LTY] and [LMFDB], while the explicit eta expression fixes the numerical input without a database query. See [Stein] for the standard Hecke and Sturm framework.

### 3.1 Exact folding formula

Put t_0=1/sqrt(60), q_0=exp(-2pi/sqrt(60)), and

\[
 S_{j,k}=\sum_{n\ge1}[q^n]P_j\,{q_0^n\over n^k}.
\]

Eta inversion [DLMF, 23.18] gives

\[
 P_j(i/(60t))=c_j(\sqrt{60}t)^3P_{\sigma(j)}(it),
\]

with sigma=(3,4,1,2) and c=(3,3,1/3,1/3). Split 4pi^2 int_0^infinity t P_j(it)dt at t_0 and substitute t=1/(60u) in the lower half. Integration of exp(-2pi n t) and t exp(-2pi n t) on [t_0,infinity) gives

\[
 L(P_j,2)=S_{j,2}+{2\pi\over\sqrt{60}}
 (S_{j,1}+c_jS_{\sigma(j),1}).
\]

Consequently L=L(P_1,2)+omega L(P_2,2)+3epsilon L(P_3,2)+3epsilon conjugate(omega)L(P_4,2). In particular this computation does not insert the conjectured Bessel value to recover L.

### 3.2 Eta coefficient and arithmetic tails

For a fixed exponent vector, let R=3/4 and

\[
 M=3^{\lceil J\rceil},\qquad
 J=\sum_{e_d>0}{e_dR^d\over1-R^d}
   +\sum_{e_d<0}{(-e_d)R^d\over(1-R^d)^2}.
\]

Bounding log(1+x) by x and -log(1-x) by x/(1-x) bounds the eta product on |q|=R by exp(J)<=M. Its positive integral q-shift has modulus at most 1 there. Cauchy's estimate and the certified inequality q_0<1/2 then give, after n=N,

\[
 \left|S_{j,k}-\sum_{n=1}^N[q^n]P_j\,{q_0^n\over n^k}\right|
 \le {3M\over(N+1)^k}(2/3)^{N+1}.
\]

Pi is enclosed using alternating Machin series pi=16atan(1/5)-4atan(1/239), with complete alternating remainders. Positive square roots use integer-square-root lower and upper bounds. A positive exponential Taylor series with a geometric tail encloses exp(2pi/sqrt(60)); its reciprocal encloses q_0. All subsequent real operations round outward on an integer grid.

### 3.3 Independent binary modular check

The second modular program uses binary fixed-point arithmetic, Euler expansion in q^2, Cauchy radius 2/3, and q_0<9/20. The resulting complete tail for S_k has factor

\[
 {40M\over13(N+1)^k}(27/40)^{N+1}.
\]

For an independently written fold, put u=eta_2^3 eta_30^3/(eta_6^3 eta_10^3), so uP_1=P_2, uP_3=P_4 and u is invariant under positive-axis Fricke inversion. Let S_(+,j,k), S_(-,j,k) use coefficients of u^j(P_1+3P_3), u^j(P_1-3P_3). Set

\[
 U_j=4S_{+,j,1}+{\sqrt{60}\over\pi}S_{+,j,2},\qquad
 V_j={\sqrt{60}\over\pi}S_{-,j,2},
\]
\[
 U=U_0+{r-1\over2}U_1,\qquad V=V_0-(r+2)V_1.
\]

These follow from the direct fold by adding and subtracting the paired Mellin values. Exact coefficient algebra gives

\[
 L={\pi\over2\sqrt{60}}[(1+\epsilon)U+(1-\epsilon)V].
\]

This second description supplies the modular cross-check; it is not used to initialize the direct Mellin computation. The old binary engine also outputs j=2,3 quantities, retained only as reproducibility controls and not needed for this paper's conjectures.

## 4. The exact determinant factor and odd-moment normalization

This derivation specializes public results of Chuang [Chu] and makes no Bessel-to-modular conjecture. Use zero-based row and column indices. Define

\[
 X=\begin{pmatrix}A_0&A_1&A_2\\\pi B_0&\pi B_1&\pi B_2\\
 \pi^2C_0&\pi^2C_1&\pi^2C_2\end{pmatrix},\quad
 R=\operatorname{diag}(-1,i,1),\quad T=\operatorname{diag}(32,8,2).
\]

Proposition 24 of [Chu] gives their period matrix P=RXT. For k=5, the middle and full pairings coincide. Let D be their Poincare pairing matrix. Proposition 13 says its first row is (0,0,8/15). It is symmetric, so the last row of D^(-1) is (15/8,0,0). Accordingly J=T D^(-1) T has last row (120,0,0), while (J^(-1))_(0,0)=0.

The Betti formula in [Chu, Proposition 18], with the row phases R removed, is

\[
 H=\begin{pmatrix}1/4&0&1/8\\0&-1/10&0\\1/8&0&3/40\end{pmatrix}.
\]

Proposition 32 of [Chu] therefore becomes X J X^t=64pi^6 H. Corollary 39 at matrix size three gives

\[
 \det X={2\pi^9\over75\sqrt{15}}.
\]

First, invert the pairing relation to obtain

\[
 X^{-1}={JX^tH^{-1}\over64\pi^6},\qquad
 (X^{-1})_{2,1}=-{75B_0\over4\pi^5}.
\]

The cofactor formula gives the same entry as -pi^2(A_0C_1-A_1C_0)/det X. Equating the two expressions proves

\[
 \boxed{D_{5,\mathrm{even}}={\pi^2\over2\sqrt{15}}D_{5,\mathrm{odd}}.}
\]

Second, (J^(-1))_(0,0)=0 in X^t H^(-1)X=64pi^6 J^(-1) yields

\[
 12A_0^2-40\pi^2A_0C_0+40\pi^4C_0^2-5\pi^2B_0^2=0.
\]

Since B_0>0 by its integral definition, the paper verifier reconstructs it as

\[
 B_0={1\over\pi}\sqrt{{12A_0^2-40\pi^2A_0C_0+40\pi^4C_0^2\over5}}.
\]

The radicand is enclosed away from zero before taking the root. This supplies the odd comparison without needing the unspecified entries of D. The raw source engines additionally output a historical reconstruction of the complete B row using a fuller normalized pairing matrix. In the present paper that full row is only a consistency control; `check_paper.py` computes its own positive conic root from the independent A_0,C_0 intervals. The even determinant is calculated from A_0,C_0,A_1,C_1, not from the conjectured odd target.

## 5. Exact algebraic compatibility, executions, and independence

The modular targets satisfy the same conic identically. With a=(3-r)/4, b=20+2r, and k=4r-10, the independently folded U,V coordinates satisfy

\[
 A_*={\pi^3\over60}(bU^2-2kaUV-ba^2V^2),\quad
 C_*={\pi\over5}(U^2+aUV-a^2V^2),\quad
 B_*={\pi^2\over\sqrt{15}}(U^2+a^2V^2).
\]

The standard-library program `audit_algebra.py` verifies these three changes of coordinates and their conic compatibility as exact polynomial identities in Q(sqrt(5)). This checks the normalization; it does not prove source and target periods coincide.

The freshly executed settings were:

| Calculation | Arithmetic grid | Source cutoff | Initial H cutoff | Eta cutoff |
|---|---|---:|---:|---:|
| Decimal source | 1100 decimal places | 1000 | 2300 | — |
| Binary source | 4200 binary places | 1050 | 2600 | — |
| Direct four-eta Mellin | 750 decimal places | — | — | 2400 |
| Independent eta fold | 3000 binary places | — | — | 2600 |

`paper_check.json` records all rational endpoints, not just rounded displays. Every comparison interval contains zero. Its actual strict discrepancy bounds are:

| Comparison | Decimal source | Binary source |
|---|---:|---:|
| A_0-A_* | 10^(-358) | 10^(-376) |
| C_0-C_* | 10^(-360) | 10^(-377) |
| D_(5,odd)-B_* | 10^(-358) | 10^(-376) |
| D_(5,even)-E_* | 10^(-352) | 10^(-369) |

The article deliberately quotes the weaker bound 10^(-352) for both determinants. The new conic-based odd calculation is tighter than the historical full-row cross-product calculation. All nine binary source intervals nest in the decimal source intervals. The two modular values differ by less than 10^(-422) in both real and imaginary parts.

The two source implementations share the explicitly stated analytic identities, but use different recurrences and arithmetic backends; their Domb and initial-H algorithms also differ. The two modular computations share eta inversion, but use different coefficient organizations, Cauchy radii, and arithmetic backends. This is implementation independence, not independence of every mathematical input. Conic and determinant identities also constrain the reported tests, so the multiple comparisons must not be counted as independent transcendence statements.

`verify.py` copies only the programs into an empty temporary workspace and regenerates every JSON output. Stored outputs are used only afterward for exact comparison. The default run neither changes the distributed records nor uses the network. `--write` explicitly regenerates reference records. The files in `logs/` distinguish the current publication replay from the earlier archive replay. No external human review or formal verification is asserted.

## References and precise dependencies

[LTY] Z. L. Lim, F.-T. Tu, and J.-D. Yu, *Galois representations and modularity for twisted moments of Kloosterman sums*, arXiv:2609.13576v1 (11 September 2026), Theorem 1.2.2(iii), Appendix A.3, equations (42) and (44). https://arxiv.org/pdf/2609.13576v1

[Chu] P.-H. Chuang (with an appendix in collaboration with J.-D. Yu), *On the periods of twisted moments of the Kloosterman connection*, Ramanujan Journal 65 (2024), 1227–1277; arXiv:2306.15216v2. Propositions 13, 18, 24, 32, and Corollary 39 are used in Section 4. https://arxiv.org/html/2306.15216v2

[BBBG] D. H. Bailey, J. M. Borwein, D. Broadhurst, and M. L. Glasser, *Elliptic integral evaluations of Bessel moments and applications*, Journal of Physics A 41 (2008), 205203; arXiv:0801.0891v2. The triple moments, their recurrence, and the Domb quartic moments are classical inputs, not claimed as discoveries of this note. https://arxiv.org/html/0801.0891v2

[Zhou] Y. Zhou, *Hilbert transforms and sum rules of Bessel moments*, Ramanujan Journal 48 (2019), 159–172; arXiv:1706.01068v2, especially equation (3.1). https://arxiv.org/html/1706.01068v2

[DLMF] NIST Digital Library of Mathematical Functions, sections 10.25 and 10.32 (modified Bessel functions and integral representations), 11.2 and 11.5 (modified Struve equation, series, and integral representations), and 23.18 (Dedekind eta transformations). https://dlmf.nist.gov/10.32 ; https://dlmf.nist.gov/11.2 ; https://dlmf.nist.gov/11.5 ; https://dlmf.nist.gov/23.18

[Stein] W. A. Stein, *Modular Forms, a Computational Approach*, Graduate Studies in Mathematics 79, American Mathematical Society (2007). Standard Hecke operators and Sturm bounds. Author's text: https://wstein.org/books/modform/modform/

[LMFDB] The LMFDB Collaboration, *The L-functions and Modular Forms Database*, orbit 60.3.b.a. https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/60/3/b/a/

The supplied files reproduce computations, not copies of those publications. The general Hilbert-transform, eta modularity, and period-pairing theorems are external analytic inputs; finite program assertions do not replace them.
