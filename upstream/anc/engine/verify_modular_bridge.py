#!/usr/bin/env python3
"""Exact audit of the two proof normalizations and the Chuang cofactor.

No period, numerical value, or database record is read. Standard library only.
The general modularity, Sturm and period-pairing results are analytic inputs.
"""
from __future__ import annotations
from fractions import Fraction as Q
from pathlib import Path
from math import factorial, comb, gcd
import json
from certify_lvalue import EXPONENTS, coefficients_log
from verify_norm import FTERMS, add, mul, scale, conj, ONE

Z=(Q(0),)*4
R=(Q(0),Q(1),Q(0),Q(0))
I=(Q(0),Q(0),Q(1),Q(0))
def k(a=0,b=0,c=0,d=0):return tuple(map(Q,(a,b,c,d)))
def sub(a,b):return add(a,scale(b,-1))
def power(a,n):
    v=ONE
    for _ in range(n):v=mul(v,a)
    return v
def sigma(a):return a[0],-a[1],-a[2],a[3]
def mm(A,B):
    return [[sum((A[i][h]*B[h][j] for h in range(len(B))),Q(0)) for j in range(len(B[0]))] for i in range(len(A))]
def invert(A):
    n=len(A);a=[list(map(Q,row))+[Q(int(i==j)) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        p=next(i for i in range(j,n) if a[i][j]);a[j],a[p]=a[p],a[j]
        den=a[j][j];a[j]=[x/den for x in a[j]]
        for i in range(n):
            if i!=j:
                den=a[i][j];a[i]=[x-den*y for x,y in zip(a[i],a[j])]
    return [r[n:] for r in a]

def main():
    N=2916
    P=[coefficients_log(e,N) for e in EXPONENTS]
    F=[P[0][n]-3*P[2][n] for n in range(N+1)]
    carrier=[0]*(37)
    for c,e in FTERMS:
        a=coefficients_log(e,36)
        for n in range(37):carrier[n]+=c*a[n]
    assert carrier==F[:37] # Both sides modular; Sturm bound36.
    for n in range(37):assert F[81*n]-2*F[9*n]+81*F[n]==0
    weights=(ONE,k(0,Q(1,5),0,Q(-2,5)),k(1,Q(4,5),-2,Q(2,5)),k(0,1,2,0))
    coeff=[tuple(sum(weights[j][h]*P[j][n] for j in range(4)) for h in range(4)) for n in range(N+1)]
    # The eight-eta/U3 specification in the norm proof is the same four-eta f.
    cs=(k(Q(5,12),Q(1,20),1,Q(19,60)),
        k(Q(-7,36),Q(-13,540),Q(29,108),Q(11,90)),
        k(Q(-1,12),Q(-1,20),0,Q(1,60)),
        k(Q(-1,36),Q(-7,540),Q(-1,108),Q(-1,90)))
    for n in range(37):
        assert tuple(sum(cs[j][h]*F[3**j*n] for j in range(4)) for h in range(4))==coeff[n]
    # The Hecke eigenidentities needed for both proofs, including orthogonality.
    for p in (3,5,7,11):
        chi=0 if 60%p==0 else (1 if pow((-15)%p,(p-1)//2,p)==1 else -1)
        for n in range(1,37):
            lhs=coeff[p*n]
            if n%p==0:lhs=add(lhs,scale(coeff[n//p],chi*p*p))
            assert lhs==mul(coeff[p],coeff[n])
    assert coeff[3]==k(0,1,-2,0)
    pairs={(a,b) for a,b in [(coeff[7],coeff[11]),(conj(coeff[7]),conj(coeff[11])),
                            (sigma(coeff[7]),sigma(coeff[11])),(conj(sigma(coeff[7])),conj(sigma(coeff[11])))]}
    assert len(pairs)==4
    # The depleted twist: an equality at level300, whose Sturm bound is180.
    a5=sigma(coeff[5])
    for n in range(181):
        chi5=0 if n%5==0 else (1 if n%5 in (1,4) else -1)
        rhs=sigma(coeff[n])
        if n%5==0:rhs=sub(rhs,mul(a5,sigma(coeff[n//5])))
        assert scale(coeff[n],chi5)==rhs
    lam=coeff[3]
    poly=scale(add(sub(k(51),scale(lam,13)),sub(scale(power(lam,2),3),power(lam,3))),Q(1,30))
    kap=k(Q(-2,5),Q(2,5),Q(-3,10),Q(3,10))
    assert mul(R,mul(sub(ONE,scale(a5,Q(1,25))),kap))==poly
    # Cross-conjugate rewriting of the individual-period candidate.
    assert scale(mul(sub(scale(R,7),scale(I,5)),kap),10)==k(125,-13,125,-41)
    # Euler factor (Rankin product in the norm proof): rational function
    # identity, separately for chi=+1 and -1, verified as finite polynomials.
    # For +1, (T/(1-aT+p2T2))'=(1-p2T2)/(1-aT+p2T2)^2.
    # For -1, even part of 1/(1-aT-p2T2) has numerator1-p2T2.
    # These follow by one differentiation and common-denominator addition;
    # verify_phase_and_euler.py supplies the independent symbolic derivation.
    beta=mul(sub(ONE,scale(conj(coeff[3]),Q(1,9))),sub(ONE,scale(conj(coeff[5]),Q(1,25))))
    delta=k(Q(25,288),Q(2,288),Q(5,288),Q(-1,288))
    assert scale(beta,Q(75,128))==mul(conj(delta),k(Q(15,2),Q(-3,2)))
    # Chuang Prop18: Euler polynomials at0, with E0=1 and
    # 2 E_n=-sum_{j<n} binom(n,j)E_j.
    E=[Q(1)]
    for n in range(1,6):E.append(-sum(Q(comb(n,j))*E[j] for j in range(n))/2)
    B=[[Q((-1)**(a+1)*factorial(5-a)*factorial(5-b),2*factorial(5)*factorial(5-a-b))*E[5-a-b] for a in range(3)] for b in range(3)]
    expected=[[Q(1,4),Q(0),Q(-1,8)],[Q(0),Q(1,10),Q(0)],[Q(-1,8),Q(0),Q(3,40)]]
    assert B==expected
    phases=(k(-1),I,ONE)
    H=[[Q(0) for j in range(3)] for i in range(3)]
    for i in range(3):
        for j in range(3):
            phase=mul(phases[i],phases[j])
            if B[i][j]:
                assert phase[1:]==(0,0,0)
                H[i][j]=B[i][j]/phase[0]
    Hi=invert(H)
    assert Hi==[[Q(24),Q(0),Q(-40)],[Q(0),Q(-10),Q(0)],[Q(-40),Q(0),Q(80)]]
    # Prop13: firstrowD=(0,0,8/15). The vector(15/8,0,0) is
    # the last row of D^{-1}, independently of the other entries.
    for d11,d12,d22 in ((0,0,0),(1,0,0),(0,1,0),(0,0,1)):
        D=[[Q(0),Q(0),Q(8,15)],[Q(0),Q(8,15),Q(d12)],[Q(8,15),Q(d12),Q(d22)]]
        assert mm([[Q(15,8),Q(0),Q(0)]],D)==[[Q(0),Q(0),Q(1)]]
    assert 2*Q(15,8)*32==120
    assert 120*Hi[1][1]/64==Q(-75,4)
    # detM3/pi^6=2/(3^(3/2)*5^(5/2))=2/(75sqrt15).
    assert Q(4,3**3*5**5)==Q(4,75**2*15)
    # Equating inverse entries gives De/B=pi^2/(2sqrt15).
    assert Q(75,4)*Q(2,75)==Q(1,2)
    # Distinguish the character in the ordinary norm Rankin series from
    # chi_-4 Sym^2. At p13, chi_-15=-1 but chi_-4=+1.
    assert pow((-15)%13,6,13)==12 and pow((-4)%13,6,13)==1
    out={'status':'PASS','carrier_bridge_through_q':36,'norm_and_numerical_form_identical':True,
         'cyclic_polynomial':'U^4-2U^2+81','Hecke_eigenchecks':144,
         'orthogonal_Galois_eigenlines':4,'depleted_twist_Sturm_bound':180,
         'twist_factor_and_cross_conjugate_rewriting':'exact',
         'Chuang_Betti_matrix':[[str(z) for z in row] for row in B],
         'dephased_H_inverse':[[str(z) for z in row] for row in Hi],
         'cofactor_coefficient':'-75/4','determinant_relation':'D_even=pi^2*D_odd/(2sqrt15)',
         'twisted_symmetric_square_value_evaluated':False,'Bessel_equalities_proved':False}
    (Path(__file__).resolve().parents[1]/'outputs'/'modular_bridge.json').write_text(json.dumps(out,indent=2)+'\n')
    print('EXACT BRIDGE PASS: same f; norm constant; depleted twist; Chuang cofactor.')
    print('Twisted symmetric-square critical value NOT evaluated.')
if __name__=='__main__':main()
