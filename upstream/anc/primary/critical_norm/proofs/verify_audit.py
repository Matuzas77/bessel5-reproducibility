"""Supplementary exact checks for the v4 critical-norm proof.
No period values are used. Finite checks supplement, not replace, the
all-index analytic arguments in the full proof.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd
import json, sys, sympy as sp
ROOT=Path(__file__).resolve().parents[1]
from kernel_coefficients import tau, chi4, chiD
from verify_exact import predicted_S_ratio

def main():
    records=[]
    # Direct integer cyclotomic reductions at higher bad-prime powers and
    # at discriminants zero, principal, and nonprincipal.
    cases=[(1,8,240),(1,62,480),(1,1,540),(7,5,1500),
           (3,2,2940),(1,2,420),(9,6,840),(2,3,600),
           (1,2,240),(5,2,480),(1,7,540),(19,1,1500),
           (1,5,2940),(1,3,420),(3,2,840),(11,2,600)]
    x=sp.Symbol('x')
    for n,l,c in cases:
        counts=[0]*c
        for r in range(c): counts[(r*r+l*r+n)%c]+=chi4(r)
        nonzero=[(j,a) for j,a in enumerate(counts) if a]
        raw=[0]*c
        for d in range(c):
            if gcd(d,c)!=1:continue
            char=chiD(60,d)
            for j,a in nonzero:raw[(d*j)%c]+=char*a
        ratio=predicted_S_ratio(n,l,c)
        for d in range(60):raw[(c//60)*d]-=ratio*chiD(60,d)
        p=sp.Poly.from_dict({(j,):a for j,a in enumerate(raw) if a},x,domain=sp.ZZ)
        assert p.rem(sp.Poly(sp.cyclotomic_poly(c,x),x)).is_zero
        records.append({'n':n,'l':l,'c':c,'Delta':l*l-4*n,'ratio_to_sqrt60':ratio})
    # Delta=0 is NOT passed through a fundamental-discriminant routine.
    # For n=s^2,l=2s with s odd, L2(u)=2 chi4(s)t(1-2t)/(1-t),
    # t=2^(-1-4u). For s even it is identically zero.
    boundary_checks=0
    for s in range(-6,7):
        for a in range(9):
            expected=0
            if s%2 and a>=2 and a%2==0:
                b=a//2
                expected=chi4(s)*(1 if b==1 else -1)*2**(b+1)
            assert tau(s*s,2*s,2,a)==expected,(s,a)
            boundary_checks+=1
    t=sp.Symbol('t')
    boundary=2*t*(1-2*t)/(1-t)
    assert boundary.subs(t,F(1,2))==0
    assert sp.diff(boundary,t).subs(t,F(1,2))==-4
    # Exact gamma-function coefficient comparison for Weber--Schafheitlin.
    u=sp.Symbol('u',real=True)
    rhs=sp.gamma(2-2*u)*sp.sin(sp.pi*u)/8
    dup=2**(-2-2*u)*sp.sqrt(sp.pi)*sp.gamma(F(3,2)-u)/sp.gamma(u)
    ratio=sp.combsimp(rhs/dup)
    assert sp.simplify(sp.expand_func(ratio)-1)==0,ratio
    # Formal FE and spectral-factor consistency, independent of decimals.
    r=sp.sqrt(5);I=sp.I
    delta=(7-I*r)*(7-I)/576
    R3=(1+(1+4*I*r)/27)*(1+(15+20*I)/125)
    assert sp.simplify(128*sp.conjugate(delta)/(75*R3)-F(1,8))==0
    out={'status':'PASS','higher_prime_power_cyclotomic_checks':records,
         'Delta_zero_character_sum_checks':boundary_checks,
         'Delta_zero_L2_factor':'2*chi4(s)*t*(1-2*t)/(1-t), t=2^(-1-4u), s odd; 0 for s even',
         'Delta_zero_derivative':'8*chi4(s)*log(2)',
         'Weber_Schafheitlin_coefficients':'exact gamma duplication/reflection check PASS',
         'L3_norm_scalar':'1/8',
         'scope':'finite analytic-normalization controls; no Bessel equality or formalized proof claimed'}
    target=ROOT/'outputs'/'audit.json';target.parent.mkdir(exist_ok=True)
    target.write_text(json.dumps(out,indent=2)+'\n')
    print('PASS 16 higher-prime-power cyclotomic checks; 117 Delta=0 checks; gamma and spectral scalars')
if __name__=='__main__': main()
