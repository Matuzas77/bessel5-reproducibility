#!/usr/bin/env python3
"""Exact phase, norm-scaling and Euler-factor identities.

These are algebraic identities in freely varying real X,Y. They do not
identify any Bessel integral with a modular value.
"""
from pathlib import Path
import json
import sympy as s


def main():
    r=s.sqrt(5); ii=s.I
    X,Y,Z=s.symbols('X Y Z', real=True)
    b0=15*s.sqrt(3)*(r-1)/2
    norm0=3*(5-r)/2
    alpha=((125-13*r)+ii*(125-41*r))/4
    Cnum=(42+6*r)*X+(33*r-69)*Y
    Abar=((125-13*r)*X+(41*r-125)*Y)/4 # A*/pi
    Cbar=Cnum/4 # pi*C*
    B=b0*Z
    zeta=(6+4*r+ii*(8-3*r))/15
    eps=(5+4*r+ii*(2*r-10))/15
    checks={}
    def check(name,z):
        assert s.simplify(s.expand_complex(z))==0,name
        checks[name]=True
    check('zeta_is_unit',zeta*s.conjugate(zeta)-1)
    check('epsilon_is_unit',eps*s.conjugate(eps)-1)
    check('norm_to_B_coefficient',b0/norm0-s.sqrt(15))
    check('phase_zeta_epsilon',zeta-eps*(2+ii)/r)
    check('phase_times_a3',zeta*(r-2*ii)/3-(4-3*ii)/5)
    lhs=2*Cbar/s.sqrt(3)+2*ii*s.sqrt(s.Rational(3,5))*(Abar-s.Rational(5,3)*Cbar)
    check('complex_period_identity',lhs-b0*zeta*(X-ii*Y))
    check('coefficient_gamma',b0*zeta-(s.sqrt(3)*(7+r)+ii*s.sqrt(3)*(11*r-23)/2))
    conic=12*Abar**2-40*Abar*Cbar+40*Cbar**2-5*B**2
    check('conic_coefficient',conic-s.Rational(3375,2)*(3-r)*(X**2+Y**2-Z**2))
    # Injectivity of the two real linear coordinates used in the phase.
    det=s.det(s.Matrix([[s.diff(s.re(lhs),X),s.diff(s.re(lhs),Y)],
                       [s.diff(s.im(lhs),X),s.diff(s.im(lhs),Y)]]))
    assert s.simplify(det)!=0
    checks['phase_real_linear_map_invertible']=True
    # Two possible unramified characters: recurrence b_j=1+chi+...+chi^j.
    a,p,T=s.symbols('a p T')
    for ch in (1,-1):
        A=1/(1-a*T+ch*p**2*T**2)
        direct=s.diff(T*A,T) if ch==1 else (A+A.subs(T,-T))/2
        proposed=(1-p**2*T**2)/((1-a*T+ch*p**2*T**2)*(1-ch*a*T+ch*p**2*T**2))
        assert s.cancel(direct-proposed)==0
        checks[f'Euler_factor_chi_{ch}']=True
    # The three omitted Euler factors of zeta(2).
    zz=s.Rational(1,6)*(1-s.Rational(1,4))*(1-s.Rational(1,9))*(1-s.Rational(1,25))
    assert zz==s.Rational(8,75)
    checks['zeta_imprimitive_constant']=True
    beta=(1-(r+2*ii)/9)*(1-r*(1-2*ii)/25)
    delta=(25+2*r+ii*(5-r))/288
    check('Rankin_norm_constant',s.Rational(75,128)*beta-s.conjugate(delta)*norm0)
    # Exact Maass-raising Stokes identity, using dbar(y)=i/2.
    ff,gg,gd,yy,pi=s.symbols('f gbar dbar_gbar y pi', nonzero=True)
    integrand=ff*yy*(-gd/(2*pi*ii)-gg/(4*pi*yy))
    derivative=-ff*(yy*gd+ii*gg/2)/(2*pi*ii)
    assert s.simplify(integrand-derivative)==0
    checks['raising_Stokes_identity']=True
    # Weight-one Eisenstein series equals an average of two binary theta series.
    # Both theta series have discriminant -15 and level 15; Sturm bound is 2.
    nmax=40
    theta=[0]*(nmax+1)
    for x in range(-10,11):
        for y in range(-10,11):
            for value in (x*x+x*y+4*y*y,2*x*x+x*y+2*y*y):
                if value<=nmax:theta[value]+=1
    def ch(n):return int(s.kronecker_symbol(-15,n))
    eisen=[2]+[2*sum(ch(d) for d in s.divisors(n)) for n in range(1,nmax+1)]
    assert theta==eisen
    checks['binary_theta_identity_q0_through40']=True
    out={'status':'PASS','checks':checks,'check_count':len(checks),
         'phase':'(6+4sqrt(5)+i(8-3sqrt(5)))/15',
         'scope':'Exact algebraic equivalence of the candidate triple to a norm-scale and a phase comparison. Neither Bessel comparison is proved.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
