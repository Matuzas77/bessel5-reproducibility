#!/usr/bin/env python3
"""Rational, directed hypergeometric evaluation of A_j,C_j and pairing-derived B_j.
No Bessel, gamma, eta, or stored numerical period is evaluated or read.
The interval backend is attributed to the retained full-period-bridge release.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,comb
from fixed_interval import Context,ceildiv

ROOT=Path(__file__).resolve().parents[1]

def add_positive_tail(x,bound):
    ctx=x.ctx; up=ceildiv(bound.numerator*ctx.scale,bound.denominator)
    return ctx.from_bounds(x.lo,x.hi+up)

def hypergeom(ctx,N,den):
    t=ctx.ball(1); s=t
    for n in range(N):
        d=(n+1)**3 if den==1 else (n+1)*(n+2)**2
        t=t*ctx.ball(Q((2*n+1)**3,32*d));s=s+t
    # Next term <= t/4, and every subsequent ratio <=1/4.
    return ctx.from_bounds(s.lo,s.hi+ceildiv(t.hi,3))

def normalized_triples(ctx,N,HN):
    pi=ctx.pi();sq3=ctx.ball(3).sqrt()
    H0=hypergeom(ctx,HN,1); H1=hypergeom(ctx,HN,2)
    h=[sq3*pi**3*H0/8,sq3*pi**3*H1/1152]
    for n in range(N-1):
        aa=Q(20*n*n+40*n+21,18*(n+2)**2)
        bb=Q((2*n+1)**4,144*(n+1)**2*(n+2)**2)
        h.append(h[-1]*aa-h[-2]*bb)
    assert all(t.lo>0 for t in h), 'Insufficient arithmetic precision for minimal solution'
    return h,{'H0':H0.exact_json(),'H1':H1.exact_json()}

def C_tail(N,j):
    n=N+1;m=n+j
    # K0(t)^3 <= (pi/(2t))^(3/2)exp(-3t).
    # First omitted positive term <= pi^2 sqrt6/4 * rational.
    rational=Q(factorial(2*n)*factorial(4*m-2),
               4**(n+2*m-1)*factorial(n)**4*factorial(2*m-1)*9**m)
    # For j<=2,n>=20 the bound-term ratios are <1/2.
    assert 2*n+1>0 and n>=20 and 0<=j<=2
    assert Q((2*n+1)*(16*(n+j)**2-1),72*(n+1)**3)<Q(1,2)
    # 2*(pi^2 sqrt6/4)<25/2.
    return Q(25,2)*rational

def domb(N):
    D=[1,4]
    for n in range(2,N+1):
        num=2*(2*n-1)*(5*n*n-5*n+2)*D[-1]-64*(n-1)**3*D[-2]
        assert num % n**3==0
        D.append(num//n**3)
    return D[:N+1]

def deriv_factor(n,j):
    p=1
    for k in range(j):p*=(2*n+1-2*k)**2
    return p

def cross(x,y):
    return [x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0]]

def G_times(v):
    return [v[0]/75-Q(121,75)*v[1]+v[2],-Q(121,75)*v[0]+v[1],v[0]]

def main():
    p=argparse.ArgumentParser();p.add_argument('--digits',type=int,default=900)
    p.add_argument('--N',type=int,default=850);p.add_argument('--hyper-N',type=int,default=1800)
    p.add_argument('--out',default='source_coarse.json');ar=p.parse_args()
    if ar.N<20 or ar.digits<50:raise ValueError('N>=20 and digits>=50 required')
    start=time.monotonic();ctx=Context(ar.digits);pi=ctx.pi();sq15=ctx.ball(15).sqrt()
    h,hinfo=normalized_triples(ctx,ar.N+2,ar.hyper_N)
    C=[];tails=[];V=[]
    Ds=domb(ar.N)
    for j in range(3):
        z=ctx.ball(0); vv=ctx.ball(0)
        for n in range(ar.N+1):
            fac=1
            for k in range(1,j+1):fac*=(n+k)**2
            z=z+h[n+j]*(4**j*comb(2*n,n)*fac)
            den=4**n*(2*n+1)**2*comb(2*n,n)**2
            vv=vv+ctx.ball(Q(Ds[n]*deriv_factor(n,j),den))
        tail=C_tail(ar.N,j);C.append(add_positive_tail(z,tail));tails.append(str(tail))
        # V_j=(D_c^j StruveCompanion)(1), NOT the uncorrected raw weighted Struve.
        vt=Q((2*ar.N+3)**(2*j),4**(ar.N+1))
        V.append(add_positive_tail(pi*vv/8,vt))
    A=[pi**2*(2*C[j]-V[j]) for j in range(3)]
    D=[pi**2*c for c in C]
    e=[-2*sq15/pi**3*n for n in cross(G_times(A),G_times(D))]
    B=[x/pi for x in e]
    vals={'A':A,'B_pairing_derived':B,'C':C,'differentiated_Struve':V}
    assert all(x.lo>0 for row in [A,B,C] for x in row)
    def disp(x,n=70):
        s=10**n; lo=x.lo*s//ctx.scale;hi=ceildiv(x.hi*s,ctx.scale)
        def f(k):
            sign='-' if k<0 else '';a,b=divmod(abs(k),s);return f'{sign}{a}.{b:0{n}d}'
        return [f(lo),f(hi)]
    out={'settings':vars(ar),'periods':{k:[x.exact_json() for x in v] for k,v in vals.items()},
         'displays':{k:[disp(x) for x in v] for k,v in vals.items()},
         'hypergeometric_initial_intervals':hinfo,'C_tail_rational_bounds':tails,
         'all_errors_enclosed':True,'special_functions_evaluated':[],
         'A_and_C_use_pairing':False,'B_uses_inherited_pairing_and_orientation':True,
         'Bessel_eta_identity_proved':False}
    (ROOT/'outputs'/ar.out).write_text(json.dumps(out,indent=2)+'\n')
    print('SOURCE CERTIFICATE COMPLETE',ar.out,'seconds',round(time.monotonic()-start,2),flush=True)
    for k,row in vals.items():
        for j,x in enumerate(row):
            print(k,j,disp(x,55),'width_decimal_order',len(str(x.hi-x.lo))-ar.digits,flush=True)
if __name__=='__main__':main()
