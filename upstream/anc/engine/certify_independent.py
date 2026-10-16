#!/usr/bin/env python3
"""Independent binary-directed source certificate.
No imports from primary implementation/backend; standard library only.
Uses direct central-binomial hypergeometric terms, the T_n/(2n)! recurrence,
explicit Domb binomial sums, and conic/linear reconstruction for the B row.
"""
from __future__ import annotations
import argparse,json,time,sys
from fractions import Fraction
from math import comb,factorial,isqrt
from pathlib import Path
sys.set_int_max_str_digits(100000)
ROOT=Path(__file__).resolve().parents[1]
S=1

def up(a,b):
    if b<=0: raise ValueError('positive divisor required')
    return -((-a)//b)

def exact(q=0):
    q=Fraction(q);return (q.numerator*S//q.denominator,up(q.numerator*S,q.denominator))
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def scale(a,q):
    q=Fraction(q);a=a if q>=0 else neg(a);q=abs(q)
    return(a[0]*q.numerator//q.denominator,up(a[1]*q.numerator,q.denominator))
def mul(a,b):
    v=[x*y for x in a for y in b];return(min(v)//S,up(max(v),S))
def div(a,b):
    if b[0]<=0<=b[1]:raise ArithmeticError('division through zero')
    if b[1]<0:return div(neg(a),neg(b))
    return mul(a,(S*S//b[1],up(S*S,b[0])))
def sqrt(a):
    if a[0]<0:raise ArithmeticError('sqrt interval not nonnegative')
    lo=isqrt(a[0]*S);hi=isqrt(a[1]*S)
    if hi*hi<a[1]*S:hi+=1
    return(lo,hi)
def power(a,n):
    z=exact(1)
    for _ in range(n):z=mul(z,a)
    return z

def atan_inv(q):
    total=exact(0);qp=q;n=0
    while True:
        term=exact(Fraction(1,(2*n+1)*qp))
        total=add(total,term if n%2==0 else neg(term))
        n+=1;qp*=q*q
        den=(2*n+1)*qp
        if den>S:
            # Alternating-series error has the sign of the next term.
            e=up(S,den)
            return (total[0],total[1]+e) if n%2==0 else(total[0]-e,total[1])

def hypergeom_direct(M):
    total0=exact(0);total1=exact(0);d=1
    for n in range(M+1):
        c=comb(2*n,n);q=Fraction(c**3,d)
        total0=add(total0,exact(q));total1=add(total1,exact(q/(n+1)**2));d*=256
    tail0=exact(q/3)[1];tail1=exact(q/(3*(M+1)**2))[1]
    return (total0[0],total0[1]+tail0),(total1[0],total1[1]+tail1)

def domb_explicit(N):
    central=[comb(2*k,k) for k in range(N+1)];D=[]
    for n in range(N+1):
        v=0;bk=1
        for k in range(n//2+1):
            t=bk*bk*central[k]*central[n-k]
            v+=t if 2*k==n else 2*t
            bk=bk*(n-k)//(k+1)
        D.append(v)
    return D

def C_tail_by_ratio(N,j):
    n=N+1;m=n+j
    # Same proved analytic K bound; separately written exact first term and ratio.
    f=Fraction(comb(2*n,n)*factorial(4*m-2),
               4**(n+2*m-1)*factorial(n)**2*factorial(2*m-1)*9**m)
    # pi^2*sqrt(6)/4 <25/4 and a geometric ratio at most1/2.
    return Fraction(25,2)*f

def Gv(a):
    return [add(sub(scale(a[0],Fraction(1,75)),scale(a[1],Fraction(121,75))),a[2]),
            add(scale(a[0],Fraction(-121,75)),a[1]),a[0]]

def encode(a,places):
    D=10**places;l=a[0]*D//S;h=up(a[1]*D,S)
    return {'lower_integer':str(l),'upper_integer':str(h),'denominator_power10':places,'width_integer':str(h-l)}

def main():
    global S
    ap=argparse.ArgumentParser();ap.add_argument('--bits',type=int,default=4200)
    ap.add_argument('--N',type=int,default=1050);ap.add_argument('--hyper-N',type=int,default=2600)
    ap.add_argument('--out',default='source_independent.json');args=ap.parse_args()
    if args.N<20 or args.bits<100:raise ValueError('N>=20 and bits>=100 required')
    S=2**args.bits;start=time.monotonic()
    pi=sub(scale(atan_inv(5),16),scale(atan_inv(239),4));pi2=mul(pi,pi);pi3=mul(pi2,pi)
    sq3=sqrt(exact(3));H0,H1=hypergeom_direct(args.hyper_N)
    common=mul(sq3,pi3)
    # r_n=T_n/(2n)!, distinct from the primary normalization.
    r=[scale(mul(common,H0),Fraction(1,8)),scale(mul(common,H1),Fraction(1,576))]
    for n in range(args.N+1):
        a=Fraction(2*(20*n*n+40*n+21),9*(2*n+4)*(2*n+3))
        b=Fraction((2*n+1)**3,9*(2*n+2)*(2*n+3)*(2*n+4))
        r.append(sub(scale(r[-1],a),scale(r[-2],b)))
    if any(x[0]<=0 for x in r):raise ArithmeticError('precision not sufficient')
    Ds=domb_explicit(args.N);C=[];W=[]
    for j in range(3):
        cv=exact(0);wv=exact(0)
        for n in range(args.N+1):
            cc=comb(2*n,n);fac=1
            for k in range(1,2*j+1):fac*=2*n+k
            cv=add(cv,scale(r[n+j],Fraction(cc*cc*fac,4**n)))
            pp=1
            for k in range(j):pp*=(2*n+1-2*k)**2
            wv=add(wv,exact(Fraction(Ds[n]*pp,4**n*(2*n+1)**2*cc**2)))
        tail=exact(C_tail_by_ratio(args.N,j))[1]
        C.append((cv[0],cv[1]+tail))
        wv=scale(mul(pi,wv),Fraction(1,8))
        wt=exact(Fraction((2*args.N+3)**(2*j),4**(args.N+1)))[1]
        W.append((wv[0],wv[1]+wt))
    A=[mul(pi2,sub(scale(C[j],2),W[j])) for j in range(3)]
    D=[mul(pi2,c)for c in C]
    # Independent B recovery: conic at the first column, then two linear equations.
    rad=scale(add(sub(scale(mul(A[0],A[0]),12),scale(mul(A[0],D[0]),40)),scale(mul(D[0],D[0]),40)),Fraction(1,5))
    E0=sqrt(rad);GA,GD=Gv(A),Gv(D)
    den=sub(mul(GA[1],GD[2]),mul(GA[2],GD[1]))
    E1=div(mul(E0,sub(mul(GA[2],GD[0]),mul(GA[0],GD[2]))),den)
    E2=div(mul(E0,sub(mul(GA[0],GD[1]),mul(GA[1],GD[0]))),den)
    B=[div(x,pi)for x in [E0,E1,E2]]
    if any(x[0]<=0 for x in A+B+C):raise ArithmeticError('positivity not established')
    places=args.bits*3//10-20
    out={'settings':vars(args),'binary_denominator_exponent':args.bits,
         'periods':{k:[encode(x,places)for x in row]for k,row in [('A',A),('B_pairing_derived',B),('C',C),('differentiated_Struve',W)]},
         'independent_of_primary_code_and_backend':True,'Domb_algorithm':'explicit binomial sum at every index',
         'hypergeom_algorithm':'closed central-binomial terms',
         'recurrence_normalization':'T_n/(2n)!','B_reconstruction':'positive conic then two linear orthogonality equations',
         'all_errors_enclosed':True,'Bessel_eta_identity_proved':False}
    (ROOT/'outputs'/args.out).write_text(json.dumps(out,indent=2)+'\n')
    print('INDEPENDENT BINARY SOURCE CERTIFICATE PASSED',args.out,'seconds',round(time.monotonic()-start,3))
if __name__=='__main__':main()
