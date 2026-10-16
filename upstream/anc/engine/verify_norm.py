#!/usr/bin/env python3
"""Exact independent verification of the level-60 norm certificate.

Standard library only. The eta series are constructed by logarithmic
 differentiation, not by the product multiplication used to find the witness.
No Bessel or L-value data are read. q^0 through q^60 are checked; q^36 is
already the weight-three level-60 Sturm bound.
"""
from fractions import Fraction as Q
from math import gcd
from pathlib import Path
import json
import time

D=(1,2,3,4,5,6,10,12,15,20,30,60)
# Explicit eta functions needed by the witness (rows of the independently
# modular eta table, omitting zero witness coefficients).
EXPS={
0:(3,0,0,0,0,0,0,0,3,0,0,0),
1:(2,0,0,1,0,0,0,0,2,0,0,1),
2:(2,0,0,0,0,1,1,0,2,0,0,0),
3:(1,0,0,2,0,0,0,0,1,0,0,2),
4:(1,0,0,1,0,1,1,0,1,0,0,1),
5:(1,0,0,0,0,2,2,0,1,0,0,0),
6:(0,2,0,0,0,0,0,1,0,1,2,0),
7:(0,1,0,0,0,0,0,2,0,2,1,0),
8:(0,0,3,0,3,0,0,0,0,0,0,0),
9:(0,0,2,0,2,0,0,1,0,1,0,0),
12:(-1,3,3,0,0,0,1,0,0,0,0,0),
13:(-1,2,3,0,0,1,0,0,0,0,1,0),
14:(-1,1,4,0,1,1,0,0,0,0,0,0),
15:(-1,1,3,0,0,1,0,1,0,1,0,0)}
RCO={0:Q(47,32256),1:Q(-9883,48384),2:Q(11779,48384),
3:Q(-1367,3024),4:Q(18059,24192),5:Q(-505,1344),
6:Q(-329,1728),7:Q(-137,288),8:Q(47,1792),9:Q(83,672),
12:Q(-2,315),13:Q(-61,1890),14:Q(-1,1260),15:Q(2,21)}
FTERMS=[(-8,{12:3,20:3}),(8,{4:1,6:2,10:2,60:1}),
(-4,{4:2,6:1,10:1,60:2}),(-6,{3:1,5:1,12:2,20:2}),
(3,{3:2,5:2,12:1,20:1}),(1,{3:3,5:3}),
(-3,{2:1,3:2,5:2,30:1}),(2,{2:2,3:1,5:1,30:2})]

def chi(n):
    if gcd(n,15)>1:return 0
    return 1 if n%15 in (1,2,4,8) else -1

def mobius(d):
    sign=1
    for p in (2,3,5):
        if d%p==0:
            if d%(p*p)==0:return 0
            sign=-sign
    return sign

def sigma_table(N):
    a=[0]*(N+1)
    for d in range(1,N+1):
        for n in range(d,N+1,d):a[n]+=d
    return a

def eta(rs,N,sigma):
    shift=sum(d*r for d,r in zip(D,rs))//24
    if shift<0:raise ValueError('negative eta q-order')
    L=N-shift
    if L<0:return [0]*(N+1)
    logs=[0]*(L+1)
    for d,r in zip(D,rs):
        if r:
            for j in range(1,L//d+1):logs[d*j]-=r*d*sigma[j]
    a=[1]+[0]*L
    for n in range(1,L+1):
        z=sum(logs[j]*a[n-j] for j in range(1,n+1))
        assert z%n==0, ('nonintegral eta coefficient',n)
        a[n]=z//n
    return [0]*shift+a

def modular_eta_check(rs):
    assert sum(rs)==6
    assert sum(d*r for d,r in zip(D,rs))%24==0
    assert sum((60//d)*r for d,r in zip(D,rs))%24==0
    # Character squareclass is 15, including negative exponents.
    parity=[]
    for p in (2,3,5):
        parity.append(sum(r*valuation(d,p) for d,r in zip(D,rs))%2)
    assert parity==[0,1,1]
    orders=[]
    for c in D:
        o=Q(60,24)*sum(Q(gcd(c,d)**2*r,gcd(c,60//c)*c*d) for d,r in zip(D,rs))
        assert o>0
        orders.append(str(o))
    return orders

def valuation(n,p):
    k=0
    while n%p==0:n//=p;k+=1
    return k

# Biquadratic arithmetic in basis 1,sqrt(5),i,i*sqrt(5).
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def mul(a,b):
    z=[Q(0)]*4
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            z[i^j]+=x*y*(5 if i&j&1 else 1)*(-1 if i&j&2 else 1)
    return tuple(z)
def scale(a,c):return tuple(x*c for x in a)
def conj(a):return (a[0],a[1],-a[2],-a[3])
ONE=(Q(1),Q(0),Q(0),Q(0))

def run():
    start=time.monotonic();B=60;N=49*B;sigma=sigma_table(N)
    R=[Q(0)]*(N+1); cusp={}
    for k,rs in EXPS.items():
        cusp[str(k)]=modular_eta_check(rs)
        a=eta(rs,N,sigma)
        c=RCO[k]
        R=[x+c*y for x,y in zip(R,a)]
    # E_3(1,chi_-15) = -8 + sum_{n>=1} sum_{d|n}chi(d)d^2 q^n.
    e=[Q(0)]*(N+1)
    # L(chi,-2)/2 from the cubic Bernoulli polynomial, not stored -8.
    B3=lambda x:x*x*x-Q(3,2)*x*x+Q(1,2)*x
    e[0]=-Q(15**2,6)*sum(chi(a)*B3(Q(a,15)) for a in range(1,16))
    assert e[0]==-8
    for d in range(1,N+1):
        z=chi(d)*d*d
        if z:
            for n in range(d,N+1,d):e[n]+=z
    for n in range(N+1):
        if n%2==0:R[n]+=Q(1,132608)*e[n//2]
        if n%4==0:R[n]-=Q(1,16576)*e[n//4]
    # Rebuild F independently by the logarithmic recurrence.
    F=[0]*(max(27*B,729)+1)
    for c,ds in FTERMS:
        rs=tuple(ds.get(d,0) for d in D);modular_eta_check(rs)
        a=eta(rs,len(F)-1,sigma)
        F=[x+c*y for x,y in zip(F,a)]
    coeffs=(Q(5,48),Q(-35,432),Q(-5,144),Q(-5,432))
    Hf=[sum(coeffs[j]*F[3**j*n] for j in range(4)) for n in range(B+1)]
    g=[Q(0)]*(B+1)
    g[0]=-sum(Q(chi(a)*a,15) for a in range(1,16))/2
    assert g[0]==1
    for d in range(1,B+1):
        for n in range(d,B+1,d):g[n]+=chi(d)
    E=[Q(0)]*(B+1)
    divisors30=[1,2,3,5,6,10,15,30]
    for d in divisors30:
        c=Q(25*mobius(d),16*d*d);E[0]+=c
        for n in range(1,B//(60//d)+1):E[n*(60//d)]-=24*c*sigma[n]
    assert E[0]==1
    H=[sum(g[j]*E[n-j] for j in range(n+1))-Q(n,12)*g[n] for n in range(B+1)]
    residual=[]
    for n in range(B+1):
        T2=R[49*n]-49*(1+int(n%7==0))*R[n]
        if n%49==0:T2+=2401*R[n//49]
        residual.append(H[n]-Hf[n]-T2-64*R[n])
    assert not any(residual)
    hk=[sum(coeffs[j]*F[3**(j+k)] for j in range(4)) for k in range(4)]
    assert hk==[Q(25,72),Q(5,18),Q(5,72),Q(5,9)]
    lam=(Q(0),Q(1),Q(-2),Q(0));lam2=mul(lam,lam);lam3=mul(lam2,lam)
    numerator=add(add(scale(ONE,hk[3]),scale(lam,hk[2])),add(scale(add(lam2,scale(ONE,-2)),hk[1]),scale(add(lam3,scale(lam,-2)),hk[0])))
    denominator=add(scale(lam3,4),scale(lam,-4))
    delta=(Q(25,288),Q(2,288),Q(5,288),Q(-1,288))
    assert mul(delta,denominator)==numerator
    beta=mul(add(ONE,scale((Q(0),Q(1),Q(2),Q(0)),Q(-1,9))),add(ONE,scale((Q(0),Q(1),Q(0),Q(-2)),Q(-1,25))))
    normconstant=(Q(15,2),Q(-3,2),Q(0),Q(0))
    assert scale(beta,Q(75,128))==mul(conj(delta),normconstant)
    # Exact absence of holomorphic-projection normalization ambiguity.
    assert Q(25,16)*sum(Q(mobius(d),d*d) for d in divisors30)==1
    assert Q(25,16)*sum(Q(mobius(d),d*60) for d in divisors30)==Q(1,144)
    result={'status':'PASS','method':'independent integer logarithmic eta recurrence and Fraction arithmetic','eta_witness_count':len(EXPS),'eta_carrier_count':len(FTERMS),'cusp_orders':cusp,'sturm_bound':36,'coefficient_checks_including_constant':B+1,'last_verified_q_power':B,'required_eta_coefficients':N,'witness_identity':'H-Hf=(T7^2+64)R','Hf_coefficients':list(map(str,coeffs)),'h_k':list(map(str,hk)),'delta_basis_1_sqrt5_i_i_sqrt5':list(map(str,delta)),'beta_basis':list(map(str,beta)),'norm_constant_basis':list(map(str,normconstant)),'arithmetic':'exact fractions'}
    path=Path(__file__).resolve().parents[1]/'outputs'/'norm_exact.json';path.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: exact norm witness q^0..q^60; eta modularity; eigenline projection; norm constant')
    print('Raw Petersson convention: no index normalization')
if __name__=='__main__':run()
