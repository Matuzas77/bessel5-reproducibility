#!/usr/bin/env python3
"""Independent binary-directed eta calculation in Q=q^2.
Imports the independently written binary arithmetic, not the retained primary
interval backend or eta code. A different Cauchy radius (2/3) is used.
The modular folding theorem and original exponents are inherited inputs.
"""
from __future__ import annotations
import argparse,json,time,hashlib
from math import ceil
from fractions import Fraction as R
from pathlib import Path
import certify_independent as b
ROOT=Path(__file__).resolve().parents[1]

def exp_positive(x):
    assert 0<x[0]<=x[1]<b.S
    total=b.exact(1);term=b.exact(1);n=0
    while True:
        n+=1;term=b.scale(b.mul(term,x),R(1,n));total=b.add(total,term)
        nxt=b.scale(b.mul(term,x),R(1,n+1))
        if nxt[1]<=2:return (total[0],total[1]+2*nxt[1])
        if n>10000:raise ArithmeticError('exponential stopping limit')

def euler_even(E,N):
    shift=sum(2*d*e for d,e in E.items())//24
    assert shift>0 and sum(2*d*e for d,e in E.items())==shift*24
    L=(N-shift)//2;cs=[0]*(L+1);cs[0]=1
    for d,e in E.items():
        for k in range(d,L+1,d):
            if e>0:
                for _ in range(e):
                    for n in range(L,k-1,-1):cs[n]-=cs[n-k]
            elif e<0:
                for _ in range(-e):
                    for n in range(k,L+1):cs[n]+=cs[n-k]
    full=[0]*(N+1)
    for n,v in enumerate(cs):full[shift+2*n]=v
    return full

def M_bound(E):
    r=R(2,3);x=0
    for d,e in E.items():
        v=r**(2*d)
        x+=e*v/(1-v) if e>0 else -e*v/(1-v)**2
    return 3**ceil(x)

def dot(a,c):
    total=b.exact(0)
    for x,y in zip(a,c):total=b.add(total,b.mul(x,y))
    return total

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--N',type=int,default=2600)
    ap.add_argument('--bits',type=int,default=3000);ap.add_argument('--out',default='eta_independent.json');args=ap.parse_args()
    b.S=2**args.bits;start=time.monotonic()
    pi=b.sub(b.scale(b.atan_inv(5),16),b.scale(b.atan_inv(239),4))
    r=b.sqrt(b.exact(5));r60=b.sqrt(b.exact(60));alpha=b.div(r60,pi)
    q=b.div(b.exact(1),exp_positive(b.div(b.scale(pi,2),r60)))
    assert q[1]*20<9*b.S
    qp=[b.exact(1)]
    for n in range(args.N):qp.append(b.mul(qp[-1],q))
    Uv=[];Vv=[];allarrays=[];majorants=[]
    for j in range(4):
        E1={1:1+3*j,3:2-3*j,5:4-3*j,15:-1+3*j}
        E3={1:-1+3*j,3:4-3*j,5:2-3*j,15:1+3*j}
        aa=euler_even(E1,args.N);bb=euler_even(E3,args.N);allarrays.extend([aa,bb])
        cc=[x+3*y for x,y in zip(aa,bb)];dd=[x-3*y for x,y in zip(aa,bb)]
        M=M_bound(E1)+3*M_bound(E3);majorants.append(M)
        tail=R(40*M,13)*R(27,40)**(args.N+1)
        sums=[]
        for arr in [cc,dd]:
            s1=b.exact(0);s2=b.exact(0)
            for n in range(1,args.N+1):
                if arr[n]:
                    s1=b.add(s1,b.scale(qp[n],R(arr[n],n)))
                    s2=b.add(s2,b.scale(qp[n],R(arr[n],n*n)))
            e1=b.exact(tail/(args.N+1))[1];e2=b.exact(tail/(args.N+1)**2)[1]
            sums.append(((s1[0]-e1,s1[1]+e1),(s2[0]-e2,s2[1]+e2)))
        Uv.append(b.add(b.scale(sums[0][0],4),b.mul(alpha,sums[0][1])))
        Vv.append(b.mul(alpha,sums[1][1]))
    U=b.add(Uv[0],b.mul(b.scale(b.sub(r,b.exact(1)),R(1,2)),Uv[1]))
    V=b.sub(Vv[0],b.mul(b.add(r,b.exact(2)),Vv[1]))
    cu=[b.scale(r,R(-9,20)),b.scale(r,R(31,20)),b.sub(b.exact(R(1,4)),b.scale(r,R(33,10))),b.scale(r,R(-9,20))]
    cv=[b.scale(r,R(-13,20)),b.scale(b.add(b.scale(r,5),b.exact(7)),R(-1,4)),
        b.scale(b.add(b.exact(105),b.scale(r,61)),R(1,40)),b.scale(b.add(b.exact(5),b.scale(r,3)),R(9,40))]
    B=dot(cu,Uv);D=dot(cv,Vv);places=args.bits*3//10-20
    # Compare digest through the existing primary fine cutoff, without reading its output.
    prefix=[]
    for arr in allarrays:prefix.extend(arr[:2401])
    out={'settings':vars(args),'quantities':{k:b.encode(v,places)for k,v in [('U',U),('V',V),('companion_U',B),('companion_V_without_i',D)]},
         'moments_U':[b.encode(v,places)for v in Uv],'moments_V':[b.encode(v,places)for v in Vv],
         'Cauchy_radius':'2/3','q_upper_bound':'9/20','geometric_ratio':'27/40','component_majorants':majorants,
         'coefficient_prefix_2400_digest':hashlib.sha256(json.dumps(prefix).encode()).hexdigest(),
         'coefficient_prefix_entries':len(prefix),'all_period_errors_enclosed':True,
         'independent_of_primary_code_and_backend':True,'eta_identity_proved':False}
    (ROOT/'outputs'/args.out).write_text(json.dumps(out,indent=2)+'\n')
    print('INDEPENDENT ETA CERTIFICATE PASSED',args.out,'seconds',round(time.monotonic()-start,3))
if __name__=='__main__':main()
