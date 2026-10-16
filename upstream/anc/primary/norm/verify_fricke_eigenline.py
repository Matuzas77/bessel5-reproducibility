#!/usr/bin/env python3
"""Remove the database/newform-classification dependency from Fricke proportionality.

Integer/Fraction q arithmetic only. Together with the Riemann--Roch dimension
calculation in the report, this proves the required good-Hecke eigenspace at
level 60 is one dimensional. The exact Fricke scalar is checked separately by
retained Dedekind-eta transformation code.
"""
from fractions import Fraction as Q
from math import gcd
from pathlib import Path
import json,time
from verify_norm_independent import eta,sigma_table,modular_eta_check,EXPS,FTERMS,D,chi

ALL=dict(EXPS)
ALL.update({10:(0,0,0,3,0,0,0,0,0,0,0,3),
11:(0,0,0,0,0,0,0,3,0,3,0,0),
16:(-1,0,3,1,0,0,0,0,0,3,0,0),
17:(2,-1,0,0,0,2,0,0,2,0,1,0)})

def ch(d,n):
    if d==1:return 1
    if d==-15:return chi(n)
    if d==-3:return 0 if n%3==0 else (1 if n%3==1 else -1)
    if d==5:return 0 if n%5==0 else (1 if n%5 in (1,4) else -1)
    raise ValueError(d)

def e3(phi,psi,dil,N):
    a=[Q(0)]*(N+1)
    if phi==1:
        cond=abs(psi)
        B3=lambda x:x*x*x-Q(3,2)*x*x+Q(1,2)*x
        a[0]=-Q(cond*cond,6)*sum(ch(psi,j)*B3(Q(j,cond)) for j in range(1,cond+1))
    for d in range(1,N//dil+1):
        z=ch(psi,d)*d*d
        if z:
            for k in range(1,N//(d*dil)+1):a[d*k*dil]+=z*ch(phi,k)
    return a

def transpose(A):return list(map(list,zip(*A)))
def rref(A):
    A=[[Q(v) for v in row] for row in A];m=len(A);n=len(A[0]);piv=[];r=0
    for c in range(n):
        p=next((i for i in range(r,m) if A[i][c]),None)
        if p is None:continue
        A[r],A[p]=A[p],A[r];q=A[r][c];A[r]=[v/q for v in A[r]]
        for i in range(m):
            if i!=r and A[i][c]:
                q=A[i][c];A[i]=[x-q*y for x,y in zip(A[i],A[r])]
        piv.append(c);r+=1
        if r==m:break
    return A,piv

def mm(A,B):
    BT=transpose(B)
    return [[sum((x*y for x,y in zip(row,col)),Q(0)) for col in BT] for row in A]
def ident(n):return [[Q(int(i==j)) for j in range(n)] for i in range(n)]
def mv(A,v):return [sum((x*y for x,y in zip(row,v)),Q(0)) for row in A]

def run():
    start=time.monotonic();sig=sigma_table(2916)
    forms=[eta(ALL[k],252,sig) for k in range(18)]
    for k in range(18):modular_eta_check(ALL[k])
    for ph,ps in [(1,-15),(-15,1),(-3,5),(5,-3)]:
        for dil in (1,2,4):forms.append(e3(ph,ps,dil,252))
    M=[[b[n] for b in forms] for n in range(37)]
    _,rows=rref(transpose(M));assert len(rows)==30
    aug=[M[i]+ident(30)[j] for j,i in enumerate(rows)]
    red,piv=rref(aug);assert piv==list(range(30));inv=[row[30:] for row in red]
    cols=[]
    for b in forms:
        tar=[b[7*n]-(49*b[n//7] if n%7==0 else 0) for n in range(37)]
        col=mv(inv,[tar[i] for i in rows]);assert mv(M,col)==tar;cols.append(col)
    T=transpose(cols);P=mm(T,T)
    for i in range(30):P[i][i]+=64
    _,piv=rref(P);assert len(piv)==26
    F=[0]*2917
    for c,ds in FTERMS:
        v=eta(tuple(ds.get(d,0) for d in D),2916,sig)
        F=[x+c*y for x,y in zip(F,v)]
    assert all(F[81*n]-2*F[9*n]+81*F[n]==0 for n in range(37))
    assert all(F[7*n]-(49*F[n//7] if n%7==0 else 0)==Q(2,9)*(7*F[3*n]+F[27*n]) for n in range(37))
    assert all(F[11*n]-(121*F[n//11] if n%11==0 else 0)==F[n]-F[9*n] for n in range(37))
    _,pivW=rref([[F[3**j*n] for j in range(4)] for n in range(1,37)])
    assert len(pivW)==4
    # Torsion-free Gamma0(60), 12 singular cusps, genus 7.
    # Check singularity of the character for each denominator and each numerator.
    cusp_count=0
    for c in D:
        h=gcd(c,60//c)
        cusp_count+=sum(gcd(a,h)==1 for a in range(1,h+1))
        w=60//gcd(c*c,60)
        for a in range(1,2*c+1):
            if gcd(a,c)==1:
                dd=1+a*c*w;assert gcd(dd,60)==1 and chi(dd)==1
    assert cusp_count==12
    index=60*(Q(3,2)*Q(4,3)*Q(6,5));assert index==144
    assert not any((x*x+1)%60==0 for x in range(60))
    assert not any((x*x+x+1)%60==0 for x in range(60))
    genus=1+index/12-Q(cusp_count,2);assert genus==7
    weight3_degree=3*index/12
    dimM=weight3_degree+1-genus;assert dimM==30
    out={'status':'PASS','coefficient_arithmetic':'integer eta logarithmic recurrence and Fraction row reduction','index':int(index),'cusps':cusp_count,'all_cusps_character_singular':True,'genus':int(genus),'weight_three_line_bundle_degree':int(weight3_degree),'full_modular_space_dimension':int(dimM),'independent_modular_functions':30,'T7_squared_plus64_rank':26,'kernel_dimension':4,'U3_cyclic_dimension':4,'U3_polynomial':'U^4-2U^2+81','T7_on_cyclic_space':'(2/9)*(7*U+U^3)','T11_on_cyclic_space':'1-U^2','conclusion':'Each joint (T7,T11) eigenspace belonging to the four conjugates of f is one dimensional in the full level-60 space. Fricke conjugation is therefore proportional to coefficient conjugation, without a database newform assumption.','elapsed_seconds':round(time.monotonic()-start,3)}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':run()
