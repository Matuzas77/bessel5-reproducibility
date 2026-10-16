"""Exact first Fourier coefficient of f|_3 W_60 via Dedekind eta transformations.
Cyclotomic arithmetic is over Q[z]/(z^216-z^108+1), z=exp(2pi i/648).
No decimal period or numerical phase enters this calculation.
"""
from fractions import Fraction as Q
from math import gcd
from pathlib import Path
import json,time
import sympy as sp
from recompute_modular_value import TERMS
BASE=Path(__file__).resolve().parent;ORDER=648

def saw(q):return Q(0) if q.denominator==1 else q-q.numerator//q.denominator-Q(1,2)
def dedekind(d,c):return sum((saw(Q(n,c))*saw(Q(d*n,c)) for n in range(1,c)),Q(0))
def sqrt_factor(q):
    # all eta terms here have only a possible sqrt(3) factor after W60 normalization
    factors=dict(sp.factorint(q.numerator))
    for p,e in sp.factorint(q.denominator).items():factors[p]=factors.get(p,0)-e
    coef=Q(1);rad=1
    for p,e in factors.items():
        k=e//2;coef*=Q(int(p))**int(k)
        if e%2:rad*=int(p)
    assert rad in (1,3),(q,rad)
    return coef,rad

def reduce_poly(p):
    v=[Q(0)]*ORDER
    for k,c in p.items():v[k%ORDER]+=c
    for k in range(ORDER-1,215,-1):
        v[k-108]+=v[k];v[k-216]-=v[k];v[k]=Q(0)
    return {k:c for k,c in enumerate(v[:216]) if c}
def add(p,q):
    out=dict(p)
    for k,c in q.items():out[k]=out.get(k,Q(0))+c
    return {k:c for k,c in out.items() if c}
def shift(p,n,coef=Q(1)):return reduce_poly({(k+n)%ORDER:coef*c for k,c in p.items()})

def first_U_fricke(j):
    m=3**j;ans={}
    for b in range(m):
        for coefficient,exps in TERMS:
            params=[];initial=Q(0);phase=Q(0);root=Q(60)**3
            for d,r in exps.items():
                g=gcd(d*b,m);a=d*b//g;c=m//g;d0=pow(a,-1,c) if c>1 else 0
                alpha=Q(60*g*g,m*d);beta=-Q(d0,c)
                initial+=r*alpha/24
                phase+=r*(Q(a,12*c)-dedekind(d0,c)-Q(1,4))
                root*=Q(g,d)**r
                params.append((alpha,beta,r))
            if initial>1:continue
            st=initial*ORDER;assert st.denominator==1
            ph=phase*ORDER/2;assert ph.denominator==1
            states={(int(st),0):1}
            for alpha,beta,r in params:
                stride=alpha*ORDER;assert stride.denominator==1
                stride=int(stride);terms=[(0,0,1)];n=1
                while stride*n*(3*n-1)//2<=ORDER-int(st):
                    for k in [n*(3*n-1)//2,n*(3*n+1)//2]:
                        if stride*k<=ORDER-int(st):
                            z=beta*k*ORDER;assert z.denominator==1
                            terms.append((stride*k,int(z)%ORDER,(-1)**n))
                    n+=1
                for _ in range(r):
                    nxt={}
                    for (qpower,zpower),c0 in states.items():
                        for qp,zp,sg in terms:
                            if qpower+qp<=ORDER:
                                key=(qpower+qp,(zpower+zp)%ORDER)
                                nxt[key]=nxt.get(key,0)+c0*sg
                    states={k:c for k,c in nxt.items() if c}
            v={zp:Q(c0) for (qp,zp),c0 in states.items() if qp==ORDER}
            factor,rad=sqrt_factor(root);v=shift(v,int(ph),factor*coefficient/m)
            if rad==3:v=add(shift(v,54),shift(v,-54))
            ans=add(ans,v)
    return reduce_poly(ans)

def main():
    start=time.time();wj=[first_U_fricke(j) for j in range(4)]
    cs=[(Q(5,12),Q(1,20),Q(1),Q(19,60)),(-Q(7,36),-Q(13,540),Q(29,108),Q(11,90)),(-Q(1,12),-Q(1,20),Q(0),Q(1,60)),(-Q(1,36),-Q(7,540),-Q(1,108),-Q(1,90))]
    rational={};root5={}
    for c,w in zip(cs,wj):
        rational=add(rational,add(shift(w,0,c[0]),shift(w,162,c[2])))
        root5=add(root5,add(shift(w,0,c[1]),shift(w,162,c[3])))
    rational=reduce_poly(rational);root5=reduce_poly(root5)
    # lambda=i epsilon = (10-2sqrt5 + i(5+4sqrt5))/15
    expect0=add({0:Q(2,3)},{162:Q(1,3)})
    expect5=add({0:Q(-2,15)},{162:Q(4,15)})
    assert rational==expect0,(rational,expect0)
    assert root5==expect5,(root5,expect5)
    out={'status':'PASS','U3_orbit_Frcke_q1_coefficients':[{str(k):str(c) for k,c in w.items()} for w in wj],'f_Frcke_q1_rational':{str(k):str(c) for k,c in rational.items()},'f_Frcke_q1_sqrt5':{str(k):str(c) for k,c in root5.items()},'lambda':'(10-2sqrt(5)+i*(5+4sqrt(5)))/15','epsilon':'(5+4sqrt(5)+i*(2sqrt(5)-10))/15','method':'Exact Dedekind eta transformation and cyclotomic polynomial reduction','elapsed':time.time()-start,'scope':'Exact q coefficient of the Fricke transform. The current package proves proportionality by the full-space T7,T11 eigenline calculation, so no database/newform-classification assumption is needed. No Bessel identity inferred.'}
    (BASE/'verify_fricke_exact.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
