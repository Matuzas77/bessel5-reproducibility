from pathlib import Path
ROOT=Path(__file__).resolve().parent
from fractions import Fraction as F
from math import gcd
from functools import lru_cache
import sympy as sp
import json, time
N=60
units=[i for i in range(N) if gcd(i,N)==1]
def chi(a):
    a%=15
    if gcd(a,15)!=1: return 0
    return (1 if a%3==1 else -1)*(1 if a%5 in (1,4) else -1)
raw=[(c,d) for c in range(N) for d in range(N) if gcd(gcd(c,d),N)==1]
can={}
for v in raw:
    if v in can: continue
    orb=[((t*v[0]%N,t*v[1]%N),t) for t in units]
    rep=min(x for x,t in orb)
    for t in units: can[(t*rep[0]%N,t*rep[1]%N)]=(rep,chi(t))
reps=sorted(set(rep for rep,f in can.values()))
ix={v:i for i,v in enumerate(reps)}
D=2*len(reps)

def add(v,i,c,d,a=1):
    rep,s=can[(c%N,d%N)]
    k=2*ix[rep]+i
    v[k]=v.get(k,0)+s*a
    if v[k]==0: del v[k]
    return v

def right(i,c,d,h):
    a,b,cc,dd=h
    out={}
    if i==0: vals=((1,cc),(0,dd))
    else: vals=((1,a),(0,b))
    for ii,v in vals:
        if v: add(out,ii,c*a+d*cc,c*b+d*dd,v)
    return out

def sumv(*vs):
    out={}
    for v in vs:
        for k,a in v.items():
            out[k]=out.get(k,0)+a
            if out[k]==0: del out[k]
    return out

def vec(v): return [v.get(i,0) for i in range(D)]
S=(0,-1,1,0); U=(0,-1,1,-1); U2=(-1,1,-1,0)
rels=[]
for (c,d) in reps:
 for i in (0,1):
    a=add({},i,c,d)
    rels.extend([sumv(a,right(i,c,d,S)),sumv(a,right(i,c,d,U),right(i,c,d,U2))])
print('D,relations',D,len(rels),flush=True)
t=time.time()
RR,pivs=sp.polys.matrices.DomainMatrix.from_Matrix(sp.Matrix([vec(v) for v in rels])).convert_to(sp.QQ).rref()
RR=RR.to_Matrix()
free=[i for i in range(D) if i not in pivs]
qD=len(free)
print('quotient',qD,'rref seconds',time.time()-t,flush=True)
red=sp.zeros(D,qD)
for j,k in enumerate(free): red[k,j]=1
for j,k in enumerate(pivs):
 for h,f in enumerate(free): red[k,h]=-RR[j,f]
assert sp.Matrix([vec(v) for v in rels])*red==sp.zeros(len(rels),qD)
print('free',free,flush=True)

def reduce(v):
    out=sp.zeros(1,qD)
    for k,a in v.items(): out+=a*red[k,:]
    return out

# elementary continued-fraction modular-symbol decomposition.
# Polynomial AX+BY. For a unimodular g=(a,b;c,d), g.P evaluated at tau is P(d tau-b,-c tau+a).
# A path {alpha,beta} with polynomial P is [-g^-1 P,g] from infinity to zero, as g{0,infty}={b/d,a/c}.

def xgcd(a,b):
    if b==0: return (abs(a),1 if a>=0 else -1,0)
    g,x,y=xgcd(b,a%b); return g,y,x-(a//b)*y

def lift(c,d):
    if gcd(c,d)!=1: raise ValueError((c,d))
    gg,x,y=xgcd(d,c) # d*x+c*y=1
    return (x,-y,c,d)

def mul(g,h):
 a,b,c,d=g;e,f,gc,hc=h
 return (a*e+b*gc,a*f+b*hc,c*e+d*gc,c*f+d*hc)
def mob(g,x):
 a,b,c,d=g
 if x is None:
    return None if c==0 else F(a,c)
 bot=c*x+d
 return None if bot==0 else (a*x+b)/bot

def leftpoly(g,P):
 a,b,c,d=g; A,B=P
 return (A*d-B*c,-A*b+B*a)
def inversepoly(g,P):
 a,b,c,d=g; A,B=P
 # g^-1 action for det1
 return (A*a+B*c,A*b+B*d)

@lru_cache(None)
def cf_edges(x):
 # edges from infinity to rational x
 if x is None:return ()
 num,den=x.numerator,x.denominator
 cf=[]
 while den:
    a=num//den;cf.append(a);num,den=den,num-a*den
 pm2,pm1=0,1; qm2,qm1=1,0
 edges=[]
 for a in cf:
    p=a*pm1+pm2; q=a*qm1+qm2
    # orient g{0,infty} from previous convergent to current
    # g = (p, previous_p; q, previous_q), scale first col to det1
    det=p*qm1-pm1*q
    assert abs(det)==1
    gg=(det*p,pm1,det*q,qm1)
    edges.append(gg)
    pm2,pm1=pm1,p; qm2,qm1=qm1,q
 return tuple(edges)

def from_inf(x,P):
    out={}
    for g in cf_edges(x):
       A,B=inversepoly(g,P)
       _,_,c,d=g
       add(out,1,c,d,A);add(out,0,c,d,B)
    return out

def symb(P,alpha,beta):
    neg={k:-a for k,a in from_inf(alpha,P).items()}
    return sumv(neg,from_inf(beta,P))

def maninleft(g,i,c,d):
    gg=lift(c,d)
    P=leftpoly(gg,(1,0) if i else (0,1))
    P=leftpoly(g,P)
    alpha=mob(mul(g,gg),F(0));beta=mob(mul(g,gg),None)
    return symb(P,alpha,beta)

# test every generator on identity left action
for k in range(D):
 c,d=reps[k//2]; i=k%2
 assert reduce(maninleft((1,0,0,1),i,c,d))==red[k,:],(k,c,d)
print('CF identity checks pass',flush=True)

def Hecke(p):
    mats=[((1,r,0,p),1) for r in range(p)]
    if N%p: mats.append(((p,0,0,1),chi(p)))
    T=sp.zeros(qD,qD)
    for row,k in enumerate(free):
      c,d=reps[k//2];i=k%2
      z={}
      for mat,ch in mats:
        v=maninleft(mat,i,c,d)
        z=sumv(z,{j:ch*a for j,a in v.items()})
      T[row,:]=reduce(z)
    # reconstruct full images independently to check annihilation of relations
    full=sp.zeros(D,qD)
    for k in range(D):
      c,d=reps[k//2];i=k%2;z={}
      for mat,ch in mats:
        vv=maninleft(mat,i,c,d)
        z=sumv(z,{j:ch*a for j,a in vv.items()})
      full[k,:]=reduce(z)
    defect=full-red*T
    print('hecke',p,'well-defined',defect==sp.zeros(D,qD),flush=True)
    assert defect==sp.zeros(D,qD)
    return T

if __name__=='__main__':
    t=time.time()
    Ts={}
    for p in [3,7,11]:
      T=Hecke(p);Ts[str(p)]=T
      print('T',p,'charpoly',sp.factor(T.charpoly().as_expr()),'sec',time.time()-t,flush=True)
    for p,T in Ts.items():
      for pp,TT in Ts.items(): assert T*TT==TT*T
    # standard winding symbol (X times path 0 to infinity)
    wind=reduce(symb((1,0),F(0),None))
    # quadratic additive twist, conductor5
    twist=reduce(sumv(*[symb((1,-F(a,5)),F(a,5),None) if a in (1,4) else {j:-v for j,v in symb((1,-F(a,5)),F(a,5),None).items()} for a in range(1,5)]))
    out={'reps':reps,'free':free,'reduction':[[str(x) for x in row] for row in red.tolist()], 'T':{p:[[str(x) for x in row] for row in T.tolist()] for p,T in Ts.items()},'winding':[str(x) for x in wind],'twist5':[str(x) for x in twist]}
    open(str(ROOT/'results/modsym_data.json'),'w').write(json.dumps(out))
