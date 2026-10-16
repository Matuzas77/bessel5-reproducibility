#!/usr/bin/env python3
"""Independent verification: character Manin relations, ceiling CF, exact fractions.
Does not import SymPy or the primary modular-symbol matrices. Reads only the
explicit rational witness as data, and reconstructs its identity from paths.
"""
from fractions import Fraction as Q
from math import gcd, ceil
from functools import lru_cache
from pathlib import Path
import json,time
ROOT=Path(__file__).resolve().parent
LEVEL=60
units=[a for a in range(LEVEL)if gcd(a,LEVEL)==1]
def chi(a):
    if gcd(a,15)!=1:return 0
    return (1 if a%3==1 else -1)*(1 if a%5 in (1,4)else -1)
lookup={}
for c in range(LEVEL):
 for d in range(LEVEL):
  if gcd(gcd(c,d),LEVEL)!=1 or (c,d)in lookup:continue
  rep=max((c*u%LEVEL,d*u%LEVEL)for u in units)
  for u in units:lookup[(rep[0]*u%LEVEL,rep[1]*u%LEVEL)]=(rep,chi(u))
reps=sorted(set(v[0]for v in lookup.values()),reverse=True)
idx={v:k for k,v in enumerate(reps)}

def add(v,k,a):
 if not a:return v
 v[k]=v.get(k,Q(0))+a
 if not v[k]:del v[k]
 return v

def merge(v,w,scale=1):
 v=v.copy()
 for k,a in w.items():add(v,k,scale*a)
 return v

def gen(i,c,d,a=1):
 rep,e=lookup[(c%60,d%60)]
 return {2*idx[rep]+i:Q(a)*e} if a else {}

relations=[]
# i=0 represents Y, i=1 represents X. Relations are explicitly expanded.
for c,d in reps:
 relations.append(merge(gen(0,c,d),gen(1,d,-c)))
 relations.append(merge(gen(1,c,d),gen(0,d,-c),-1))
 # U=(0,-1;1,-1); U^2=(-1,1;-1,0)
 v=gen(0,c,d);v=merge(v,gen(1,d,-c-d));v=merge(v,gen(0,d,-c-d),-1);v=merge(v,gen(1,-c-d,c),-1);relations.append(v)
 v=gen(1,c,d);v=merge(v,gen(0,d,-c-d),-1);v=merge(v,gen(1,-c-d,c),-1);v=merge(v,gen(0,-c-d,c));relations.append(v)

pivots={}
for rr in reversed(relations):
 row=rr.copy()
 while row:
  k=max(row)
  if k in pivots:row=merge(row,pivots[k],-row[k])
  else:
   val=row[k];pivots[k]={j:a/val for j,a in row.items()};break
free=[j for j in range(288)if j not in pivots]
assert len(free)==48

def reduce(v):
 out=v.copy()
 for k in sorted(pivots,reverse=True):
  if k in out:out=merge(out,pivots[k],-out[k])
 return out
assert all(not reduce(v)for v in relations)

@lru_cache(None)
def inverse_lift(c,d):
 if c==0:c=60
 while gcd(c,d)!=1:d+=60
 # Find a,b with ad-bc=1 by iterative Euclid.
 rr,ss=d,c;u,v=1,0;w,z=0,1
 while ss:
  q=rr//ss;rr,ss=ss,rr-q*ss;u,v=v,u-q*v;w,z=z,w-q*z
 assert rr==1
 return (u,-w,c,d)

def compose(A,B):
 a,b,c,d=A;e,f,g,h=B
 return (a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h)
def image(A,x):
 a,b,c,d=A
 if x is None:return None if c==0 else Q(a,c)
 den=c*x+d
 return (a*x+b)/den if den else None

def act_poly(A,P):
 a,b,c,d=A;u,v=P
 return (d*u-c*v,-b*u+a*v)

@lru_cache(None)
def negative_cf_edges(x):
 if x is None:return ()
 n,d=x.numerator,x.denominator
 aa=[]
 while d:
  a=-((-n)//d);aa.append(a);n,d=d,a*d-n
 p0,p1=0,1;q0,q1=-1,0;out=[]
 for a in aa:
  p=a*p1-p0;q=a*q1-q0
  det=p*q1-p1*q
  assert det in(-1,1)
  out.append((det*p,p1,det*q,q1))
  p0,p1=p1,p;q0,q1=q1,q
 assert Q(p1,q1)==x
 return tuple(out)

def inf_path(x,P):
 out={};u,v=P
 for a,b,c,d in negative_cf_edges(x):
  # inverse action P(aX+bY,cX+dY)
  out=merge(out,gen(1,c,d,a*u+c*v))
  out=merge(out,gen(0,c,d,b*u+d*v))
 return out

def path(P,a,b):return merge(inf_path(b,P),inf_path(a,P),-1)

@lru_cache(None)
def single_hecke(k,p):
 c,d=reps[k//2];i=k%2
 base=inverse_lift(c,d)
 P=act_poly(base,(i,1-i))
 out={}
 mat=[((1,j,0,p),1)for j in range(p)]
 if 60%p:mat.append(((p,0,0,1),chi(p)))
 for B,ch in mat:
  C=compose(B,base)
  raw=path(act_poly(B,P),image(C,Q(0)),image(C,None))
  out=merge(out,raw,ch)
 return reduce(out)

def hecke(v,p):
 out={}
 for k,a in reduce(v).items():out=merge(out,single_hecke(k,p),a)
 return out

# Alternative CF reconstructs the identity on all 288 original generators.
for k in range(288):
 c,d=reps[k//2];i=k%2;g=inverse_lift(c,d)
 assert reduce(path(act_poly(g,(i,1-i)),image(g,Q(0)),image(g,None)))==reduce({k:Q(1)})

# Explicit rational witness supplied in the primary artifact, interpreted as raw
# Manin symbols, not as rows of its matrix representation.
dat=json.loads((ROOT/'twist_witness.json').read_text())
R={}
for row in dat['witness']:
 R=merge(R,gen(row['monomial_index'],row['c'],row['d'],Q(row['coefficient'])))
w=reduce(path((Q(1),Q(0)),Q(0),None))
t={}
for a in range(1,5):t=merge(t,path((Q(1),-Q(a,5)),Q(a,5),None),1 if a in(1,4)else -1)
t=reduce(t)
poly={};cur=w
for a in map(Q,dat['twist_polynomial']):
 poly=merge(poly,cur,a);cur=hecke(cur,3)
rhs=merge(hecke(hecke(R,7),7),reduce(R),64)
residual=merge(merge(t,poly,-1),rhs,-1)
assert not residual,residual
# Verify U3,T7 actually descend: reconstruct every raw-generator image and test
# the complete relation set; no primary reduction matrix is imported.
for p in (3,7):
 for rel in relations:
  v={}
  for k,a in rel.items():v=merge(v,single_hecke(k,p),a)
  assert not v
print('PASS independent ceiling-CF witness, 288 generator reconstructions, all 576 relations for U3 and T7; quotient dimension48.')
json.dump({'status':'pass','quotient_dimension':48,'manin_generators':288,'relations':576,'continued_fraction':'ceiling negative CF','witness_residual':{},'uses_primary_matrices':False,'Bessel_identity_proved':False},open(ROOT.parent/'outputs'/'twist_symbols_exact.json','w'),indent=2)
