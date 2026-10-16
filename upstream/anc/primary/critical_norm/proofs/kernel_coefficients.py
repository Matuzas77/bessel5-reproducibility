"""Exact critical Poincare-kernel coefficients; analytic proof is in REPORT.md."""
from fractions import Fraction as F
from math import isqrt
from functools import lru_cache
import sympy as sp

def chi4(n): return 0 if n%2==0 else (1 if n%4==1 else -1)
def chip(n,p):
 n%=p
 return 0 if not n else (1 if pow(n,(p-1)//2,p)==1 else -1)
def vp(n,p):
 n=abs(n);v=0
 if not n:return 0
 while n%p==0:n//=p;v+=1
 return v
@lru_cache(None)
def tau_cached(n,l,p,a):
 pa=p**a;mod=pa*(4 if p==2 else p);ans=0
 for r in range(mod):
  Q=r*r+l*r+n
  if Q%pa:continue
  ans+=chi4(r)*chi4(Q//pa) if p==2 else chip(Q//pa,p)
 return ans
def tau(n,l,p,a):
 mod=p**a*(4 if p==2 else p)
 return tau_cached(n%mod,l%mod,p,a)
def bad(n,l):
 D=l*l-4*n;parts=[]
 for p in (2,3,5):
  vals=[tau(n,l,p,a) for a in range(vp(D,p)+(2 if p==2 else 1))]
  parts.append(sum((F(v*((-1)**a if p==5 else 1),p**a) for a,v in enumerate(vals)),F(0)))
 return parts
@lru_cache(None)
def fundamental(Delta):
 sf=1;f=1
 for p,e in sp.factorint(abs(Delta)).items():sf*=int(p)**(e%2);f*=int(p)**(e//2)
 if Delta<0:sf=-sf
 if sf%4==1:D=sf;F0=f
 else:D=4*sf;F0=f//2;assert f%2==0
 assert D*F0**2==Delta
 return D,F0
@lru_cache(None)
def hfund(D):
 assert D<0 and D%4 in(0,1)
 h=0
 for a in range(1,isqrt((-D)//3)+1):
  for b in range(-a,a+1):
   if (b*b-D)%(4*a):continue
   c=(b*b-D)//(4*a)
   if c<a or ((abs(b)==a or a==c) and b<0):continue
   h+=1
 return F(h,3 if D==-3 else 2 if D==-4 else 1)
def chiD(D,p):return int(sp.kronecker_symbol(D,p))
def ep(D,p):
 v=vp(D,p);t=F(chiD(60,p),p)
 E=(1-t)*sum((p**(a//2)*t**a for a in range(v+1)),F(0))
 if v%2==0 and chip(D//p**v,p)==1:E+=2*p**(v//2)*t**(v+1)
 return E

def kernel_coefficient(n,sign=1):
 main=F(isqrt(n)*chi4(isqrt(n)) if isqrt(n)**2==n else 0);corr=F(0)
 for l in range(1,isqrt(4*n-1)+1):
  delta=l*l-4*n;lb=bad(n,l)
  if 0 in lb:continue
  D,F0=fundamental(60*delta)
  z=hfund(D)*F0*l*F(1,600)
  for b in lb:z*=b
  for p in(2,3,5):z*=F(1-F(chiD(D,p),p),1-F(1,p*p))
  for p in sp.factorint(abs(delta)):
   p=int(p)
   if p in(2,3,5):continue
   z*=ep(delta,p)/F(1-F(1,p*p),1-F(chiD(D,p),p))
  corr+=z
 return main+sign*corr
if __name__=='__main__':
 import json,time,argparse
 from pathlib import Path
 ap=argparse.ArgumentParser();ap.add_argument('--N',type=int,default=60);args=ap.parse_args()
 st=time.time();v=[F(0)]+[kernel_coefficient(n) for n in range(1,args.N+1)]
 out={'N':args.N,'coefficients':[str(x) for x in v],'status':'EXACT_COEFFICIENTS_ANALYTIC_PROOF_IN_REPORT'}
 p=Path(__file__).resolve().parents[1]/'outputs'/f'kernel_{args.N}.json';p.write_text(json.dumps(out,indent=2)+'\n')
 print('Elapsed',time.time()-st,'first40',v[:40])
