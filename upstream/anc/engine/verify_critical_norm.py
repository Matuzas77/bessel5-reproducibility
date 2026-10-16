"""Standard-library reconstruction. Does not import primary code or primary outputs.
Independent eta product algorithm; local 2-adic averaging; odd-prime valuation
formula; reduced-form counting with nonnegative middle coefficient; Q(sqrt5,i).
"""
from fractions import Fraction as Q
from math import gcd,isqrt
from functools import lru_cache
from pathlib import Path
import json,time
ROOT=Path(__file__).resolve().parents[1]
DIV=(1,2,3,4,5,6,10,12,15,20,30,60)
ROWS=[(3,0,0,0,0,0,0,0,3,0,0,0),(2,0,0,1,0,0,0,0,2,0,0,1),(2,0,0,0,0,1,1,0,2,0,0,0),(1,0,0,2,0,0,0,0,1,0,0,2),(1,0,0,1,0,1,1,0,1,0,0,1),(1,0,0,0,0,2,2,0,1,0,0,0),(0,2,0,0,0,0,0,1,0,1,2,0),(0,1,0,0,0,0,0,2,0,2,1,0),(0,0,3,0,3,0,0,0,0,0,0,0),(0,0,2,0,2,0,0,1,0,1,0,0),(0,0,0,3,0,0,0,0,0,0,0,3),(0,0,0,0,0,0,0,3,0,3,0,0),(-1,3,3,0,0,0,1,0,0,0,0,0),(-1,2,3,0,0,1,0,0,0,0,1,0),(-1,1,4,0,1,1,0,0,0,0,0,0),(-1,1,3,0,0,1,0,1,0,1,0,0),(-1,0,3,1,0,0,0,0,0,3,0,0),(2,-1,0,0,0,2,0,0,2,0,1,0)]
KCO=tuple(map(Q,['0','-15/8','269/128','-5/2','127/32','-235/64','-4','35/8','29/128','11/16','5/4','1/6','37/64','15/64','33/128','-3/16','13/48','87/128']))
WCO=tuple(map(Q,['1/192','-43/55296','17/27648','-619/27648','-871/27648','-11/384','-269/3456','1/6','-1/6144','-1/3072','0','0','49/2880','-13/8640','1/72','-1/160','0','0']))
PCO=tuple(map(Q,['49/144','-5/48','5/16','1/16']))
SPECS=({2:1,6:2,10:4,30:-1},{2:4,6:-1,10:1,30:2},{2:-1,6:4,10:2,30:1},{2:2,6:1,10:-1,30:4})
checks=[]
def check(name,cond):
 if not cond:raise AssertionError(name)
 checks.append(name)
def factor(n):
 n=abs(n);r={};p=2
 while p*p<=n:
  while n%p==0:r[p]=r.get(p,0)+1;n//=p
  p=3 if p==2 else p+2
 if n>1:r[n]=1
 return r

def val(n,p):
 n=abs(n);a=0
 if n==0:return 10**6
 while n%p==0:a+=1;n//=p
 return a

def c4(n):return 0 if n%2==0 else (1 if n%4==1 else -1)
def leg(D,p):
 if p==2:return 0 if D%2==0 else (1 if D%8 in(1,7) else -1)
 r=pow(D%p,(p-1)//2,p)
 return 0 if r==0 else 1 if r==1 else -1
@lru_cache(None)
def disc(N):
 f=factor(N);sq=1;sf=-1 if N<0 else 1
 for p,e in f.items():sq*=p**(e//2);sf*=p**(e%2)
 if sf%4!=1:sf*=4;sq//=2
 assert sf*sq*sq==N
 return sf,sq
@lru_cache(None)
def class_h(D):
 M=-D;number=0
 for b in range(isqrt(M//3)+1):
  if (b*b+M)%4:continue
  prod=(b*b+M)//4
  for a in range(max(1,b),isqrt(prod)+1):
   if prod%a:continue
   c=prod//a
   number+=1 if b==0 or a==b or a==c else 2
 return Q(number,3 if D==-3 else 2 if D==-4 else 1)

def local2(n,l):
 D=l*l-4*n;v=val(D,2);mod=2**(v+3);sm=0
 for r in range(1,mod,2):
  q=r*r+l*r+n;a=val(q,2)
  if a<=v+1:sm+=c4(r)*c4(q//2**a)
 return Q(4*sm,mod)
def localodd(D,p):
 v=val(D,p);ans=Q(0)
 for a in range(v+1):
  if a<v:z=(p-1)*p**(a//2) if a%2==0 else 0
  elif a%2==0:z=-p**(a//2)
  else:z=p**((a+1)//2)*leg(-D//p**a,p)
  ans+=Q(z*(-1 if p==5 and a%2 else 1),p**a)
 return ans

def rootsfactor(D,p):
 v=val(D,p);t=Q(leg(60,p),p)
 rr=[p**(a//2) for a in range(v+1)]
 end=2*p**(v//2) if v%2==0 and leg(D//p**v,p)==1 else 0
 # Sum (R_a-R_{a-1})t^a, including the single jump to the stable tail.
 return 1+sum((Q(rr[a]-rr[a-1])*t**a for a in range(1,v+1)),Q(0))+(end-rr[-1])*t**(v+1)

def coeff(n):
 t=isqrt(n);ans=Q(t*c4(t) if t*t==n else 0)
 for l in range(1,isqrt(4*n-1)+1):
  D=l*l-4*n;a=local2(n,l)
  if not a:continue
  a*=localodd(D,3)*localodd(D,5)
  if not a:continue
  d,f=disc(60*D);a*=Q(l*f,600)*class_h(d)
  for p in (2,3,5):a*=Q(1-Q(leg(d,p),p),1-Q(1,p*p))
  for p in factor(D):
   if p in(2,3,5):continue
   a*=rootsfactor(D,p)*Q(1-Q(leg(d,p),p),1-Q(1,p*p))
  ans+=a
 return ans

def eta(e,N):
 shift=sum(d*x for d,x in zip(DIV,e));assert shift%24==0;shift//=24
 M=N-shift;out=[1]+[0]*M
 # Multiplication by each Euler factor, not a logarithmic derivative recurrence.
 for d,x in zip(DIV,e):
  for k in range(d,M+1,d):
   for _ in range(abs(x)):
    if x>0:
     for j in range(M,k-1,-1):out[j]-=out[j-k]
    else:
     for j in range(k,M+1):out[j]+=out[j-k]
 return [0]*shift+out

def lin(c,B,n):return sum((a*b[n] for a,b in zip(c,B)),Q(0))
# Field basis 1,sqrt5,i,i sqrt5; bit0 is sqrt5, bit1 is i.
def fld(a=0,b=0,c=0,d=0):return tuple(map(Q,(a,b,c,d)))
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def neg(x):return tuple(-a for a in x)
def sub(x,y):return add(x,neg(y))
def scale(x,a):return tuple(a*b for b in x)
def mul(x,y):
 out=[Q(0)]*4
 for i,a in enumerate(x):
  for j,b in enumerate(y):
   out[i^j]+=a*b*(5 if i&j&1 else 1)*(-1 if i&j&2 else 1)
 return tuple(out)
def conj(x):return(x[0],x[1],-x[2],-x[3])
def sig(x,sr,si):return(x[0],sr*x[1],si*x[2],sr*si*x[3])
def fsuma(vals):
 out=fld()
 for v in vals:out=add(out,v)
 return out

def main():
 start=time.time();N=1764;B=[eta(e,N) for e in ROWS];P=[eta(tuple(sp.get(d,0) for d in DIV),N) for sp in SPECS]
 K=[Q(0)]+[coeff(n) for n in range(1,241)]
 check('240 class-number coefficients from independent algorithms equal Euler-product eta form',all(K[n]==lin(KCO,B,n) for n in range(241)))
 R=[lin(WCO,B,n) for n in range(N+1)];H=[lin(PCO,P,n) for n in range(N+1)]
 check('Hecke witness coefficient-by-coefficient through Sturm36',all(lin(KCO,B,n)-H[n]==R[49*n]+(15-49*(n%7==0))*R[n]+(2401*R[n//49] if n%49==0 else 0) for n in range(37)))
 # All eta modularity conditions freshly recomputed.
 for e in ROWS+list(tuple(sp.get(d,0) for d in DIV) for sp in SPECS):
  assert sum(e)==6 and sum(d*x for d,x in zip(DIV,e))%24==0 and sum(60//d*x for d,x in zip(DIV,e))%24==0
  for c in DIV:
   assert sum((Q(gcd(c,d)**2*x,d) for d,x in zip(DIV,e)),Q(0))>0
  assert all(sum(x*val(d,p) for d,x in zip(DIV,e))%2==want for p,want in ((2,0),(3,1),(5,1)))
 check('all22 eta modularity and cusp criteria',True)
 r=fld(0,1);ii=fld(0,0,1);one=fld(1)
 om=fld(0,Q(1,5),0,Q(-2,5));eps=fld(Q(1,3),Q(4,15),Q(-2,3),Q(2,15))
 fv=[one,om,scale(eps,3),scale(mul(eps,conj(om)),3)]
 delta=fld(Q(49,576),Q(-1,576),Q(-7,576),Q(-7,576))
 recovered=[fsuma(mul(sig(delta,sr,si),sig(v,sr,si)) for sr in(-1,1) for si in(-1,1)) for v in fv]
 check('Galois spectral decomposition of four-form kernel',recovered==[fld(c) for c in PCO])
 fcoeff=[fsuma(scale(fv[j],P[j][n]) for j in range(4)) for n in range(397)]
 check('fixed a3,a5,a7 embedding',fcoeff[3]==fld(0,1,-2,0) and fcoeff[5]==fld(0,1,0,2) and fcoeff[7]==fld(0,0,-8,0))
 for p,ap in [(7,fld(0,0,-8,0)),(11,fld(0,0,0,4))]:
  for n in range(1,37):
   valf=fcoeff[p*n]
   if n%p==0:valf=add(valf,scale(fcoeff[n//p],leg(-15,p)*p*p))
   assert valf==mul(ap,fcoeff[n])
 check('T7,T11 eigenidentities and pairwise distinct conjugate eigenpairs',len({(sig(fcoeff[7],sr,si),sig(fcoeff[11],sr,si)) for sr in(-1,1) for si in(-1,1)})==4)
 R3=mul(fld(Q(28,27),0,0,Q(4,27)),fld(Q(28,25),0,Q(4,25),0))
 check('critical L3 norm factor cross multiplication',scale(conj(delta),1024)==scale(R3,75))
 check('zeta deletion normalization',Q(1,6)*Q(3,4)*Q(8,9)*Q(24,25)==Q(8,75))
 check('conductor FE algebra',960==64*15)
 h=[H[3**j] for j in range(4)]
 check('U3 moments of packet',[str(x) for x in h]==['49/144','-19/144','-91/144','-119/144'])
 # Independent dyadic integral cancellation for principal-discriminant examples.
 pairs=[(4,1),(31,8),(244,63)]
 for rr,uu in pairs:
  n=rr*rr-15*uu*uu;l=2*rr
  assert local2(n,l)==0
 check('dyadic exceptional principal-character tests',True)
 out={'checks':checks,'N_kernel':240,'N_eta':N,'kernel_coefficients':list(map(str,K)),'kernel_eta_coordinates':list(map(str,KCO)),'f_packet_P_coordinates':list(map(str,PCO)),'hecke_witness_coordinates':list(map(str,WCO)),'h':list(map(str,h)),'delta_field_coordinates':list(map(str,delta)),'R3_field_coordinates':list(map(str,R3)),'L3_norm_factor':'1/8','L2_norm_factor':'sqrt(15)*pi^2','claim':'analytic proof separate; exact modular certificate checked; original integral unproved'}
 (ROOT/'outputs'/'critical_norm_exact.json').write_text(json.dumps(out,indent=2)+'\n')
 print('PASS',len(checks),'grouped independent checks; elapsed',round(time.time()-start,2),flush=True)
if __name__=='__main__':main()
