"""Independent q-expansion reconstruction, not using archived L values.
Candidate eta carrier and bad coefficients are audited against each other below.
Numerical integrations are exponentially convergent diagnostics, not interval proofs.
"""
from pathlib import Path
import json,math,time
import mpmath as mp
import sympy as sp

TERMS=[(-8,{12:3,20:3}),(8,{4:1,6:2,10:2,60:1}),(-4,{4:2,6:1,10:1,60:2}),(-6,{3:1,5:1,12:2,20:2}),(3,{3:2,5:2,12:1,20:1}),(1,{3:3,5:3}),(-3,{2:1,3:2,5:2,30:1}),(2,{2:2,3:1,5:1,30:2})]
def eta_coeff(N):
 out=[0]*(N+1)
 for c,exps in TERMS:
  shift=sum(d*r for d,r in exps.items())//24
  a=[0]*(N+1-shift);a[0]=1
  # Repeated in-place multiplication by (1-q^(dk)); exact integers.
  for d,r in exps.items():
   for k in range(1,(N-shift)//d+1):
    m=d*k
    for _ in range(r):
     for j in range(len(a)-1,m-1,-1):a[j]-=a[j-m]
  for n,an in enumerate(a):out[n+shift]+=c*an
 return out
# Q(sqrt(5), i): basis (1,sqrt(5),i,i sqrt(5))
z=(0,0,0,0);one=(1,0,0,0)
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def sc(a,c):return tuple(c*x for x in a)
def mul(a,b):
 c=[0]*4
 for j in range(4):
  for k in range(4):
   c[j^k]+=a[j]*b[k]*(5 if j&k&1 else 1)*(-1 if j&k&2 else 1)
 return tuple(c)
def ell(a):return a[0]-2*a[1]+a[2]-a[3]
def numeric(a):return mp.mpc(a[0]+mp.sqrt(5)*a[1],a[2]+mp.sqrt(5)*a[3])
def coefficients(F):
 N=len(F)-1;coef=[z]*(N+1);coef[1]=one
 spf=list(range(N+1))
 for p in range(2,N+1):
  if spf[p]!=p:continue
  for k in range(p*p,N+1,p):
   if spf[k]==k:spf[k]=p
 apr={2:z,3:(0,1,-2,0),5:(0,1,0,2)}
 for p in range(7,N+1):
  if spf[p]!=p:continue
  c3=int(sp.kronecker_symbol(-3,p));c5=int(sp.kronecker_symbol(5,p));f=F[p]
  if c3==1:ap=(f,0,0,0) if c5==1 else (0,0,f,0)
  elif c5==1:ap=(0,0,0,-f)
  else:
   assert f%2==0
   ap=(0,-f//2,0,0)
  apr[p]=ap
 for n in range(2,N+1):
  p=spf[n];m=n//p
  if m%p:coef[n]=mul(apr[p],coef[m])
  else:coef[n]=add(mul(apr[p],coef[m]),sc(coef[m//p],-(int(sp.kronecker_symbol(-15,p)) if 60%p else 0)*p*p))
 return coef

def main():
 mp.mp.dps=160;N=1100
 t=time.time();F=eta_coeff(N);aa=coefficients(F)
 bad=[(n,F[n],ell(aa[n])) for n in range(1,N+1) if F[n]!=ell(aa[n])]
 print('carrier checks:',N,'mismatches',len(bad),bad[:12],flush=True)
 a=[numeric(x) for x in aa]
 y0=1/mp.sqrt(60)
 def f(y):return mp.fsum(a[n]*mp.exp(-2*mp.pi*n*y) for n in range(1,N+1))
 eps=f(y0)/mp.conj(f(y0))
 print('eps',mp.nstr(eps,100),flush=True)
 errs=[]
 for fac in [mp.mpf('0.7'),mp.mpf('1.25'),mp.mpf('1.7')]:
  y=fac*y0
  err=abs(f(1/(60*y))-eps*60**mp.mpf('1.5')*y**3*mp.conj(f(y)))
  errs.append(mp.nstr(err,8));print('Fricke',fac,mp.nstr(err,8),flush=True)
 h=2*mp.pi/mp.sqrt(60)
 L2=mp.fsum(a[n]*mp.exp(-h*n)*(1+h*n)/n**2+eps*h*mp.conj(a[n])*mp.exp(-h*n)/n for n in range(1,N+1))
 print('L2',mp.nstr(L2,155),flush=True)
 alpha=(125-13*mp.sqrt(5)+1j*(125-41*mp.sqrt(5)))/4
 candidate=mp.pi*mp.re(alpha*L2**2)
 print('Ccandidate',mp.nstr(candidate,155),flush=True)
 data={'scope':'Independent exact q-series reconstruction; numerical L-value not interval-certified; eta carrier identity and eigenform identification need source audit.','N':N,'coefficient_mismatches':bad[:20],'first_coefficients':[list(x) for x in aa[:41]],'Fricke_errors':errs,'epsilon':{'re':mp.nstr(eps.real,160),'im':mp.nstr(eps.imag,160)},'L2':{'re':mp.nstr(L2.real,160),'im':mp.nstr(L2.imag,160)},'candidate':mp.nstr(candidate,160),'elapsed':time.time()-t}
 Path(__file__).with_suffix('.json').write_text(json.dumps(data,indent=2))
 print('elapsed',time.time()-t,flush=True)
if __name__=='__main__':main()
