"""Exact finite certificates. No Bessel or modular period value is an input."""
from pathlib import Path
from fractions import Fraction as F
from math import gcd,isqrt
import json,time
import sympy as s
from kernel_coefficients import kernel_coefficient,chi4,chip,chiD,vp,tau,hfund
from modular_basis import DS,EXPS,eta_series
ROOT=Path(__file__).resolve().parents[1]
KCO=[0,F(-15,8),F(269,128),F(-5,2),F(127,32),F(-235,64),-4,F(35,8),F(29,128),F(11,16),F(5,4),F(1,6),F(37,64),F(15,64),F(33,128),F(-3,16),F(13,48),F(87,128)]
WCO=[F(x) for x in ['1/192', '-43/55296', '17/27648', '-619/27648', '-871/27648', '-11/384', '-269/3456', '1/6', '-1/6144', '-1/3072', '0', '0', '49/2880', '-13/8640', '1/72', '-1/160', '0', '0']]
PCO=[F(49,144),F(-5,48),F(5,16),F(1,16)]
PSPECS=[{2:1,6:2,10:4,30:-1},{2:4,6:-1,10:1,30:2},{2:-1,6:4,10:2,30:1},{2:2,6:1,10:-1,30:4}]
checks=[]
def check(name,cond):
 if not cond:raise AssertionError(name)
 checks.append(name)
def lin(co,B,n):return sum((F(c)*b[n] for c,b in zip(co,B)),F(0))
def hecke(a,p,n):return a[p*n]+(chiD(-15,p)*p*p*a[n//p] if n%p==0 else 0)
def oddtau_formula(n,l,p,a):
 D=l*l-4*n;v=vp(D,p)
 if a>v:return 0
 if a<v:return (p-1)*p**(a//2) if a%2==0 else 0
 if v%2==0:return -p**(v//2)
 return p**((v+1)//2)*chip(-D//p**v,p)
def roots(n,l,p,a):
 if a==0:return 1
 m=p**a;return sum((r*r+l*r+n)%m==0 for r in range(m))
def predicted_S_ratio(n,l,c):
 k=c//60;bad=1;out=1
 for p in (2,3,5):
  a=0
  while k%p==0:a+=1;k//=p;bad*=p
  out*=tau(n,l,p,a)*((-1)**a if p==5 else 1)
 out*=bad
 for p,a in s.factorint(k).items():
  p=int(p);out*=chiD(60,p)**a*p**a*(roots(n,l,p,a)-roots(n,l,p,a-1))
 return out

def main():
 st=time.time();N=1764
 cusp=[]
 for e in EXPS+[[sp.get(d,0) for d in DS] for sp in PSPECS]:
  assert sum(e)==6 and sum(d*v for d,v in zip(DS,e))%24==0 and sum(60//d*v for d,v in zip(DS,e))%24==0
  vals=[F(60,24*gcd(c,60//c)*c)*sum((F(gcd(c,d)**2*v,d) for d,v in zip(DS,e)),F(0)) for c in DS]
  assert min(vals)>0
  char=[sum(v*vp(d,p) for d,v in zip(DS,e))%2 for p in (2,3,5)]
  assert char==[0,1,1]
  cusp.append([str(v) for v in vals])
 check('22 eta quotients: weight, both congruences, all cusp orders, character',True)
 B=[eta_series(e,N) for e in EXPS];P=[eta_series([sp.get(d,0) for d in DS],N) for sp in PSPECS]
 K=[F(0)]+[kernel_coefficient(n) for n in range(1,201)]
 check('critical class-number coefficients equal explicit eta form through q^200',all(K[n]==lin(KCO,B,n) for n in range(201)))
 H=[lin(PCO,P,n) for n in range(N+1)];R=[lin(WCO,B,n) for n in range(N+1)]
 check('direct Hecke witness through Sturm bound36',all(lin(KCO,B,n)-H[n]==R[49*n]+(15-49*(n%7==0))*R[n]+(2401*R[n//49] if n%49==0 else 0) for n in range(37)))
 # Reconstruct the good Hecke matrices of the four-form space, using the Sturm injection.
 A=s.Matrix([[p[n] for p in P] for n in range(37)])
 good={}
 for p in (7,11):
  M=s.Matrix.hstack(*(A.gauss_jordan_solve(s.Matrix([hecke(col,p,n) for n in range(37)]))[0] for col in P))
  good[str(p)]=[[str(v) for v in M.row(i)] for i in range(4)]
  check(f'T_{p} full four-form modular identity',A*M==s.Matrix([[hecke(col,p,n) for col in P] for n in range(37)]))
 U=s.Matrix.hstack(*(A.gauss_jordan_solve(s.Matrix([col[3*n] for n in range(37)]))[0] for col in P))
 check('U3 polynomial on four-form space',U**4-2*U**2+81*s.eye(4)==s.zeros(4))
 r=s.sqrt(5);I=s.I;om=(1-2*I)/r;eps=(5+4*r+I*(2*r-10))/15
 fv=s.Matrix([1,om,3*eps,3*eps*s.conjugate(om)]);lam=r-2*I
 check('fixed f eigenvalues U3,T7,T11',all(s.simplify(v)==0 for v in (U-lam*s.eye(4))*fv) and all(s.simplify(v)==0 for v in (s.Matrix(good['7'])+8*I*s.eye(4))*fv) and all(s.simplify(v)==0 for v in (s.Matrix(good['11'])-4*I*r*s.eye(4))*fv))
 h=[H[3**j] for j in range(4)]
 delta=s.simplify(s.expand_complex((h[3]+lam*h[2]+(lam**2-2)*h[1]+(lam**3-2*lam)*h[0])/(4*lam**3-4*lam)))
 check('exact f spectral coefficient',s.simplify(delta-(7-I*r)*(7-I)/576)==0)
 badR=(1+(1+4*I*r)/27)*(1+(15+20*I)/125)
 check('pruned zeta factor',F(1,6)*F(3,4)*F(8,9)*F(24,25)==F(8,75))
 check('L3 to Petersson norm constant',s.simplify(128*s.conjugate(delta)/(75*badR))==s.Rational(1,8))
 # Formal local identities, including unequal roots.
 aa,bb,z=s.symbols('aa bb z');even=(1+aa*bb*z)/((1-aa**2*z)*(1-bb**2*z))
 sym=1/((1-aa**2*z)*(1-aa*bb*z)*(1-bb**2*z))
 check('square-index Euler identity',s.cancel(even/sym-(1-(aa*bb*z)**2))==0)
 check('local inverse factors at3 and5',s.expand((1+(1-4*I*r)*z)*(1+(1+4*I*r)*z))==1+2*z+81*z*z and s.expand((1-(-15+20*I)*z)*(1-(-15-20*I)*z))==1+30*z+625*z*z)
 check('functional-equation ratio at2 and3',s.simplify(2*s.sqrt(960)/(2*s.pi**s.Rational(3,2))/s.sqrt(s.pi)-8*s.sqrt(15)/s.pi**2)==0)
 # Tests of all-index local lemmas, including exceptional square discriminants.
 oddtests=0;twotests=0
 for n in range(1,21):
  for l in range(0,13):
   D=l*l-4*n
   if D==0:continue
   for p in (3,5):
    for a in range(vp(D,p)+3):
     assert tau(n,l,p,a)==oddtau_formula(n,l,p,a);oddtests+=1
   for a in range(vp(D,2)+2,vp(D,2)+5):
    assert tau(n,l,2,a)==0;twotests+=1
 check('odd local valuation formula grid',True);check('2-adic support cutoff grid',True)
 primitivepairs=0
 for u in range(8):
  for v in range(8):
   if u%2==0 and v%2==0:continue
   Q=v*v-15*u*u
   assert (Q//2**vp(Q,2))%4==1;primitivepairs+=1
 check('primitive norm15 odd part is1mod4; parity classes',True)
 # Exact cyclotomic Gauss/CRT controls, before using any L or class-number formula.
 X=s.Symbol('X');gauss=[]
 for n,l in [(1,1),(7,5),(3,2),(1,62)]:
  for c in (60,120,180,300,780):
   v=[0]*c
   us=[(d,chiD(60,d)) for d in range(c) if gcd(d,c)==1]
   for rr in range(c):
    cr=chi4(rr)
    if not cr:continue
    q=(rr*rr+l*rr+n)%c
    for dd,cd in us:v[dd*q%c]+=cr*cd
   rat=predicted_S_ratio(n,l,c)
   for dd in range(60):v[(c//60)*dd]-=rat*chiD(60,dd)
   poly=s.Poly.from_dict({(i,):x for i,x in enumerate(v) if x},X,domain=s.ZZ)
   assert poly.rem(s.Poly(s.cyclotomic_poly(c,X),X)).is_zero
   gauss.append({'n':n,'l':l,'c':c,'ratio_to_sqrt60':rat})
 check('20 cyclotomic direct sums equal local CRT formula',True)
 out={'checks':checks,'N_kernel':200,'N_eta':N,'eta_divisors':DS,'eta_exponents':EXPS,'cusp_orders':cusp,'kernel_coefficients':list(map(str,K)),'kernel_eta_coordinates':list(map(str,KCO)),'f_packet_P_coordinates':list(map(str,PCO)),'hecke_witness_coordinates':list(map(str,WCO)),'good_hecke_matrices':good,'U3':[[str(v) for v in U.row(i)] for i in range(4)],'h':[str(x) for x in h],'delta_f':str(delta),'L3_norm_factor':'1/8','L2_norm_factor':'sqrt(15)*pi^2','local_test_counts':{'odd':oddtests,'2_support':twotests,'primitive_norm_pairs':primitivepairs},'cyclotomic_tests':gauss,'claim':'modular_norm_identity_proved_in_written_argument; original_Bessel_eta_identity_unproved'}
 (ROOT/'outputs'/'exact.json').write_text(json.dumps(out,indent=2)+'\n')
 print('PASS',len(checks),'grouped checks; elapsed',round(time.time()-st,2),flush=True)
if __name__=='__main__':main()
