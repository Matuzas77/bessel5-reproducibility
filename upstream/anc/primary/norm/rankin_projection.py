#!/usr/bin/env python3
"""Primary constructive Rankin projection, exact SymPy arithmetic.
Independent verification is in verify_norm_independent.py.
No Bessel data or numerical L-values are used.
"""
from pathlib import Path
import sys,json,time
from fractions import Fraction as Q
from math import gcd
import sympy as s
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'input'))
from recompute_modular_value import eta_coeff
D=json.loads((HERE/'eta_basis.json').read_text())['D']
exps=json.loads((HERE/'eta_basis.json').read_text())['exponents']

def chi(d,n):return int(s.kronecker_symbol(d,n)) if d!=1 else 1

def etaq(rs,N):
 shift=sum(d*r for d,r in zip(D,rs))//24;a=[0]*(N+1);a[shift]=1
 for d,r in zip(D,rs):
  for k in range(d,N+1,d):
   for _ in range(abs(r)):
    if r>0:
     for n in range(N,k-1,-1):a[n]-=a[n-k]
    else:
     for n in range(k,N+1):a[n]+=a[n-k]
 return a

def generalized_Lminus2(dis):
 N=abs(dis);B=sum(s.Rational(chi(dis,a))*s.bernoulli(3,s.Rational(a,N)) for a in range(1,N+1))*N**2
 return -B/3

def e3(phi,psi,dil,N):
 a=[s.Rational(0)]*(N+1)
 if phi==1:a[0]=generalized_Lminus2(psi)/2
 for d in range(1,N//dil+1):
  z=chi(psi,d)*d*d
  if z:
   for k in range(1,N//(dil*d)+1):a[dil*d*k]+=z*chi(phi,k)
 return a

def g1(N):
 g=[s.Rational(0)]*(N+1)
 g[0]=-sum(chi(-15,a)*s.Rational(a,15) for a in range(1,16))/2
 for d in range(1,N+1):
  x=chi(-15,d)
  if x:
   for n in range(d,N+1,d):g[n]+=x
 return g

def holproj(N):
 E=[s.Rational(0)]*(N+1)
 for d in s.divisors(30):
  fac=s.Rational(25,16)*s.mobius(d)/d**2;E[0]+=fac;scale=60//d
  for k in range(1,N//scale+1):
   for n in range(k,N//scale+1,k): E[n*scale]-=24*k*fac
 g=g1(N);H=[sum(g[j]*E[n-j] for j in range(n+1))-s.Rational(n,12)*g[n] for n in range(N+1)]
 return g,E,H

def run():
 t=time.time(); N=252
 bases=[etaq(e,N) for e in exps]
 pairs=[(1,-15),(-15,1),(-3,5),(5,-3)]
 for ph,ps in pairs:
  for d in [1,2,4]:bases.append(e3(ph,ps,d,N))
 M=s.Matrix([[b[n] for b in bases] for n in range(37)])
 piv=M.T.rref()[1];inds=list(piv)
 print('M rank',len(inds),'secs',time.time()-t,flush=True)
 Inv=M[inds,:].inv()
 def solve(b):
  c=Inv*s.Matrix([b[i] for i in inds]); assert M*c==s.Matrix(b[:37]); return c
 cols=[]
 for b in bases:
  vec=[b[n*7]+(-49*b[n//7] if n%7==0 else 0) for n in range(37)]
  cols.append(solve(vec))
 T=s.Matrix.hstack(*cols);x=s.symbols('x');cp=s.factor(T.charpoly(x).as_expr())
 print('T7 charpoly',cp,flush=True)
 P=T*T+64*s.eye(30); ker=P.nullspace(); im=P.columnspace()
 print('ker',len(ker),'im',len(im),flush=True)
 ff=eta_coeff(972)
 W=s.Matrix.hstack(*[solve([ff[3**j*n] for n in range(37)]) for j in range(4)])
 assert P*W==s.zeros(30,4)
 image_pivots=list(P.rref()[1]); V=W.row_join(P[:,image_pivots]);assert V.det()!=0
 g,E,H=holproj(36);hv=solve(H);sol=V.inv()*hv;hcoords=sol[:4,0]
 print('H coords in F,UF..:',list(hcoords),flush=True)
 hh=[sum(hcoords[j]*ff[3**(j+k)] for j in range(4)) for k in range(4)]
 r=s.sqrt(5);lam=r-2*s.I
 delta=s.simplify(s.expand_complex((hh[3]+lam*hh[2]+(lam**2-2)*hh[1]+(lam**3-2*lam)*hh[0])/(4*lam**3-4*lam)))
 eta=(1-(r+2*s.I)/9)*(1-r*(1-2*s.I)/25)
 const=s.simplify(s.expand_complex(s.Rational(75,128)*eta/s.conjugate(delta)))
 assert s.simplify(const-3*(5-r)/2)==0
 print('g0',g[0]);print('delta=',delta);print('norm/|L|2 =',const,'/pi4')
 Bconst=s.simplify((15*s.sqrt(3)*(r-1)/2)/const)
 print('B candidate / (pi4 norm)=',Bconst)
 def sm(M):return [[str(z) for z in row] for row in M.tolist()]
 out={'status':'PASS','T7_characteristic_polynomial':str(cp),'M_rank':len(inds),'kernel_dimension':len(ker),'complement_dimension':len(im),'H_cyclic_coordinates':[str(z) for z in hcoords], 'image_pivots':image_pivots, 'witness_coefficients':[str(sol[4+image_pivots.index(i),0]) if i in image_pivots else '0' for i in range(30)],'hk':[str(z) for z in hh],'delta':str(delta),'bad_factor':str(s.expand(eta)),'norm_over_Lnorm_pi_minus4':str(const),'B_over_pi4_norm':str(Bconst),'E1_constant':str(g[0]),'sturm_pivots':inds,'M':sm(M),'T7':sm(T),'H_coeffs':[str(z) for z in H],'elapsed':time.time()-t}
 (HERE/'rankin_projection.json').write_text(json.dumps(out,indent=2)+'\n')
 print('elapsed',time.time()-t,flush=True)
if __name__=='__main__':run()
