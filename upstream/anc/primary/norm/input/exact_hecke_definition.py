"""Define the normalized cusp form from an eta carrier and finite Hecke algebra.
No Bessel data or numerical L values enter this code.
"""
from pathlib import Path
from math import isqrt,gcd
import json,time
import sympy as s
from recompute_modular_value import TERMS,eta_coeff,coefficients

def sparse_eta_coeff(N):
 ans=[0]*(N+1)
 for c,exps in TERMS:
  shift=sum(d*r for d,r in exps.items())//24
  nmax=N-shift;v=[0]*(nmax+1);v[0]=1
  for d,r in exps.items():
   sparse=[(0,1)];j=1
   while d*j*(3*j-1)//2<=nmax:
    for k in [d*j*(3*j-1)//2,d*j*(3*j+1)//2]:
     if k<=nmax:sparse.append((k,(-1)**j))
    j+=1
   for _ in range(r):
    nxt=[0]*(nmax+1)
    for offset,sg in sparse:
     for k in range(nmax-offset+1):nxt[k+offset]+=sg*v[k]
    v=nxt
  for n,x in enumerate(v):ans[n+shift]+=c*x
 return ans

def chi60(n):return int(s.kronecker_symbol(-15,n)) if gcd(n,60)==1 else 0

def derive():
 t=time.time();N=2916;F=sparse_eta_coeff(N)
 assert F[:181]==eta_coeff(180)
 W=s.Matrix(36,4,lambda m,j:F[3**j*(m+1)])
 assert W.rank()==4
 # Exact eta modularity checks
 ds=s.divisors(60)
 for _,ex in TERMS:
  assert sum(ex.values())==6
  assert sum(d*r for d,r in ex.items())%24==0
  assert sum((60//d)*r for d,r in ex.items())%24==0
  prod=s.prod(d**r for d,r in ex.items());assert s.sqrt(s.Rational(prod,15)).is_Integer
  orders=[s.Rational(60,24)*sum(s.Rational(gcd(c,d)**2*r,gcd(c,60//c)*c*d) for d,r in ex.items()) for c in ds]
  assert min(orders)>0
 assert all(F[81*n]-2*F[9*n]+81*F[n]==0 for n in range(1,37))
 matrices={}
 # Hecke algebra is spanned by T_1,...,T_36 (Sturm-duality).
 # T_n commutes with U_3. Thus T_n F in the U_3 cyclic span proves stability.
 for n in range(1,37):
  target=s.Matrix([sum(chi60(d)*d*d*F[m*n//(d*d)] for d in s.divisors(gcd(m,n))) for m in range(1,37)])
  col=W.gauss_jordan_solve(target)[0]
  # U3 matrix in cyclic basis F,U3F,U3²F,U3³F
  U=s.Matrix([[0,0,0,-81],[1,0,0,0],[0,1,0,2],[0,0,1,0]])
  T=s.Matrix.hstack(*[U**j*col for j in range(4)])
  matrices[n]=T
 sqrt5=s.sqrt(5);a3=sqrt5-2*s.I
 c=(U-a3*s.eye(4)).nullspace()[0]
 norm=(s.Matrix([[F[3**j] for j in range(4)]])*c)[0]
 c=c.applyfunc(lambda x:s.expand(s.radsimp(x/norm)))
 print('f=sum c_j U3^j F',list(c),flush=True)
 eig=[]
 for n,T in matrices.items():
  Tc=T*c
  av=s.expand(s.radsimp(Tc[0]/c[0]))
  assert all(s.simplify(x)==0 for x in Tc-av*c)
  eig.append(av)
 # identify all four rational components and form coefficients through 108.
 ca=coefficients(F[:109]);evalf=lambda a:a[0]+sqrt5*a[1]+s.I*a[2]+s.I*sqrt5*a[3]
 assert all(s.simplify(sum(c[j]*F[3**j*n] for j in range(4))-evalf(ca[n]))==0 for n in range(1,109))
 print('modularity, all Sturm Hecke generators, 108 q coefficients PASS',flush=True)
 out={'status':'PASS','sturm_bound':36,'cyclic_rank':4,'U3_minimal_polynomial':'X^4-2X^2+81','projection_coefficients':[str(x) for x in c],'first_36_Hecke_eigenvalues':[str(x) for x in eig],'first_108_coefficients_match':True,'scope':'Eta modularity and finite Hecke stability prove a normalized simultaneous eigenform. The identification of the label uses the external database record. No Bessel-period equality is proved.','elapsed':time.time()-t}
 Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))
 print('elapsed',time.time()-t,flush=True)
 return c
if __name__=='__main__':derive()
