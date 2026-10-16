from pathlib import Path
ROOT=Path(__file__).resolve().parent
from fractions import Fraction as F
from math import gcd
import json
# Arithmetic in Q(sqrt(5),i), coordinates 1,r,i,ir.
def K(a=0,b=0,c=0,d=0): return tuple(F(x) for x in (a,b,c,d))
def plus(x,y):return tuple(a+b for a,b in zip(x,y))
def neg(x):return tuple(-a for a in x)
def scale(x,a):return tuple(a*b for b in x)
def times(x,y):
 out=[F(0)]*4
 for a,v in enumerate(x):
  for b,w in enumerate(y):
   ir=(a//2)+(b//2);rr=(a%2)+(b%2)
   out[(ir%2)*2+rr%2]+=v*w*(-1)**(ir//2)*5**(rr//2)
 return tuple(out)
def sigma(x):return (x[0],-x[1],-x[2],x[3])
Z=K();ONE=K(1);R=K(0,1);II=K(0,0,1)
ome=K(0,F(1,5),0,F(-2,5))
eps=K(F(1,3),F(4,15),F(-2,3),F(2,15))
omebar=K(0,F(1,5),0,F(2,5))
coeff=(ONE,ome,scale(eps,3),scale(times(eps,omebar),3))
a3=K(0,1,-2);a7=K(0,0,-8);a5=K(0,1,0,2)
def chi15(n):
 n%=15
 return 0 if gcd(n,15)>1 else (1 if n%3==1 else -1)*(1 if n%5 in(1,4)else -1)
def chi5(n):return 0 if n%5==0 else (1 if n%5 in(1,4)else -1)
exps=[{2:1,6:2,10:4,30:-1},{2:4,6:-1,10:1,30:2},{2:-1,6:4,10:2,30:1},{2:2,6:1,10:-1,30:4}]

def eta_product(ex,M):
 shift=sum(d*a for d,a in ex.items())//24
 assert sum(d*a for d,a in ex.items())%24==0
 L=[0]*(M+1)
 for d,r in ex.items():
  for k in range(1,M//d+1):
   # q d/dq log product = -sum r d k q^(dk)/(1-q^(dk))
   for n in range(d*k,M+1,d*k): L[n]-=r*d*k
 b=[0]*(M+1);b[0]=1
 for n in range(1,M-shift+1):
  z=sum(L[j]*b[n-j] for j in range(1,n+1));assert z%n==0
  b[n]=z//n
 return [0]*shift+b[:M+1-shift]

M=252
P=[eta_product(ex,M) for ex in exps]
f=[tuple(sum(F(P[j][n])*coeff[j][k]for j in range(4))for k in range(4))for n in range(M+1)]
assert f[1]==ONE and f[3]==a3 and f[5]==a5 and f[7]==a7
for n in range(1,37):
 assert f[3*n]==times(a3,f[n]),('U3',n)
 got=plus(f[7*n],scale(f[n//7],chi15(7)*49)) if n%7==0 else f[7*n]
 assert got==times(a7,f[n]),('T7',n)
 assert f[5*n]==times(a5,f[n]),('U5',n)
for n in range(1,181):
 got=plus(sigma(f[n]),neg(times(sigma(a5),sigma(f[n//5])))) if n%5==0 else sigma(f[n])
 assert scale(f[n],chi5(n))==got,('twist',n)
# cusp orders for four eta products
divs=[d for d in range(1,61) if 60%d==0]
cusps=[]
character_products=[]
for ex in exps:
 prod=F(1)
 for d,rr in ex.items():prod*=F(d)**rr
 expected_root=(40,40,360,360)[len(character_products)]
 assert prod==15*expected_root**2
 character_products.append(str(prod))
 assert sum(ex.values())==6
 assert sum(d*r for d,r in ex.items())%24==0
 assert sum((60//d)*r for d,r in ex.items())%24==0
 orders=[sum(F(60*gcd(c,d)**2*r,24*gcd(c*c,60)*d)for d,r in ex.items())for c in divs]
 assert all(o>0 for o in orders)
 cusps.append([str(o)for o in orders])
kap=K(F(-4,10),F(4,10),F(-3,10),F(3,10))
p=[F(17,10),F(-13,30),F(1,10),F(-1,30)]
value=Z; power=ONE
for a in p:value=plus(value,scale(power,a));power=times(power,a3)
expected=times(R,times(plus(ONE,scale(sigma(a5),F(-1,25))),kap))
assert value==expected
out={'coefficients_checked':M,'eigenform_Sturm_bound':36,'twist_level':300,'twist_Sturm_bound':180,'four_product_character_products':character_products,'four_product_cusp_denominators':divs,'four_product_cusp_orders':cusps,'fixed_coefficients':{str(n):list(map(str,f[n]))for n in range(M+1)},'critical_ratio_kappa':list(map(str,kap)),'twisted_modular_symbol_eigenvalue':list(map(str,value)),'status':'pass','Bessel_identity_proved':False}
open(str(ROOT/'results/modular_coefficients.json'),'w').write(json.dumps(out,indent=2))
print('PASS eta cusp conditions, U3,T7,U5 eigenidentities, twist to level300 through q^180, exact critical-value factor.')
