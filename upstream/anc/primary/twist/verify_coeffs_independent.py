#!/usr/bin/env python3
"""Independent coefficients: pentagonal eta products and binary theta enumeration.
No symbolic algebra package and no primary coefficient result file is used.
"""
from fractions import Fraction as F
from math import isqrt,gcd
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
N=252

def conv(a,b,n=N):
 c=[0]*(n+1)
 for i,v in enumerate(a):
  if v:
   for j,w in enumerate(b[:n+1-i]):
    if w:c[i+j]+=v*w
 return c

def pentagonal(d):
 a=[0]*(N+1)
 # |j| <= N+1 is deliberately generous; exact cutoff below.
 for j in range(-N-1,N+2):
  k=d*j*(3*j-1)//2
  if k<=N:a[k]=(-1)**j if j>=0 else (-1)**(-j)
 return a

e2=pentagonal(2);e10=pentagonal(10)
h0=conv(conv(e2,e2),conv(e10,e10))
h=[0]+h0[:-1]
th1=[0]*(N+1);th2=[0]*(N+1)
R=2*isqrt(N)+2
# Both Q1,Q2 >= (x*x+y*y)/2; this box contains every relevant pair.
for x in range(-R,R+1):
 for y in range(-R,R+1):
  a=x*x+x*y+4*y*y;b=2*x*x+x*y+2*y*y
  if a<=N:th1[a]+=1
  if b<=N:th2[b]+=1
g=[(a+b)//2 for a,b in zip(th1,th2)];d=[(a-b)//2 for a,b in zip(th1,th2)]
assert all((a+b)%2==0 for a,b in zip(th1,th2))
g2=[g[n//2] if n%2==0 else 0 for n in range(N+1)]
d2=[d[n//2] if n%2==0 else 0 for n in range(N+1)]
h3=[h[n//3] if n%3==0 else 0 for n in range(N+1)]
Ps=[conv(h,g2),conv(h,d2),conv(h3,g2),conv(h3,d2)]

# Field coordinates are (a+br)+(c+dr)i. Multiplication is via quadratic pairs,
# rather than the primary four-basis convolution.
def qadd(x,y):return(x[0]+y[0],x[1]+y[1])
def qmul(x,y):return(x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def qneg(x):return(-x[0],-x[1])
def C(a=0,b=0,c=0,d=0):return((F(a),F(b)),(F(c),F(d)))
def add(x,y):return(qadd(x[0],y[0]),qadd(x[1],y[1]))
def mul(x,y):return(qadd(qmul(x[0],y[0]),qneg(qmul(x[1],y[1]))),qadd(qmul(x[0],y[1]),qmul(x[1],y[0])))
def smul(a,x):return((a*x[0][0],a*x[0][1]),(a*x[1][0],a*x[1][1]))
def sig(x):return((x[0][0],-x[0][1]),(-x[1][0],x[1][1]))
zero=C();one=C(1);rr=C(0,1)
ome=C(0,F(1,5),0,F(-2,5));ep=C(F(1,3),F(4,15),F(-2,3),F(2,15))
cs=[one,ome,smul(3,ep),C(0,1,2)]
f=[]
for n in range(N+1):
 v=zero
 for p,c in zip(Ps,cs):v=add(v,smul(p[n],c))
 f.append(v)
lam=C(0,1,-2);a5=C(0,1,0,2);a7=C(0,0,-8)
assert f[1]==one and f[3]==lam and f[5]==a5 and f[7]==a7
for n in range(1,37):
 assert f[3*n]==mul(lam,f[n])
 val=add(f[7*n],smul(-49,f[n//7])) if n%7==0 else f[7*n]
 assert val==mul(a7,f[n])
 assert f[5*n]==mul(a5,f[n])
for n in range(1,181):
 chi=0 if n%5==0 else (1 if n%5 in(1,4)else -1)
 rhs=add(sig(f[n]),smul(-1,mul(sig(a5),sig(f[n//5])))) if n%5==0 else sig(f[n])
 assert smul(chi,f[n])==rhs

kap=C(F(-2,5),F(2,5),F(-3,10),F(3,10))
real_rows=[add(sig(c),smul(-1,mul(kap,c)))for c in cs]
det=qadd(qmul(real_rows[1][0],real_rows[3][1]),qneg(qmul(real_rows[3][0],real_rows[1][1])))
assert det==(F(4),F(-4,5))
# l2=Ay*x+By*z; l4=Aw*x+Bw*z
Ay=C(F(3,4),F(1,4));By=C(F(-9,4),F(-9,4))
Aw=C(F(1,4),F(-1,4));Bw=C(F(5,4),F(1,4))
for j,coords in enumerate(([one,Ay,zero,Aw],[zero,By,one,Bw])):
 val=zero
 for c,p in zip(real_rows,coords):val=add(val,mul(c,p))
 assert val==zero

Lx=add(add(cs[0],mul(cs[1],Ay)),mul(cs[3],Aw))
Lz=add(add(cs[2],mul(cs[1],By)),mul(cs[3],Bw))
gamma=C(F(125,4),F(-13,4),F(125,4),F(-41,4))
actual=[mul(gamma,mul(Lx,Lx))[0],mul(gamma,smul(2,mul(Lx,Lz)))[0],mul(gamma,mul(Lz,Lz))[0]]
assert actual==[(F(25),F(-25)),(F(900),F(0)),(F(-675),F(-675))]
p=[F(17,10),F(-13,30),F(1,10),F(-1,30)];v=zero;power=one
for a in p:v=add(v,smul(a,power));power=mul(power,lam)
assert v==mul(rr,mul(add(one,smul(F(-1,25),sig(a5))),kap))
print('PASS independent theta/pentagonal eigenform and twist identities, determinant, two-period plane and candidate quadratic.')
json.dump({'status':'pass','coefficient_method':'pentagonal eta and binary theta enumeration','max_coefficient':N,'twist_bound':180,'period_plane_verified':True,'uses_primary_coefficients':False,'Bessel_identity_proved':False},open(ROOT/'results/independent_coefficients.json','w'),indent=2)
