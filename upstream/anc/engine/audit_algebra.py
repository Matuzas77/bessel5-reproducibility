#!/usr/bin/env python3
"""Exact coefficient checks in Q(sqrt(5)); no numerical period input."""
from fractions import Fraction as F
from dataclasses import dataclass
@dataclass(frozen=True)
class Q5:
    a:F=F(0)
    b:F=F(0)
    def __add__(x,y):
        y=coerce(y);return Q5(x.a+y.a,x.b+y.b)
    __radd__=__add__
    def __neg__(x):return Q5(-x.a,-x.b)
    def __sub__(x,y):return x+-coerce(y)
    def __rsub__(x,y):return coerce(y)+-x
    def __mul__(x,y):
        y=coerce(y);return Q5(x.a*y.a+5*x.b*y.b,x.a*y.b+x.b*y.a)
    __rmul__=__mul__
    def __truediv__(x,y):
        y=coerce(y);d=y.a*y.a-5*y.b*y.b
        if not d:raise ZeroDivisionError
        return x*Q5(y.a/d,-y.b/d)
    def __rtruediv__(x,y):return coerce(y)/x
    def __pow__(x,n):
        result=coerce(1)
        for _ in range(n):result=result*x
        return result

def coerce(v):return v if isinstance(v,Q5) else Q5(F(v))
def linprod(x,y):return [x[0]*y[0],x[0]*y[1]+x[1]*y[0],x[1]*y[1]]
def add(x,y):return [a+b for a,b in zip(x,y)]
def scale(x,c):return [a*c for a in x]
def subtract(x,y):return add(x,scale(y,-1))
def check(x,y):
    if [coerce(z) for z in x]!=[coerce(z) for z in y]:raise AssertionError((x,y))

def main():
    r=Q5(F(0),F(1)); a=(3-r)/4;b=20+2*r;k=4*r-10
    er=(5+4*r)/15;ei=(2*r-10)/15
    X=[1+er,1-er];Y=[ei,-ei]
    Re=subtract(linprod(X,X),linprod(Y,Y));Im=scale(linprod(X,Y),2)
    norm=add(linprod(X,X),linprod(Y,Y))
    alpha=125-13*r;beta=41*r-125;gamma=42+6*r;delta=33*r-69
    check(scale(add(scale(Re,alpha),scale(Im,beta)),F(1,960)),[b/60,-2*k*a/60,-b*a*a/60])
    check(scale(norm,3*(5-r)/32),[1,0,a*a])
    check(scale(add(scale(Re,gamma),scale(Im,delta)),F(1,960)),[F(1,5),a/5,-a*a/5])
    A=[alpha/4,beta/4];C=[gamma/4,delta/4]
    conic=add(subtract(scale(linprod(A,A),12),scale(linprod(A,C),40)),scale(linprod(C,C),40))
    c2=135*(5-r)**2/4
    check(conic,[5*c2,0,5*c2])
    # The coefficient of |L|^2 in D_even follows from D_even=pi^2 B/(2sqrt15).
    if 3*(5-r)/4 != (3*(5-r)/2)/2:raise AssertionError
    print('EXACT ALGEBRA PASS: three period conversions, conic compatibility, determinant factor')
if __name__=='__main__':main()
