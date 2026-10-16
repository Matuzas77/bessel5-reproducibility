"""Directed fixed-point intervals, standard library only.

Every endpoint is an integer divided by Context.scale. Transcendental functions
use convergent real series and explicit remainders. This is not a claim of a
formal proof-assistant implementation. The interval algorithm is documented in
REPORT.md and checked by a separate rational-series implementation.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt


def ceildiv(n: int, d: int) -> int:
    if d <= 0:
        raise ValueError('positive divisor required')
    return -((-n)//d)

class Context:
    def __init__(self, digits: int = 120):
        if digits < 30:
            raise ValueError('use at least 30 fixed-point digits')
        self.digits=digits; self.scale=10**digits
        self._pi=None; self._log2=None
    def ball(self, value=0):
        if isinstance(value, Ball):
            if value.ctx is not self: raise ValueError('different contexts')
            return value
        f=Fraction(value)
        return Ball(self,(f.numerator*self.scale)//f.denominator,
                    ceildiv(f.numerator*self.scale,f.denominator))
    def from_bounds(self, lo,hi):
        return Ball(self,lo,hi)
    def pi(self):
        if self._pi is None:
            self._pi=16*self._atan_inverse(5)-4*self._atan_inverse(239)
        return self._pi
    def _atan_inverse(self,q):
        x=self.ball(Fraction(1,q)); term=x; x2=x*x; total=self.ball(0)
        n=0
        while True:
            add=term/(2*n+1)
            total=total+add if n%2==0 else total-add
            term=term*x2; n+=1
            nxt=term/(2*n+1)
            if nxt.abs_upper() <= 2:
                return total.widen(nxt.abs_upper())
            if n>self.digits*4:raise ArithmeticError('atan series limit')
    def log2(self):
        if self._log2 is None:self._log2=self.ball(2)._log_near()
        return self._log2

@dataclass(frozen=True)
class Ball:
    ctx: Context
    lo: int
    hi: int
    def __post_init__(self):
        if self.lo>self.hi:raise ValueError('reversed interval')
    def _b(self,other):return self.ctx.ball(other)
    def __add__(self,other):
        b=self._b(other);return Ball(self.ctx,self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return Ball(self.ctx,-self.hi,-self.lo)
    def __sub__(self,other):return self+-self._b(other)
    def __rsub__(self,other):return self._b(other)+-self
    def __mul__(self,other):
        b=self._b(other); ss=self.ctx.scale
        v=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi]
        return Ball(self.ctx,min(v)//ss,ceildiv(max(v),ss))
    __rmul__=__mul__
    def inv(self):
        if self.lo<=0<=self.hi:raise ZeroDivisionError('interval contains zero')
        s2=self.ctx.scale**2
        if self.lo>0:return Ball(self.ctx,s2//self.hi,ceildiv(s2,self.lo))
        return -(-self).inv()
    def __truediv__(self,other):return self*self._b(other).inv()
    def __rtruediv__(self,other):return self._b(other)*self.inv()
    def __pow__(self,n):
        if not isinstance(n,int):raise TypeError('integer powers only')
        if n<0:return self.inv()**(-n)
        z=self.ctx.ball(1);b=self
        while n:
            if n&1:z=z*b
            b=b*b;n//=2
        return z
    def sqrt(self):
        if self.lo<0:raise ValueError('sqrt of negative interval')
        s=self.ctx.scale
        l=isqrt(self.lo*s);v=self.hi*s;u=isqrt(v)
        if u*u<v:u+=1
        return Ball(self.ctx,l,u)
    def abs_upper(self):return max(abs(self.lo),abs(self.hi))
    def widen(self,n):return Ball(self.ctx,self.lo-n,self.hi+n)
    def _log_near(self):
        if self.lo<=0:raise ValueError('log domain')
        z=(self-1)/(self+1);z2=z*z;term=z;total=self.ctx.ball(0);n=0
        # Valid for z of either sign; bound the remainder by its absolute value.
        q=Ball(self.ctx,0,z2.abs_upper())
        if q.hi>=self.ctx.scale:raise ValueError('log series too far from 1')
        while True:
            total=total+term/(2*n+1)
            term=term*z2;n+=1
            tail=2*Ball(self.ctx,0,term.abs_upper())/((2*n+1)*(1-q))
            if tail.hi<=4:return (2*total).widen(tail.hi)
            if n>self.ctx.digits*8:raise ArithmeticError('log series limit')
    def log(self):
        if self.lo<=0:raise ValueError('log domain')
        s=self.ctx.scale;middle=(self.lo+self.hi)//2
        k=0
        while middle>=2*s:middle//=2;k+=1
        while middle<s:middle*=2;k-=1
        y=self/(2**k) if k>=0 else self*(2**(-k))
        out=y._log_near()
        return out+k*self.ctx.log2() if k else out
    def dilog_positive(self):
        # Used only for 0<=z<=~1/2; exact geometric remainder.
        if self.lo<0 or self.hi*4>3*self.ctx.scale:
            raise ValueError('positive Li2 argument outside [0,3/4]')
        total=self.ctx.ball(0);power=self;n=1
        while True:
            total=total+power/(n*n)
            power=power*self;n+=1
            tail=power/(n*n*(1-self))
            if tail.hi<=4:return Ball(self.ctx,total.lo,total.hi+tail.hi)
            if n>self.ctx.digits*16:raise ArithmeticError('Li2 series limit')
    def li2_negative(self):
        """Returns Li_2(-self), for self>=0, real branch."""
        if self.lo<0:raise ValueError('negative Li2 magnitude')
        if self.hi==0:return self.ctx.ball(0)
        if self.lo>self.ctx.scale:
            v=self.inv()
            return -self.ctx.pi()**2/6-self.log()**2/2-v.li2_negative()
        y=self/(1+self)
        return -y.dilog_positive()-(1+self).log()**2/2
    def atanh(self):
        if self.lo<=-self.ctx.scale or self.hi>=self.ctx.scale:
            raise ValueError('atanh domain')
        return ((1+self).log()-(1-self).log())/2
    def contains(self,x):
        q=Fraction(x)*self.ctx.scale
        return self.lo<=q<=self.hi
    def width_fraction(self):return Fraction(self.hi-self.lo,self.ctx.scale)
    def exact_json(self):
        return {'lower_integer':str(self.lo),'upper_integer':str(self.hi),
                'denominator_power10':self.ctx.digits,
                'width_integer':str(self.hi-self.lo)}
