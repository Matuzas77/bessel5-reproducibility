"""Eta cusp forms used in the exact finite certificate. Modularity is proved in REPORT.md."""
from fractions import Fraction
from math import isqrt
from pathlib import Path
import json, sys
import sympy as s
sys.path.insert(0,str(Path(__file__).parent))
from kernel_coefficients import kernel_coefficient, chi4
DS=(1,2,3,4,5,6,10,12,15,20,30,60)
EXPS=[
[3,0,0,0,0,0,0,0,3,0,0,0],
[2,0,0,1,0,0,0,0,2,0,0,1],
[2,0,0,0,0,1,1,0,2,0,0,0],
[1,0,0,2,0,0,0,0,1,0,0,2],
[1,0,0,1,0,1,1,0,1,0,0,1],
[1,0,0,0,0,2,2,0,1,0,0,0],
[0,2,0,0,0,0,0,1,0,1,2,0],
[0,1,0,0,0,0,0,2,0,2,1,0],
[0,0,3,0,3,0,0,0,0,0,0,0],
[0,0,2,0,2,0,0,1,0,1,0,0],
[0,0,0,3,0,0,0,0,0,0,0,3],
[0,0,0,0,0,0,0,3,0,3,0,0],
[-1,3,3,0,0,0,1,0,0,0,0,0],
[-1,2,3,0,0,1,0,0,0,0,1,0],
[-1,1,4,0,1,1,0,0,0,0,0,0],
[-1,1,3,0,0,1,0,1,0,1,0,0],
[-1,0,3,1,0,0,0,0,0,3,0,0],
[2,-1,0,0,0,2,0,0,2,0,1,0]]

def eta_series(e,N):
 sh=sum(d*v for d,v in zip(DS,e));assert sh%24==0
 sh//=24;assert sh>=0
 M=N-sh;b=[0]*(M+1)
 for d,v in zip(DS,e):
  if not v:continue
  for m in range(d,M+1,d):
   for k in range(m,M+1,m):b[k]-=v*m
 a=[1]+[0]*M
 for n in range(1,M+1):
  v=sum(b[k]*a[n-k] for k in range(1,n+1));assert v%n==0
  a[n]=v//n
 return [0]*sh+a

