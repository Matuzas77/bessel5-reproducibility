from pathlib import Path
ROOT=Path(__file__).resolve().parent
import sympy as s,json
D=json.load(open(str(ROOT/'results/modsym_data.json')))
T={p:s.Matrix([[s.Rational(x) for x in row]for row in v])for p,v in D['T'].items()}
w=s.Matrix([[s.Rational(x)for x in D['winding']]])
t=s.Matrix([[s.Rational(x)for x in D['twist5']]])
U=T['3'];G=T['7']**2+64*s.eye(48)
pol=s.Rational(17,10)*s.eye(48)-s.Rational(13,30)*U+s.Rational(1,10)*U**2-s.Rational(1,30)*U**3
q=t-w*pol
sol=s.linsolve((G.T,q.T));x=list(sol)[0]
syms=set().union(*(a.free_symbols for a in x))
x=s.Matrix([[a.subs({t:0 for t in syms})for a in x]])
assert x*G==q
print('witness',[(k,str(a))for k,a in enumerate(x) if a]);print('count',sum(a!=0 for a in x))
print('manin symbols',[(D['reps'][D['free'][k]//2], 'X' if D['free'][k]%2 else 'Y',str(a))for k,a in enumerate(x)if a])
D['twist_witness']=[str(a)for a in x]
D['twist_polynomial']=['17/10','-13/30','1/10','-1/30']
open(str(ROOT/'results/modsym_data.json'),'w').write(json.dumps(D))
