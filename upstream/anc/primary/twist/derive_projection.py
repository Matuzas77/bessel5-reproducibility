from pathlib import Path
ROOT=Path(__file__).resolve().parent
import sympy as s,json,time
D=json.load(open(str(ROOT/'results/modsym_data.json')))
T={p:s.Matrix([[s.Rational(x) for x in row]for row in v])for p,v in D['T'].items()}
w=s.Matrix([[s.Rational(x)for x in D['winding']]])
t=s.Matrix([[s.Rational(x)for x in D['twist5']]])
B=s.Matrix.hstack(*(T['7']**2+64*s.eye(48)).nullspace())
print('kernel',B.shape,flush=True)
rows=(B.T.rref()[1]);print('indep rows',rows)
Br=B[list(rows),:]
U=Br.inv()*(T['3']*B)[list(rows),:]
assert T['3']*B==B*U
T7=Br.inv()*(T['7']*B)[list(rows),:]
print('Uchar',s.factor(U.charpoly().as_expr()),flush=True)
print('T7relation',T7-s.Rational(2,9)*(7*U+U**3)==s.zeros(8),flush=True)
w8=w*B;t8=t*B
print('w8',w8,'t8',t8,flush=True)
# solve t8=w8*R(U) for polynomial degree<=3
M=s.Matrix.vstack(w8,w8*U,w8*U**2,w8*U**3).T
print('ranks',M.rank(),M.row_join(t8.T).rank())
sol=s.linsolve((M,t8.T));print('sol',sol)
r=s.sqrt(5);i=s.I
lamb=r-2*i
kap=(r-1)*(4+3*i)/10
a5sig=-r+2*i*r
c=s.expand(r*(1-a5sig/25)*kap)
print('predicted scalar',s.simplify(c))
if sol is not s.EmptySet:
 coeff=list(sol)[0];pred=sum(a*lamb**k for k,a in enumerate(coeff))
 print('calculated scalar',s.simplify(s.expand(pred)),'difference',s.simplify(pred-c))
print('all U eigen proof',flush=True)
# t8-w8poly equality polynomial valid all8 vectors -> no algebraic field numerics
out={'kernel_basis':[[str(x)for x in row]for row in B.tolist()], 'independent_rows':list(rows),'U':[[str(x)for x in row]for row in U.tolist()], 'w8':[str(x)for x in w8], 't8':[str(x)for x in t8], 'polynomial':[str(x)for x in list(sol)[0]] if sol is not s.EmptySet else None,'comparison_scalar':str(c)}
open(str(ROOT/'results/modsym_projection.json'),'w').write(json.dumps(out,indent=2))
