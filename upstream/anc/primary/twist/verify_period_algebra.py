from pathlib import Path
import sympy as s,json
ROOT=Path(__file__).resolve().parent
r=s.sqrt(5);I=s.I
x,y,z,w=s.symbols('l1 l2 l3 l4',real=True)
om=(1-2*I)/r;ep=(5+4*r+I*(2*r-10))/15
coeff=[1,om,3*ep,3*ep*s.conjugate(om)]
def sig(c):
    return s.expand(s.conjugate(s.expand(c)).subs(r,-r))
L=sum(c*v for c,v in zip(coeff,(x,y,z,w)))
Ls=sum(sig(c)*v for c,v in zip(coeff,(x,y,z,w)))
kap=(r-1)*(4+3*I)/10
q=s.expand(Ls-kap*L)
real=s.simplify(s.re(q));imag=s.simplify(s.im(q))
M=s.Matrix([[s.diff(eq,v)for v in (y,w)]for eq in (real,imag)])
assert s.simplify(M.det())!=0
sol=s.solve([real,imag],(y,w),simplify=True)
expected={y:((3+r)*x-9*(1+r)*z)/4,w:((1-r)*x+(5+r)*z)/4}
assert all(s.simplify(sol[k]-v)==0 for k,v in expected.items())
gamma=(125-13*r+I*(125-41*r))/4
Q=s.expand(s.re(gamma*L**2))
q0=0;q1=0
for mon,c in s.Poly(Q,x,y,z,w).terms():
 c=s.expand(s.simplify(c));aa=c.coeff(r,0);bb=c.coeff(r,1)
 assert s.simplify(c-aa-r*bb)==0
 term=x**mon[0]*y**mon[1]*z**mon[2]*w**mon[3]
 q0+=aa*term;q1+=bb*term
assert s.simplify(q1.subs(expected))==0
red=s.expand(s.simplify(Q.subs(expected)))
assert s.simplify(red-25*((1-r)*x*x+36*x*z-27*(1+r)*z*z))==0
Lred=s.simplify(s.expand(L.subs(expected)))
assert s.simplify(Lred-(2*r*(1-2*I)*x/5+(8*r/5+I*(5+9*r/5))*z))==0
# The conjugate value in the other embedding; coefficient conjugation is not a
# Galois operation on arbitrary transcendental L-values.
Lr=s.expand(L).subs(r,-r)
assert s.simplify((Lr-s.conjugate(kap)*s.conjugate(L)).subs(expected))==0
# Companion candidates, exact simplification on the same modular plane.
B=s.expand(15*(r-1)*L*s.conjugate(L)/2) # B*/sqrt3
C=s.expand(((42+6*r)*s.re(L**2)+(33*r-69)*s.im(L**2))/4) # pi C*
assert s.simplify(B.subs(expected)-30*((r-1)*x*x-8*x*z+9*(1+r)*z*z))==0
assert s.simplify(C.subs(expected)-30*((1-r)*x*x+18*x*z-9*(1+r)*z*z))==0
# Candidate conic, not an evaluation of any Bessel period.
BB=s.simplify(B.subs(expected));CC=s.simplify(C.subs(expected))
assert s.simplify(12*red**2-40*red*CC+40*CC**2-15*BB**2)==0
# Clear rational quadratic relation, with no radicals in its coefficients.
Qnull=s.expand(4*q1)
mat=s.hessian(Qnull,(x,y,z,w))/2
assert mat.det()!=0
print('PASS exact twist-to-period plane, nonzero 2x2 determinant, rational quadratic cancellation, all three candidate forms, candidate conic.')
json.dump({'status':'pass','linear_relation_determinant':str(s.simplify(M.det())), 'linear_relations':{str(k):str(v)for k,v in expected.items()},'rational_null_quadric':str(Qnull),'null_quadric_matrix':[[str(x)for x in rr]for rr in mat.tolist()],'candidate_A_over_pi':str(red),'candidate_B_over_sqrt3':str(BB),'candidate_pi_C':str(CC),'L_in_two_periods':str(s.expand(Lred)),'Bessel_identity_proved':False},open(ROOT/'results/period_algebra.json','w'),indent=2)
