#!/usr/bin/env python3
"""Compare the paper's formulas using freshly generated rational intervals."""
from __future__ import annotations
import json
from pathlib import Path
from fractions import Fraction
from fixed_interval import Context, ceildiv
ROOT=Path(__file__).resolve().parents[1]

def load(name):
    return json.loads((ROOT/'outputs'/name).read_text())

def import_interval(ctx,record):
    p=record['denominator_power10']; lo=int(record['lower_integer']); hi=int(record['upper_integer'])
    if p<=ctx.digits:
        scale=10**(ctx.digits-p); return ctx.from_bounds(lo*scale,hi*scale)
    scale=10**(p-ctx.digits); return ctx.from_bounds(lo//scale,ceildiv(hi,scale))

def bound(v):
    m=v.abs_upper()
    if not m: raise AssertionError('An unexpectedly exact zero is not a numerical certificate')
    exponent=len(str(m))-v.ctx.digits
    if not m<10**(v.ctx.digits+exponent): raise AssertionError('Invalid strict bound')
    return exponent

def decimal_interval(v,places=45):
    factor=10**places; scale=v.ctx.scale
    def fmt(n):
        sign='-' if n<0 else ''; integer,remainder=divmod(abs(n),factor)
        return f'{sign}{integer}.{remainder:0{places}d}'
    return [fmt(v.lo*factor//scale),fmt(ceildiv(v.hi*factor,scale))]

def main():
    ctx=Context(1100); pi=ctx.pi(); r=ctx.ball(5).sqrt(); r15=ctx.ball(15).sqrt()
    primary=load('source_primary.json'); binary=load('source_binary.json')
    direct=load('modular_direct.json'); eta=load('eta_binary.json')
    x=import_interval(ctx,direct['L_real']); y=import_interval(ctx,direct['L_imaginary'])
    real_square=x*x-y*y; imag_square=2*x*y; norm=x*x+y*y
    targets={
       'A0':pi/4*((125-13*r)*real_square+(41*r-125)*imag_square),
       'C0':((42+6*r)*real_square+(33*r-69)*imag_square)/(4*pi),
       'D5_odd':3*r15*(5-r)/2*norm,
       'D5_even':3*pi*pi*(5-r)/4*norm,
    }
    # Separate modular implementation: binary arithmetic, q^2 Euler expansion,
    # Cauchy radius 2/3, and the independently normalized U,V fold.
    U=import_interval(ctx,eta['quantities']['U']); V=import_interval(ctx,eta['quantities']['V'])
    er=(5+4*r)/15; ei=(2*r-10)/15; factor=pi/(2*ctx.ball(60).sqrt())
    xb=factor*((1+er)*U+(1-er)*V); yb=factor*ei*(U-V)
    for name,v in [('L_real',x-xb),('L_imaginary',y-yb)]:
        if not v.contains(0) or bound(v)>-400:
            raise AssertionError('Modular cross-check failed: '+name)
    rows={}; nesting=[]
    for label,source in [('decimal',primary),('binary',binary)]:
        A=[import_interval(ctx,z) for z in source['periods']['A']]
        C=[import_interval(ctx,z) for z in source['periods']['C']]
        B=[import_interval(ctx,z) for z in source['periods']['B_pairing_derived']]
        scaled_C=pi*pi*C[0]
        odd=((12*A[0]**2-40*A[0]*scaled_C+40*scaled_C**2)/5).sqrt()/pi
        if not (odd-B[0]).contains(0):
            raise AssertionError('Conic and full-matrix odd reconstructions disagree')
        actual={'A0':A[0],'C0':C[0],'D5_odd':odd,'D5_even':A[0]*C[1]-A[1]*C[0]}
        errors={name:actual[name]-target for name,target in targets.items()}
        for name,v in errors.items():
            if not v.contains(0) or bound(v)>-300:
                raise AssertionError('The conjecture failed the stated check: '+label+' '+name)
        determinant_identity=actual['D5_even']-pi*pi/(2*r15)*odd
        if not determinant_identity.contains(0):
            raise AssertionError('Even/odd determinant normalization failed')
        rows[label]={
          'actual_intervals':{k:v.exact_json() for k,v in actual.items()},
          'actual_decimal_intervals':{k:decimal_interval(v) for k,v in actual.items()},
          'residuals':{k:v.exact_json() for k,v in errors.items()},
          'strict_absolute_bound_exponents':{k:bound(v) for k,v in errors.items()},
          'determinant_relation_residual':determinant_identity.exact_json(),
        }
    for sector in ['A','C','B_pairing_derived']:
        for j in range(3):
            p=import_interval(ctx,primary['periods'][sector][j]); b=import_interval(ctx,binary['periods'][sector][j])
            if not p.lo<=b.lo<=b.hi<=p.hi:
                raise AssertionError('Independent source interval did not nest')
            nesting.append(f'{sector}[{j}]')
    result={
       'status':'PASS',
       'L_real_decimal_interval':decimal_interval(x),
       'L_imaginary_decimal_interval':decimal_interval(y),
       'target_intervals':{k:v.exact_json() for k,v in targets.items()},
       'source_comparisons':rows,
       'modular_cross_check':{'real_bound_exponent':bound(x-xb),'imaginary_bound_exponent':bound(y-yb)},
       'source_intervals_nesting':nesting,
       'normalization':'L(f,s)=sum a_n n^(-s); the square is taken before Re and Im',
       'B_row_uses_established_Bessel_pairing':True,
       'paper_odd_value_reconstructed_from_positive_conic_root':True,
       'legacy_full_B_row_is_only_a_consistency_control':True,
       'A_C_rows_use_pairing':False,
       'D5_even_is_computed_as_actual_A_C_minor':True,
       'all_stated_tail_and_rounding_errors_enclosed':True,
       'conjectures_proved':False,
    }
    (ROOT/'outputs'/'paper_check.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PUBLICATION CHECK PASS')
    for label,row in rows.items(): print(label,row['strict_absolute_bound_exponents'])
    print('Modular comparison',result['modular_cross_check'])
    print('CONJECTURES_PROVED=False')
if __name__=='__main__': main()
