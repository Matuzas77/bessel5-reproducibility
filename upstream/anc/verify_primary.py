#!/usr/bin/env python3
"""Optional second exact reconstruction using SymPy (no Bessel quadrature).

Copy the primary code into a temporary tree, rebuild the norm projection and
Manin-symbol matrices, and verify both raw witnesses. No saved periods or
matrices initialize these computations. Requires SymPy and mpmath.
"""
from pathlib import Path
import subprocess,sys,tempfile,shutil,time,json
ROOT=Path(__file__).resolve().parent
JOBS=[
 ('norm','input/exact_hecke_definition.py'),
 ('norm','verify_fricke_eigenline.py'),
 ('norm','input/verify_fricke_exact.py'),
 ('norm','rankin_projection.py'),
 ('norm','verify_norm_independent.py'),
 ('norm','verify_phase_and_euler.py'),
 ('twist','derive_symbols.py'),
 ('twist','build_witness.py'),
 ('twist','derive_projection.py'),
 ('twist','verify_symbols_independent.py'),
 ('twist','verify_modular_coeffs.py'),
 ('twist','verify_coeffs_independent.py'),
 ('twist','verify_period_algebra.py'),
 ('critical_norm','proofs/verify_exact.py'),
 ('critical_norm','proofs/verify_audit.py'),
]
def main():
    if not __debug__:raise RuntimeError('Run without -O: assertions are part of the verifier')
    t=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='bessel5_exact_') as d:
        work=Path(d)/'primary';shutil.copytree(ROOT/'primary',work)
        (work/'twist/results').mkdir(exist_ok=True)
        (work/'critical_norm/outputs').mkdir(exist_ok=True)
        for group,name in JOBS:
            cwd=work/group
            p=subprocess.run([sys.executable,'-B',str(cwd/name)],cwd=cwd,
                             capture_output=True,text=True,timeout=240)
            print('RUN',group,name,flush=True)
            if p.returncode:
                print(p.stdout,p.stderr);raise RuntimeError('Failed: '+name)
            print('PASS',group,name,flush=True)
        for group,name in [('norm','rankin_projection.json'),('norm','verify_norm_independent.json'),
                           ('twist','results/independent_symbols.json')]:
            print(group,name,json.loads((work/group/name).read_text()).get('status'))
        for name, stored in [('exact.json','critical_norm_primary.json'),('audit.json','critical_norm_audit.json')]:
            assert json.loads((work/'critical_norm/outputs'/name).read_text())==json.loads((ROOT/'outputs'/stored).read_text()), name
            print('EXACT RECORD MATCH:',stored)
    print('FINAL PASS: 15 exact reconstruction jobs; no Bessel equality asserted. Analytic proofs are supplied separately.')
    print('Elapsed seconds:',round(time.monotonic()-t,2))
if __name__=='__main__':main()
