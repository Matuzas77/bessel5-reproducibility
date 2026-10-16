#!/usr/bin/env python3
"""Fresh non-destructive replay of numerical claims and exact proof witnesses.

Python >=3.11, standard library only. Default: recompute in an empty temporary
output directory, compare every JSON record to the delivered record, and discard
the temporary copy. --write records a new reference run. This is not a proof of
any exact conjectural identity.
"""
from __future__ import annotations
import argparse, hashlib, json, shutil, subprocess, sys, tempfile, time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
JOBS=[
 ('audit_algebra.py',[],None),
 ('verify_norm.py',[],'norm_exact.json'),
 ('verify_twist_symbols.py',[],'twist_symbols_exact.json'),
 ('verify_modular_bridge.py',[],'modular_bridge.json'),
 ('verify_critical_norm.py',[],'critical_norm_exact.json'),
 ('certify_source.py',['--digits','1100','--N','1000','--hyper-N','2300','--out','source_primary.json'],'source_primary.json'),
 ('certify_independent.py',['--out','source_binary.json'],'source_binary.json'),
 ('certify_lvalue.py',[],'modular_direct.json'),
 ('certify_eta_independent.py',['--out','eta_binary.json'],'eta_binary.json'),
 ('check_paper.py',[],'paper_check.json'),
]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true',help='replace recorded JSON outputs with fresh calculations')
    args=parser.parse_args()
    if sys.version_info<(3,11):raise RuntimeError('Python 3.11 or later is required')
    if not __debug__:raise RuntimeError('Run without -O: assertions are part of the verifier')
    started=time.monotonic();log=[]
    def say(text):print(text,flush=True);log.append(text)
    say('Python '+sys.version.split()[0]+'; standard-library-only replay')
    with tempfile.TemporaryDirectory(prefix='bessel5_paper_') as temp:
        work=Path(temp);shutil.copytree(ROOT/'engine',work/'engine',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        (work/'outputs').mkdir()
        for script,arguments,output in JOBS:
            command=[sys.executable,'-B',str(work/'engine'/script),*arguments]
            task=subprocess.run(command,cwd=work,text=True,capture_output=True,timeout=240)
            say(task.stdout.rstrip())
            if task.returncode:
                say(task.stderr);raise RuntimeError('Failed: '+script)
            if output:
                fresh=work/'outputs'/output
                if args.write:
                    (ROOT/'outputs').mkdir(exist_ok=True);shutil.copyfile(fresh,ROOT/'outputs'/output)
                else:
                    saved=ROOT/'outputs'/output
                    if json.loads(fresh.read_text())!=json.loads(saved.read_text()):
                        raise AssertionError('Fresh result differs from saved record: '+output)
                    say('EXACT RECORD MATCH: '+output)
    say('EXACT MODULAR WITNESSES: norm, twist, same-form bridge, Chuang cofactor and critical-norm witness PASS.')
    say('FINAL PASS: all four checked comparisons exceed 300 decimal places.')
    say('BESSEL_CONJECTURES_PROVED=False. Eq.(3.5) has an analytic proof plus these finite witnesses; code does not formalize analysis.')
    say('Elapsed seconds: '+str(round(time.monotonic()-started,2)))
    if args.write:
        (ROOT/'logs').mkdir(exist_ok=True)
        (ROOT/'logs'/'reference_run.txt').write_text('\n'.join(log)+'\n')
if __name__=='__main__':main()
