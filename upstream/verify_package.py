#!/usr/bin/env python3
"""Audit release integrity and compare the two exact critical-norm certificates.
This is a lightweight check of stored data, not an analytic proof or a full replay.
For fresh computations run python -B anc/verify.py and anc/verify_primary.py.
"""
from pathlib import Path
import hashlib, json, re, sys

def main():
    root=Path(__file__).resolve().parent
    if not __debug__:
        raise RuntimeError('Run without -O; assertions are part of this audit.')
    manifest=root/'SHA256SUMS'
    if not manifest.is_file():
        raise FileNotFoundError('Missing SHA256SUMS')
    count=0
    for line in manifest.read_text().splitlines():
        if not line.strip(): continue
        digest, relative=line.split('  ',1)
        target=(root/relative).resolve()
        if root not in target.parents:
            raise ValueError('Unsafe manifest path: '+relative)
        actual=hashlib.sha256(target.read_bytes()).hexdigest()
        if actual!=digest: raise AssertionError('Hash mismatch: '+relative)
        count+=1
    print('PASS SHA-256 manifest:',count,'files')
    out=root/'anc/outputs'
    primary=json.loads((out/'critical_norm_primary.json').read_text())
    independent=json.loads((out/'critical_norm_exact.json').read_text())
    assert primary['kernel_coefficients']==independent['kernel_coefficients'][:201]
    for key in ['kernel_eta_coordinates','f_packet_P_coordinates','hecke_witness_coordinates','h','L3_norm_factor','L2_norm_factor']:
        assert primary[key]==independent[key],key
    assert primary['L3_norm_factor']=='1/8'
    assert primary['L2_norm_factor']=='sqrt(15)*pi^2'
    print('PASS 201 common class-number coefficients and six exact fields')
    audit=json.loads((out/'critical_norm_audit.json').read_text())
    assert audit['status']=='PASS' and audit['Delta_zero_character_sum_checks']==117
    assert len(audit['higher_prime_power_cyclotomic_checks'])==16
    assert sum(x['ratio_to_sqrt60']!=0 for x in audit['higher_prime_power_cyclotomic_checks'])==8
    print('PASS stored extra audit: 16 cyclotomic sums, including 8 nonzero, and 117 boundary cases')
    tex=(root/'main.tex').read_text()
    equations={key:body for key,body in re.findall(r'\\begin\{equation\}\\label\{([^}]+)\}(.*?)\\end\{equation\}',tex,re.S)}
    assert len(equations)==9
    assert r'\stackrel?' not in equations['eq:factorization']
    for key in ['eq:main','eq:Dnorm','eq:LTY']:
        assert r'\stackrel?' in equations[key],key
    paper=json.loads((out/'paper_check.json').read_text())
    assert paper['status']=='PASS' and paper['conjectures_proved'] is False
    print('PASS article claim boundaries: equation (3.5) exact; Bessel equalities conjectural')
    print('FINAL PASS. Stored checks do not formalize the analytic continuation proof.')
if __name__=='__main__': main()
