#!/usr/bin/env python3
"""BSE application checks for V2DLAVG; accepts optional baseline and MPI launcher.

python3 tests/test_v2dlavg_app.py /path/to/wtb.x --baseline /path/to/old-wtb.x
Add --mpiexec /path/to/mpiexec to check two-rank assembly as well.
"""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
from test_nonorth_optics_app import fixture, rows, close, conductivity, sheet


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable',type=Path)
    parser.add_argument('--baseline',type=Path)
    parser.add_argument('--mpiexec',type=Path)
    args=parser.parse_args()
    exe=args.executable.resolve(strict=True)
    root=Path(tempfile.mkdtemp(prefix='wtb-v2dlavg-app-'))
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='2')
    results={}

    def run(name,potential,temperature=0,height=20,kernel='fast',binary=exe,ranks=1,ediel=5.1622):
        path=root/name;path.mkdir();(path/'out').mkdir()
        fixture(path,0.,0.)
        text=(path/'tb.dat').read_text().splitlines();text[5]=f'0 0 {height}'
        (path/'tb.dat').write_text('\n'.join(text)+'\n')
        config={'DFT':'"S"','SYSDIM':'"2D"','PARAMS_FILE':'"tb.dat"','PARAMS_FORMAT':'"text"',
            'BSE_CENTER_FILE':'"position.dat"','OUTPUT':'"./out/"','CALC_DATA':'"./out/"',
            'NTHREADS':2,'NGX':8,'NGY':8,'NGZ':1,'NOCPF':1,'NBANDSV':1,'NBANDSC':1,
            'SHIFT_1':.017,'SHIFT_2':.019,'SIGMA':.02,'CSHIFT':.08,'RNMD':'F',
            'IPA':'T','BSE':'T','SPEC':'T','CPOL':'T','DTDIAG':'T','BSE_ALGO':'cheev',
            'COULOMB_POT':potential,'EDIEL':ediel,'ENSPECI':.1,'ENSPECF':8.,'NESPEC':101,
            'TEMP':temperature,'TA':'"VE"','ST':.2,'PHAVG':.005,
            'THMAX':height*1.e-10,'NI':1.33,'NS':1.6}
        inp=''.join(f'{k}= {v}\n' for k,v in config.items());(path/'input.dat').write_text(inp)
        command=[str(binary)] if ranks==1 else [str(args.mpiexec.resolve()),'-n',str(ranks),str(binary)]
        proc=subprocess.run(command,cwd=path,env=dict(env,WTB_KERNEL=kernel),input=inp,
            capture_output=True,text=True,timeout=120)
        (path/'run.log').write_text(proc.stdout+proc.stderr)
        assert proc.returncode==0,path/'run.log'
        for p in ['ipa','bse']:
            conductivity(path/'out',p);sheet(path/'out',p,1.33,1.6,height*1.e-10)
        results[name]={'temperature_K':temperature,'potential':potential,'height_A':height,
            'ranks':ranks,'kernel':kernel,'exit_code':proc.returncode}
        return path

    old=run('v2dl','V2DL')
    new=run('v2dlavg','V2DLAVG')
    # Single transition per k: the new q=0 contribution adds v0 to every
    # diagonal of this zero-temperature Hamiltonian, without changing off-diagonals.
    hx=math.pi/8;hy=math.pi/(4*8)
    average=(hx*math.asinh(hy/hx)+hy*math.asinh(hx/hy))/(hx*hy)
    v0=-90.4756*average/(4*64*5.1622)
    eold=rows(old/'out/bse_oscf.dat');enew=rows(new/'out/bse_oscf.dat')
    worst=max(abs(b[0]-a[0]-v0) for a,b in zip(eold,enew))
    assert len(eold)==len(enew)==64 and worst<1.5e-4,(v0,worst)
    results['q0_shift']={'predicted_eV':v0,'maximum_residual_eV':worst}
    close(rows(old/'out/ipa_oscf.dat'),rows(new/'out/ipa_oscf.dat'),'IPA independent of Coulomb choice')

    for kernel in ['old','check']:
        p=run('v2dlavg-'+kernel,'V2DLAVG',kernel=kernel)
        close([[r[0]] for r in rows(p/'out/bse_oscf.dat')],[[r[0]] for r in enew],kernel+' spectrum')
        close(rows(p/'out/bse_diel_xx.dat'),rows(new/'out/bse_diel_xx.dat'),kernel+' dielectric',rtol=2.e-4)
    tall=run('v2dlavg-height40','V2DLAVG',height=40)
    close(rows(tall/'out/bse_oscf.dat'),enew,'BSE vacuum independence',rtol=2.e-4)
    close(rows(tall/'out/bse_absorptance.dat'),rows(new/'out/bse_absorptance.dat'),
        'sheet absorptance vacuum independence',rtol=2.e-4)
    run('v2dlavg-temperature','V2DLAVG',temperature=300)
    if args.baseline:
        baseline=args.baseline.resolve(strict=True)
        for potential in ['V2DL','V2DTAVG']:
            before=run('before-'+potential,potential,binary=baseline)
            after=old if potential=='V2DL' else run('after-'+potential,potential)
            for f in ['bse_oscf.dat','bse_diel_xx.dat','ipa_opt_cond_real.dat','bse_absorptance.dat']:
                close(rows(before/'out'/f),rows(after/'out'/f),'legacy '+potential+' '+f,rtol=2.e-4)
    if args.mpiexec:
        p=run('v2dlavg-mpi2','V2DLAVG',ranks=2)
        close([[r[0]] for r in rows(p/'out/bse_oscf.dat')],[[r[0]] for r in enew],
            'MPI eigenvalues',atol=1.5e-4)
        close(rows(p/'out/bse_diel_xx.dat'),rows(new/'out/bse_diel_xx.dat'),
            'MPI dielectric',atol=2.e-5,rtol=3.e-4)
    (root/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    print(f'v2dlavg_app: PASS; artifacts: {root}')


if __name__=='__main__':main()
