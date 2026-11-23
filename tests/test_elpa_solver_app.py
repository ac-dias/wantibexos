#!/usr/bin/env python3
"""Application regression of the ELPA BSE solvers (BSE_ALGO= elpa1, elpa2, BSE_NSTATES=).

Usage: python3 tests/test_elpa_solver_app.py /absolute/path/to/wtb.x --mpiexec /path/to/mpiexec [--ranks 4]

The executable must be built with -DMPI -DELPA (makefiles/makefile-*+elpa).  The generated two-orbital
DFT=S insulator of test_nonorth_optics_app.py is run on a 24 x 24 mesh (576 exciton states):
  * BSE_ALGO= cheev on one process is the reference;
  * elpa (= elpa1), elpa1 and elpa2 on 2 and --ranks MPI ranks must reproduce its exciton energies and
    optical strengths, for all states and for BSE_NSTATES= 1, 50 and more than the dimension, at TEMP= 0 and 80;
  * BSE_NSTATES= with a ScaLAPACK solver must be refused.
An ELPA library whose 2-stage solver returns wrong eigenvalues makes the elpa2 runs stop with the message of
bse_elpa_verify (or, if that check is passed, fail the comparison here).  No material data are needed.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_nonorth_optics_app import fixture, rows

MESH = 24
DIM = MESH*MESH


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--mpiexec', type=Path, required=True)
    parser.add_argument('--ranks', type=int, default=4)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    mpiexec = args.mpiexec.resolve(strict=True)
    root = Path(tempfile.mkdtemp(prefix='wtb-elpa-app-'))
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='2')
    runs = 0

    def run(name, algo, ranks, nstates=None, temp=0, expect_ok=True):
        nonlocal runs
        path = root/name
        path.mkdir(); (path/'out').mkdir()
        fixture(path, 0., .6)
        config = {'DFT': '"S"', 'PARAMS_FILE': '"tb.dat"', 'PARAMS_FORMAT': '"text"',
                  'OUTPUT': '"./out/"', 'CALC_DATA': '"./out/"', 'NTHREADS': 2,
                  'NGX': MESH, 'NGY': MESH, 'NGZ': 1, 'NOCPF': 1, 'NBANDSC': 1, 'NBANDSV': 1,
                  'SIGMA': 0., 'CSHIFT': .03, 'FERMI_SHIFT': 0., 'IPA': 'F', 'OPT_BZ': 'F',
                  'BSE': 'T', 'BSE_ALGO': f'"{algo}"', 'DTDIAG': 'T', 'CPOL': 'T', 'RNMD': 'F', 'SPEC': 'T',
                  'ENSPECI': 1., 'ENSPECF': 4., 'NESPEC': 9, 'THMAX': 1.e-8, 'NI': 1.33, 'NS': 1.6,
                  'TEMP': temp, 'TA': '"VE"', 'ST': .2, 'PHAVG': .005}
        if nstates is not None:
            config['BSE_NSTATES'] = nstates
        text = ''.join(f'{k}= {v}\n' for k, v in config.items())
        (path/'input.dat').write_text(text)
        command = [str(executable)] if ranks == 1 else [str(mpiexec), '-n', str(ranks), str(executable)]
        proc = subprocess.run(command, cwd=path, env=env, capture_output=True, text=True, input=text, timeout=600)
        (path/'run.log').write_text(proc.stdout+proc.stderr)
        runs += 1
        if not expect_ok:
            assert proc.returncode != 0, (name, 'should have been refused', path/'run.log')
            return proc.stdout+proc.stderr
        assert proc.returncode == 0, (name, path/'run.log')
        return rows(path/'out'/'bse_oscf.dat'), (path/'out'/'log_bse_optics.dat').read_text()

    def close(table, ref, label, energy_tol=2.e-4, strength_tol=1.e-3):
        scale = max(abs(r[j]) for r in ref for j in range(1, 7))
        for i, (a, b) in enumerate(zip(table, ref)):
            assert abs(a[0]-b[0]) <= energy_tol, (label, 'state', i+1, 'energy', a[0], b[0])
            for j in range(1, 7):
                assert abs(a[j]-b[j]) <= strength_tol*scale, (label, 'state', i+1, 'column', j, a[j], b[j])

    worst = 0.
    for temp in (0, 80):
        reference, _ = run(f'ref-T{temp}', 'cheev', 1, temp=temp)
        assert len(reference) == DIM
        for ranks in sorted({2, args.ranks}):
            for algo in ('elpa', 'elpa1', 'elpa2'):
                if temp and algo == 'elpa':
                    continue
                table, log = run(f'{algo}-r{ranks}-T{temp}', algo, ranks, temp=temp)
                assert len(table) == DIM, (algo, ranks, 'row count')
                close(table, reference, f'{algo} {ranks} ranks T={temp}')
                worst = max(worst, max(abs(a[0]-b[0]) for a, b in zip(table, reference)))
                for nstates in (1, 50, DIM+100):
                    if algo == 'elpa':
                        continue
                    table, log = run(f'{algo}-r{ranks}-T{temp}-n{nstates}', algo, ranks, nstates=nstates, temp=temp)
                    expected = min(nstates, DIM)
                    assert len(table) == expected, (algo, ranks, nstates, 'row count', len(table))
                    close(table, reference[:expected], f'{algo} {ranks} ranks T={temp} BSE_NSTATES={nstates}')
                    assert ('lowest of the' in log) == (nstates < DIM), (algo, ranks, nstates, 'log_bse_optics.dat')
    message = run('refused-cheev-nstates', 'cheev', args.ranks, nstates=10, expect_ok=False)
    assert 'BSE_NSTATES= needs BSE_ALGO= elpa1 or elpa2' in message, message
    print(f'elpa_solver_app: PASS ({runs} runs; largest energy difference to cheev {worst:.1e} eV); artifacts: {root}')


if __name__ == '__main__':
    main()
