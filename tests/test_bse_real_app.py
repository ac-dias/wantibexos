#!/usr/bin/env python3
"""Application regression of the real time-reversal basis of the BSE Hamiltonian (BSE_REAL= T).

Usage: python3 tests/test_bse_real_app.py /absolute/path/to/wtb.x [--mpiexec /path/to/mpiexec [--ranks 4]] [--elpa]

The generated two-orbital DFT=S insulator of test_nonorth_optics_app.py is run on k meshes with one (odd mesh) and
four (even mesh) self-conjugate k points, and on a non-square mesh:
  * BSE_REAL= F with BSE_ALGO= cheev on one process is the reference;
  * BSE_REAL= T with cheev, cheevd, cheevr and cheevx on one process, and, with --mpiexec, ScaLAPACK (PSSYEV) on 2 and
    --ranks MPI ranks (and ELPA 1-stage with --elpa, for an executable built with -DELPA) must give the same exciton
    energies, and the same optical strengths summed over the degenerate levels (inside a degenerate level the
    eigenvectors, so the strength of a single state, are free);
  * the exciton wavefunction of the ground state (BSE_WF= T) must have the same moduli and k-space density;
  * the inputs that BSE_REAL= T cannot take (BSE_HAM_SAVE, TEMP > 0, a shifted mesh, a complex H(R)) must be refused.
No material data are needed.
"""
import argparse
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_nonorth_optics_app import fixture, rows


def clusters(table, groups=None, tol=2.e-5):
    """summed strengths of each group of degenerate levels: the groups of consecutive levels closer than tol, or
    those given (the groups of the reference, applied to the table of the same rows)"""
    if groups is None:
        groups = []
        for r in table:
            if groups and r[0]-groups[-1][1] <= tol:
                groups[-1][0] += 1
            else:
                groups.append([1, r[0]])
            groups[-1][1] = r[0]
        groups = [g[0] for g in groups]
    out, i = [], 0
    for size in groups:
        rowsum = [sum(table[i+k][1+j] for k in range(size)) for j in range(6)]
        out.append((table[i][0], size, rowsum))
        i += size
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('--mpiexec', type=Path)
    parser.add_argument('--ranks', type=int, default=4)
    parser.add_argument('--elpa', action='store_true')
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    mpiexec = args.mpiexec.resolve(strict=True) if args.mpiexec else None
    root = Path(tempfile.mkdtemp(prefix='wtb-real-app-'))
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='2')
    runs = 0

    def run(name, mesh, algo, ranks, real, temp=0, expect_ok=True, wavefunction=False, extra=None, shift=None, fix=None):
        nonlocal runs
        path = root/name
        path.mkdir(); (path/'out').mkdir()
        fixture(path, 0., .6)
        if fix:
            fix(path)
        config = {'DFT': '"S"', 'PARAMS_FILE': '"tb.dat"', 'PARAMS_FORMAT': '"text"',
                  'OUTPUT': '"./out/"', 'CALC_DATA': '"./out/"', 'NTHREADS': 2,
                  'NGX': mesh[0], 'NGY': mesh[1], 'NGZ': 1, 'NOCPF': 1, 'NBANDSC': 1, 'NBANDSV': 1,
                  'SIGMA': 0., 'CSHIFT': .03, 'FERMI_SHIFT': 0., 'IPA': 'F', 'OPT_BZ': 'F',
                  'BSE': 'T', 'BSE_ALGO': f'"{algo}"', 'BSE_REAL': 'T' if real else 'F',
                  'DTDIAG': 'T', 'CPOL': 'T', 'RNMD': 'F', 'SPEC': 'T',
                  'ENSPECI': 1., 'ENSPECF': 4., 'NESPEC': 9, 'THMAX': 1.e-8, 'NI': 1.33, 'NS': 1.6,
                  'TEMP': temp, 'TA': '"VE"', 'ST': .2, 'PHAVG': .005}
        if shift:
            config.update({'SHIFT_1': shift, 'SHIFT_2': shift})
        if wavefunction:
            config.update({'BSE_WF': 'T', 'EXC_WF_I': 1, 'EXC_WF_F': 1})
        if extra:
            config.update(extra)
        text = ''.join(f'{k}= {v}\n' for k, v in config.items())
        (path/'input.dat').write_text(text)
        command = [str(executable)] if ranks == 1 else [str(mpiexec), '-n', str(ranks), str(executable)]
        proc = subprocess.run(command, cwd=path, env=env, capture_output=True, text=True, input=text, timeout=900)
        (path/'run.log').write_text(proc.stdout+proc.stderr)
        runs += 1
        if not expect_ok:
            assert proc.returncode != 0, (name, 'should have been refused', path/'run.log')
            return proc.stdout+proc.stderr
        assert proc.returncode == 0, (name, path/'run.log')
        return rows(path/'out'/'bse_oscf.dat'), (path/'out'/'log_bse_optics.dat').read_text()

    def wavefunction(name):
        return rows(root/name/'out'/'exc_wf_1.dat'), rows(root/name/'out'/'exc_den_1.dat')

    def compare(table, reference, label, energy_tol=2.e-4, strength_tol=2.e-3):
        assert len(table) == len(reference), (label, 'row count')
        for i, (a, b) in enumerate(zip(table, reference)):
            assert abs(a[0]-b[0]) <= energy_tol, (label, 'state', i+1, 'energy', a[0], b[0])
        ca = clusters(reference)
        cb = clusters(table, [c[1] for c in ca])
        scale = max(abs(v) for c in ca for v in c[2])
        for i, (a, b) in enumerate(zip(ca, cb)):
            for j in range(6):
                assert abs(a[2][j]-b[2][j]) <= strength_tol*scale, (label, 'level', i+1, 'column', j+1, a[2][j], b[2][j])

    worst = 0.
    for mesh in ((24, 24), (9, 9), (12, 10)):
        tag = f'{mesh[0]}x{mesh[1]}'
        dim = mesh[0]*mesh[1]
        reference, _ = run(f'ref-{tag}', mesh, 'cheev', 1, False)
        assert len(reference) == dim
        for algo in ('cheev', 'cheevd', 'cheevr', 'cheevx'):
            table, log = run(f'real-{algo}-{tag}', mesh, algo, 1, True)
            compare(table, reference, f'BSE_REAL {algo} {tag}')
            worst = max(worst, max(abs(a[0]-b[0]) for a, b in zip(table, reference)))
            assert 'BSE_REAL: time-reversal real basis' in log, (algo, tag, 'log_bse_optics.dat')
            assert 'SYMMETRY OF THE SYSTEM' in log and 'Point group:' in log, (algo, tag, 'symmetry section of the log')
        if mpiexec:
            for ranks in sorted({2, args.ranks}):
                table, log = run(f'real-scalapack-r{ranks}-{tag}', mesh, 'cheev', ranks, True)
                compare(table, reference, f'BSE_REAL ScaLAPACK {ranks} ranks {tag}')
                worst = max(worst, max(abs(a[0]-b[0]) for a, b in zip(table, reference)))
                if args.elpa:
                    table, log = run(f'real-elpa-r{ranks}-{tag}', mesh, 'elpa1', ranks, True)
                    compare(table, reference, f'BSE_REAL ELPA {ranks} ranks {tag}')
                    worst = max(worst, max(abs(a[0]-b[0]) for a, b in zip(table, reference)))

    # exciton wavefunction of the ground state
    mesh = (24, 24)
    run('wf-ref', mesh, 'cheev', 1, False, wavefunction=True)
    wf_ref, den_ref = wavefunction('wf-ref')
    for name, algo, ranks in [('wf-real-cheevd', 'cheevd', 1)] + ([('wf-real-scalapack', 'cheev', args.ranks)] if mpiexec else []):
        run(name, mesh, algo, ranks, True, wavefunction=True)
        wf, den = wavefunction(name)
        assert len(wf) == len(wf_ref) and len(den) == len(den_ref), (name, 'row count')
        for i, (a, b) in enumerate(zip(wf, wf_ref)):
            assert (a[0], a[4], a[5]) == (b[0], b[4], b[5]), (name, 'row', i+1, 'integer columns')
            assert abs(math.hypot(a[6], a[7])-math.hypot(b[6], b[7])) <= 3.e-4, (name, 'row', i+1, 'modulus')
        for i, (a, b) in enumerate(zip(den, den_ref)):
            assert abs(a[3]-b[3]) <= 3.e-4, (name, 'k point', i+1, 'density')

    # refused inputs
    message = run('refused-temp', (12, 12), 'cheev', 1, True, temp=80, expect_ok=False)
    assert 'BSE_REAL= T is not implemented for TEMP > 0' in message, message
    message = run('refused-save', (12, 12), 'cheev', 1, True, expect_ok=False, extra={'BSE_HAM_SAVE': 'T'})
    assert 'does not support BSE_HAM_SAVE' in message, message
    message = run('refused-shift', (12, 12), 'cheev', 1, True, expect_ok=False, shift=0.3)
    assert 'symmetric under k -> -k' in message, message

    def complex_hopping(path):
        text = (path/'tb.dat').read_text().splitlines()
        for i, line in enumerate(text):
            parts = line.split()
            if len(parts) == 8 and parts[0] == '1' and parts[3] == '1' and parts[4] == '2':
                parts[6] = '0.05'      # imaginary part of one hopping
                text[i] = ' '.join(parts)
        (path/'tb.dat').write_text('\n'.join(text)+'\n')
    message = run('refused-complex', (12, 12), 'cheev', 1, True, expect_ok=False, fix=complex_hopping)
    assert 'needs a real H(R)' in message, message

    print(f'bse_real_app: PASS ({runs} runs; largest energy difference to the complex solver {worst:.1e} eV); artifacts: {root}')


if __name__ == '__main__':
    main()
