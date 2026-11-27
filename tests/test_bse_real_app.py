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
  * cheevd (SSYEVD) on a 54x54 mesh, N = 2916, must agree with cheev: its workspace 1+6N+2N^2 is above 2^24 numbers, where a
    single precision workspace query rounds an odd size down (MKL, N >= 2895);
  * the log names the routine that was called (SSYEV, SSYEVD, SSYEVR, SSYEVX for BSE_REAL= T);
  * the inputs that BSE_REAL= T cannot take (BSE_HAM_SAVE, TEMP > 0, a shifted mesh, a complex H(R)) must be refused.
The same comparisons are made for a generated spinor model with spin-orbit coupling (T^2 = -1, the Kramers real basis), in
the spin up / spin down layout of DFT=S and in the interleaved layout of DFT=W, on meshes with 4, 1 and 0 self-conjugate k
points; a time-reversal-breaking term, and an odd number of bands where a Kramers doublet would be cut, must be refused.
No material data are needed.
"""
import argparse
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_nonorth_optics_app import fixture, rows


def fixture_soc(path, fermi, scissor, zeeman=0., dft='S'):
    """Spinor version of the two-orbital fixture: 2 orbitals x 2 spins, systype SOC.

    DFT=S: tb.dat, the overlap S0 x 1, spin up orbitals first, then spin down (as siesta2wtb.py writes them).
    DFT=W: hr.dat, orthogonal basis, orbitals interleaved (2m-1 up, 2m down, as Wannier90 writes spinors).
    H(R) = A(R) x 1 + i sum_a b_a(R) x sigma_a with real A, b_a and b_a(-R) = -b_a(R)^T, which is exactly the form that
    commutes with T = i sigma_y K (T^2 = -1).  A complex H(R) in the usual sense: the time-reversal representation is the
    spin flip, not the identity.  zeeman adds a real sigma_z term that breaks time reversal."""
    a_r = [[[-1., .35], [.35, 1.6]], [[0., .4], [.1, 0.]], [[0., .1], [.4, 0.]]]
    s_r = [[[1., .08], [.08, 1.]], [[.02, .07], [.01, .015]], [[.02, .01], [.07, .015]]]
    # b_a(R = 0) antisymmetric; b_a(+1) = B_a, b_a(-1) = -B_a^T
    b0 = {'x': [[0., .12], [-.12, 0.]], 'y': [[0., .08], [-.08, 0.]], 'z': [[0., .15], [-.15, 0.]]}
    b1 = {'x': [[.05, .03], [-.02, .04]], 'y': [[-.04, .02], [.06, .01]], 'z': [[.03, -.05], [.02, -.06]]}
    bm1 = {c: [[-b1[c][j][i] for j in range(2)] for i in range(2)] for c in b1}
    b = [b0, b1, bm1]
    tau = [[0., 0., 0.], [.32, -.12, .08]]
    norb = 4

    def split(index):
        """spin (0 up, 1 down) and spatial orbital of the file index"""
        return divmod(index, 2) if dft == 'S' else divmod(index, 2)[::-1]

    if dft == 'S':
        lines = ['SOC', str(scissor), str(fermi), '1 0 0', '0 4 0', '0 0 4', str(norb), '3', '# x y z a b h ih s # #']
    else:
        lines = ['SOC', str(scissor), str(fermi), '1 0 0', '0 4 0', '0 0 4', '# hr', str(norb), '3', '1 1 1']
    for ir, r in enumerate([0, 1, -1]):
        for col in range(norb):
            for row in range(norb):
                (sr_, o), (sc_, p) = split(row), split(col)
                re = a_r[ir][o][p] if sr_ == sc_ else 0.
                im = 0.
                bx, by, bz = b[ir]['x'][o][p], b[ir]['y'][o][p], b[ir]['z'][o][p]
                # i (bx sigma_x + by sigma_y + bz sigma_z) in the up/down blocks
                if (sr_, sc_) == (0, 0):
                    im += bz
                    re += zeeman*(1. if o == p and ir == 0 else 0.)
                elif (sr_, sc_) == (1, 1):
                    im -= bz
                    re -= zeeman*(1. if o == p and ir == 0 else 0.)
                elif (sr_, sc_) == (0, 1):
                    im += bx
                    re += by
                else:
                    im += bx
                    re -= by
                if dft == 'S':
                    ov = s_r[ir][o][p] if sr_ == sc_ else 0.
                    lines.append(f'{r} 0 0 {row+1} {col+1} {re} {im} {ov}')
                else:
                    lines.append(f'{r} 0 0 {row+1} {col+1} {re} {im}')
    if dft == 'S':
        (path/'tb.dat').write_text('\n'.join(lines)+'\n')
        basis = ['# index species x y z']
        for row in range(norb):
            t = tau[split(row)[1]]
            basis.append(f'{row+1} X {t[0]} {t[1]} {t[2]}')
        (path/'basis_set-soc').write_text('\n'.join(basis)+'\n')
    else:
        (path/'hr.dat').write_text('\n'.join(lines)+'\n')


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

    def run(name, mesh, algo, ranks, real, temp=0, expect_ok=True, wavefunction=False, extra=None, shift=None, fix=None,
            soc=None, zeeman=0., nbands=None):
        nonlocal runs
        path = root/name
        path.mkdir(); (path/'out').mkdir()
        if soc:
            fixture_soc(path, 0., .6, zeeman=zeeman, dft=soc)
        else:
            fixture(path, 0., .6)
        if fix:
            fix(path)
        nv, nc = nbands or ((2, 2) if soc else (1, 1))
        config = {'DFT': f'"{soc or "S"}"', 'PARAMS_FILE': '"hr.dat"' if soc == 'W' else '"tb.dat"', 'PARAMS_FORMAT': '"text"',
                  'OUTPUT': '"./out/"', 'CALC_DATA': '"./out/"', 'NTHREADS': 2,
                  'NGX': mesh[0], 'NGY': mesh[1], 'NGZ': 1, 'NOCPF': 2 if soc else 1, 'NBANDSC': nc, 'NBANDSV': nv,
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

    def compare(table, reference, label, energy_tol=2.e-4, strength_tol=2.e-3, cluster_tol=2.e-5):
        assert len(table) == len(reference), (label, 'row count')
        for i, (a, b) in enumerate(zip(table, reference)):
            assert abs(a[0]-b[0]) <= energy_tol, (label, 'state', i+1, 'energy', a[0], b[0])
        ca = clusters(reference, tol=cluster_tol)
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
            # the keyword is the name of the complex driver: the real path calls the real symmetric one
            routine = {'cheev': 'SSYEV', 'cheevd': 'SSYEVD', 'cheevr': 'SSYEVR', 'cheevx': 'SSYEVX'}[algo]
            assert re.search(rf'Eigensolver: {routine} \(real symmetric', log), (algo, tag, 'solver named in the log')
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

    # a matrix above N = 2895: SSYEVD needs 1+6N+2N^2 > 2^24 numbers, an odd count that a single precision workspace query
    # rounds down by one (MKL), and the routine then refuses LWORK.  54x54 mesh, N = 2916, against SSYEV
    mesh = (54, 54)
    assert 1+6*mesh[0]*mesh[1]+2*(mesh[0]*mesh[1])**2 > 2**24
    reference, _ = run('real-cheev-54x54', mesh, 'cheev', 1, True)
    table, log = run('real-cheevd-54x54', mesh, 'cheevd', 1, True)
    compare(table, reference, 'BSE_REAL cheevd 54x54 (workspace above 2^24 numbers)')
    worst = max(worst, max(abs(a[0]-b[0]) for a, b in zip(table, reference)))

    # spin-orbit coupling: time reversal T = i sigma_y K, T^2 = -1 (Kramers real basis).  The levels come in pairs a few
    # 1e-4 eV apart, so the strengths are summed over the levels closer than 1e-3 eV
    soc_tol = 1.e-3
    soc_note = 'spin-orbit coupling, time reversal T = i sigma_y K (T^2 = -1)'
    for layout in ('S', 'W'):
        for mesh in ((12, 12), (9, 9), (12, 10)) if layout == 'S' else ((12, 12), (9, 9)):
            tag = f'{layout}-{mesh[0]}x{mesh[1]}'
            reference, _ = run(f'soc-ref-{tag}', mesh, 'cheev', 1, False, soc=layout)
            assert len(reference) == 4*mesh[0]*mesh[1]
            for algo in ('cheev', 'cheevd', 'cheevr', 'cheevx') if layout == 'S' else ('cheevd',):
                table, log = run(f'soc-real-{algo}-{tag}', mesh, algo, 1, True, soc=layout)
                compare(table, reference, f'BSE_REAL SOC {algo} {tag}', cluster_tol=soc_tol)
                worst = max(worst, max(abs(a[0]-b[0]) for a, b in zip(table, reference)))
                assert soc_note in log, (algo, tag, 'log_bse_optics.dat')
            if mpiexec and (layout == 'S' or mesh == (12, 12)):
                for ranks in sorted({2, args.ranks}) if layout == 'S' and mesh != (12, 10) else [args.ranks]:
                    table, log = run(f'soc-real-scalapack-r{ranks}-{tag}', mesh, 'cheev', ranks, True, soc=layout)
                    compare(table, reference, f'BSE_REAL SOC ScaLAPACK {ranks} ranks {tag}', cluster_tol=soc_tol)
                    worst = max(worst, max(abs(a[0]-b[0]) for a, b in zip(table, reference)))
                    if args.elpa and layout == 'S':
                        table, log = run(f'soc-real-elpa-r{ranks}-{tag}', mesh, 'elpa1', ranks, True, soc=layout)
                        compare(table, reference, f'BSE_REAL SOC ELPA {ranks} ranks {tag}', cluster_tol=soc_tol)
                        worst = max(worst, max(abs(a[0]-b[0]) for a, b in zip(table, reference)))
    # half-step mesh: no self-conjugate k point, so no Kramers pairs at all, and an odd number of bands is allowed
    for nbands in ((2, 2), (1, 1)):
        reference, _ = run(f'soc-ref-half-{nbands[0]}', (12, 12), 'cheev', 1, False, shift=1/24, soc='S', nbands=nbands)
        table, log = run(f'soc-real-half-{nbands[0]}', (12, 12), 'cheev', 1, True, shift=1/24, soc='S', nbands=nbands)
        compare(table, reference, f'BSE_REAL SOC half-step mesh, {nbands[0]} bands', cluster_tol=soc_tol)
        assert ' 72 pairs of k points (k, -k) and 0 self-conjugate k points' in log, log[:2000]
    # exciton wavefunction: the ground state is a pair, so the densities are summed over the two states
    def pair_density(name):
        wf, den = {}, {}
        for n in (1, 2):
            for r in rows(root/name/'out'/f'exc_wf_{n}.dat'):
                key = (r[0], r[4], r[5])
                wf[key] = wf.get(key, 0.)+r[6]**2+r[7]**2
            for i, r in enumerate(rows(root/name/'out'/f'exc_den_{n}.dat')):
                den[i] = den.get(i, 0.)+r[3]
        return wf, den
    run('soc-wf-ref', (12, 12), 'cheev', 1, False, wavefunction=True, extra={'EXC_WF_F': 2}, soc='S')
    wf_ref, den_ref = pair_density('soc-wf-ref')
    for name, algo, ranks in [('soc-wf-real-cheevd', 'cheevd', 1)] + ([('soc-wf-real-scalapack', 'cheev', args.ranks)] if mpiexec else []):
        run(name, (12, 12), algo, ranks, True, wavefunction=True, extra={'EXC_WF_F': 2}, soc='S')
        wf, den = pair_density(name)
        assert wf.keys() == wf_ref.keys() and den.keys() == den_ref.keys(), (name, 'rows')
        for key in wf:
            assert abs(wf[key]-wf_ref[key]) <= 3.e-4, (name, key, 'summed |psi|^2', wf[key], wf_ref[key])
        for key in den:
            assert abs(den[key]-den_ref[key]) <= 3.e-4, (name, 'k point', key+1, 'density')
    # refusals: a term that breaks time reversal (a real sigma_z coupling), and Kramers doublets cut by the window
    message = run('soc-refused-zeeman', (12, 12), 'cheev', 1, True, expect_ok=False, soc='S', zeeman=.3)
    assert 'needs a time-reversal invariant tight-binding model' in message, message
    message = run('soc-refused-odd-bands', (12, 12), 'cheev', 1, True, expect_ok=False, soc='S', nbands=(1, 1))
    assert 'needs NBANDSV and NBANDSC even' in message, message

    # refused inputs
    message = run('refused-temp', (12, 12), 'cheev', 1, True, temp=80, expect_ok=False)
    assert 'BSE_REAL= T is not implemented for TEMP > 0' in message, message
    message = run('refused-save', (12, 12), 'cheev', 1, True, expect_ok=False, extra={'BSE_HAM_SAVE': 'T'})
    assert 'does not support BSE_HAM_SAVE' in message, message
    message = run('refused-shift', (12, 12), 'cheev', 1, True, expect_ok=False, shift=0.3)
    assert 'symmetric under k -> -k' in message, message
    # the mesh is tested first, before the tight-binding file is read: a missing file does not matter
    message = run('refused-shift-early', (12, 12), 'cheev', 1, True, expect_ok=False, shift=1/36,
                  extra={'PARAMS_FILE': '"missing-file.dat"'})
    assert 'symmetric under k -> -k' in message and 'SHIFT_1= 0.027778 with NGX= 12' in message, message
    assert 'multiple of 1/(2 NGX) = 0.041667' in message and 'missing-file' not in message, message
    # a half-step shift keeps the mesh symmetric (no self-conjugate k point)
    reference, _ = run('ref-halfstep', (12, 12), 'cheev', 1, False, shift=1/24)
    table, log = run('real-halfstep', (12, 12), 'cheev', 1, True, shift=1/24)
    compare(table, reference, 'BSE_REAL half-step shifted mesh')
    assert ' 72 pairs of k points (k, -k) and 0 self-conjugate k points' in log, log[:2000]

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
