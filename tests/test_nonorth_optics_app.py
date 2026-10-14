#!/usr/bin/env python3
"""Application regression on a generated two-orbital DFT=S insulator.

Usage: python3 tests/test_nonorth_optics_app.py /absolute/path/to/wtb.x
No material data or third-party Python packages are needed. Runs use a fresh
system temporary directory, printed on completion and retained for inspection.
"""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


def fixture(path, fermi, scissor):
    h = [[[-1., .35], [.35, 1.6]], [[0., .4], [.1, 0.]], [[0., .1], [.4, 0.]]]
    s = [[[1., .08], [.08, 1.]], [[.02, .07], [.01, .015]], [[.02, .01], [.07, .015]]]
    tau = [[0., 0., 0.], [.32, -.12, .08]]
    lines = ['NP', str(scissor), str(fermi), '1 0 0', '0 4 0', '0 0 4', '2', '3', '# x y z a b h ih s # #']
    positions = ['synthetic position matrix', '2', '3']
    for ir, r in enumerate([0, 1, -1]):
        for b in range(2):
            for a in range(2):
                lines.append(f'{r} 0 0 {a+1} {b+1} {h[ir][a][b]} 0 {s[ir][a][b]}')
                xyz = []
                for axis in range(3):
                    dip = (0.12*(axis+1) if a != b else (0.03 if a == 0 else -0.02)) if ir == 0 else 0.
                    val = dip + 0.5*(tau[a][axis]+tau[b][axis]+(r if axis == 0 else 0.))*s[ir][a][b]
                    xyz.extend([str(val), '0'])
                positions.append(f'{r} 0 0 {a+1} {b+1} ' + ' '.join(xyz))
    (path/'tb.dat').write_text('\n'.join(lines)+'\n')
    (path/'position.dat').write_text('\n'.join(positions)+'\n')
    (path/'basis_set-NP').write_text('# index species x y z\n1 X 0 0 0\n2 X .32 -.12 .08\n')


def rows(path):
    result = [[float(v) for v in line.split()] for line in path.read_text().splitlines()
              if line.strip() and not line.lstrip().startswith('#')]
    assert result and all(math.isfinite(v) for row in result for v in row), path
    return result


def close(a, b, label, atol=3.e-6, rtol=3.e-5):
    assert len(a) == len(b), label
    error = 0.
    for x, y in zip(a, b):
        assert len(x) == len(y), label
        for u, v in zip(x, y):
            assert math.isclose(u, v, abs_tol=atol, rel_tol=rtol), (label, u, v)
            error = max(error, abs(u-v))
    return error


def conductivity(path, prefix):
    hbar = 6.582119569e-16
    for suffix, diel_stem, channels in [
            ('', 'diel_', [('xx', 1.), ('yy', 1.), ('zz', 1.), ('xy', 0.), ('xz', 0.), ('yz', 0.)]),
            ('-pol', 'diel-pol_', [('x', 1.), ('y', 1.), ('z', 1.), ('sp', 1.), ('sm', 1.)])]:
        real = rows(path/f'{prefix}_opt_cond_real{suffix}.dat')
        imag = rows(path/f'{prefix}_opt_cond_imag{suffix}.dat')
        assert len(real) == len(imag), (prefix, suffix, 'conductivity row count')
        assert all(len(r) == len(channels)+1 for r in real+imag), (prefix, suffix, 'conductivity columns')
        for j, (component, delta) in enumerate(channels, start=1):
            label = f'{prefix} conductivity {component}'
            diel = rows(path/f'{prefix}_{diel_stem}{component}.dat')
            assert len(diel) == len(real), (label, 'dielectric row count')
            assert all(len(d) == 3 for d in diel), (label, 'dielectric columns')
            for d, re, im in zip(diel, real, imag):
                close([[d[0], d[0]]], [[re[0], im[0]]], label+' energy alignment')
                factor = d[0]/(4*math.pi*hbar)
                close([[re[j]/factor, im[j]/factor]],
                      [[d[2], -(d[1]-delta)]], label)


def sheet(path, prefix, ni, ns, thmax):
    # Normal incidence on a sheet of conductivity sigma between indices ni and ns.
    clight = 3.e8
    for suffix, channels in [('', ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']), ('-pol', ['x', 'y', 'z', 'sp', 'sm'])]:
        real = rows(path/f'{prefix}_opt_cond_real{suffix}.dat')
        imag = rows(path/f'{prefix}_opt_cond_imag{suffix}.dat')
        tr, rf, ab = (rows(path/f'{prefix}_{name}{suffix}.dat') for name in ['transmittance', 'reflectance', 'absorptance'])
        assert len(tr) == len(rf) == len(ab) == len(real), (prefix, suffix, 'sheet row count')
        assert all(len(r) == len(channels)+1 for r in tr+rf+ab), (prefix, suffix, 'sheet columns')
        for j, component in enumerate(channels, start=1):
            label = f'{prefix} sheet {component}'
            for re, im, t, r, a in zip(real, imag, tr, rf, ab):
                close([[t[0], r[0], a[0]]], [[re[0]]*3], label+' energy alignment')
                y = 4*math.pi*thmax*complex(re[j], im[j])/clight
                d = abs(ni+ns+y)**2
                close([[t[j], r[j], a[j]]], [[4*ni*ns/d, abs(ni-ns-y)**2/d, 4*ni*y.real/d]], label)
                # Diagonal and circular channels are passive; cross terms are not observables.
                if component not in ('xy', 'xz', 'yz'):
                    assert a[j] > -3.e-6, (label, a[j])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    root = Path(tempfile.mkdtemp(prefix='wtb-nonorth-app-'))
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='2')
    runs = 0
    results = {}

    def run(mode, dipoles, eta, variant, fermi=0., shift=0., scissor=.6, cshift=.03):
        nonlocal runs
        path = root/f'{mode}-d{int(dipoles)}-s{eta}-{variant}'
        path.mkdir(); (path/'out').mkdir()
        fixture(path, fermi, scissor)
        config = {'DFT': '"S"', 'PARAMS_FILE': '"tb.dat"', 'PARAMS_FORMAT': '"text"',
                  'OUTPUT': '"./out/"', 'CALC_DATA': '"./out/"', 'NTHREADS': 2,
                  'NGX': 3, 'NGY': 1, 'NGZ': 1, 'NOCPF': 1, 'NBANDSC': 1, 'NBANDSV': 1,
                  'SIGMA': eta, 'CSHIFT': cshift, 'FERMI_SHIFT': shift,
                  'IPA': 'T' if mode == 'ipa' else 'F', 'OPT_BZ': 'T' if mode == 'bz' else 'F',
                  'BSE': 'T' if mode.startswith('bse') else 'F', 'BSE_ALGO': 'cheev',
                  'DTDIAG': 'T', 'CPOL': 'T', 'RNMD': 'F', 'SPEC': 'F' if mode == 'bz' else 'T',
                  'ENSPECI': 1., 'ENSPECF': 4., 'NESPEC': 9, 'THMAX': 1.e-8, 'NI': 1.33, 'NS': 1.6,
                  'TEMP': 80 if mode == 'bse-temp' else 0, 'TA': '"VE"', 'ST': .2, 'PHAVG': .005}
        if dipoles:
            config['BSE_CENTER_FILE'] = '"position.dat"'
        (path/'input.dat').write_text(''.join(f'{k}= {v}\n' for k, v in config.items()))
        proc = subprocess.run([str(executable)], cwd=path, env=env, capture_output=True, text=True, input=(path/'input.dat').read_text(), timeout=60)
        (path/'run.log').write_text(proc.stdout+proc.stderr)
        assert proc.returncode == 0, path/'run.log'
        out = path/'out'
        names = ['bz_act_'+c+'.dat' for c in ['x', 'y', 'z', 'sp', 'sm']] if mode == 'bz' else ['ipa_oscf.dat', 'ipa_oscf-pol.dat']
        if mode.startswith('bse'):
            names += ['bse_oscf.dat', 'bse_oscf-pol.dat']
        data = {name: rows(out/name) for name in names}
        if mode != 'bz':
            conductivity(out, 'ipa')
            sheet(out, 'ipa', 1.33, 1.6, 1.e-8)
            if mode.startswith('bse'):
                conductivity(out, 'bse')
                sheet(out, 'bse', 1.33, 1.6, 1.e-8)
        runs += 1
        return data

    for mode in ['ipa', 'bz', 'bse', 'bse-temp']:
        worst = 0.
        for dipoles in [False, True]:
            for eta in [0., .08]:
                base = run(mode, dipoles, eta, 'base')
                for variant, fermi, shift in [('header', .4, 0.), ('shift', 0., .4), ('both', .15, .25)]:
                    shifted = run(mode, dipoles, eta, variant, fermi, shift)
                    for name in base:
                        worst = max(worst, close(base[name], shifted[name], mode+' '+name))
                # Spectral broadening must leave the transition matrix elements unchanged.
                broadened = run(mode, dipoles, eta, 'cshift', cshift=.17)
                for name in base:
                    close(base[name], broadened[name], 'CSHIFT independence '+mode+' '+name)
                if mode != 'bz':
                    unscissored = run(mode, dipoles, eta, 'scissor', scissor=0.)
                    for name in ['ipa_oscf.dat', 'ipa_oscf-pol.dat']:
                        close([r[1:] for r in base[name]], [r[1:] for r in unscissored[name]], 'scissor independence')
                        close([[r[0]-.6] for r in base[name]], [[r[0]] for r in unscissored[name]], 'scissor transition shift')
        results[mode] = {'max_reference_shift_difference': worst}
    # With one transition at each k, OPT_BZ and IPA must report the same strengths.
    ipa = rows(root/'ipa-d1-s0.08-base/out/ipa_oscf.dat')
    for j, component in enumerate(['x', 'y', 'z'], start=1):
        bz = rows(root/f'bz-d1-s0.08-base/out/bz_act_{component}.dat')
        # These routines enumerate equivalent reciprocal-grid images in different orders.
        close(sorted([[r[j]] for r in ipa]), sorted([[r[3]] for r in bz]), 'IPA vs OPT_BZ '+component)
    report = {'runs': runs, 'results': results, 'scope': 'synthetic two-orbital insulator; no material validation'}
    (root/'results.json').write_text(json.dumps(report, indent=2)+'\n')
    print(f'nonorth_optics_app: PASS ({runs} runs); artifacts: {root}')


if __name__ == '__main__':
    main()
