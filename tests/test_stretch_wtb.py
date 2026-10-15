#!/usr/bin/env python3
"""utils/stretch_wtb.py on a generated two-orbital non-orthogonal model (needs numpy only, no wtb.x).

Usage: python3 tests/test_stretch_wtb.py

For alpha and E0 the stretched file must hold H' = alpha H + (1 - alpha) E0 S and the unchanged S, and its band energies at several k must be
E0 + alpha (E - E0) with the same eigenvectors (|<c'|S|c>| = 1), in the binary and the text layout and across them; the Fermi level of the
header goes to E0 + alpha (E_F - E0), the scissor stays; alpha = 1 gives the same binary file; alpha <= 0 and an output that is the input are
refused.
"""
import contextlib
import io
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'utils'))
import stretch_wtb
import wtb_tbfile as tbf

H = [[[-1., .35], [.35, 1.6]], [[0., .4], [.1, 0.]], [[0., .1], [.4, 0.]]]        # images R = 0, +1, -1 along x
S = [[[1., .08], [.08, 1.]], [[.02, .07], [.01, .015]], [[.02, .01], [.07, .015]]]
RVEC = [[0., 0., 0.], [1., 0., 0.], [-1., 0., 0.]]
LATTICE = [[1., 0., 0.], [0., 4., 0.], [0., 0., 4.]]
FERMI, SCISSOR = -0.4, 0.3
failures = []


def check(ok, what):
    print(('ok   ' if ok else 'FAIL ') + what)
    if not ok:
        failures.append(what)


def write_model(path):
    w = tbf.BinaryWriter(str(path), 'NP', SCISSOR, FERMI, LATTICE, RVEC, 2)
    for h, s in zip(H, S):
        w.add_block(np.array(h), np.zeros((2, 2)), np.array(s))
    w.close()


def bands(tb, kx):
    """generalized eigenvalues of H(k) c = E S(k) c and the eigenvectors, k along x"""
    ph = np.exp(1j * kx * tb['rvec'][:, 0])
    Hk, Sk = np.tensordot(ph, tb['H'], axes=(0, 0)), np.tensordot(ph, tb['S'], axes=(0, 0))
    L = np.linalg.cholesky(Sk)
    Li = np.linalg.inv(L)
    e, u = np.linalg.eigh(Li @ Hk @ Li.conj().T)
    return e, Li.conj().T @ u, Sk


def run(src, dst, *args):
    with contextlib.redirect_stdout(io.StringIO()):
        return stretch_wtb.main([str(src), str(dst), *args])


with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    binary, text = tmp / 'tb.bin', tmp / 'tb.dat'
    write_model(binary)
    tbf.binary_to_text(str(binary), str(text))
    orig = tbf.read_tb(str(binary))

    for alpha, e0 in ((1.3, None), (0.8, -3.0), (1.4, 0.5)):
        e0v = orig['fermi'] if e0 is None else e0
        extra = [] if e0 is None else ['--e0', str(e0)]
        outs = {'binary': tmp / 'a.bin', 'text': tmp / 'a.dat', 'text->binary': tmp / 'b.bin', 'binary->text': tmp / 'b.dat'}
        run(binary, outs['binary'], '--alpha', str(alpha), *extra)
        run(text, outs['text'], '--alpha', str(alpha), *extra)
        run(text, outs['text->binary'], '--alpha', str(alpha), '--format', 'binary', *extra)
        run(binary, outs['binary->text'], '--alpha', str(alpha), '--format', 'text', *extra)
        for name, path in outs.items():
            tb = tbf.read_tb(str(path))
            ref = alpha * orig['H'] + (1 - alpha) * e0v * orig['S']
            tag = 'alpha %.1f, E0 %s, %s' % (alpha, 'E_F' if e0 is None else e0, name)
            check(np.abs(tb['H'] - ref).max() < 2e-6, tag + ': H = alpha H + (1-alpha) E0 S')
            check(np.abs(tb['S'] - orig['S']).max() < 1e-6, tag + ': S unchanged')
            check(abs(tb['fermi'] - (e0v + alpha * (orig['fermi'] - e0v))) < 1e-6 and abs(tb['scissor'] - SCISSOR) < 1e-6,
                  tag + ': Fermi level to E0 + alpha (E_F - E0), scissor as it was')
            check(np.abs(tb['rvec'] - orig['rvec']).max() < 1e-6 and np.abs(tb['lattice'] - orig['lattice']).max() < 1e-6,
                  tag + ': translations and lattice unchanged')
        tb = tbf.read_tb(str(outs['binary']))
        worst, overlap = 0.0, 1.0
        for kx in (0.0, 0.7, 1.9, 3.1):
            e0k, c0, s0 = bands(orig, kx)
            e1k, c1, s1 = bands(tb, kx)
            worst = max(worst, np.abs(e1k - (e0v + alpha * (e0k - e0v))).max())
            overlap = min(overlap, np.abs(np.abs(np.einsum('ia,ij,jb->ab', c1.conj(), s1, c0)).diagonal()).min())
        check(worst < 1e-5 and overlap > 1 - 1e-5,
              'alpha %.1f, E0 %s: band energies E0 + alpha (E - E0) to %.1e eV, eigenvectors unchanged (overlap %.7f)' % (
                  alpha, 'E_F' if e0 is None else e0, worst, overlap))

    run(binary, tmp / 'one.bin', '--alpha', '1')
    check((tmp / 'one.bin').read_bytes() == binary.read_bytes(), 'alpha = 1 writes the same binary file')

    for args, what in (((str(binary), str(tmp / 'z.bin'), '--alpha', '0'), 'alpha = 0 refused'),
                       ((str(binary), str(tmp / 'z.bin'), '--alpha', '-1'), 'alpha < 0 refused'),
                       ((str(binary), str(binary), '--alpha', '1.2'), 'output = input refused')):
        try:
            with contextlib.redirect_stderr(io.StringIO()):
                stretch_wtb.main(list(args))
            check(False, what)
        except SystemExit as e:
            check(e.code not in (None, 0), what)

print('\n%s' % ('all passed' if not failures else '%d FAILED' % len(failures)))
sys.exit(1 if failures else 0)
