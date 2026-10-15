#!/usr/bin/env python3
"""Stretch the bands of a tight-binding file of WanTiBEXOS (DFT= "S": H(R), S(R) of a non-orthogonal basis) by a factor
alpha about the energy E0, before wtb.x reads it.

  python3 stretch_wtb.py tb-NP.bin tb-NP-s1.2.bin --alpha 1.2                  # about the Fermi level of the header
  python3 stretch_wtb.py tb-NP.dat tb-NP-s1.2.dat --alpha 1.2 --e0 -4.5         # about E0 = -4.5 eV (text layout)
  python3 stretch_wtb.py tb-NP.dat tb-NP-s1.2.bin --alpha 1.2 --format binary   # and into the other layout

Every block of the file becomes

  H'(R) = alpha H(R) + (1 - alpha) E0 S(R),        S'(R) = S(R),

so that H'(k) c = [E0 + alpha (E - E0)] S(k) c for every eigenpair of H(k) c = E S(k) c: all eigenvalues move to
E0 + alpha (E - E0), the eigenvectors are those of the original file, and the gap and the widths of all the bands are
multiplied by alpha (alpha > 1 widens them, 0 < alpha < 1 narrows them). The matrices of an orthogonal basis (S = delta_R0
times the identity) are stretched in the same way. The Fermi level of the header goes to E0 + alpha (E_F - E0), the same
place in the stretched spectrum (it stays where it is for the default E0 = E_F, with the same occupations, whatever
alpha); the scissor of line 2 is not changed, and wtb.x adds it to the conduction bands of the stretched model.

What does not change: the eigenvectors, so the Berry curvature and the interband dipoles. The optical vertex of wtb.x,
<c|dH/dk - (Ec+Ev)/2 dS/dk|v> / (Ec-Ev+i SIGMA), has dH/dk of the stretched H and Ec, Ev of the stretched bands: the factor
alpha cancels for SIGMA = 0. For SIGMA > 0 its regularization is 1/alpha as large; multiply SIGMA and CSHIFT by alpha to keep
them as they were. The position matrix (BSE_CENTER_FILE) and the form-factor file (BSE_FF_FILE) belong to the basis, not to
H, and stay as they are. The kernel of the BSE (COULOMB_POT, EDIEL) is not stretched.

One alpha for all the bands: a stretch of the valence bands different from that of the conduction bands is not a linear
operation on H(R) (it needs the projectors on the bands, which depend on k) and is not done here.

The output has the layout of the input unless --format says otherwise. The values are those of the input in double precision,
written in single precision in the binary layout (the precision wtb.x keeps) and with 12 digits in the text layout. Numpy
is needed (as for wtb_tbfile.py, which this script uses to read and write the file).
"""
import argparse
import os
import sys

import numpy as np

import wtb_tbfile as tbf


def stretch_block(reH, imH, S, alpha, e0):
    """H' = alpha H + (1 - alpha) E0 S on the elements of one lattice vector (Re H, Im H, S as read), double precision"""
    S = np.asarray(S, dtype=np.float64)
    return alpha * np.asarray(reH, dtype=np.float64) + (1.0 - alpha) * e0 * S, alpha * np.asarray(imH, dtype=np.float64)


def stretch_fermi(fermi, alpha, e0):
    return e0 + alpha * (fermi - e0)


def _images_with_translation(fname, binary):
    """(header, generator of (i, j, Re H, Im H, S, R) per lattice vector, R the Cartesian translation)"""
    if binary:
        hd, images = tbf.iter_binary(fname)
        return hd, ((i, j, re, im, s, hd['rvec'][iR]) for iR, (i, j, re, im, s) in enumerate(images))
    return tbf.iter_text(fname)


def _stretched(images, alpha, e0):
    for i, j, re, im, s, r in images:
        re2, im2 = stretch_block(re, im, s, alpha, e0)
        yield i, j, re2, im2, s, r


def _text_header_lines(fname):
    with open(fname) as f:
        return [f.readline() for _ in range(9)]


def _write_binary(fname, hd, fermi, images):
    if 'rvec' in hd:                                      # binary input: every translation is in the header, stream
        w = tbf.BinaryWriter(fname, hd['systype'], hd['scissor'], fermi, hd['lattice'], hd['rvec'], hd['nbasis'])
        for i, j, re, im, s, r in images:
            w.add(i, j, re, im, s)
    else:                                                 # the header of the binary file needs every translation first
        rv, kept = [], []
        for i, j, re, im, s, r in images:
            keep = (re != 0) | (im != 0) | (s != 0)
            kept.append((i[keep], j[keep], re[keep], im[keep], s[keep]))
            rv.append(r)
        w = tbf.BinaryWriter(fname, hd['systype'], hd['scissor'], fermi, hd['lattice'], np.array(rv), hd['nbasis'])
        for c in kept:
            w.add(*c)
    w.close()
    return w.stored


def _write_text(fname, hd, fermi, images, raw_header):
    n = hd['nbasis']
    row, col = np.indices((n, n))
    row, col = row.T.ravel() + 1, col.T.ravel() + 1       # text order: column outer, row (the first index) inner
    with open(fname, 'w') as f:
        if raw_header is not None:                        # text input: its header as it is, but for the Fermi level
            raw_header[2] = raw_header[2] if fermi == hd['fermi'] else '%.12e\n' % fermi
            f.writelines(raw_header)
        else:
            f.write('%s\n%.12e\n%.12e\n' % (hd['systype'], hd['scissor'], fermi))
            for r in hd['lattice']:
                f.write('%.12e  %.12e  %.12e\n' % tuple(r))
            f.write('%d\n%d\n#rcell x   rcell y   rcell z   i   j   ReH   ImH   S\n' % (n, hd['nvec']))
        for i, j, re, im, s, r in images:
            blk = np.zeros((3, n, n))                     # [component, row, column]
            blk[0][i - 1, j - 1], blk[1][i - 1, j - 1], blk[2][i - 1, j - 1] = re, im, s
            table = np.column_stack([np.tile(r, (n * n, 1)), row, col, blk[0].T.ravel(), blk[1].T.ravel(), blk[2].T.ravel()])
            np.savetxt(f, table, fmt=['%.12e'] * 3 + ['%d'] * 2 + ['%.12e'] * 3, delimiter='   ')


def main(argv):
    ap = argparse.ArgumentParser(description='Stretch the bands of a WanTiBEXOS tight-binding file (DFT= "S") by alpha about E0: '
                                 'H\' = alpha H + (1 - alpha) E0 S.')
    ap.add_argument('input', help='tight-binding file, binary or text layout (detected)')
    ap.add_argument('output', help='stretched tight-binding file')
    ap.add_argument('--alpha', type=float, required=True, help='stretch factor of the band energies about E0 (> 0)')
    ap.add_argument('--e0', type=float, default=None, help='energy of the fixed point, eV, in the energy reference of the file '
                    '(default: the Fermi level of the header)')
    ap.add_argument('--format', choices=('binary', 'text'), default=None, help='layout of the output (default: that of the input)')
    args = ap.parse_args(argv)

    if not args.alpha > 0.0:
        ap.error('--alpha must be positive (alpha <= 0 inverts or collapses the bands)')
    if os.path.exists(args.output) and os.path.samefile(args.input, args.output):
        ap.error('the output file is the input file')

    binary = tbf.is_binary(args.input)
    hd, images = _images_with_translation(args.input, binary)
    e0 = hd['fermi'] if args.e0 is None else args.e0
    fermi = stretch_fermi(hd['fermi'], args.alpha, e0)
    out_binary = binary if args.format is None else (args.format == 'binary')
    images = _stretched(images, args.alpha, e0)

    if out_binary:
        stored = _write_binary(args.output, hd, fermi, images)
        layout = 'binary, %d stored elements' % stored
    else:
        _write_text(args.output, hd, fermi, images, None if binary else _text_header_lines(args.input))
        layout = 'text'
    print('%s: %s, N = %d, %d lattice vectors; alpha = %g about E0 = %.6f eV (%s); Fermi level %.6f -> %.6f eV; scissor %.6f eV '
          'as it was' % (args.output, layout, hd['nbasis'], hd['nvec'], args.alpha, e0,
                         'Fermi level of the header' if args.e0 is None else '--e0', hd['fermi'], fermi, hd['scissor']))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
