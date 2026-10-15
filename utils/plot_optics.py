#!/usr/bin/env python3
"""Optical conductivity and absorptance of a WanTiBEXOS calculation of a two-dimensional sheet, IPA and BSE.

  python3 plot_optics.py --run "stretched 1.3=../03-bexos/out" --run "DFT=../03-bexos-dft/out" --input ../03-bexos/input.dat

One or more --run LABEL=DIRECTORY (up to four; the directory is the OUTPUT of the BeXos run, which holds ipa_*.dat and bse_*.dat).
Figure (optics.png, or -o): rows IPA and BSE, columns Re sigma / sigma_0 and the sheet absorptance (%), with the in-plane average of
the xx and yy components. A table of the same curves is written next to it (optics.csv) and the position and height of the main
peak (0.5-7 eV) of Re sigma / sigma_0 is printed for every curve.

  * BeXos writes the optical conductivity *_opt_cond_real.dat in s^-1 (Gaussian units) for the slab of thickness THMAX;
    the sheet conductivity is sigma_sheet = 4 pi eps0 THMAX sigma[s^-1] (S) and sigma_0 = e^2/(4 hbar) = 6.0853e-5 S.
    THMAX is taken from --input (the BeXos input file) or --thmax, 2.0e-9 m if neither is given: it must be the THMAX of the run.
  * The sheet absorptance is the *_absorptance.dat of BeXos (1 - T - R of the sheet between NI and NS, as THMAX, NI and NS of the
    input define them), in-plane average of xx and yy.

Needs numpy and matplotlib only.
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

EPS0 = 8.8541878128e-12          # F/m
SIGMA0 = 6.085337018e-5          # e^2 / (4 hbar), S
COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#e87ba4']      # blue, orange, aqua, magenta
SURFACE, INK, INK2, MUTED = '#fcfcfb', '#0b0b0b', '#52514e', '#8a897f'

plt.rcParams.update({
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
    'text.color': INK, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
    'axes.edgecolor': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
    'axes.grid': True, 'grid.color': '#e4e3dd', 'grid.linewidth': 0.6, 'axes.axisbelow': True,
    'font.size': 10, 'axes.titlesize': 11, 'axes.titleweight': 'normal', 'legend.frameon': False,
})


def thmax_from_input(path):
    for line in open(path):
        m = re.match(r'\s*THMAX\s*=\s*([-+0-9.eEdD]+)', line)
        if m:
            return float(m.group(1).replace('d', 'e').replace('D', 'e'))
    raise SystemExit('THMAX= not found in %s' % path)


def table(path):
    if not Path(path).exists():
        raise SystemExit('%s not found: is this the OUTPUT directory of a BeXos run with SPEC= T and %s?' % (
            path, 'IPA= T' if 'ipa_' in Path(path).name else 'BSE= T'))
    return np.loadtxt(path)


def curves(folder, prefix, thmax):
    """energy (eV), Re sigma_sheet / sigma_0 and absorptance (%), in-plane average of xx and yy"""
    sr = table(Path(folder) / ('%s_opt_cond_real.dat' % prefix))
    ab = table(Path(folder) / ('%s_absorptance.dat' % prefix))
    e = sr[:, 0]
    sigma = 0.5 * (sr[:, 1] + sr[:, 2]) * 4 * np.pi * EPS0 * thmax / SIGMA0
    absorptance = 100 * 0.5 * (ab[:, 1] + ab[:, 2])
    return e, sigma, absorptance


def peak(e, y, lo=0.5, hi=7.0):
    m = (e >= lo) & (e <= hi)
    i = np.flatnonzero(m)[np.argmax(y[m])]
    return e[i], y[i]


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--run', action='append', required=True, metavar='LABEL=DIR', help='BeXos OUTPUT directory, with a label')
    ap.add_argument('--input', help='the BeXos input file (for THMAX)')
    ap.add_argument('--thmax', type=float, help='THMAX of the run, metres (default: from --input, else 2.0e-9)')
    ap.add_argument('--emax', type=float, default=10.0, help='upper end of the energy axis, eV (default 10)')
    ap.add_argument('-o', '--output', default='optics.png', help='figure (default optics.png); the table gets the same name with .csv')
    args = ap.parse_args(argv)
    if len(args.run) > len(COLORS):
        ap.error('at most %d runs' % len(COLORS))
    thmax = args.thmax or (thmax_from_input(args.input) if args.input else 2.0e-9)
    print('THMAX = %g m' % thmax)

    runs = []
    for item in args.run:
        label, _, folder = item.rpartition('=')
        runs.append((label or Path(folder).parent.name, folder))

    fig, axs = plt.subplots(2, 2, figsize=(11, 7.4), sharex=True)
    columns, header = [], ['energy_eV']
    for r, (prefix, name) in enumerate((('ipa', 'IPA'), ('bse', 'BSE'))):
        for (label, folder), color in zip(runs, COLORS):
            e, s, a = curves(folder, prefix, thmax)
            if not columns:
                columns.append(e)
            columns += [s, a]
            header += ['%s_%s_ReSigma_over_sigma0' % (name, label.replace(' ', '_')), '%s_%s_absorptance_percent' % (name, label.replace(' ', '_'))]
            for c, y in enumerate((s, a)):
                axs[r, c].plot(e, y, color=color, lw=1.6, label=label)
            pe, py = peak(e, s)
            print('%s, %-24s peak of Re sigma / sigma_0 in 0.5-7 eV: %.3f eV (%.2f); absorptance there %.2f %%' % (
                name, label, pe, py, a[np.argmin(abs(e - pe))]))
        axs[r, 0].set_ylabel('Re $\\sigma$ / $\\sigma_0$,  $\\sigma_0 = e^2/4\\hbar$')
        axs[r, 1].set_ylabel('Sheet absorptance (%)')
        axs[r, 0].set_title('%s optical conductivity' % name)
        axs[r, 1].set_title('%s absorptance' % name)
    for ax in axs.ravel():
        ax.set_xlim(0, args.emax)
        ax.set_ylim(bottom=0)
        visible = [ln.get_ydata()[ln.get_xdata() <= args.emax].max() for ln in ax.get_lines()]
        ax.set_ylim(0, 1.12 * max(visible))
    for ax in axs[1]:
        ax.set_xlabel('Photon energy (eV)')
    h, lab = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc='upper center', ncol=len(runs), labelcolor=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(args.output, dpi=170)
    csv = str(Path(args.output).with_suffix('.csv'))
    np.savetxt(csv, np.column_stack(columns), delimiter=',', header=','.join(header), comments='', fmt='%.8g')
    print('wrote %s and %s' % (args.output, csv))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
