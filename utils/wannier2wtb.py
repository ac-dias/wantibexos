#!/usr/bin/env python3
"""
wannier2wtb.py: a Wannier90 model as WanTiBEXOS input (DFT= "W").

    python3 wannier2wtb.py seedname [--efermi E]
                           [--formfactor --mesh NGX NGY NGZ [--ff-ecut 100]
                            [--ff-radius 5.0] [--yes]]

Reads seedname.wout (the lattice vectors, the final centres and spreads) and
seedname_hr.dat, and writes seedname-NP.dat: the header of the WanTiBEXOS
tight-binding file (NP, scissor 0, the Fermi level --efermi, a1-a3 in A)
followed by the hr.dat. The Wannier90 position matrix seedname_r.dat
(write_rmn = .true.) is BSE_CENTER_FILE as it is. Unpolarized models only;
the hr.dat must not need use_ws_distance (WanTiBEXOS does not apply
seedname_wsvec.dat).

--formfactor also writes seedname_ff.bin, the form factors
<m,0|exp(iQ.r)|n,R> of the Wannier functions (utils/wtb_formfactor.py), for
the direct term of the BSE on the k mesh --mesh and the exchange term up to
--ff-ecut eV. They are integrated over the Wannier functions Wannier90 plots
on a real-space grid: seedname_00001.xsf, ... from wannier_plot = .true.,
wannier_plot_format = xcrysden, wannier_plot_mode = crystal and a
wannier_plot_supercell that holds every function, on the UNK files of
pw2wannier90 with write_unk = .true. and reduce_unk = .false. (the grid of
the density, on which the sums are exact). Wannier90 divides each function
by its phase where it is largest and writes the real part, so the functions
are normalised on the grid and their phases are recovered from
seedname_r.dat: the position matrix of the plotted functions against
Wannier90's own, pair by pair along the largest elements. Real Wannier
functions only (no spin-orbit; the Im/Re ratio Wannier90 prints is small).
The radius of function m is --ff-radius times the square root of its spread;
the radii set the lattice vectors of the file, whose size is printed (and
above 1 GB confirmed) before the plots are read; --yes skips the question.
"""
import argparse
import os
import re
import sys

import numpy as np


def read_wout(path):
    """lattice (rows a_i, A), final centres (A) and spreads (A^2), k grid"""
    txt = open(path).read()
    lat = np.zeros((3, 3))
    for i, x, y, z in re.findall(r"a_([123])\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)", txt)[:3]:
        lat[int(i) - 1] = float(x), float(y), float(z)
    final = txt[txt.rindex("Final State"):] if "Final State" in txt else txt
    rows = re.findall(r"WF centre and spread\s+(\d+)\s+\(\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+)\s*\)\s+([\d.]+)",
                      final)
    nw = max(int(r[0]) for r in rows)
    centres, spreads = np.zeros((nw, 3)), np.zeros(nw)
    for r in rows[:nw]:
        centres[int(r[0]) - 1] = float(r[1]), float(r[2]), float(r[3])
        spreads[int(r[0]) - 1] = float(r[4])
    g = re.search(r"Grid size\s*=\s*(\d+)\s*x\s*(\d+)\s*x\s*(\d+)", txt)
    nk = tuple(int(v) for v in g.groups()) if g else None
    return lat, centres, spreads, nk


def read_hr(path):
    """num_wann and the lattice vectors of a Wannier90 hr.dat"""
    L = open(path).read().split("\n")
    nw, nR = int(L[1]), int(L[2])
    i = 3 + (nR + 14) // 15
    R = np.array([[int(v) for v in l.split()[:3]] for l in L[i:i + nR * nw * nw:nw * nw]])
    return nw, R


def read_rdat(path):
    """{(R1, R2, R3, m, n): <m,0|r|n,R> (3,) complex} of a Wannier90 _r.dat"""
    L = open(path).read().split("\n")
    nw, nR = int(L[1]), int(L[2])
    out = {}
    for l in L[3:3 + nR * nw * nw]:
        s = l.split()
        v = [float(x) for x in s[5:11]]
        out[tuple(int(x) for x in s[:5])] = np.array([v[0] + 1j * v[1], v[2] + 1j * v[3], v[4] + 1j * v[5]])
    return out


def read_xsf(path):
    """grid dimensions, origin (A), spanning vectors (A) and values (x fastest)"""
    L = open(path).read().split("\n")
    i = next(j for j, l in enumerate(L) if "BEGIN_DATAGRID_3D" in l)
    n = [int(v) for v in L[i + 1].split()]
    origin = np.array([float(v) for v in L[i + 2].split()])
    span = np.array([[float(v) for v in L[i + 3 + k].split()] for k in range(3)])
    vals = []
    for l in L[i + 6:]:
        if "END_DATAGRID" in l:
            break
        vals.extend(l.split())
    data = np.array(vals[:n[0] * n[1] * n[2]], dtype=float)
    return n, origin, span, data.reshape(n[2], n[1], n[0]).transpose(2, 1, 0)


def write_params(path, lat, efermi, hr):
    with open(path, "w") as f:
        f.write("NP\n{:.6f}\n{:.6f}\n".format(0.0, efermi))
        f.write("".join("{:16.10f} {:16.10f} {:16.10f}\n".format(*v) for v in lat))
        f.write(open(hr).read())


def plotted_groups(wff, seed, lat, centres, radius, nk):
    """the plotted Wannier functions as window groups on their common grid"""
    groups, grid, norms = [], None, []
    for m in range(len(centres)):
        n, origin, span, data = read_xsf("{}_{:05d}.xsf".format(seed, m + 1))
        # the spanning vectors run from the first grid point to the last: n - 1 steps
        M = np.rint(np.linalg.norm(lat, axis=1) * (np.array(n) - 1)
                    / np.linalg.norm(span, axis=1)).astype(int)
        S = np.zeros(3, dtype=int)
        for i in range(3):
            if n[i] % M[i] == 0:                 # Wannier90: S cells of M points
                S[i] = n[i] // M[i]
            elif (n[i] - 1) % M[i] == 0:         # the first point repeated at the end
                S[i] = (n[i] - 1) // M[i]
                data = np.take(data, range(n[i] - 1), axis=i)
            else:
                raise SystemExit("{}: the grid is not commensurate with the lattice".format(seed))
        if grid is None:
            grid = wff.Grid(lat, M)
        elif np.any(grid.M != M):
            raise SystemExit("the plots of {} do not share one grid".format(seed))
        f0 = origin @ np.linalg.inv(lat) * M
        lo = np.rint(f0).astype(int)
        if np.abs(f0 - lo).max() > 1e-3:
            raise SystemExit("{}: the origin of the plot is not a grid point".format(seed))
        # the window is one period along a_i when it spans the k grid there
        per = [S[i] == (nk[i] if nk else 0) for i in range(3)]
        norm = float(np.sum(data ** 2)) * grid.dV
        norms.append(norm)
        groups.append(wff.Group([m], centres[m], radius[m], lo=lo,
                                values=(data / np.sqrt(norm))[None].astype(np.float32), periodic=per))
    return groups, grid, np.array(norms)


def pair_dipole(wff, grid, A, B, R):
    """<m,0| r |n,R> of two plotted functions (window groups)"""
    rc = np.asarray(R) @ grid.lat
    lo, shape = wff._pair_region(grid, A, B, B.centre + rc)
    wa, wb = A.on(lo, shape)[0], B.on(lo - np.asarray(R) * grid.M, shape)[0]
    r = grid.points(lo, shape)
    return np.einsum("abc,abcx->x", np.conj(wa) * wb, r) * grid.dV


def fix_phases(wff, grid, groups, rdat):
    """p_n with w_n = p_n w~_n (w~ the plot), from r_mn(R) = conj(p_m) p_n r~_mn(R),
    along the largest off-diagonal elements of seedname_r.dat"""
    nw = len(groups)
    cand = sorted(((float(np.linalg.norm(v)), k) for k, v in rdat.items()
                   if (k[3] != k[4] or any(k[:3])) and abs(k[3] - k[4]) + sum(map(abs, k[:3])) > 0),
                  reverse=True)
    p = {0: 1.0 + 0j}
    checks = []
    while len(p) < nw:
        for size, k in cand:
            m, n = k[3] - 1, k[4] - 1
            if (m in p) == (n in p) or size < 1e-6:
                continue
            R = np.array(k[:3])
            rt = pair_dipole(wff, grid, groups[m], groups[n], R)
            r = rdat[k]
            c = np.vdot(rt, r) / np.vdot(rt, rt)        # r = c rt, c = conj(p_m) p_n
            if m in p:
                p[n] = p[m] * c / abs(c)
            else:
                p[m] = p[n] * np.conj(c) / abs(c)
            checks.append((k, abs(c)))
            break
        else:
            raise SystemExit("the phases of functions {} are not fixed by seedname_r.dat".format(
                sorted(set(range(nw)) - set(p))))
    ph = np.array([p[m] for m in range(nw)])
    for m, g in enumerate(groups):
        g.values = (g.values * ph[m]).astype(np.complex64)
    return ph, checks


def main():
    ap = argparse.ArgumentParser(description="Wannier90 model -> WanTiBEXOS input (DFT= \"W\")")
    ap.add_argument("seedname")
    ap.add_argument("--efermi", type=float, default=0.0, help="Fermi level written in the file (eV)")
    ap.add_argument("--formfactor", action="store_true",
                    help="also write seedname_ff.bin, the form factors of the Wannier functions")
    ap.add_argument("--mesh", type=int, nargs=3, default=None, help="the k mesh of the BSE (NGX NGY NGZ)")
    ap.add_argument("--ff-ecut", type=float, default=100.0, help="exchange G up to this energy (eV)")
    ap.add_argument("--ff-radius", type=float, default=5.0, help="radius / sqrt(spread) of each function")
    ap.add_argument("--yes", action="store_true", help="write the form factors even above 1 GB")
    a = ap.parse_args()
    if a.formfactor and a.mesh is None:
        ap.error("--formfactor needs --mesh NGX NGY NGZ, the k-point mesh of the BSE")
    seed = a.seedname
    lat, centres, spreads, nk = read_wout(seed + ".wout")
    nw, Rhr = read_hr(seed + "_hr.dat")
    if nw != len(centres):
        raise SystemExit("{}.wout has {} functions, {}_hr.dat {}".format(seed, len(centres), seed, nw))
    if a.formfactor:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import wtb_formfactor as wff
        radius = a.ff_radius * np.sqrt(spreads)
        R = wff.lattice_vectors(lat, centres, radius, periodic=[bool(np.any(Rhr[:, i])) for i in range(3)])
        Qd, Qe = wff.direct_q(lat, a.mesh), wff.exchange_g(lat, a.ff_ecut)
        ffname = seed + "_ff.bin"
        wff.confirm_size(ffname, nw, len(R), len(Qd) + len(Qe), a.yes, len(Qd), len(Qe))
    if os.path.exists(seed + "_wsvec.dat"):
        # after each "R1 R2 R3 m n" line: the number of images, then their shifts
        rows = [l.split() for l in open(seed + "_wsvec.dat").read().split("\n")[1:] if l.strip()]
        if any((len(r) == 1 and r[0] != "1") or (len(r) == 3 and r != ["0", "0", "0"]) for r in rows):
            print("warning: {}_wsvec.dat has images other than R itself (use_ws_distance); "
                  "WanTiBEXOS does not apply them".format(seed))
    write_params(seed + "-NP.dat", lat, a.efermi, seed + "_hr.dat")
    print("wrote {0}-NP.dat (DFT= \"W\", PARAMS_FILE); BSE_CENTER_FILE: {0}_r.dat".format(seed))
    if not a.formfactor:
        return 0

    if not os.path.exists(seed + "_r.dat"):
        raise SystemExit("--formfactor needs {}_r.dat (write_rmn = .true.) for the phases".format(seed))
    groups, grid, norms = plotted_groups(wff, seed, lat, centres, radius, nk)
    ph, checks = fix_phases(wff, grid, groups, read_rdat(seed + "_r.dat"))
    kept = []
    for g in groups:                     # the norm within the radius
        lo, shape = grid.box(g.centre, g.radius, maxwidth=g.period)
        pts = grid.points(lo, shape) - g.centre
        w = g.on(lo, shape)[0]
        kept.append(float(np.sum(np.abs(w) ** 2 * (np.sum(pts * pts, axis=-1) <= g.radius ** 2))) * grid.dV)
    wr = wff.FFWriter(ffname, lat, centres, R, np.vstack([Qd, Qe]), len(Qd), a.mesh, a.ff_ecut)
    wff.compute(wr, grid, groups, R, [(Qd, 0), (Qe, len(Qd))])
    wr.close()
    ff = wff.read_ff(ffname)
    iq0 = int(np.argmin(np.abs(Qd).sum(axis=1)))
    S = np.asarray(ff["F"][iq0], dtype=complex)
    dS = max(float(np.abs(S[iR] - (np.eye(nw) if not np.any(Ri) else 0.0)).max()) for iR, Ri in enumerate(R))
    rdat = read_rdat(seed + "_r.dat")
    dr = max(float(np.abs(pair_dipole(wff, grid, groups[k[3] - 1], groups[k[4] - 1], np.array(k[:3]))
                          - rdat[k]).max()) for k, _ in checks)
    print("form factors: {} written (grid {}x{}x{} per cell); plots normalised by {:.4g}-{:.4g}; "
          "norm within the radii {:.5f}-{:.5f}".format(ffname, *grid.M, 1 / np.sqrt(norms.max()),
                                                        1 / np.sqrt(norms.min()), min(kept), max(kept)))
    print("  phases of the plots: {}; max |r_mn(R) - Wannier90's| on the {} pairs that fixed them: {:.1e} A".format(
        " ".join("{:+.0f}".format(np.angle(p, deg=True)) for p in ph), len(checks), dr))
    print("  max |F(R;0) - delta| {:.1e}".format(dS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
