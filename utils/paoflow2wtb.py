#!/usr/bin/env python3
"""
paoflow2wtb.py: the PAOFLOW tight-binding Hamiltonian of a Quantum ESPRESSO
run, written as WanTiBEXOS input.

    python3 paoflow2wtb.py prefix.save [--configuration minimal] [--basispath DIR]
                           [--pthr 0.95] [--shift auto] [--shift-type 1]
                           [--efermi E] [--seedname paoflow] [--npool 1]
                           [--formfactor --mesh NGX NGY NGZ [--ff-ecut 100]
                            [--ff-spacing A] [--ff-tail 1e-6] [--yes]]
    mpirun -np N python3 paoflow2wtb.py ...      (PAOFLOW's own MPI parallelism)

prefix.save is the save directory of a pw.x run on a Monkhorst-Pack grid.
PAOFLOW projects the Kohn-Sham states on its pseudo-atomic orbitals
(PAOFLOW.projections: --configuration, --basispath), keeps the states that
project well (PAOFLOW.projectability: --pthr, and --shift, the energy eta in
eV above the Fermi level from which the Hamiltonian is filtered, or auto)
and builds H(k) on the grid in the Loewdin-orthonormalised basis
(PAOFLOW.pao_hamiltonian: --shift-type). The script turns H(k) into H(R)
and writes the same model for both tight-binding readers of WanTiBEXOS, in
the current directory:

    DFT= "W"   <seedname>-NP.dat  WanTiBEXOS header + Wannier90 hr.dat
               <seedname>_r.dat   Wannier90 position matrix with the orbital
                                  centres (the atomic positions) on the R = 0
                                  diagonal and nothing else: BSE_CENTER_FILE
    DFT= "S"   tb-NP.dat          the layout of siesta2wtb.py, S(R) = 1
               basis_set-NP       the orbital centres
    <seedname>-info.txt           the orbitals and the checks below
    <seedname>_ff.bin             with --formfactor: the form factors
                                  <m,0|exp(iQ.r)|n,R> of the Loewdin orbitals
                                  (utils/wtb_formfactor.py), for the direct
                                  term of the BSE on the k mesh --mesh and the
                                  exchange term up to --ff-ecut eV

The Loewdin orbitals are the real-space form of PAOFLOW's orthonormalised
atomic Bloch sums at the k-points of its grid (calc_atwfc_k, ortho_atwfc_k,
on the wavefunction cutoff of the pw.x run): |m,0> = N^-1/2 sum_k |m k>, one
FFT per orbital over the supercell of the grid, on enough points per cell
for the whole G sphere of a pair and the largest Q (--ff-spacing can refine
it). The radius of an orbital is where the Loewdin orbital itself has less
than --ff-tail of its norm outside, on the supercell where it is built: the
pairs closer than the sum of their radii make the lattice vectors of the
file, each integrated over the box around both spheres. The Loewdin orbitals
reach further than the pseudo-atomic orbitals they are made of (for the
box-state bases standard and extended, up to 11% of the norm lies outside
the radius of the pseudo-atomic orbital), so these give only a first, smaller
estimate of the file size, printed (and above 1 GB confirmed) before PAOFLOW
starts; the size from the Loewdin radii is printed, and above 1 GB confirmed
again when it grew, before the integrals. --yes skips both questions.

PAOFLOW's Bloch sums carry the structure factor exp(-i(k+G).tau_m) of QE's
atomic wavefunctions, so H_mn(k) = sum_R exp(+ik.R) H_mn(R) with
H_mn(R) = <m,0|H|n,R>, and H(R) = (1/N) sum_k exp(-ik.R) H(k) on the grid of
N k-points, modulo the grid's supercell (PAOFLOW's own HRs = ifftn(Hks) is H
at -R). Each H_mn(R) is put at the supercell images R + T with the shortest
|R + T + tau_n - tau_m|, split equally between ties, as Wannier90's
use_ws_distance does for each pair: H(k) on the grid does not change, and
between the grid points it is the minimal-image interpolation.

Energies are PAOFLOW's: eV, zero at the Fermi level of the DFT run. The
Fermi level written in the files is mid-gap on the grid for an insulator
and 0 otherwise (--efermi overrides it); the scissor line is 0. Only
unpolarized runs (NP). Tested with PAOFLOW 3.0.0 and Quantum ESPRESSO 7.6.
"""
import argparse
import os
import sys

import numpy as np

BOHR_A = 0.529177210903
TWO_PI = 2.0 * np.pi


def ws_fold(HRbox, nk, lat, tau, tol=1e-5):
    """Split every H_mn(R) over the minimal images of R + tau_n - tau_m.

    HRbox: (nawf, nawf, nk1, nk2, nk3), R in the box [0, nk); lat: rows a1,
    a2, a3 (A); tau: (nawf, 3) orbital centres (A). Returns the lattice
    vectors (nR, 3), H(R) (nR, nawf, nawf) and the image count of every
    (m, n, R) of the box.
    """
    nk = np.asarray(nk)
    nawf = HRbox.shape[0]
    rng = [range(-2, 3) if n > 1 else range(0, 1) for n in nk]
    T = np.array([(a, b, c) for a in rng[0] for b in rng[1] for c in rng[2]])
    dtau = tau[None, :, :] - tau[:, None, :]           # tau_n - tau_m, (m, n, 3)
    out = {}
    nimg = np.zeros((nawf, nawf) + tuple(nk), dtype=int)
    for r in np.ndindex(*nk):
        Rc = np.array(r)[None, :] + T * nk[None, :]    # candidate images
        vec = (Rc @ lat)[:, None, None, :] + dtau[None, :, :, :]
        d = np.linalg.norm(vec, axis=-1)               # (nT, m, n)
        ties = d <= d.min(axis=0)[None] + tol
        cnt = ties.sum(axis=0)
        nimg[(slice(None), slice(None)) + r] = cnt
        h = HRbox[(slice(None), slice(None)) + r]
        for it in np.nonzero(ties.any(axis=(1, 2)))[0]:
            key = tuple(int(v) for v in Rc[it])
            out.setdefault(key, np.zeros((nawf, nawf), dtype=complex))
            out[key] += ties[it] / cnt * h
    keys = sorted(out, key=lambda R: (abs(R[0]) + abs(R[1]) + abs(R[2]), R))
    return np.array(keys, dtype=int), np.array([out[k] for k in keys]), nimg


def hr_to_hk(R, HR, kcrys):
    """H(k) = sum_R exp(+2 pi i k.R) H(R) at crystal k-points, (nk, nawf, nawf)"""
    ph = np.exp(1j * TWO_PI * (np.asarray(kcrys) @ R.T))
    return np.einsum("kr,rmn->kmn", ph, HR)


def _lat_lines(lat):
    return "".join("{:16.10f} {:16.10f} {:16.10f}\n".format(*v) for v in lat)


def write_w90(path, R, HR, lat, efermi, header):
    """DFT= "W": NP, scissor, Fermi level, a1-a3 (A), then a Wannier90
    hr.dat (row index fastest, as Wannier90 writes it), degeneracies 1"""
    nR, nawf = len(R), HR.shape[1]
    with open(path, "w") as f:
        f.write("NP\n{:.6f}\n{:.6f}\n".format(0.0, efermi))
        f.write(_lat_lines(lat))
        f.write(header + "\n")
        f.write("{:12d}\n{:12d}\n".format(nawf, nR))
        for i in range(0, nR, 15):
            f.write("".join("{:5d}".format(1) for _ in range(min(15, nR - i))) + "\n")
        for r, h in zip(R, HR):
            for n in range(nawf):
                for m in range(nawf):
                    f.write("{:5d}{:5d}{:5d}{:5d}{:5d}{:16.9f}{:16.9f}\n".format(
                        r[0], r[1], r[2], m + 1, n + 1, h[m, n].real, h[m, n].imag))


def write_rdat(path, tau, header):
    """Wannier90 position matrix: the R = 0 block, the centres (A) on the
    diagonal; WanTiBEXOS takes the centres from it"""
    nawf = len(tau)
    with open(path, "w") as f:
        f.write(header + "\n")
        f.write("{:12d}\n{:12d}\n".format(nawf, 1))
        for n in range(nawf):
            for m in range(nawf):
                x = tau[m] if m == n else (0.0, 0.0, 0.0)
                f.write("{:5d}{:5d}{:5d}{:5d}{:5d}".format(0, 0, 0, m + 1, n + 1) +
                        "".join("{:16.10f}{:16.10f}".format(v, 0.0) for v in x) + "\n")


def write_tbnp(path, R, HR, lat, efermi):
    """DFT= "S": the tb-NP.dat of siesta2wtb.py (R Cartesian in A, row
    index fastest), S = 1 at R = 0"""
    nR, nawf = len(R), HR.shape[1]
    with open(path, "w") as f:
        f.write("NP\n{:.6f}\n{:.6f}\n".format(0.0, efermi))
        f.write(_lat_lines(lat))
        f.write("{}\n{}\n".format(nawf, nR))
        f.write("#rcell x   rcell y   rcell z   i   j   ReH   ImH   S\n")
        for r, rc, h in zip(R, R @ lat, HR):
            zero = not np.any(r)
            for n in range(nawf):
                for m in range(nawf):
                    s = 1.0 if (zero and m == n) else 0.0
                    f.write("{:.10f}   {:.10f}   {:.10f}   {}   {}   {:.9f}   {:.9f}   {:.9f}\n".format(
                        rc[0], rc[1], rc[2], m + 1, n + 1, h[m, n].real, h[m, n].imag, s))


def write_basis_set(path, orbitals):
    """basis_set-NP: index, species, centre (A), l, m, spin"""
    with open(path, "w") as f:
        f.write("bindex aspecie ax ay az l m spin\n")
        for i, o in enumerate(orbitals):
            f.write("{}   {}   {:.10f}   {:.10f}   {:.10f}   {}   {}   0\n".format(
                i + 1, o["species"], o["tau"][0], o["tau"][1], o["tau"][2], o["l"], o["m"]))


def gsphere(kbohr, bg, ecut):
    """Miller indices (3, npw) of |k+G|^2 <= ecut (Ry, bohr), bg columns b_i"""
    nmax = int(np.ceil(np.sqrt(ecut) / np.min(np.linalg.norm(bg, axis=0)))) + 3
    r = np.arange(-nmax, nmax + 1)
    h = np.array(np.meshgrid(r, r, r, indexing="ij")).reshape(3, -1)
    kg = h.T @ bg.T + kbohr
    return h[:, np.sum(kg * kg, axis=1) <= ecut]


def basis_records(pf, a):
    """the basis records PAOFLOW.projections will build, without projecting"""
    from PAOFLOW.projection.do_atwfc_proj import build_aewfc_basis, build_pswfc_basis_all
    arry, attr = pf.data_controller.data_dicts()
    if a.configuration in (None, "minimal"):
        return build_pswfc_basis_all(pf.data_controller)[0]
    from PAOFLOW.inputs.basis_presets import resolve_configuration
    if a.basispath is not None:
        attr["basispath"] = os.path.join(a.basispath, "")
    arry["configuration"] = resolve_configuration(pf.data_controller, a.configuration)
    return build_aewfc_basis(pf.data_controller)[0]


def orbital_radius(rec, tail):
    """radius (A) outside which the radial function r R(r) of the record has
    less than tail of its norm"""
    r, w = np.asarray(rec["r"], float), np.asarray(rec["wfc"], float)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] ** 2 + w[:-1] ** 2) * np.diff(r))])
    return float(r[np.argmax(1.0 - cum / cum[-1] < tail)]) * BOHR_A


def loewdin_radii(full, grid, nk, centre, tail, width=0.01):
    """radius (A) of each orbital of full (the orbitals of one centre on the
    supercell of the k grid, where they are periodic) outside which it has less
    than tail of its norm: |w|^2 against the distance to the centre, the minimum
    image in the supercell, in bins of width"""
    D = np.array(full.shape[1:])
    fc = np.asarray(centre, float) @ np.linalg.inv(grid.lat)       # crystal coordinates
    f = []
    for i in range(3):
        x = np.arange(D[i]) / grid.M[i] - fc[i]
        f.append(x - nk[i] * np.round(x / nk[i]))
    x23 = f[1][:, None, None] * grid.lat[1] + f[2][None, :, None] * grid.lat[2]
    nb = int(0.5 * float(np.sum(nk * np.linalg.norm(grid.lat, axis=1))) / width) + 2
    hist = np.zeros((full.shape[0], nb))
    for i1 in range(D[0]):
        d = np.linalg.norm(x23 + f[0][i1] * grid.lat[0], axis=-1).ravel()
        b = np.minimum((d / width).astype(int), nb - 1)
        for j in range(full.shape[0]):
            hist[j] += np.bincount(b, weights=np.abs(full[j, i1]).ravel() ** 2, minlength=nb)
    out = np.cumsum(hist[:, ::-1], axis=1)[:, ::-1]        # the norm from each bin outwards
    rad = np.empty(full.shape[0])
    for j in range(full.shape[0]):
        small = out[j] < tail * out[j, 0]
        rad[j] = width * (np.argmax(small) if small.any() else nb)
    return rad


def formfactor_setup(pf, a):
    """rank 0: the sets of the form-factor file and a first estimate of its
    size, with the radii of the pseudo-atomic orbitals, before projecting"""
    import wtb_formfactor as wff
    arry, attr = pf.data_controller.data_dicts()
    lat = np.asarray(arry["a_vectors"]) * attr["alat"] * BOHR_A
    recs = basis_records(pf, a)
    tau = np.array([np.asarray(b["tau"], float) * BOHR_A for b in recs])
    rad = np.array([orbital_radius(b, a.ff_tail) for b in recs])
    # no lattice vectors along a direction PAOFLOW's grid does not sample
    nk = (int(attr["nk1"]), int(attr["nk2"]), int(attr["nk3"]))
    R = wff.lattice_vectors(lat, tau, rad, periodic=[n > 1 for n in nk])
    Qd, Qe = wff.direct_q(lat, a.mesh), wff.exchange_g(lat, a.ff_ecut)
    name = a.seedname + "_ff.bin"
    wff.confirm_size(name, len(recs), len(R), len(Qd) + len(Qe), a.yes, len(Qd), len(Qe),
                     hint=", or raise --ff-tail")
    return {"name": name, "lat": lat, "tau": tau, "rad": rad, "R": R, "Qd": Qd, "Qe": Qe}


def formfactors(pf, a, ff):
    """rank 0: the Loewdin orbitals in real space and their form factors"""
    import wtb_formfactor as wff
    from PAOFLOW.projection.do_atwfc_proj import calc_atwfc_k, ortho_atwfc_k
    arry, attr = pf.data_controller.data_dicts()
    basis = arry["basis"]
    lat, tau, rad = ff["lat"], ff["tau"], ff["rad"]
    nk = np.array([int(attr["nk1"]), int(attr["nk2"]), int(attr["nk3"])])
    bg = np.asarray(arry["b_vectors"]).T * 2.0 * np.pi / attr["alat"]      # columns b_i, 1/bohr
    ecut = float(attr["ecutwfc"])
    # enough points for the G spheres of both orbitals of a pair and the largest
    # Q (then the sums have no aliasing); --ff-spacing can only refine that
    hmax = np.abs(gsphere(np.zeros(3), bg, ecut)).max(axis=1) + 1
    qmax = np.ceil(np.abs(np.vstack([ff["Qd"], ff["Qe"]])).max(axis=0)).astype(int)
    grid = wff.Grid.with_spacing(lat, a.ff_spacing or 1e9, minimum=2 * hmax + qmax + 1)
    D = grid.M * nk
    chis = []
    for idx in np.ndindex(*nk):
        kb = (np.array(idx, float) / nk) @ bg.T
        mill = gsphere(kb, bg, ecut)
        gk = {"xk": kb, "igwx": mill.shape[1], "mill": mill, "bg": bg, "gamma_only": False}
        chis.append((tuple((np.array(idx)[:, None] + nk[:, None] * mill) % D[:, None]),
                     ortho_atwfc_k(calc_atwfc_k(basis, gk))))
    norm = float(np.prod(D)) / (float(np.prod(nk)) * np.sqrt(abs(np.linalg.det(lat))))
    # atoms: orbitals with one centre share a box
    centres = []
    for m, t in enumerate(tau):
        for c in centres:
            if np.linalg.norm(c["tau"] - t) < 1e-6:
                c["orbs"].append(m)
                break
        else:
            centres.append({"tau": t, "orbs": [m]})
    # each orbital on one whole period of the grid (the supercell of PAOFLOW's
    # k grid): its tail is kept wherever the other orbital of a pair is large
    for c in centres:
        c["full"] = np.zeros((len(c["orbs"]),) + tuple(D), dtype=np.complex64)
    for m in range(len(tau)):
        C = np.zeros(tuple(D), dtype=complex)
        for K, chi in chis:
            C[K] = chi[m]
        w = np.fft.ifftn(C) * norm
        for c in centres:
            if m in c["orbs"]:
                c["full"][c["orbs"].index(m)] = w
    # the radii of the Loewdin orbitals themselves, and with them the lattice
    # vectors, the boxes and the size of the file
    rad = np.array(rad, dtype=float)
    for c in centres:
        rad[c["orbs"]] = loewdin_radii(c["full"], grid, nk, c["tau"], a.ff_tail)
        c["rad"] = float(rad[c["orbs"]].max())
    R = wff.lattice_vectors(lat, tau, rad, periodic=[n > 1 for n in nk])
    print("form factors: radii of the Loewdin orbitals {:.2f}-{:.2f} A (pseudo-atomic orbitals "
          "{:.2f}-{:.2f} A), {} lattice vectors (first estimate {})".format(
              rad.min(), rad.max(), ff["rad"].min(), ff["rad"].max(), len(R), len(ff["R"])), flush=True)
    wff.confirm_size(ff["name"], len(tau), len(R), len(ff["Qd"]) + len(ff["Qe"]),
                     a.yes or len(R) <= len(ff["R"]), len(ff["Qd"]), len(ff["Qe"]),
                     hint=", or raise --ff-tail")
    ff["R"] = R
    kept = np.zeros(len(tau))
    for c in centres:
        lo, shape = grid.box(c["tau"], c["rad"], maxwidth=D)
        ix = np.ix_(*[(lo[i] + np.arange(shape[i])) % D[i] for i in range(3)])
        for j, m in enumerate(c["orbs"]):
            kept[m] = float(np.sum(np.abs(c["full"][j][ix]) ** 2)) * grid.dV
    groups = [wff.Group(c["orbs"], c["tau"], c["rad"], full=c["full"]) for c in centres]
    Qc = np.vstack([ff["Qd"], ff["Qe"]])
    wr = wff.FFWriter(ff["name"], lat, tau, R, Qc, len(ff["Qd"]), a.mesh, a.ff_ecut)
    wff.compute(wr, grid, groups, R, [(ff["Qd"], 0), (ff["Qe"], len(ff["Qd"]))])
    wr.close()
    # F(R; 0) is the overlap of orthonormal orbitals
    F = wff.read_ff(ff["name"])
    iq0 = int(np.argmin(np.abs(ff["Qd"]).sum(axis=1)))
    S = np.asarray(F["F"][iq0], dtype=complex)
    dS = max(float(np.abs(S[iR] - (np.eye(len(tau)) if not np.any(R) else 0.0)).max())
             for iR, R in enumerate(ff["R"]))
    return ("form factors: {} written (BSE_FF_FILE; grid {}x{}x{} per cell); norm within the radii "
            "{:.6f}-{:.6f}; max |F(R;0) - delta| {:.1e}").format(
                ff["name"], *grid.M, kept.min(), kept.max(), dS)


def convert(pf, a):
    """rank 0: H(k) of the PAOFLOW run -> H(R) -> the files"""
    arry, attr = pf.data_controller.data_dicts()
    if int(attr["nspin"]) != 1 or attr.get("dftSO", False):
        raise SystemExit("paoflow2wtb.py: only unpolarized runs (nspin = 1, no spin-orbit)")
    lat = np.asarray(arry["a_vectors"]) * attr["alat"] * BOHR_A
    nk = (int(attr["nk1"]), int(attr["nk2"]), int(attr["nk3"]))
    Hk = np.asarray(arry["Hks"])[..., 0]                  # (nawf, nawf, nk1, nk2, nk3)
    nawf = Hk.shape[0]
    orbs = [{"species": b["atom"], "label": b["label"], "l": int(b["l"]), "m": int(b["m"]),
             "tau": np.asarray(b["tau"], dtype=float) * BOHR_A} for b in arry["basis"]]
    tau = np.array([o["tau"] for o in orbs])
    box = np.fft.fftn(Hk, axes=(2, 3, 4)) / float(np.prod(nk))
    R, HR, nimg = ws_fold(box, nk, lat, tau)

    # checks: H(k) on the grid, and the bands at the k-points of the DFT run
    kg = np.array([np.array(i) / np.array(nk) for i in np.ndindex(*nk)])
    Hg = np.transpose(Hk, (2, 3, 4, 0, 1)).reshape(-1, nawf, nawf)
    dH = float(np.abs(hr_to_hk(R, HR, kg) - Hg).max())
    Eg = np.linalg.eigvalsh(Hg)
    kq = np.asarray(arry["kpnts"]) @ np.linalg.inv(np.asarray(arry["b_vectors"]))
    Eq = np.linalg.eigvalsh(hr_to_hk(R, HR, kq))
    Edft = np.asarray(arry["my_eigsmat"])[:, :, 0].T[:, :nawf]
    eta = float(attr["shift"])
    nocc = int(round(float(attr["nelec"]))) // 2
    dE = float(np.abs(Eq[:, :nocc] - Edft[:, :nocc]).max()) if 0 < nocc <= nawf else float("nan")
    gap = None
    if 0 < nocc < nawf and Eg[:, nocc].min() > Eg[:, nocc - 1].max() + 1e-3:
        gap = (float(Eg[:, nocc - 1].max()), float(Eg[:, nocc].min()))
    if a.efermi is not None:
        ef = a.efermi
    else:
        ef = 0.5 * (gap[0] + gap[1]) if gap else 0.0

    try:
        from importlib.metadata import version
        pver = version("PAOFLOW")
    except Exception:
        pver = "?"
    head = ("paoflow2wtb.py: PAOFLOW {}, basis {}, pthr {}, shift {:.4f} eV, {}x{}x{} grid, "
            "minimal images").format(pver, a.configuration, attr["pthr"], eta, *nk)
    write_w90(a.seedname + "-NP.dat", R, HR, lat, ef, head)
    write_rdat(a.seedname + "_r.dat", tau, head)
    write_tbnp("tb-NP.dat", R, HR, lat, ef)
    write_basis_set("basis_set-NP", orbs)
    lines = [head, "save directory: {}".format(os.path.abspath(a.savedir)),
             "orbitals {}, lattice vectors {}, pairs split between images {}".format(
                 nawf, len(R), int((nimg > 1).sum())),
             "bands below the shift (bnd) {}, occupied bands {}".format(int(attr["bnd"]), nocc),
             "max |H(k) - PAOFLOW H(k)| on the grid: {:.2e} eV".format(dH),
             "max |E - E_DFT| of the occupied bands at the {} DFT k-points: {:.4f} eV".format(len(kq), dE),
             ("grid VBM {:.4f} eV, CBM {:.4f} eV".format(*gap) if gap else "no gap on the grid")
             + ", Fermi level in the files {:.6f} eV".format(ef),
             "", "# index species shell l m x y z (A)"]
    lines += ["{:4d} {:3s} {:4s} {} {:2d} {:14.8f} {:14.8f} {:14.8f}".format(
        i + 1, o["species"], o["label"], o["l"], o["m"], *o["tau"]) for i, o in enumerate(orbs)]
    open(a.seedname + "-info.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines[:7]))
    print("wrote {0}-NP.dat {0}_r.dat (DFT= \"W\"), tb-NP.dat basis_set-NP (DFT= \"S\"), "
          "{0}-info.txt".format(a.seedname))


def main():
    ap = argparse.ArgumentParser(description="PAOFLOW tight-binding Hamiltonian -> WanTiBEXOS input")
    ap.add_argument("savedir", help="the pw.x save directory (prefix.save)")
    ap.add_argument("--configuration", default="minimal",
                    help="PAOFLOW basis: minimal (the UPF wavefunctions), standard, extended")
    ap.add_argument("--basispath", default=None, help="PAOFLOW's BASIS directory (standard, extended)")
    ap.add_argument("--pthr", type=float, default=0.95, help="projectability threshold")
    ap.add_argument("--shift", default="auto", help="eta in eV, or auto")
    ap.add_argument("--shift-type", type=int, default=1, help="PAOFLOW's shift_type")
    ap.add_argument("--efermi", type=float, default=None, help="Fermi level written in the files (eV)")
    ap.add_argument("--seedname", default="paoflow", help="name of the DFT= \"W\" files")
    ap.add_argument("--npool", type=int, default=1)
    ap.add_argument("--formfactor", action="store_true",
                    help="also write <seedname>_ff.bin, the form factors of the orbitals")
    ap.add_argument("--mesh", type=int, nargs=3, default=None, help="the k mesh of the BSE (NGX NGY NGZ)")
    ap.add_argument("--ff-ecut", type=float, default=100.0, help="exchange G up to this energy (eV)")
    ap.add_argument("--ff-spacing", type=float, default=None,
                    help="grid spacing of the integrals (A); default: the grid of the cutoff")
    ap.add_argument("--ff-tail", type=float, default=1e-6,
                    help="norm left outside the radius of each Loewdin orbital (the box-state bases "
                         "standard and extended need a larger one: see README)")
    ap.add_argument("--yes", action="store_true", help="write the form factors even above 1 GB")
    a = ap.parse_args()
    if a.formfactor and a.mesh is None:
        ap.error("--formfactor needs --mesh NGX NGY NGZ, the k-point mesh of the BSE")

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    try:
        from PAOFLOW import PAOFLOW
        pf = PAOFLOW.PAOFLOW(savedir=a.savedir, outputdir="paoflow2wtb.out", npool=a.npool,
                             smearing=None, verbose=False)
        ff = None
        if a.formfactor:
            ok = 0
            if comm.Get_rank() == 0:
                ff = formfactor_setup(pf, a)
                ok = 1
            if not comm.bcast(ok, root=0):
                raise SystemExit(1)
        pf.projections(configuration=a.configuration,
                       basispath=None if a.basispath is None else os.path.join(a.basispath, ""))
        pf.projectability(pthr=a.pthr, shift=a.shift if a.shift == "auto" else float(a.shift))
        pf.pao_hamiltonian(shift_type=a.shift_type)
        if comm.Get_rank() == 0:
            convert(pf, a)
            if ff is not None:
                print(formfactors(pf, a, ff), flush=True)
        comm.Barrier()
    except BaseException as e:
        # an exception on one rank would leave the others waiting forever
        if comm.Get_size() > 1:
            sys.stderr.write("paoflow2wtb.py: {}\n".format(e))
            comm.Abort(1)
        raise
    return 0


if __name__ == "__main__":
    sys.exit(main())
