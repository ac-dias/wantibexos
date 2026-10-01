"""
wtb_formfactor.py: orbital form factors of a tight-binding basis for
WanTiBEXOS, shared by siesta2wtb.py, paoflow2wtb.py and wannier2wtb.py
(option --formfactor).

The form factor of orbitals m and n is

    F_mn(R; Q) = <m,0| exp(iQ.r) |n,R> = int phi_m(r)^* exp(iQ.r) phi_n(r - R) d^3r,

phi_n(r - R) the orbital n of the cell at the lattice vector R and r measured
from the origin of the lattice (the frame of the orbital centres). With it the
Bloch states |c k> = sum_m c_m(k) |k,m> of the tight-binding model have the
exact pair densities

    <c k| exp(iQ.r) |c' k'> = sum_mn c_m(k)^* c'_n(k') sum_R exp(ik'.R) F_mn(R; Q),

Q = k - k' + G, instead of the point-centre form sum_m c_m^* c'_m exp(iQ.t_m)
(F_mn(R; Q) -> delta_mn delta_R0 exp(iQ.t_m) for point charges at t_m;
F_mn(R; 0) is the overlap S_mn(R)). Two sets of Q are written:

  * direct term: every mesh vector q = (i/N1, j/N2, l/N3) of the BSE mesh
    (--mesh, NGX NGY NGZ) as its shortest images q - G, all of them on the zone
    boundary, as bse_q_images of WanTiBEXOS chooses them (ties within 1e-4);
  * exchange term (local fields): the reciprocal lattice vectors G != 0 with
    hbar^2 G^2/2m <= --ff-ecut (eV).

The integrals are sums over a real-space grid aligned with the lattice vectors
(M_i points along a_i), where exp(iQ.r) separates in the crystal coordinates
of r and Q: three successive one-dimensional sums per pair of orbital groups.
Orbitals that vanish beyond their radius (SIESTA's) are integrated over the
overlap of their spheres; orbitals with tails (PAOFLOW's Loewdin orbitals, the
Wannier functions) over the box around both spheres, so that the tail of each
is kept where the other is large. Lattice vectors R are kept when some pair
of orbitals is closer than the sum of their radii; the file size is known
from that before anything is computed, and above 1 GB the scripts ask before
going on (or stop when they cannot ask), unless --yes is given.

File <name>_ff.bin, little-endian stream (Fortran access='stream'):

    character(8)  'WTBFF001'
    integer(4)    norb, nR, nQ, ndir, nexc, mesh(3)
    real(8)       ecut (eV)
    real(8)       lat(3,3)        lat(:,i) = a_i (A)
    real(8)       tau(3,norb)     orbital centres (A)
    integer(4)    R(3,nR)         lattice vectors, units of a_i
    real(8)       Qc(3,nQ)        Q in units of b_i: the ndir direct ones, then
                                  the nexc exchange ones
    complex(4)    F(norb,norb,nR,nQ)   F(m,n,iR,iQ) = <m,0|exp(iQ.r)|n,R>
"""
import sys

import numpy as np

MAGIC = b"WTBFF001"
LIMIT = 1.0e9                     # bytes: above it the scripts ask before writing
HB2M2 = 3.80998212                # hbar^2/2m_e, eV A^2
TOL_IMG = 1.0e-4                  # bse_q_images


def recip(lat):
    """rows b_i with a_i.b_j = 2 pi delta_ij"""
    return 2.0 * np.pi * np.linalg.inv(np.asarray(lat, dtype=float)).T


def direct_q(lat, mesh):
    """crystal coordinates of the shortest images q - G of every vector q of
    the mesh, all images of equal length (the choice of bse_q_images)"""
    B = recip(lat)
    mesh = np.asarray(mesh, dtype=int)
    cand = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], float)
    out = []
    for idx in np.ndindex(*mesh):
        d = np.array(idx, float) / mesh
        g = np.rint(d) + cand
        qc = d[None, :] - g
        ql = np.linalg.norm(qc @ B, axis=1)
        out.extend(qc[ql <= ql.min() * (1.0 + TOL_IMG) + 1.0e-6])
    return np.array(out)


def exchange_g(lat, ecut):
    """crystal coordinates (integers) of the G != 0 with hbar^2 G^2/2m <= ecut,
    by length"""
    if ecut <= 0:
        return np.zeros((0, 3))
    B = recip(lat)
    gmax = np.sqrt(ecut / HB2M2)
    n = [int(np.ceil(gmax / np.linalg.norm(b))) + 1 for b in B]
    h = np.array(np.meshgrid(*[np.arange(-m, m + 1) for m in n], indexing="ij")).reshape(3, -1).T
    g2 = np.sum((h @ B) ** 2, axis=1)
    keep = (g2 > 1e-12) & (HB2M2 * g2 <= ecut * (1 + 1e-12))
    h, g2 = h[keep], g2[keep]
    return h[np.lexsort((h[:, 2], h[:, 1], h[:, 0], np.round(g2, 9)))].astype(float)


def lattice_vectors(lat, tau, radius, periodic=(True, True, True)):
    """R (units of a_i) for which some pair of orbitals m (cell 0), n (cell R)
    is closer than radius_m + radius_n, by length"""
    lat = np.asarray(lat, float)
    tau = np.asarray(tau, float)
    radius = np.asarray(radius, float)
    reach = 2.0 * radius.max() + (np.linalg.norm(np.ptp(tau, axis=0)) if len(tau) > 1 else 0.0)
    # a box of lattice vectors large enough for that reach
    recl = np.linalg.norm(recip(lat), axis=1) / (2.0 * np.pi)       # 1/(plane spacing)
    nmax = [int(np.ceil(reach * recl[i])) + 1 if periodic[i] else 0 for i in range(3)]
    R = np.array(np.meshgrid(*[np.arange(-m, m + 1) for m in nmax], indexing="ij")).reshape(3, -1).T
    # per species-free test: the shortest orbital distance for each R
    d = (R @ lat)[:, None, None, :] + tau[None, None, :, :] - tau[None, :, None, :]
    ok = (np.linalg.norm(d, axis=-1) < radius[None, :, None] + radius[None, None, :]).any(axis=(1, 2))
    R = R[ok]
    return R[np.lexsort((R[:, 2], R[:, 1], R[:, 0], np.abs(R).sum(axis=1)))]


def header_bytes(norb, nR, nQ):
    return 8 + 8 * 4 + 8 + 9 * 8 + 3 * 8 * norb + 3 * 4 * nR + 3 * 8 * nQ


def file_size(norb, nR, nQ):
    return header_bytes(norb, nR, nQ) + 8 * norb * norb * nR * nQ


def confirm_size(path, norb, nR, nQ, yes, ndir=None, nexc=None):
    """print the size the form-factor file will have; above LIMIT ask (or stop)"""
    size = file_size(norb, nR, nQ)
    parts = "" if ndir is None else " ({} direct, {} exchange)".format(ndir, nexc)
    print("form factors: {} will take {:.3f} GB: {} orbitals, {} lattice vectors, {} Q vectors{}".format(
        path, size / 1e9, norb, nR, nQ, parts), flush=True)
    if size <= LIMIT or yes:
        return True
    msg = "{} would take {:.2f} GB, more than 1 GB".format(path, size / 1e9)
    if sys.stdin is not None and sys.stdin.isatty():
        ans = input(msg + ". Write it? [y/N] ")
        if ans.strip().lower() in ("y", "yes"):
            return True
        raise SystemExit("stopped before computing the form factors")
    raise SystemExit(msg + "; run again with --yes to write it anyway (or reduce --mesh or --ff-ecut)")


class FFWriter:
    """the form-factor file: header, then F on a memory map filled pair by pair"""

    def __init__(self, path, lat, tau, R, Qc, ndir, mesh, ecut):
        norb, nR, nQ = len(tau), len(R), len(Qc)
        self.shape = (nQ, nR, norb, norb)            # C order: m fastest = F(m,n,iR,iQ)
        with open(path, "wb") as f:
            f.write(MAGIC)
            np.array([norb, nR, nQ, ndir, nQ - ndir, *mesh], dtype="<i4").tofile(f)
            np.array([ecut], dtype="<f8").tofile(f)
            np.asarray(lat, dtype="<f8").tofile(f)
            np.asarray(tau, dtype="<f8").tofile(f)
            np.asarray(R, dtype="<i4").tofile(f)
            np.asarray(Qc, dtype="<f8").tofile(f)
            self.offset = f.tell()
            f.truncate(self.offset + 8 * int(np.prod(self.shape)))
        self.F = np.memmap(path, dtype="<c8", mode="r+", offset=self.offset, shape=self.shape)

    def close(self):
        self.F.flush()
        del self.F


def read_ff(path):
    """the arrays of a form-factor file; F as a read-only memory map (nQ, nR, n, m)"""
    with open(path, "rb") as f:
        if f.read(8) != MAGIC:
            raise ValueError("{} is not a WanTiBEXOS form-factor file".format(path))
        norb, nR, nQ, ndir, nexc, m1, m2, m3 = np.fromfile(f, dtype="<i4", count=8)
        ecut = float(np.fromfile(f, dtype="<f8", count=1)[0])
        lat = np.fromfile(f, dtype="<f8", count=9).reshape(3, 3)
        tau = np.fromfile(f, dtype="<f8", count=3 * norb).reshape(norb, 3)
        R = np.fromfile(f, dtype="<i4", count=3 * nR).reshape(nR, 3)
        Qc = np.fromfile(f, dtype="<f8", count=3 * nQ).reshape(nQ, 3)
        off = f.tell()
    F = np.memmap(path, dtype="<c8", mode="r", offset=off, shape=(nQ, nR, norb, norb))
    return {"lat": lat, "tau": tau, "R": R, "Qc": Qc, "ndir": int(ndir), "nexc": int(nexc),
            "mesh": (int(m1), int(m2), int(m3)), "ecut": ecut, "F": F}


class Grid:
    """real-space grid of M_i points along each lattice vector a_i"""

    def __init__(self, lat, M):
        self.lat = np.asarray(lat, float)
        self.M = np.asarray(M, dtype=int)
        self.dV = abs(np.linalg.det(self.lat)) / float(np.prod(self.M))

    @classmethod
    def with_spacing(cls, lat, spacing, minimum=(1, 1, 1)):
        lat = np.asarray(lat, float)
        M = [max(int(np.ceil(np.linalg.norm(a) / spacing)), m) for a, m in zip(lat, minimum)]
        return cls(lat, M)

    def box(self, centre, radius, maxwidth=None):
        """the index box (lo, shape) of grid points of the sphere around centre;
        maxwidth: at most that many points along a_i (one period of an
        orbital that is periodic along a_i), centred on the centre"""
        B = recip(self.lat) / (2.0 * np.pi)                   # frac = r @ B.T
        fc = np.asarray(centre, float) @ B.T
        half = radius * np.linalg.norm(B, axis=1)             # half-width in fractions of a_i
        lo = np.floor((fc - half) * self.M).astype(int)
        hi = np.ceil((fc + half) * self.M).astype(int) + 1
        if maxwidth is not None:
            for i in range(3):
                if maxwidth[i] and hi[i] - lo[i] > maxwidth[i]:
                    lo[i] = int(np.round(fc[i] * self.M[i])) - int(maxwidth[i]) // 2
                    hi[i] = lo[i] + int(maxwidth[i])
        return lo, hi - lo

    def points(self, lo, shape):
        """Cartesian points of the box, shape (L1, L2, L3, 3)"""
        f = [(lo[i] + np.arange(shape[i])) / self.M[i] for i in range(3)]
        F = np.stack(np.meshgrid(*f, indexing="ij"), axis=-1)
        return F @ self.lat


class Transform:
    """sum_x rho(x) exp(2 pi i Qc.x) dV on the grid, for one set of Q, by one-
    dimensional sums in x3, x2 and x1 (x = n/M, crystal coordinates)"""

    def __init__(self, grid, Qc):
        self.g = grid
        self.Qc = np.asarray(Qc, float)
        self.u3, self.i3 = np.unique(np.round(self.Qc[:, 2], 12), return_inverse=True)
        self.u2, self.i2 = np.unique(np.round(self.Qc[:, 1], 12), return_inverse=True)

    def __call__(self, rho, lo):
        """rho (P, L1, L2, L3) on the box at grid index lo -> (P, nQ)"""
        if len(self.Qc) == 0:
            return np.zeros((rho.shape[0], 0), complex)
        M = self.g.M
        L = rho.shape[1:]
        x = [(lo[i] + np.arange(L[i])) / M[i] for i in range(3)]
        E3 = np.exp(2j * np.pi * np.outer(x[2], self.u3))
        T = np.tensordot(rho, E3, axes=([3], [0]))                       # (P, L1, L2, n3)
        E2 = np.exp(2j * np.pi * np.outer(x[1], self.u2))
        V = np.einsum("pabc,bd->padc", T, E2, optimize=True)             # (P, L1, n2, n3)
        E1 = np.exp(2j * np.pi * np.outer(x[0], self.Qc[:, 0]))          # (L1, nQ)
        return np.einsum("paq,aq->pq", V[:, :, self.i2, self.i3], E1, optimize=True) * self.g.dV


class Group:
    """orbitals with one centre and one radius, in one of three forms:

      confined  zero outside their box: values (n, L1, L2, L3) at grid index lo
      periodic  full (n, D1, D2, D3), one whole period D of the grid (a
                supercell): the orbitals of a k grid of that size
      window    values (n, L1, L2, L3) at lo, known only there (a plot), except
                along the axes flagged in periodic, where the window is one
                period of the orbitals

    Periodic and window orbitals have tails beyond the radius: a pair of them
    is integrated over the box around both spheres (see compute)."""

    def __init__(self, index, centre, radius, lo=None, values=None, full=None, periodic=None):
        self.index = list(index)
        self.centre = np.asarray(centre, float)
        self.radius = float(radius)
        self.lo = None if lo is None else np.asarray(lo, dtype=int)
        self.values = values
        self.full = full
        self.periodic = None if periodic is None else [bool(p) for p in periodic]
        if full is not None:
            self.period = np.array(full.shape[1:])
        elif periodic is not None:
            self.period = np.array([values.shape[1 + i] if periodic[i] else 0 for i in range(3)])
        else:
            self.period = None

    @property
    def extended(self):
        return self.full is not None or self.periodic is not None

    def on(self, lo, shape):
        """the values on the box (lo, shape): periodic, or zero outside a window"""
        if self.full is not None:
            ix = [(lo[i] + np.arange(shape[i])) % self.period[i] for i in range(3)]
            return self.full[:, ix[0][:, None, None], ix[1][None, :, None], ix[2][None, None, :]]
        L = self.values.shape[1:]
        ix, ok = [], []
        for i in range(3):
            j = lo[i] - self.lo[i] + np.arange(shape[i])
            if self.periodic[i]:
                j %= L[i]
            ok.append((j >= 0) & (j < L[i]))
            ix.append(np.clip(j, 0, L[i] - 1))
        out = self.values[:, ix[0][:, None, None], ix[1][None, :, None], ix[2][None, None, :]]
        return out * (ok[0][:, None, None] & ok[1][None, :, None] & ok[2][None, None, :])


def _pair_region(grid, A, B, cb):
    """the box around both spheres (A at its centre, B at cb) for periodic
    groups, at most one period D along each axis (centred between them)"""
    loA, shA = grid.box(A.centre, A.radius)
    loB, shB = grid.box(cb, B.radius)
    lo = np.minimum(loA, loB)
    hi = np.maximum(loA + shA, loB + shB)
    for i in range(3):
        D = [p[i] for p in (A.period, B.period) if p[i] > 0]
        D = min(D) if D else 0
        if D and hi[i] - lo[i] > D:
            mid = (lo[i] + hi[i]) // 2
            lo[i], hi[i] = mid - D // 2, mid - D // 2 + D
    return lo, hi - lo


def compute(writer, grid, groups, R, Qsets, offsets=(0,)):
    """F for every pair of groups (m in A at cell 0, n in B at cell R) into the
    writer's map; Qsets: [(Qc block, first column)] covering the file's Q list;
    offsets: the same F is written at orbital indices + each offset (the spin
    blocks of a spin-diagonal basis). Confined groups are integrated over the
    intersection of their boxes, periodic ones over the box around both
    spheres, so that each orbital's tail is kept where the other is large"""
    lat = grid.lat
    trans = [(Transform(grid, q), c0) for q, c0 in Qsets]
    Rc = np.asarray(R) @ lat
    npair = 0
    for iR, (Ri, rc) in enumerate(zip(R, Rc)):
        shift = np.asarray(Ri) * grid.M
        for A in groups:
            for Bg in groups:
                if np.linalg.norm(Bg.centre + rc - A.centre) >= A.radius + Bg.radius:
                    continue
                if A.extended and Bg.extended:
                    lo, shape = _pair_region(grid, A, Bg, Bg.centre + rc)
                    wA, wB = A.on(lo, shape), Bg.on(lo - shift, shape)
                else:
                    loB = Bg.lo + shift
                    lo = np.maximum(A.lo, loB)
                    hi = np.minimum(A.lo + np.array(A.values.shape[1:]), loB + np.array(Bg.values.shape[1:]))
                    if np.any(hi <= lo):
                        continue
                    sa = tuple(slice(lo[i] - A.lo[i], hi[i] - A.lo[i]) for i in range(3))
                    sb = tuple(slice(lo[i] - loB[i], hi[i] - loB[i]) for i in range(3))
                    wA, wB = A.values[(slice(None),) + sa], Bg.values[(slice(None),) + sb]
                for ia, m in enumerate(A.index):
                    rho = np.conj(wA[ia])[None] * wB                                  # (nB, box)
                    for tr, c0 in trans:
                        f = tr(rho, lo).T.astype(np.complex64)                        # (nq, nB)
                        for off in offsets:
                            writer.F[c0:c0 + f.shape[0], iR, [j + off for j in Bg.index], m + off] = f
                npair += 1
    return npair


def bloch(ff, kc, iq):
    """<k+Q,m| exp(iQ.r) |k,n> = sum_R exp(ik.R) F_mn(R; Q) for the Q of column
    iq at crystal k (the ket's k): (norb, norb), for checks"""
    ph = np.exp(2j * np.pi * (ff["R"] @ np.asarray(kc, float)))
    return np.einsum("r,rnm->mn", ph, np.asarray(ff["F"][iq], dtype=complex))
