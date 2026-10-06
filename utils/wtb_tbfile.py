#!/usr/bin/env python3
"""The tight-binding file of WanTiBEXOS for DFT= "S" (H(R), S(R) of a non-orthogonal basis): the binary layout that
siesta2wtb.py writes by default and wtb.x reads by default (PARAMS_FORMAT= "binary"), a reader of both layouts, and a
converter between them.

  python3 wtb_tbfile.py info    tb-soc.bin|tb-soc.dat
  python3 wtb_tbfile.py convert tb-soc.dat tb-soc.bin     text -> binary  (the layout of the input file is detected)
  python3 wtb_tbfile.py convert tb-soc.bin tb-soc.dat     binary -> text

  import wtb_tbfile; tb = wtb_tbfile.read_tb("tb-soc.bin")      # either layout; the dictionary is described at read_tb

The text layout is the one wtb.x has always read (nine header lines, then one line "x y z i j ReH ImH S" per element: x y z
the Cartesian translation of the image in Angstrom, i the row, j the column, i running fastest). The binary layout holds the
same numbers in single precision, which is what wtb.x keeps in memory (the HSX file of SIESTA is single precision too, so
the conversion from siesta2wtb.py loses nothing), and only the non-zero elements (7 % of the rows of the 4 x 4 MoSi2N4
supercell, 21 % of the unit cell): 20 bytes each instead of ~130 bytes of text, and no parsing.

Binary layout, little endian, Fortran stream access (every number is read with one `read(unit) array`):

  offset  type        content
       0  char(8)     'WTBTB001'
       8  int32       version (1)
      12  int32       N, the number of basis functions (2 x orbitals for the spinor models)
      16  int32       nvec, the number of lattice vectors (images)
      20  float32     scissor (eV), line 2 of the text layout
      24  float32     Fermi level (eV), line 3
      28  char(4)     calculation type of line 1: 'NP', 'SP' or 'SOC', blank padded
      32  float32(3,3) the lattice vectors, three components each (Angstrom); read(unit) rtmp(3,3), rlat = transpose(rtmp)
      68  float32(nvec,3)  translation of each image, Cartesian, Angstrom; all x, then all y, then all z
   then, for each lattice vector in turn:
          int64       nnz, the number of stored elements of this image
          int32(nnz)  row index i (1-based, the first index of the text lines)
          int32(nnz)  column index j
          float32(nnz) Re H (eV)
          float32(nnz) Im H (eV)
          float32(nnz) S
  Elements that are not stored are zero (H, S) of the Hamiltonian block H_ij(R) = <i,0|H|j,R>.
"""
import struct
import sys

import numpy as np

MAGIC = b'WTBTB001'
VERSION = 1
HEADER = 68


def is_binary(fname):
    with open(fname, 'rb') as f:
        return f.read(8) == MAGIC


def _systype(text):
    return text.strip()[:4].ljust(4).encode('ascii')


class BinaryWriter(object):
    """f = BinaryWriter(fname, systype, scissor, fermi, lattice, rvec, nbasis); f.add(i, j, reH, imH, S) per lattice
    vector in the order of rvec (1-based index arrays and the values of the stored elements); f.close()"""

    def __init__(self, fname, systype, scissor, fermi, lattice, rvec, nbasis):
        rvec = np.asarray(rvec, dtype=np.float64).reshape(-1, 3)
        self.nvec, self.nbasis, self.added, self.stored = len(rvec), int(nbasis), 0, 0
        self.f = open(fname, 'wb')
        self.f.write(MAGIC + struct.pack('<iiiff', VERSION, self.nbasis, self.nvec, scissor, fermi) + _systype(systype))
        self.f.write(np.asarray(lattice, dtype='<f4').reshape(3, 3).tobytes())
        self.f.write(np.ascontiguousarray(rvec.T, dtype='<f4').tobytes())

    def add(self, i, j, reH, imH, S):
        n = len(i)
        self.f.write(struct.pack('<q', n))
        self.f.write(np.asarray(i, dtype='<i4').tobytes())
        self.f.write(np.asarray(j, dtype='<i4').tobytes())
        for a in (reH, imH, S):
            self.f.write(np.asarray(a, dtype='<f4').tobytes())
        self.added += 1
        self.stored += n

    def add_block(self, reH, imH, S):
        """one lattice vector as dense (N, N) arrays [row, column]: the elements where any of the three is non-zero"""
        reH, imH, S = (np.asarray(a, dtype=np.float32) for a in (reH, imH, S))
        col, row = np.nonzero(((reH != 0) | (imH != 0) | (S != 0)).T)         # sorted by column, then row
        self.add(row + 1, col + 1, reH[row, col], imH[row, col], S[row, col])

    def close(self):
        self.f.close()
        if self.added != self.nvec:
            raise RuntimeError('%d lattice vectors written, the header says %d' % (self.added, self.nvec))


def _read_binary_header(f):
    head = f.read(HEADER)
    if head[:8] != MAGIC:
        raise ValueError('not a WanTiBEXOS binary tight-binding file (no WTBTB001 at the start)')
    version, nbasis, nvec, scissor, fermi = struct.unpack('<iiiff', head[8:28])
    if version != VERSION:
        raise ValueError('binary tight-binding file of version %d, this reader knows version %d' % (version, VERSION))
    lattice = np.frombuffer(head[32:68], dtype='<f4').reshape(3, 3).astype(np.float64)
    rvec = np.frombuffer(f.read(12 * nvec), dtype='<f4').reshape(3, nvec).T.astype(np.float64)
    return {'systype': head[28:32].decode('ascii').strip(), 'scissor': scissor, 'fermi': fermi, 'lattice': lattice,
            'nbasis': nbasis, 'nvec': nvec, 'rvec': rvec}


def iter_binary(fname):
    """(header dictionary, generator of (i, j, reH, imH, S) per lattice vector: 1-based indices, float32 values)"""
    f = open(fname, 'rb')
    hd = _read_binary_header(f)

    def images():
        for _ in range(hd['nvec']):
            n = struct.unpack('<q', f.read(8))[0]
            i = np.frombuffer(f.read(4 * n), dtype='<i4')
            j = np.frombuffer(f.read(4 * n), dtype='<i4')
            v = [np.frombuffer(f.read(4 * n), dtype='<f4') for _ in range(3)]
            yield i, j, v[0], v[1], v[2]
        f.close()
    return hd, images()


def _read_text_header(f):
    lines = [f.readline() for _ in range(9)]
    systype = lines[0].split()[0]
    return {'systype': systype, 'scissor': float(lines[1].split()[0]), 'fermi': float(lines[2].split()[0]),
            'lattice': np.array([[float(v) for v in lines[3 + i].split()[:3]] for i in range(3)]),
            'nbasis': int(lines[6].split()[0]), 'nvec': int(lines[7].split()[0])}


def iter_text(fname, chunk=4000000):
    """(header dictionary, generator of (i, j, reH, imH, S, rvec) per lattice vector of the text layout)"""
    f = open(fname)
    hd = _read_text_header(f)
    n2 = hd['nbasis'] ** 2

    def images():
        for _ in range(hd['nvec']):
            parts = []
            left = n2
            while left:
                lines = [f.readline() for _ in range(min(chunk, left))]
                if not lines[-1]:
                    raise ValueError('%s: the file ends inside a lattice vector' % fname)
                parts.append(np.loadtxt(lines, ndmin=2))
                left -= len(lines)
            a = np.vstack(parts) if len(parts) > 1 else parts[0]
            yield a[:, 3].astype(np.int64), a[:, 4].astype(np.int64), a[:, 5], a[:, 6], a[:, 7], a[-1, :3]
        f.close()
    return hd, images()


def read_tb(fname):
    """The tight-binding file in either layout, as a dictionary: systype, scissor, fermi, lattice (rows are the lattice
    vectors, Angstrom), nbasis, nvec, rvec (nvec, 3: Cartesian translation of each image, Angstrom), H (nvec, N, N)
    complex, S (nvec, N, N), indexed [image, row, column] as the matrices H_ij(R) = <i,0|H|j,R> of the text layout."""
    if is_binary(fname):
        hd, images = iter_binary(fname)
        pairs = ((x + (None,)) for x in images)
    else:
        hd, images = iter_text(fname)
        pairs = images
    n, nvec = hd['nbasis'], hd['nvec']
    H = np.zeros((nvec, n, n), dtype=np.complex128)
    S = np.zeros((nvec, n, n))
    rvec = np.zeros((nvec, 3))
    for iR, (i, j, re, im, s, r) in enumerate(pairs):
        H[iR, i - 1, j - 1] = re.astype(np.float64) + 1j * im
        S[iR, i - 1, j - 1] = s
        if r is not None:
            rvec[iR] = r
    if 'rvec' in hd:
        rvec = hd.pop('rvec')
    hd.update(H=H, S=S, rvec=rvec)
    return hd


def text_to_binary(src, dst):
    hd, images = iter_text(src)
    rv = []
    chunks = []
    # the header of the binary file needs every translation before the first image: keep the images' data (sparse)
    for i, j, re, im, s, r in images:
        keep = (re != 0) | (im != 0) | (s != 0)
        chunks.append((i[keep], j[keep], re[keep], im[keep], s[keep]))
        rv.append(r)
    w = BinaryWriter(dst, hd['systype'], hd['scissor'], hd['fermi'], hd['lattice'], np.array(rv), hd['nbasis'])
    for c in chunks:
        w.add(*c)
    w.close()
    return w.stored, hd['nvec'] * hd['nbasis'] ** 2


def binary_to_text(src, dst):
    hd, images = iter_binary(src)
    n = hd['nbasis']
    with open(dst, 'w') as f:
        f.write('%s\n%.12e\n%s\n' % (hd['systype'], hd['scissor'], hd['fermi']))
        for row in hd['lattice']:
            f.write('%.12e  %.12e  %.12e\n' % tuple(row))
        f.write('%d\n%d\n#rcell x   rcell y   rcell z   i   j   ReH   ImH   S\n' % (n, hd['nvec']))
        for iR, (i, j, re, im, s) in enumerate(images):
            blk = np.zeros((3, n, n))                                  # [component, row, column]
            blk[0][i - 1, j - 1], blk[1][i - 1, j - 1], blk[2][i - 1, j - 1] = re, im, s
            row, col = np.indices((n, n))                              # text order: column outer, row (the first index) inner
            table = np.column_stack([np.tile(hd['rvec'][iR], (n * n, 1)), (row.T.ravel() + 1), (col.T.ravel() + 1),
                                     blk[0].T.ravel(), blk[1].T.ravel(), blk[2].T.ravel()])
            np.savetxt(f, table, fmt=['%.12e'] * 3 + ['%d'] * 2 + ['%.12e'] * 3, delimiter='   ')


def main(argv):
    if len(argv) >= 2 and argv[0] == 'info':
        if is_binary(argv[1]):
            hd, images = iter_binary(argv[1])
            nnz = sum(len(x[0]) for x in images)
            print('%s: binary, %s, N = %d, %d lattice vectors, %d stored elements of %d (%.1f %%)' % (
                argv[1], hd['systype'], hd['nbasis'], hd['nvec'], nnz, hd['nvec'] * hd['nbasis'] ** 2,
                100.0 * nnz / (hd['nvec'] * hd['nbasis'] ** 2)))
        else:
            hd, _ = iter_text(argv[1])
            print('%s: text, %s, N = %d, %d lattice vectors' % (argv[1], hd['systype'], hd['nbasis'], hd['nvec']))
        return 0
    if len(argv) == 3 and argv[0] == 'convert':
        if is_binary(argv[1]):
            binary_to_text(argv[1], argv[2])
            print('%s: text written' % argv[2])
        else:
            stored, total = text_to_binary(argv[1], argv[2])
            print('%s: binary written, %d stored elements of %d (%.1f %%)' % (argv[2], stored, total, 100.0 * stored / total))
        return 0
    print(__doc__)
    return 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
