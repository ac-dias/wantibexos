#!/usr/bin/env python3
import sisl
import numpy as np
#import matplotlib.pyplot as plt
import os
import sys

############################################################################
def calctype(data):
    
    if data == 'unpolarized':
     outdata = 'NP'
        
    elif data == 'polarized':
     outdata = 'SP'
        
    elif data == 'non-colinear':
     outdata = 'SOC'
        
    elif data == 'spin-orbit':
     outdata = 'SOC'
            
    return outdata   


def ncases(var):
    # 13 significant digits: S(k) of an LCAO basis is close to singular, and
    # rounding H and S to 6 decimals grows with energy in the eigenvalues
    return "%.12e"%var

###############################################################################
# --rmatrix: position matrix of the SIESTA basis
#
# X_ij(R) = <phi_i,0| r |phi_j,R> (Angstrom, absolute coordinates, the frame
# of basis_set-*), with phi_i,0 orbital i of the home cell and phi_j,R
# orbital j of the cell R, written as tb-*_r.dat in the layout of Wannier90's
# seedname_r.dat (which WanTiBEXOS' rmn_read parses). It is the matrix of
# the non-orthogonal SIESTA basis, not of orthonormal Wannier functions, so
# WanTiBEXOS reads it as BSE_CENTER_FILE with DFT=S for the part of the
# optical dipole that H(R), S(R) and the orbital centres miss, the dipoles
# <phi_i|r - (r_i + r_j)/2|phi_j> between basis orbitals (the kernel takes
# the orbital centres of basis_set-*).
# It needs the basis functions, here the ones sisl reads from the
# *.ion.nc/*.ion.xml files next to the fdf.
#
# Each atom pair is integrated over the region where both of its orbitals
# can be non-zero: two centres A, B in bipolar coordinates (r_a, r_b, phi)
# about the A-B axis, volume element r_a r_b/d dr_a dr_b dphi; one centre in
# spherical coordinates. The phi (and theta) integrands are trigonometric
# polynomials of low degree, integrated exactly by the uniform (Gauss) rules;
# r_a and r_b use Gauss-Legendre on pieces split at every orbital cutoff (and
# wherever the r_b range crosses one), so each piece is smooth. The same
# quadrature gives the overlap, which is checked against SIESTA's S.

def _gauss_pieces(lo, hi, cuts, n):
    """Gauss-Legendre nodes and weights on [lo, hi], split at the cuts inside"""
    edges = np.unique([lo, hi] + [c for c in cuts if lo < c < hi])
    x, w = np.polynomial.legendre.leggauss(n)
    xs, ws = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a > 1e-10:
            xs.append(0.5*(b-a)*x + 0.5*(b+a))
            ws.append(0.5*(b-a)*w)
    if not xs:
        return np.empty(0), np.empty(0)
    return np.concatenate(xs), np.concatenate(ws)


def _pair_quadrature(A, B, cut_a, cut_b, n=24, nphi=12):
    """Points and weights for integrals over the region where an orbital
    on A (cutoff radii cut_a) and one on B (cut_b) are both non-zero;
    None when they do not overlap. One centre: 2n radial nodes per piece
    (cheap, and the slowest to converge)."""
    Ra, Rb = max(cut_a), max(cut_b)
    phi = 2*np.pi*np.arange(nphi)/nphi
    d = np.linalg.norm(B - A)
    if d < 1e-8:
        r, wr = _gauss_pieces(0.0, min(Ra, Rb), list(cut_a) + list(cut_b), 2*n)
        ct, wt = np.polynomial.legendre.leggauss(8)
        st = np.sqrt(1.0 - ct**2)
        pts = A + r[:, None, None, None]*np.stack(np.broadcast_arrays(
            st[:, None]*np.cos(phi), st[:, None]*np.sin(phi),
            ct[:, None]*np.ones(nphi)), -1)[None]
        w = (wr*r*r)[:, None, None]*wt[None, :, None]*(2*np.pi/nphi)
        w = np.broadcast_to(w, pts.shape[:3])
        return pts.reshape(-1, 3), w.reshape(-1)
    if d >= Ra + Rb:
        return None
    ez = (B - A)/d
    t = np.array([1.0, 0.0, 0.0]) if abs(ez[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    ex = t - ez*(t @ ez)
    ex /= np.linalg.norm(ex)
    ey = np.cross(ez, ex)
    circle = np.cos(phi)[:, None]*ex + np.sin(phi)[:, None]*ey
    cuts = list(cut_a) + [d] + [c + d for c in cut_b] + [c - d for c in cut_b] \
        + [d - c for c in cut_b]
    ra, wa = _gauss_pieces(0.0, Ra, cuts, n)
    P, W = [], []
    for x, wx in zip(ra, wa):
        lo, hi = abs(x - d), min(x + d, Rb)
        if hi - lo < 1e-10:
            continue
        rb, wb = _gauss_pieces(lo, hi, cut_b, n)
        z = (x*x - rb*rb + d*d)/(2*d)
        rho = np.sqrt(np.clip(x*x - z*z, 0.0, None))
        P.append((A + z[:, None, None]*ez + rho[:, None, None]*circle).reshape(-1, 3))
        W.append(np.repeat(wx*x*rb*wb/d*(2*np.pi/nphi), nphi))
    return np.concatenate(P), np.concatenate(W)


def position_matrix(geom, n=24):
    """X[isc, i, j, :] = <phi_i,0| r |phi_j,R_isc> and the overlap Sq from
    the same quadrature, for every supercell isc of geom (sisl order)."""
    no, ncell = geom.no, geom.n_s
    X = np.zeros((ncell, no, no, 3))
    Sq = np.zeros((ncell, no, no))
    try:
        for atom in geom.atoms.atom:
            for o in atom.orbitals:
                o.psi(np.zeros((1, 3)))
    except Exception:
        sys.exit('--rmatrix needs the basis functions: run with the '
                 '*.ion.nc or *.ion.xml files SIESTA wrote next to the fdf')
    for isc in range(ncell):
        R = geom.lattice.sc_off[isc] @ geom.cell
        for a, atoma in enumerate(geom.atoms):
            A = geom.xyz[a]
            oa = geom.a2o(a, all=True)
            for b, atomb in enumerate(geom.atoms):
                B = geom.xyz[b] + R
                q = _pair_quadrature(A, B, [o.R for o in atoma.orbitals],
                                     [o.R for o in atomb.orbitals], n)
                if q is None:
                    continue
                pts, w = q
                fa = np.array([o.psi(pts - A) for o in atoma.orbitals])*w
                fb = np.array([o.psi(pts - B) for o in atomb.orbitals])
                ob = geom.a2o(b, all=True)
                Sq[isc, oa[:, None], ob] = fa @ fb.T
                for c in range(3):
                    X[isc, oa[:, None], ob, c] = (fa*pts[:, c]) @ fb.T
    return X, Sq


def write_rmatrix(fname, geom, X, nspin, source):
    """tb-*_r.dat in the Wannier90 seedname_r.dat layout: integer R, then
    m n Re(x) Im(x) Re(y) Im(y) Re(z) Im(z) with m (the home-cell orbital)
    running fastest; nspin = 2 repeats X in both spin blocks (r is
    spin-diagonal), in the orbital order of tb-*.dat."""
    no, ncell = geom.no, geom.n_s
    with open(fname, 'w') as f:
        print(f'siesta2wtb.py: <m,0|r|n,R> (Angstrom) of the SIESTA basis of {source}', file=f)
        print(nspin*no, file=f)
        print(ncell, file=f)
        for isc in range(ncell):
            r1, r2, r3 = geom.lattice.sc_off[isc]
            for s in range(nspin):
                for n_ in range(no):
                    for s2 in range(nspin):
                        for m in range(no):
                            x = X[isc, m, n_] if s2 == s else np.zeros(3)
                            print('%5d%5d%5d%5d%5d' % (r1, r2, r3, s2*no+m+1, s*no+n_+1)
                                  + ''.join('%16.10f%16.10f' % (v, 0.0) for v in x), file=f)

###############################################################################

#inputfdf=  './teste-honpas/mos2.fdf'

inputfdf=  sys.argv[1]  
# siesta2wtb.py file.fdf --rmatrix: also write tb-*_r.dat (see position_matrix)
rmatrix = '--rmatrix' in sys.argv[2:]
# siesta2wtb.py file.fdf --formfactor --mesh NGX NGY NGZ [--ff-ecut 100]
# [--ff-spacing 0.1] [--yes]: also write tb-*_ff.bin, the form factors
# <phi_i,0|exp(iQ.r)|phi_j,R> of the basis, for the direct term of the BSE on
# that k mesh and the exchange term up to --ff-ecut eV (utils/wtb_formfactor.py).
# The size of the file is printed before anything is written; above 1 GB the
# script asks first, unless --yes.
formfactor = '--formfactor' in sys.argv[2:]
ffyes = '--yes' in sys.argv[2:]


def _option(name, count, cast, default):
    """the count values after name on the command line, or default"""
    if name not in sys.argv[2:]:
        return default
    i = sys.argv.index(name)
    values = [cast(v) for v in sys.argv[i + 1:i + 1 + count]]
    return values if count > 1 else values[0]


ffmesh = _option('--mesh', 3, int, None)
ffecut = _option('--ff-ecut', 1, float, 100.0)
ffspacing = _option('--ff-spacing', 1, float, 0.1)
if formfactor and ffmesh is None:
    sys.exit('--formfactor needs --mesh NGX NGY NGZ, the k-point mesh of the BSE')
#fermi= sys.argv[2]
fermi= 0.00

#os.system("cp ./teste-honpas/mos2.out ./teste-honpas/run.out")
geom = sisl.get_sile(inputfdf).read_geometry()
tshs = sisl.get_sile(inputfdf).read_hamiltonian(geometry=geom) #pegar hamiltoniano do siesta

#fermi = sisl.get_sile(folder).read_fermi_level()

scs = 0.0000 #scissors operator

nbasis= tshs.no
ncell=  tshs.nsc[0]*tshs.nsc[1]*tshs.nsc[2]
sptype=str(tshs.spin)

sptype2=sptype[5:-1]

if formfactor:
 import wtb_formfactor as wff
 ffsuffix, ffnspin = {'unpolarized': ('NP', 1), 'polarized': ('sp', 2),
                      'non-colinear': ('nc', 2), 'spin-orbit': ('soc', 2)}[sptype2]
 ffname = "tb-%s_ff.bin" % ffsuffix
 fflat = np.array(tshs.geometry.cell)
 fftau, ffrad = [], []
 for ia, atom in enumerate(tshs.geometry.atoms):
  for orbital in atom.orbitals:
   fftau.append(tshs.geometry.xyz[ia])
   ffrad.append(orbital.R)
 fftau, ffrad = np.array(fftau), np.array(ffrad)
 # no lattice vectors along a direction without neighbours in SIESTA's supercell
 ffR = wff.lattice_vectors(fflat, fftau, ffrad, periodic=[n > 1 for n in tshs.geometry.nsc])
 ffQd, ffQe = wff.direct_q(fflat, ffmesh), wff.exchange_g(fflat, ffecut)
 wff.confirm_size(ffname, ffnspin*tshs.no, len(ffR), len(ffQd) + len(ffQe), ffyes, len(ffQd), len(ffQe))

if sptype2 == 'unpolarized' :

 f = open("system-info-NP.txt", "w")
 print(tshs, file=f)
 f.close()

 nbasis= tshs.no
 ncell=  tshs.nsc[0]*tshs.nsc[1]*tshs.nsc[2]
 sptype=str(tshs.spin)


 f = open("tb-NP.dat", "w")

 print(calctype(sptype2),file=f,flush=True)
 print(ncases(scs),file=f,flush=True)
 #os.system("grep 'siesta:         Fermi =' run.out | awk '{print $4;}' >>  honpas_tb-NP.dat")
 print(fermi,file=f)
 print(ncases(tshs.cell[0,0]),'',ncases(tshs.cell[0,1]),'',ncases(tshs.cell[0,2]),file=f)
 print(ncases(tshs.cell[1,0]),'',ncases(tshs.cell[1,1]),'',ncases(tshs.cell[1,2]),file=f)
 print(ncases(tshs.cell[2,0]),'',ncases(tshs.cell[2,1]),'',ncases(tshs.cell[2,2]),file=f)
 print(nbasis,file=f)
 print(ncell,file=f)
 print('#rcell x',' ','rcell y',' ','rcell z',' ','i',' ','j',' ','ReH',' ','ImH',' ','S',file=f)

 #row index j runs fastest (Wannier90 _hr order), as hamiltonian_nort_input_read expects
 for i in range(0,ncell):
  for k in range(0,nbasis):
   for j in range(0,nbasis):

    a = tshs.geometry.o2sc(k+i*nbasis)[0]
    b = tshs.geometry.o2sc(k+i*nbasis)[1]
    c = tshs.geometry.o2sc(k+i*nbasis)[2]
		
    d = tshs.H[j,k+i*nbasis]
    e = tshs.S[j,k+i*nbasis]
	
    print(ncases(a),' ',ncases(b), ' ',ncases(c), ' ',j+1,' ',k+1,' ',ncases(d),' ',ncases(0.00),' ',ncases(e),file=f)


 f.close()

 f = open("basis_set-NP", "w")

 print('bindex aspecie ax ay az l m spin',file=f)

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
    xyz = tshs.geometry.xyz[ia, :]
    for orbital in atom:
        print(io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',0,file=f)
        io += 1

 f.close()


####################################################################

if sptype2 == 'polarized' :

 f = open("system-info-sp.txt", "w")
 print(tshs, file=f)
 f.close()

 nbasis= tshs.no
 ncell=  tshs.nsc[0]*tshs.nsc[1]*tshs.nsc[2]
 sptype=str(tshs.spin)


 f = open("basis_set-sp", "w")

 print('#bindex aspecie ax ay az l m spin',file=f)

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
    xyz = tshs.geometry.xyz[ia, :]
    for orbital in atom:
        print(io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',1,file=f)
        io += 1

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
    xyz = tshs.geometry.xyz[ia, :]
    for orbital in atom:
        print(nbasis+io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',-1,file=f)
        io += 1

 f.close()


 f = open("tb-sp.dat", "w")

 print(calctype(sptype2),file=f,flush=True)
 print(ncases(scs),file=f,flush=True)
 #os.system("grep 'siesta:         Fermi =' run.out | awk '{print $4;}' >>  honpas_tb-sp.dat")
 print(fermi,file=f)
 print(ncases(tshs.cell[0,0]),'',ncases(tshs.cell[0,1]),'',ncases(tshs.cell[0,2]),file=f)
 print(ncases(tshs.cell[1,0]),'',ncases(tshs.cell[1,1]),'',ncases(tshs.cell[1,2]),file=f)
 print(ncases(tshs.cell[2,0]),'',ncases(tshs.cell[2,1]),'',ncases(tshs.cell[2,2]),file=f)
 print(2*nbasis,file=f)
 print(ncell,file=f)
 print('#rcell x',' ','rcell y',' ','rcell z',' ','i',' ','j',' ','ReH',' ','ImH',' ','S',file=f)
#np.empty aloca todos os arrays vazios
#alocando todos os arrays com 0

 S = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))
 reH = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))
 imH = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))

 a = np.zeros((ncell+1,(2*nbasis)+2))
 b = np.zeros((ncell+1,(2*nbasis)+2))
 c = np.zeros((ncell+1,(2*nbasis)+2))


 for i in range(0,ncell):
  for j in range(0,nbasis): 
   for k in range(0,nbasis):

#parte up-up
  
    reH[i+1,j+1,k+1] = tshs[j,k+i*nbasis][0]
    #imH[i+1,j+1,k+1] = tshs[j,k+i*nbasis][4]
    S[i+1,j+1,k+1] = tshs.S[j,k+i*nbasis]

    

#parte up-dn

    #reH[i+1,(j+1),(k+1)+nbasis] = tshs[j,k+i*nbasis][2]
    #imH[i+1,(j+1),(k+1)+nbasis] = tshs[j,k+i*nbasis][3]
    #S[i+1,(j+1),(k+1)+nbasis] = 0.0

#parte dn-up

    #reH[i+1,(j+1)+nbasis,(k+1)] = tshs[j,k+i*nbasis][6]
    #imH[i+1,(j+1)+nbasis,(k+1)] = tshs[j,k+i*nbasis][7]
    #S[i+1,(j+1)+nbasis,(k+1)] = 0.0

#parte dn-dn

    reH[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs[j,k+i*nbasis][1]
    #imH[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs[j,k+i*nbasis][5]
    S[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs.S[j,k+i*nbasis]
    
#escrevendo a localização dos atomos da base

 for i in range(0,ncell):
  for k in range(0,nbasis): 

   a[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[0]
   b[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[1]
   c[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[2]
	
   a[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[0]
   b[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[1]
   c[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[2]
	

 #row index j runs fastest (Wannier90 _hr order), as hamiltonian_nort_input_read expects
 for i in range(1,ncell+1):
  for k in range(1,(2*nbasis)+1):
   for j in range(1,(2*nbasis)+1):
  	
#   a = tshs.geometry.o2sc(k+i*nbasis)[0]
#   b = tshs.geometry.o2sc(k+i*nbasis)[1]
#   c = tshs.geometry.o2sc(k+i*nbasis)[2]
		
#   d = tshs.H[j,k+i*nbasis]
#   e = tshs.S[j,k+i*nbasis]
	
    print(ncases(a[i,k]),' ',ncases(b[i,k]), ' ',ncases(c[i,k]), ' ',j,' ',k,' ',ncases(reH[i,j,k]),' ',ncases(imH[i,j,k]),' ',ncases(S[i,j,k]),file=f)


 f.close()

####################################################################

if sptype2 == 'non-colinear' :

 f = open("system-info-nc.txt", "w")
 print(tshs, file=f)
 f.close()

 nbasis= tshs.no
 ncell=  tshs.nsc[0]*tshs.nsc[1]*tshs.nsc[2]
 sptype=str(tshs.spin)


 f = open("basis_set-nc", "w")

 print('#bindex aspecie ax ay az l m spin',file=f)

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
    xyz = tshs.geometry.xyz[ia, :]
    for orbital in atom:
        print(io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',1,file=f)
        io += 1

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
    xyz = tshs.geometry.xyz[ia, :]
    for orbital in atom:
        print(nbasis+io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',-1,file=f)
        io += 1

 f.close()


 f = open("tb-nc.dat", "w")

 print(calctype(sptype2),file=f,flush=True)
 print(ncases(scs),file=f,flush=True)
 #os.system("grep 'siesta:         Fermi =' run.out | awk '{print $4;}' >>  honpas_tb-nc.dat")
 print(fermi,file=f)
 print(ncases(tshs.cell[0,0]),'',ncases(tshs.cell[0,1]),'',ncases(tshs.cell[0,2]),file=f)
 print(ncases(tshs.cell[1,0]),'',ncases(tshs.cell[1,1]),'',ncases(tshs.cell[1,2]),file=f)
 print(ncases(tshs.cell[2,0]),'',ncases(tshs.cell[2,1]),'',ncases(tshs.cell[2,2]),file=f)
 print(2*nbasis,file=f)
 print(ncell,file=f)
 print('#rcell x',' ','rcell y',' ','rcell z',' ','i',' ','j',' ','ReH',' ','ImH',' ','S',file=f)

#np.empty aloca todos os arrays vazios
#alocando todos os arrays com 0

 S = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))
 reH = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))
 imH = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))

 a = np.zeros((ncell+1,(2*nbasis)+2))
 b = np.zeros((ncell+1,(2*nbasis)+2))
 c = np.zeros((ncell+1,(2*nbasis)+2))


 for i in range(0,ncell):
  for j in range(0,nbasis): 
   for k in range(0,nbasis):

#parte up-up
  
     reH[i+1,j+1,k+1] = tshs[j,k+i*nbasis][0]
    #imH[i+1,j+1,k+1] = tshs[j,k+i*nbasis][4]
     S[i+1,j+1,k+1] = tshs.S[j,k+i*nbasis]

    

#parte up-dn (sisl: H_ud = D2 + i D3, H_du = D2 - i D3)

     reH[i+1,(j+1),(k+1)+nbasis] = tshs[j,k+i*nbasis][2]
     imH[i+1,(j+1),(k+1)+nbasis] = tshs[j,k+i*nbasis][3]
     S[i+1,(j+1),(k+1)+nbasis] = 0.0

#parte dn-up

     reH[i+1,(j+1)+nbasis,(k+1)] = tshs[j,k+i*nbasis][2]
     imH[i+1,(j+1)+nbasis,(k+1)] = -tshs[j,k+i*nbasis][3]
     S[i+1,(j+1)+nbasis,(k+1)] = 0.0

#parte dn-dn

     reH[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs[j,k+i*nbasis][1]
    #imH[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs[j,k+i*nbasis][5]
     S[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs.S[j,k+i*nbasis]
    
#escrevendo a localização dos atomos da base

 for i in range(0,ncell):
  for k in range(0,nbasis): 

   a[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[0]
   b[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[1]
   c[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[2]
	
   a[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[0]
   b[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[1]
   c[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[2]
	

 #row index j runs fastest (Wannier90 _hr order), as hamiltonian_nort_input_read expects
 for i in range(1,ncell+1):
  for k in range(1,(2*nbasis)+1):
   for j in range(1,(2*nbasis)+1):
  	
#   a = tshs.geometry.o2sc(k+i*nbasis)[0]
#   b = tshs.geometry.o2sc(k+i*nbasis)[1]
#   c = tshs.geometry.o2sc(k+i*nbasis)[2]
		
#   d = tshs.H[j,k+i*nbasis]
#   e = tshs.S[j,k+i*nbasis]
	
    print(ncases(a[i,k]),' ',ncases(b[i,k]), ' ',ncases(c[i,k]), ' ',j,' ',k,' ',ncases(reH[i,j,k]),' ',ncases(imH[i,j,k]),' ',ncases(S[i,j,k]),file=f)


 f.close()

####################################################################

if sptype2 == 'spin-orbit' :

 f = open("system-info-soc.txt", "w")
 print(tshs, file=f)
 f.close()
 
 nbasis= tshs.no
 ncell=  tshs.nsc[0]*tshs.nsc[1]*tshs.nsc[2]
 sptype=str(tshs.spin)

 f = open("basis_set-soc", "w")

 print('#bindex aspecie ax ay az l m spin',file=f)

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
     xyz = tshs.geometry.xyz[ia, :]
     for orbital in atom:
         print(io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',1,file=f)
         io += 1

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
     xyz = tshs.geometry.xyz[ia, :]
     for orbital in atom:
         print(nbasis+io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',-1,file=f)
         io += 1

 f.close()


 f = open("tb-soc.dat", "w")

 print(calctype(sptype2),file=f,flush=True)
 print(ncases(scs),file=f,flush=True)
 #os.system("grep 'siesta:         Fermi =' run.out | awk '{print $4;}' >>  honpas_tb-soc.dat")
 print(fermi,file=f)
 print(ncases(tshs.cell[0,0]),'',ncases(tshs.cell[0,1]),'',ncases(tshs.cell[0,2]),file=f)
 print(ncases(tshs.cell[1,0]),'',ncases(tshs.cell[1,1]),'',ncases(tshs.cell[1,2]),file=f)
 print(ncases(tshs.cell[2,0]),'',ncases(tshs.cell[2,1]),'',ncases(tshs.cell[2,2]),file=f)
 print(2*nbasis,file=f)
 print(ncell,file=f)
 print('#rcell x',' ','rcell y',' ','rcell z',' ','i',' ','j',' ','ReH',' ','ImH',' ','S',file=f)

#np.empty aloca todos os arrays vazios
#alocando todos os arrays com 0

 S = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))
 reH = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))
 imH = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))

 a = np.zeros((ncell+1,(2*nbasis)+2))
 b = np.zeros((ncell+1,(2*nbasis)+2))
 c = np.zeros((ncell+1,(2*nbasis)+2))


 for i in range(0,ncell):
  for j in range(0,nbasis): 
   for k in range(0,nbasis):

#parte up-up
  
    reH[i+1,j+1,k+1] = tshs[j,k+i*nbasis][0]
    imH[i+1,j+1,k+1] = tshs[j,k+i*nbasis][4]
    S[i+1,j+1,k+1] = tshs.S[j,k+i*nbasis]

    

#parte up-dn

    reH[i+1,(j+1),(k+1)+nbasis] = tshs[j,k+i*nbasis][2]
    imH[i+1,(j+1),(k+1)+nbasis] = tshs[j,k+i*nbasis][3]
    S[i+1,(j+1),(k+1)+nbasis] = 0.0

#parte dn-up

    reH[i+1,(j+1)+nbasis,(k+1)] = tshs[j,k+i*nbasis][6]
    imH[i+1,(j+1)+nbasis,(k+1)] = tshs[j,k+i*nbasis][7]
    S[i+1,(j+1)+nbasis,(k+1)] = 0.0

#parte dn-dn

    reH[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs[j,k+i*nbasis][1]
    imH[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs[j,k+i*nbasis][5]
    S[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs.S[j,k+i*nbasis]
    
#escrevendo a localização dos atomos da base

 for i in range(0,ncell):
  for k in range(0,nbasis): 

   a[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[0]
   b[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[1]
   c[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[2]
	
   a[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[0]
   b[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[1]
   c[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[2]
	

 #row index j runs fastest (Wannier90 _hr order), as hamiltonian_nort_input_read expects
 for i in range(1,ncell+1):
  for k in range(1,(2*nbasis)+1):
   for j in range(1,(2*nbasis)+1):
  	
#   a = tshs.geometry.o2sc(k+i*nbasis)[0]
#   b = tshs.geometry.o2sc(k+i*nbasis)[1]
#   c = tshs.geometry.o2sc(k+i*nbasis)[2]
		
#   d = tshs.H[j,k+i*nbasis]
#   e = tshs.S[j,k+i*nbasis]
	
    print(ncases(a[i,k]),' ',(b[i,k]), ' ',ncases(c[i,k]), ' ',j,' ',k,' ',ncases(reH[i,j,k]),' ',ncases(imH[i,j,k]),' ',ncases(S[i,j,k]),file=f)


 f.close()

#os.system("rm run.out")	
####################################################################

if formfactor:

 try:
  for atom in tshs.geometry.atoms.atom:
   for orbital in atom.orbitals:
    orbital.psi(np.zeros((1, 3)))
 except Exception:
  sys.exit('--formfactor needs the basis functions: run with the '
           '*.ion.nc or *.ion.xml files SIESTA wrote next to the fdf')
 ffgrid = wff.Grid.with_spacing(fflat, ffspacing)
 ffgroups = []
 for ia, atom in enumerate(tshs.geometry.atoms):
  A = tshs.geometry.xyz[ia]
  rmax = max(orbital.R for orbital in atom.orbitals)
  lo, shape = ffgrid.box(A, rmax)
  pts = ffgrid.points(lo, shape).reshape(-1, 3) - A
  values = np.array([orbital.psi(pts) for orbital in atom.orbitals], dtype=np.float32)
  ffgroups.append(wff.Group(tshs.geometry.a2o(ia, all=True), A, rmax, lo,
                            values.reshape((len(atom.orbitals),) + tuple(shape))))
 ffwriter = wff.FFWriter(ffname, fflat, np.vstack([fftau]*ffnspin), ffR, np.vstack([ffQd, ffQe]),
                         len(ffQd), ffmesh, ffecut)
 wff.compute(ffwriter, ffgrid, ffgroups, ffR, [(ffQd, 0), (ffQe, len(ffQd))],
             offsets=[s*tshs.no for s in range(ffnspin)])
 ffwriter.close()
 # F(R; Q = 0) is the overlap: against SIESTA's, on every lattice vector of the file
 ff = wff.read_ff(ffname)
 iq0 = int(np.argmin(np.abs(ffQd).sum(axis=1)))
 Ssiesta = tshs.tocsr(tshs.S_idx).toarray().reshape(nbasis, -1, nbasis).transpose(1, 0, 2)
 sc = [tuple(int(v) for v in o) for o in tshs.geometry.lattice.sc_off]
 dS = 0.0
 for iR, R in enumerate(ffR):
  S = Ssiesta[sc.index(tuple(int(v) for v in R))] if tuple(int(v) for v in R) in sc else 0.0
  dS = max(dS, float(np.abs(np.asarray(ff["F"][iq0, iR, :nbasis, :nbasis]).T - S).max()))
 print("%s written (BSE_FF_FILE; grid %dx%dx%d per cell); max |F(R;0) - S(SIESTA)| = %.1e" % (
     ffname, *ffgrid.M, dS))

####################################################################	

if rmatrix:

 suffix, nspin = {'unpolarized': ('NP', 1), 'polarized': ('sp', 2),
                  'non-colinear': ('nc', 2), 'spin-orbit': ('soc', 2)}[sptype2]
 X, Sq = position_matrix(tshs.geometry)
 # the quadrature's overlap against SIESTA's, over the whole supercell
 Ssiesta = tshs.tocsr(tshs.S_idx).toarray().reshape(nbasis, -1, nbasis).transpose(1, 0, 2)
 dS = np.abs(Ssiesta - Sq).max()
 write_rmatrix("tb-%s_r.dat" % suffix, tshs.geometry, X, nspin, inputfdf)
 print("tb-%s_r.dat written (BSE_CENTER_FILE); max |S(quadrature) - S(SIESTA)| = %.1e" % (suffix, dS))
