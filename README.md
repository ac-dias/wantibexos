# wantibexos-dev
WanTiBEXOS DEV Repository

The online documentation is available in:
[https://wantibexos.readthedocs.io/en/latest/](https://wantibexos.readthedocs.io/)

Performance build
-----------------

The checked-in `makefile.inc` is a diagnostic configuration (`-O0` with
runtime checks).  For production runs on the MacPorts MPI/OpenBLAS setup used
by this repository, copy `makefiles/makefile-osx-mpich+openblas-release` to
`makefile.inc` before building.  It enables `-O3 -march=native` but does not
enable unsafe floating-point transformations.  Use a parallel build to reduce
compilation time, for example `make -j 8`.

The BSE dielectric post-processing routines use OpenMP and take their thread
count from the `NTHREADS=` input setting.  Choose a value that matches the CPU
cores allocated to the calculation.

For memory problems during parallel run, please export the following environment variable:
export KMP_STACKSIZE=XXXmb, being XXX the amount of virtual RAM per thread, I suggest something around 300mb, but for some situations, more could be necessary.

BSE Hamiltonian restart
-----------------------

The optical BSE solver can checkpoint the Hamiltonian before diagonalization and
reuse it in a later run.  The feature is disabled by default.  To save the
matrix, add the following to the input:

```
BSE_HAM_SAVE= T
BSE_HAM_FILE= bse_hamiltonian.bin
```

`BSE_HAM_FILE` is relative to `OUTPUT=`.  For a restart, keep all physical BSE
settings and the input Hamiltonian unchanged, set `BSE_HAM_READ= T`, and leave
`BSE_HAM_SAVE` false.  The restart still evaluates the single-particle
quantities needed for the optical spectrum, but skips construction of the BSE
Hamiltonian.  Checkpoints are native binary files and should be reused with the
same executable/platform.  Their header validates the matrix dimension and
matrix layout, but cannot detect changed physical inputs.

For distributed MPI BSE runs, the checkpoint is a single global matrix written
with MPI-IO.  It can be restarted with a different MPI rank count or
process-grid layout; each rank reads its own block-cyclic part of the matrix.

Q=0 position matrix and orbital centres
---------------------------------------

The optical BSE solvers can use the position matrix <0m|r|Rn> of the
tight-binding basis, for `DFT=W` (Wannier90's `write_rmn` output,
`seedname_r.dat`) and for `DFT=S` (`siesta2wtb.py file.fdf --rmatrix`,
`tb-NP_r.dat`), with the same keyword:

```
BSE_CENTER_FILE= seedname_r.dat
```

When this setting is absent, the legacy derivative-only optical matrix element
and, for `DFT=W`, the scalar Coulomb kernel are retained unchanged. For `DFT=W`
the code reconstructs
`A(k) = sum_R exp(i k.R) <0|r|R>` exactly at each BSE k point and uses
`dH/dk - i[A,H]` for Q=0 transitions. It also extracts the Wannier centres
`t_m=<0m|r|0m>` from the R=0 diagonal and replaces the direct-kernel vertices by

```
sum_m C1_m^* C2_m exp[+i (k1-k2).t_m]
sum_n V2_n^* V1_n exp[-i (k1-k2).t_n].
```

This is the centre-resolved G=0 approximation
`V_mn(q)=V(q) exp[i q.(t_m-t_n)]`, with q the shortest image of k1-k2 (see
BSE kernel conventions below). It changes the BSE Hamiltonian, exciton
energies, and oscillator strengths. The screened scalar Coulomb potential and
its q=0 averaging are otherwise unchanged.

`BSE_CENTER_GRAD= T` (`DFT=W`) adds the next order in q of the vertices: the first-order
term of <m,0|exp(iq.r)|n,R> expanded about the pair midpoint (t_m+t_n+R)/2,

```
i q.sum_mn C1_m^* C2_n exp[+i q.(t_m+t_n)/2] [A(kb) - diag(t)]_mn
```

added to the electron vertex, and the same with V2^*, V1, -q and
exp[-i q.(t_m+t_n)/2] to the hole vertex, with A at the mean kb = (k1+k2)/2 of
the two states. Without the midpoint phase the term would depend on the
coordinate origin and break the lattice symmetry. It uses the whole position
matrix and costs one Fourier sum of r(R) per matrix element. It is not
implemented for `DFT=S`, where the form factors (`BSE_FF`, below) give the
vertices exactly, and it cannot be combined with `BSE_FF`.

The finite-temperature optical BSE (`TEMP` > 0) uses the same vertex and
kernel, and the single-particle tools (`IPA= T`, `OPT_BZ= T`) the same optical
vertex, for `DFT=W` as for `DFT=S`. The q-path BSE (`BSE_BND= T`, at zero and
finite temperature) uses the centre phases in its direct kernel and in the
exchange vertices

```
sum_m C_m^*(k+Q) V_m(k) exp[+i Q.t_m],
```

with Q the shortest image of the exciton momentum; `BSE_CENTER_GRAD` applies
to the Q=0 kernels only. For `DFT=W` the path is intentionally rejected when the companion
`seedname_wsvec.dat` is present: Wannier90 `use_ws_distance` requires its
additional Wigner-Seitz correction, which has not yet been implemented.
`BSE_CENTER_FILE` is resolved from the run directory, and all MPI ranks must be
able to read it.

`BSE_RMAT_FILE`, the former name of this keyword (the optical vertex only,
without the Wannier centres in the `DFT=W` kernel), is no longer read: an input
that has it stops with a message to rename it `BSE_CENTER_FILE`. For `DFT=S`
the results with the new name are the same; for `DFT=W` the kernel then has
its centre phases.

For `DFT=S` (SIESTA/Honpas files from `utils/siesta2wtb.py`) the orbital
centres come from the `basis_set-*` file written next to `tb-*.dat`, and the
kernel uses them with or without `BSE_CENTER_FILE`. The optical vertex is
`<c| dH/dk - (Ec+Ev)/2 dS/dk |v>/(Ec-Ev)` in the atomic gauge, which treats
each orbital as if its charge sat on its centre. `BSE_CENTER_FILE= tb-NP_r.dat`,
the position matrix `<a,0|r|b,R>` of the basis written by
`siesta2wtb.py file.fdf --rmatrix`, adds what the orbital shapes give, the
dipoles between basis orbitals

```
<c|D|v>,    D_ab(R) = <a,0| r - (t_a+t_b+R)/2 |b,R>,
```

which makes the vertex the full interband dipole of the basis. For monolayer
h-BN they lower the strength of the gap transition at K by 11% and bring the
momentum matrix elements to within 1.5% of plane-wave (Quantum ESPRESSO)
values with the same pseudopotentials, from 8-21% too large. The kernel does
not change, and neither do exciton energies.

The scissor on the second line of the tight-binding file shifts the
conduction bands, and so every transition energy (the BSE diagonal and the
IPA spectra). The optical vertex takes the unshifted energies: the scissor
leaves the eigenvectors, and so the interband dipoles and the oscillator
strengths, as they are.
The Berry curvature (`BERRY`, `BERRY_BZ`) is a property of the eigenvectors and
does not use the scissor either.

Two-dimensional Coulomb interactions
------------------------------------

`COULOMB_POT= V2DTPV` is the two-dimensional interaction of Trolle, Pedersen
and Veniard, Sci. Rep. 7, 39844 (2017), Eq. (10), per unit area like `V2D`,
`V2DK` and `V2DRK`: the Coulomb interaction of two charges spread uniformly
across a slab of thickness d, screened by the finite-thickness model dielectric
function of their Eqs. (1) and (12),

```
W(q) = e^2 F(|q| d) / (2 eps0 A |q| eps_TPV(|q|))
F(b) = 2 (b - 1 + exp(-b)) / b^2
```

with `TPV_KAPPA`, `TPV_QTF` (1/Angstrom), `TPV_HWP` (eV), `TPV_THICKNESS` = d
(Angstrom), `TPV_ALPHA` and the half-spaces `EDIEL_T`, `EDIEL_B`; `EDIEL` other
than 1 switches the screening on (`EDIEL= 1` is the bare slab interaction
F/(2 eps0 A |q|)). eps_TPV is defined as the ratio of the bare to the screened
interaction, both averaged across the slab, so F is part of W: eps_TPV goes
back to 1 at large |q|, and without F the interaction would return to the
point-charge Coulomb there. F = 1 for d = 0. The q = 0 term is W averaged over
the k-grid cell around q = 0.

`V2DTAVG` and `V2DAVGTPV` instead use the cell-averaged truncated Coulomb
`4 pi e^2 (1 - exp(-|q| L/2)) / (V q^2)`, divided by `EDIEL` or by eps_TPV(q):
the point-charge 2D Coulomb times `2 (1 - exp(-|q| L/2)) / (|q| L)`, as if the
charges were spread over the whole cell height L instead of d. Their excitons
therefore depend on the vacuum: for monolayer h-BN (`V2DTAVG`, `EDIEL= 1.7`)
the binding energy of the lowest exciton drops from 1.38 to 0.86 eV when the
cell height grows from 12.5 to 25 Angstrom, while with `V2DTPV` it does not
change.

BSE kernel conventions
----------------------

All BSE kernels (Q=0 and finite Q, zero and finite temperature, `DFT=W` and
`DFT=S`) use the standard direct term

```
W(q) <c1|c2> <v2|v1>,    q = k1 - k2 - G,
```

with the hole overlap `<v2|v1>` (earlier versions used its complex conjugate
`<v1|v2>`), and with q the shortest image of the grid difference k1-k2: the
screened potential and the Wannier-centre phases are not periodic in q, and
the raw difference made the kernel depend on how the k-grid cell is chosen.
Where several images are equally short (the zone boundary), the kernels
average over them; the finite-Q exchange term takes the shortest images of Q
in the same way. Exciton oscillator strengths are contracted as
`|sum_j A_j^* <c|r|v>_j|^2`, consistent with these overlaps (earlier versions
used `A_j`).

Together with the Wannier-centre phases of `BSE_CENTER_FILE`, this keeps the
lattice symmetry of the kernel: in monolayer h-BN the lowest bright exciton is
the degenerate E' doublet, independent of the k-grid cell,
where the earlier kernel split it by 145 meV (36x36 grid), and the q-path
energies at the three M points agree to 0.5 meV (16 meV apart without the
centre phases). BSE Hamiltonians
(`BSE_HAM_SAVE`) and q-path checkpoints written by earlier builds hold the
earlier kernel and must be regenerated.

ELPA BSE diagonalization
-------------------------

The distributed optical-BSE solver can use ELPA instead of ScaLAPACK.  Build an
MPI+OpenMP executable with `makefiles/makefile-intel(ifx)+mkl+elpa` as the
`makefile.inc` template, after loading ELPA and Intel MPI modules that were built
against the same MPI implementation.  Set `ELPA_ROOT`, or override
`ELPA_FCFLAGS` and `ELPA_LDFLAGS` with the paths supplied by the cluster.
ELPA itself must have been built with OpenMP enabled; the code passes its
existing `nthreads` setting to ELPA as `omp_threads`.

Set `BSE_ALGO=elpa` in the input and launch with at least two MPI ranks.  ELPA
uses the existing 64-by-64 BLACS block-cyclic distribution and the ELPA 2-stage
complex-Hermitian solver.  ELPA builds both triangles of the Hermitian matrix
and needs an additional distributed local matrix for the eigenvectors while
diagonalizing, so its matrix construction work and peak memory are respectively
about twice and one extra local BSE matrix per MPI rank.  Executables built
without `-DELPA` reject `BSE_ALGO=elpa` explicitly; all other distributed
choices continue to use ScaLAPACK `PCHEEV`.

BSE q-path MPI parallelism
--------------------------

The BSE q-path calculation (`BSE_BND= T`) can distribute independent exciton
momenta over MPI ranks.  The default is:

```
BSE_KPATH_MPI= Q
```

Each rank diagonalizes complete fixed-Q BSE Hamiltonians for a cyclic subset of
the requested path, and rank zero collects the eigenvalues and writes the
q-path output.  This is therefore most useful when the number of q-path points
is at least the number of MPI ranks.  `HYBRID` is reserved for a future mode
that will split ranks into multiple intra-Q communicator groups; it is not yet
implemented and is rejected explicitly.

For restartable q-path calculations, enable per-Q checkpoints:

```
BSE_KPATH_CHECKPOINT= T
BSE_KPATH_CHECKPOINT_FILE= bse_kpath_checkpoint
```

Each completed Q point writes its eigenvalue block in the output directory as
`<prefix>_q########.bin`.  On the next run with the same option, valid
files matching the Q point and BSE dimension are reused automatically; missing
or incomplete files are recalculated.  The file is first written as a temporary
file and then atomically renamed, so an interrupted write cannot replace a
previous valid checkpoint.  Checkpoints are native binary files: reuse them
only with the same executable/platform and unchanged physical BSE inputs.  Use
a different prefix or disable the option for a fresh calculation.

BSE q-grid sampling
-------------------

`KPATH_BSE=` normally points to the legacy q-path file: an even number of
endpoints, then the number of points per segment, then the endpoint records.
That format is unchanged. The same file can now instead request a uniform
reciprocal-space grid:

```
GRID
8 8 1
0.0 0.0 0.0
```

The first line selects grid mode. The second gives `Nq1 Nq2 Nq3`; the optional
third line is a shift in units of a grid spacing and defaults to `0 0 0`.
Thus the example covers an 8-by-8-by-1 grid and includes Gamma. A shift of
`0.5 0.5 0.0` creates a half-grid-shifted mesh. The grid is enumerated in the
reciprocal-lattice basis and passed to the finite-Q BSE solver in Cartesian
coordinates. In grid mode `bands_bse.dat` has four columns: Cartesian
`qx qy qz` and exciton energy. Path mode retains its existing two-column
path-distance/energy output and `KLABELS-BSE.dat` file.

BSE memory estimates
--------------------

Before allocating a dense BSE Hamiltonian, the BSE log and standard output now
report its exact storage as a default-precision complex matrix (8 bytes per
element). Optical BSE reports both the global dense equivalent and, for a
distributed run, the rank-zero local block. Finite-Q BSE reports the full
per-active-rank matrix, because MPI-over-Q assigns complete Q sectors rather
than distributing one Hamiltonian.

For a pre-run estimate, use:

```
python3 utils/bse_memory_estimate.py input.dat --nodes 4 --threads 8
```

The script reads `NGX`, `NGY`, `NGZ`, `NBANDSC`, `NBANDSV`, `BSE_ALGO`,
`PARAMS_FILE`, (for `BSE_FF`) `BSE_FF_FILE` and `BSE_FF_EXCHANGE`, which every
rank holds while it builds its part of the Hamiltonian, and (for `BSE_BND`)
`KPATH_BSE`. `--nodes` means MPI ranks, not
physical compute nodes. It reports the Hamiltonian, explicit source-array peak
per rank, and a conservative planning value that includes source-level
OpenMP eigensystem scratch and a configurable safety factor. MPI/OpenMP,
BLAS, ELPA, allocator, and operating-system overhead are outside the source
model; add a per-thread reserve with `--external-thread-mib` when appropriate.

The Siesta/Honpas Hamiltonian extract script (siesta2wtb.py) was tested in SISL version 0.16.2, could not be work in other versions.
`siesta2wtb.py file.fdf --rmatrix` also writes `tb-*_r.dat`, the position matrix `<a,0|r|b,R>` of the basis (its `BSE_CENTER_FILE`, see above); it needs the `*.ion.nc` or `*.ion.xml` files SIESTA writes next to the fdf.
The PAOFLOW script, `paoflow2wtb.py prefix.save [--configuration minimal|standard|extended] [--basispath DIR] [--pthr 0.95] [--shift auto]`, builds the PAOFLOW tight-binding Hamiltonian of a Quantum ESPRESSO run (projections, projectability, pao_hamiltonian; the save directory of a pw.x run on a Monkhorst-Pack grid) and writes it for both readers: `paoflow-NP.dat` with `paoflow_r.dat` (the orbital centres, for `BSE_CENTER_FILE`) for `DFT=W`, `tb-NP.dat` with `basis_set-NP` for `DFT=S`. H(R) is placed at the minimal images of each orbital pair, as Wannier90's `use_ws_distance`. It was tested with PAOFLOW 3.0.0, unpolarized runs only; run it with `mpirun` to use PAOFLOW's MPI parallelism.
`wannier2wtb.py seedname [--efermi E]` writes `seedname-NP.dat`, the `DFT=W` tight-binding file (the WanTiBEXOS header followed by `seedname_hr.dat`), from a Wannier90 run (`seedname.wout`, `seedname_hr.dat`); `seedname_r.dat` is its `BSE_CENTER_FILE`.

Orbital form factors
--------------------

`siesta2wtb.py`, `paoflow2wtb.py` and `wannier2wtb.py` take `--formfactor --mesh NGX NGY NGZ` to also write
`*_ff.bin`, the form factors of the orbitals of the tight-binding basis,

```
F_mn(R; Q) = <m,0| exp(iQ.r) |n,R>,
<c k| exp(iQ.r) |c' k'> = sum_mn c_m(k)^* c'_n(k') sum_R exp(ik'.R) F_mn(R; Q),
```

with which the pair densities of the tight-binding states are exact instead of point-centred
(`sum_m c_m^* c'_m exp(iQ.t_m)`); `F(R; 0)` is the overlap. The file holds two sets of Q:

- every shortest image q - G of the vectors q of the BSE mesh `NGX NGY NGZ`, ties on the zone
  boundary included, as `bse_q_images` chooses them: the direct term;
- the G != 0 with hbar^2 G^2/2m up to `--ff-ecut` (default 100 eV): the exchange term (the local
  fields).

The size of the file is printed before anything is computed; above 1 GB the scripts ask for
confirmation, and stop when they cannot ask (a batch job), unless `--yes` is given. The layout (a
little-endian stream in single precision, for `access='stream'`) is described in
`utils/wtb_formfactor.py`, which the three scripts share. The BSE reads them with `BSE_FF` (below).

- SIESTA, `siesta2wtb.py file.fdf --formfactor --mesh ...`: the numerical orbitals of the
  `*.ion.nc` or `*.ion.xml` files, integrated where they overlap (`--ff-spacing`, default 0.15 A:
  for h-BN with a DZP basis the lowest exciton is within 0.01 meV of 0.1 A, in a fifth of the
  time, and 0.3 A moves it by 10 meV); `F(R; 0)` is checked against SIESTA's overlap.
- PAOFLOW, `paoflow2wtb.py prefix.save ... --formfactor --mesh ...`: the Loewdin orbitals of the
  model, from PAOFLOW's orthonormalised Bloch sums on its k grid; `--ff-tail` (default 1e-3) sets
  their radii from their own norm. The pseudo-atomic orbitals, whose radii are smaller, give the
  first estimate of the size printed before PAOFLOW starts; the final size follows before the
  integrals. The Loewdin orbitals of the box-state bases (`--configuration standard`,
  `extended`) are delocalized: for monolayer h-BN 1% of their norm lies beyond 6-11 A. With the
  minimal basis the default tail gives a 0.042 GB file on a 36x36 mesh, with the lowest exciton
  within 0.15 meV and the exchange term unchanged against `1e-6` (0.17 GB, nine times the time);
  the box-state bases have radii of 8.6-16 A at that tail and files of 4.6 GB (standard) and
  8.5 GB (extended), 2.1 and 3.7 GB at `1e-2`, 46 and 76 GB at `1e-6`. Radii from the
  pseudo-atomic orbitals, before, left 1-11% of their norm out and moved the lowest exciton by
  5-12 meV.
- Wannier90, `wannier2wtb.py seedname --formfactor --mesh ...`: the Wannier functions Wannier90
  plots (`wannier_plot = .true.`, `wannier_plot_format = xcrysden`, `wannier_plot_mode = crystal`,
  a `wannier_plot_supercell` that holds each function, on UNK files of pw2wannier90 with
  `write_unk = .true.` and `reduce_unk = .false.`), with their phases matched to `seedname_r.dat`;
  real Wannier functions only. `--ff-radius` (default 5) times the square root of each spread sets
  their radii; the plot supercell should reach that far. The script prints the norm of each function
  within its radius and how far `F(R; 0)` is from the identity: check both.

The optical BSE (`BSE= T`, Q=0, at zero and finite temperature) uses them with

```
BSE_FF= T                    # default F: the file is not read
BSE_FF_FILE= "hbn_ff.bin"
BSE_FF_EXCHANGE= T           # optional, default F: the exchange term
BSE_FF_ECUT= 50              # optional: the G of the exchange term up to 50 eV (default: all)
```

(quote file names that contain a `/`, as for `OUTPUT`). `BSE_FF= T` replaces the point-centre
vertices of the direct kernel (the Wannier-centre phases of `BSE_CENTER_FILE` for `DFT=W`, the
orbital centres of `basis_set-*` for `DFT=S`) by the pair densities of the orbitals,

```
<c k1| exp(iq.r) |c' k2> = sum_mn c_m(k1)^* c'_n(k2) sum_R exp(ik2.R) F_mn(R; q),
```

at each shortest image q of k1-k2, for `DFT=W` and `DFT=S` alike. The optical vertex does not
change (`BSE_CENTER_FILE` gives the position matrix). The number of orbitals and the lattice must be
those of the tight-binding file, and the BSE mesh must divide the `--mesh` of the file; the solver
checks that every q of the BSE mesh is in the file, and writes to the log the largest deviation of
the overlaps of the BSE bands, computed from F(R; 0), from the identity (for `DFT=S`, F(R; 0) is
S(R)). `BSE_FF_EXCHANGE= T` adds the exchange term with the G != 0 of the file,

```
K^x(vck, v'c'k') = spinf sum_G v(G) <c k| exp(iG.r) |v k> <v' k'| exp(-iG.r) |c' k'>,
```

v(G) = e^2/(eps0 G^2 N_k Omega) the bare Coulomb interaction (for `SYSDIM= "2D"` truncated at half the
cell height, as `V2DT`; 1D and 0D systems are not implemented), spinf = 2 for the singlets of an
unpolarized (NP) model and 1 otherwise. G = 0 is left out, so the BSE gives the macroscopic
dielectric function with local fields, as Yambo's exchange term (`BSENGexx`); without the exchange
term the excitons of an unpolarized model are the triplets. In a layer the G along the axis hardly
couple in-plane transitions (for h-BN they do not move the lowest exciton); the first shell with an
in-plane component is at 32 eV, and 60 eV give 95% of the shift of 100 eV, the default of
`--ff-ecut`. `BSE_FF_ECUT` converges the term with one file. Saved Hamiltonians
(`BSE_HAM_SAVE`) record whether the form factors and the exchange term were used. `BSE_CENTER_GRAD`
cannot be combined with `BSE_FF`, and the q-path BSE (`BSE_BND`) keeps the centre phases.

Citing
   ------

   If you would like to cite the WanTiBEXOS package in a paper or presentation, the
   following references can be used:

- BibTeX::
        
        @article{DIAS_CPC_2023,
          title = {WanTiBEXOS: a Wannier based Tight Binding code for electronic band structure, excitonic and optoelectronic properties of solids},
          journal = {Computer Physics Communications},
          pages = {108636},
          year = {2023},
          issn = {0010-4655},
          doi = {https://doi.org/10.1016/j.cpc.2022.108636},
          url = {https://www.sciencedirect.com/science/article/pii/S0010465522003551},
          author = {Alexandre C. Dias and Julian F.R.V. Silveira and Fanyao Qu},
          keywords = {Tight-Binding, Wannier functions, Excitons, Electronic and optical properties}
        }
