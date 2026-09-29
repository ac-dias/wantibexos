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

Q=0 Wannier position data
--------------------------

The zero-temperature optical BSE solver can use the Wannier90 position matrix
<0m|r|Rn> both in its vertical optical vertex and in the G=0 direct Coulomb
kernel. Enable the Wannier90 `write_rmn` output and pass its `seedname_r.dat`
file explicitly:

```
BSE_CENTER_FILE= seedname_r.dat
```

When this setting is absent, the legacy derivative-only optical matrix element
and scalar Coulomb kernel are retained unchanged. When it is present, the code
reconstructs
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

`BSE_CENTER_GRAD= T` adds the next order in q of the vertices: the first-order
term of <m,0|exp(iq.r)|n,R> expanded about the pair midpoint (t_m+t_n+R)/2,

```
i q.sum_mn C1_m^* C2_n exp[+i q.(t_m+t_n)/2] [A(kb) - diag(t)]_mn
```

added to the electron vertex, and the same with V2^*, V1, -q and
exp[-i q.(t_m+t_n)/2] to the hole vertex, with A at the mean kb = (k1+k2)/2 of
the two states. Without the midpoint phase the term would depend on the
coordinate origin and break the lattice symmetry. It uses the whole position
matrix and costs one Fourier sum of r(R) per matrix element.

The finite-temperature optical BSE (`TEMP` > 0) uses the same vertex and
kernel. The q-path BSE (`BSE_BND= T`, at zero and finite temperature) uses the
centre phases in its direct kernel and in the exchange vertices

```
sum_m C_m^*(k+Q) V_m(k) exp[+i Q.t_m],
```

with Q the shortest image of the exciton momentum; `BSE_CENTER_GRAD` applies
to the Q=0 kernels only. The centre phases require the
orthonormal Wannier representation (`DFT=W`), not the nonorthogonal `DFT=S`
path. The path is also intentionally rejected when the companion
`seedname_wsvec.dat` is present: Wannier90 `use_ws_distance` requires its
additional Wigner-Seitz correction, which has not yet been implemented.
`BSE_CENTER_FILE` is resolved from the run directory, and all MPI ranks must be
able to read it.

The former `BSE_RMAT_FILE` input name remains available as an optical-only
compatibility alias. It applies `dH/dk - i[A,H]` but deliberately retains the
legacy scalar BSE kernel, so existing inputs do not silently acquire different
exciton energies. New calculations using the centre-resolved direct kernel
should use `BSE_CENTER_FILE`.

For `DFT=S` (SIESTA/Honpas files from `utils/siesta2wtb.py`) the orbital
centres come from the `basis_set-*` file written next to `tb-*.dat`, and
`BSE_CENTER_FILE` is rejected. The optical vertex is
`<c| dH/dk - (Ec+Ev)/2 dS/dk |v>/(Ec-Ev)` in the atomic gauge, which treats
each orbital as if its charge sat on its centre. `BSE_RMAT_FILE= tb-NP_r.dat`,
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
`PARAMS_FILE`, and (for `BSE_BND`) `KPATH_BSE`. `--nodes` means MPI ranks, not
physical compute nodes. It reports the Hamiltonian, explicit source-array peak
per rank, and a conservative planning value that includes source-level
OpenMP eigensystem scratch and a configurable safety factor. MPI/OpenMP,
BLAS, ELPA, allocator, and operating-system overhead are outside the source
model; add a per-thread reserve with `--external-thread-mib` when appropriate.

The Siesta/Honpas Hamiltonian extract script (siesta2wtb.py) was tested in SISL version 0.16.2, could not be work in other versions.
`siesta2wtb.py file.fdf --rmatrix` also writes `tb-*_r.dat`, the position matrix `<a,0|r|b,R>` of the basis (see `BSE_RMAT_FILE` above); it needs the `*.ion.nc` or `*.ion.xml` files SIESTA writes next to the fdf.

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
