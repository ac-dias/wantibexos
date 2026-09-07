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

Q=0 Wannier position-matrix vertex
-----------------------------------

The zero- and finite-temperature optical BSE solver can include the Wannier90 position
matrix <0m|r|Rn> in its vertical optical vertex. Enable the Wannier90
`write_rmn` output and pass its `seedname_r.dat` file explicitly:

```
BSE_RMAT_FILE= seedname_r.dat
```

When this setting is absent, the legacy derivative-only optical matrix element
is retained unchanged. When it is present, the code reconstructs
`A(k) = sum_R exp(i k.R) <0|r|R>` exactly at each BSE k point and uses
`dH/dk - i[A,H]` for Q=0 transitions. This changes IPA and BSE oscillator
strengths but does not modify the BSE Hamiltonian or exciton energies.

This is the Q=0 `BSE= T` path; finite-Q optical vertices remain unchanged. It requires the
orthonormal Wannier representation (`DFT=W`), not the nonorthogonal `DFT=S`
path. The path is also
intentionally rejected when the companion `seedname_wsvec.dat` is present:
Wannier90 `use_ws_distance` requires its additional Wigner-Seitz correction,
which has not yet been implemented. `BSE_RMAT_FILE` is resolved from the run
directory, and all MPI ranks must be able to read it.

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
uses the existing 64-by-64 BLACS block-cyclic distribution and the ELPA 1-stage
complex-Hermitian solver.  ELPA builds both triangles of the Hermitian matrix
and needs an additional distributed local matrix for the eigenvectors while
diagonalizing, so its matrix construction work and peak memory are respectively
about twice and one extra local BSE matrix per MPI rank.  Executables built
without `-DELPA` reject `BSE_ALGO=elpa` explicitly; all other distributed
choices continue to use ScaLAPACK `PCHEEV`.

Finite-temperature optical BSE
--------------------------------

The Q=0 optical calculation (`BSE= T`) uses the same distributed matrix
assembly, ScaLAPACK/ELPA eigensolvers, and optical contractions at zero and
finite temperature. For example:

```
BSE= T
TEMP= 300
TA= FA
BSE_ALGO= elpa
```

Run ELPA with at least two MPI ranks and an ELPA-enabled executable. For
ScaLAPACK use `BSE_ALGO= cheev` with multiple MPI ranks; one rank uses LAPACK.
`TA= VE` and `TA= BE` retain their existing `ST`/`PHAVG` gap corrections;
`TA= FA` applies no gap correction. `BSE_RMAT_FILE`, tensor/circular optical
outputs, wavefunction output, and Hamiltonian restart also work at finite T.
This extension concerns optical Q=0 calculations, not distributed diagonalization
inside each finite-Q `BSE_BND` sector.

Let `F_i = f_v(i)-f_c(i)` and `D_i = E_c(i)-E_v(i)`. The solver constructs
`H_T = D + sqrt(F) K sqrt(F)` and uses optical vertices `sqrt(F) h`.
For positive occupations this is the similarity transform of `D + F K`;
it yields a Hermitian eigenproblem and the same linear-response resolvent.
Zero occupations are supported without division and produce dark, decoupled
transitions. Negative or nonfinite occupation differences are rejected.
The chemical potential remains zero in the existing single-particle energy
reference. This adds no temperature-dependent screening or phonon linewidths.

The old serial temperature routine passed an upper triangle of `D + F K` to
a Hermitian LAPACK solver even when F varied between transitions. Its old
finite-T results therefore need not match the corrected formulation.
`BSE_WF` now writes normalized eigenvectors `A` of `H_T`; the corresponding
right eigenvectors of the row-weighted kernel are `sqrt(F) A`.

New finite-T checkpoints use matrix-kind 3 and reject the old row-weighted
kind-2 files. As with other Hamiltonian restarts, keep all physical inputs,
including `TEMP`, `TA`, `ST`, and `PHAVG`, unchanged: the header does not
fingerprint physical parameters. ELPA stores both triangles, whereas ScaLAPACK
stores the upper triangle; use the matching solver layout when restarting.

Run `make -j1 test-bse-temperature` for the occupation, similarity-transform,
Hermiticity, dark-state, and optical-response tests. The end-to-end test uses a
supplied Wannier fixture with two valence and two conduction bands, three built
executables, and NumPy:

```
python3 tests/test_bse_temperature_mpi.py \
  --case /path/to/bse-test-pos \
  --serial-exe /path/to/serial/wtb.x \
  --scalapack-exe /path/to/scalapack/wtb.x \
  --elpa-exe /path/to/elpa/wtb.x \
  --output /path/to/new-test-directory
```

It checks 0 K, 300 K, a 10000 K numerical occupation-weight stress test,
VE/BE gap corrections, all optical tensor components, circular polarization,
wavefunction output, restart, and rejection of legacy finite-T checkpoints.
Spectra are broadened before comparing solvers to avoid dependence on
rotations within degenerate exciton subspaces.

Standard four-way BSE benchmark
-------------------------------

The reproducible Q=0 regression runs the original `wantibexos-dev` serial
reference followed by the current serial, two-rank ScaLAPACK, and two-rank
ELPA implementations. The three current runs use the same Wannier90
position-matrix file; the legacy reference deliberately retains its original
vertex. Run it with:

```
make benchmark-bse
```

Each invocation performs fresh isolated current/ELPA builds and creates a new
timestamped directory under `benchmark-results/`. It records executable and
input hashes, exact inputs and logs, Lorentzian-broadened absorption/oscillator
figures, energy residuals, tensor-resolved totals, and degeneracy-safe
numerical metrics. Use `make benchmark-bse BENCHMARK_ARGS='--dry-run'` to
inspect the plan. The configurable paths and output inventory are documented
in [`benchmarks/README.md`](benchmarks/README.md).

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
