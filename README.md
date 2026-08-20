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
