# wantibexos-dev
WanTiBEXOS DEV Repository

The online documentation is available in:
[https://wantibexos.readthedocs.io/en/latest/](https://wantibexos.readthedocs.io/)

For memory problems during parallel run, please export the following environment variable:
export KMP_STACKSIZE=XXXmb, being XXX the amount of virtual RAM per thread, I suggest something around 300mb, but for some situations, more could be necessary.

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
