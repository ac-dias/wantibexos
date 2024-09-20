# WanTiBEXOS
WanTiBEXOS code repository

This code, in the actual version, can only be compiled with intel fortran compiler - ifort

For ifort older versions, please change the flag -qmkl in makefile for -mkl, to compile the code.

The online documentation is available in:
[https://wantibexos.readthedocs.io/en/latest/](https://wantibexos.readthedocs.io/)

For memory problems during parallel run, please export the following environment variable:

export KMP_STACKSIZE=XXXmb, 
being XXX the amount of virtual RAM per thread, I suggest something around 300mb, but for some situations, more could be necessary.

Citing
   ------

   If you would like to cite the WanTiBEXOS package in a paper or presentation, the
   following references can be used:

- BibTeX::
        
        @article{DIAS_CPC_2022,
          title = {WanTiBEXOS: a Wannier based Tight Binding code for electronic band structure, excitonic and optoelectronic properties of solids},
          journal = {Computer Physics Communications},
          pages = {108636},
          year = {2022},
          issn = {0010-4655},
          doi = {https://doi.org/10.1016/j.cpc.2022.108636},
          url = {https://www.sciencedirect.com/science/article/pii/S0010465522003551},
          author = {Alexandre C. Dias and Julian F.R.V. Silveira and Fanyao Qu},
          keywords = {Tight-Binding, Wannier functions, Excitons, Electronic and optical properties}
        }
