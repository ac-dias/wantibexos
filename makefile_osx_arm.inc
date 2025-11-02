### If gcc and opeblas are installed via hemebrew, you can use the command below to get the library prefix

OPENBLAS_PREFIX = $(shell brew --prefix openblas)

FOR   = gfortran
OMP   = -fopenmp

### To use OpenBLAS as a linear algebra library
LIBS  = -I$(OPENBLAS_PREFIX)/include -L$(OPENBLAS_PREFIX)/lib -lopenblas

### If necessary, use Apple's Accelerate library, which is optimized for M-series ARM chips.
#LIBS = -framework Accelerate

EXTRA= -Wno-argument-mismatch -Wargument-mismatch  -Wunderflow -Werror=line-truncation -DOPENBLAS
COND  = -cpp -O3 -march=native
DIR   = ./bin
