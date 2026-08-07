#FOR=gfortran-mp-15
FOR=/opt/local/bin/mpifort-mpich-mp
OMP=-fopenmp 

ELPA_BUILD=/Users/reilya/Dropbox/pacotes/elpa-2025.06.002/build-openmp
ELPA_API_VERSION=20250131

# Use the static library from the ELPA build tree.
#LIBS=$(ELPA_BUILD)/.libs/libelpa_openmp.a -L/opt/local/lib -lscalapack -lopenblas

# ELPA module files and conditional compilation.
#EXTRA=-DELPA -DELPA_API_VERSION=$(ELPA_API_VERSION) -I$(ELPA_BUILD)/modules

LIBS= -L/opt/local/lib  -lopenblas  -lscalapack 
#LIBS = -L/opt/local/lib  -lopenblas 
DIR="./bin/"
COND= -DMPI -O0 -fbacktrace -fcheck=all -g
#COND= -O0 -fbacktrace -fcheck=all -g 
