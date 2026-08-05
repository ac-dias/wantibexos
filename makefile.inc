#FOR=gfortran-mp-15
FOR=/opt/local/bin/mpifort-mpich-mp
OMP=-fopenmp 
LIBS= -L/opt/local/lib  -lopenblas  -lscalapack 
#LIBS = -L/opt/local/lib  -lopenblas 
#-qmkl=cluster
#EXTRA= -DMKL
DIR="./bin/"
COND= -DMPI -O0 -fbacktrace -fcheck=all -g
#COND= -O0 -fbacktrace -fcheck=all -g 
