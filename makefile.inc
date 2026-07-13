FOR=/opt/local/bin/mpifort-mpich-mp
OMP=-fopenmp 
#-shared-intel 
LIBS= -L/opt/local/lib  -lopenblas  -lscalapack 
#-qmkl=cluster
#EXTRA= -DMKL
DIR="./bin/"
COND= -DMPI -O0 -fbacktrace -fcheck=all -g 
#-fpp -O2 -DMPI
