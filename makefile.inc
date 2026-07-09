FOR=mpif90
OMP=-qopenmp -shared-intel 
LIBS=-qmkl=cluster
EXTRA= -DMKL
DIR="./bin/"
COND=-fpp -O2 -DMPI
