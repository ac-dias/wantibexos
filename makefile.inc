FOR=mpif90
OMP=-qopenmp -shared-intel 
LIBS=-qmkl
EXTRA= -DMKL
DIR="./bin/"
COND=-fpp -O2
