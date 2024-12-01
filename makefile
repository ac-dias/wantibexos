#This part is for ifort
FOR=ifort
OMP=-qopenmp
LIBS=-qmkl
EXTRA= -shared-intel
COND=-fpp
#Uncomment to add optimizations (-O2 is ifort's default behavior)
#FFLAGS= -O2 -xHOST 

#Uncomment this part to use GNU fortran compiler
#FOR=gfortran
#OMP=-fopenmp
#LIBS= -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl  -m64  -I"${MKLROOT}/include"
#COND=-x f95-cpp-input
#FFLAGS=-O2 -march=native


DIR="./bin/"

all :	main pp

main :
	$(FOR) $(FFLAGS) $(COND) wtb_main.f90 -o $(DIR)wtb.x $(LIBS) $(OMP)
	$(FOR) $(FFLAGS) $(COND) wtb_main.f90 -o $(DIR)wtbf.x $(LIBS) $(OMP) -DFAST1

pp :
	$(FOR) $(FFLAGS) ./utils/nc_nv_finder.f90 -o $(DIR)nc_nv_finder.x $(OMP) $(LIBS)
	$(FOR) $(FFLAGS) ./utils/param_gen.f90  -o $(DIR)param_gen.x
	$(FOR) $(FFLAGS) ./utils/param_gen_vasp.f90  -o $(DIR)param_gen_vasp.x
	$(FOR) $(FFLAGS) ./utils/dirgapf.f90  -o $(DIR)dirgapf.x

test :
	$(FOR) $(FFLAGS) $(COND) wtb_main.f90 -o $(DIR)wtb.x $(LIBS) $(OMP)
