include makefile.inc

SUB_OBJ=\
	./build/subroutines/*.o \


SUB_PROGRAM=\
	./build/subprograms/*.o \

																			

all : main pp
	
pp :
	$(FOR) ./utils/nc_nv_finder.F90 -o $(DIR)nc_nv_finder.x $(OMP) $(LIBS)
	$(FOR) ./utils/param_gen.F90  -o $(DIR)param_gen.x
	$(FOR) ./utils/param_gen_vasp.F90  -o $(DIR)param_gen_vasp.x
	$(FOR) ./utils/absorbance.F90  -o $(DIR)absorbance.x	
	$(FOR) ./utils/slme/pce-code.f90 ./utils/slme/pce-subs.f90  -o $(DIR)pce.x
	$(FOR) ./utils/huckel2wtb/src/overlaps_jc.f90 ./utils/huckel2wtb/src/diagonalize.f90 ./utils/huckel2wtb/src/Huckel_TB.f90 -o $(DIR)huckel2wtb.x $(LIBS)
	cp ./utils/*.py  $(DIR)
	rm ./*.mod


main :  subprograms
	$(FOR) ./src/wtb_main.F90 $(SUB_OBJ) $(SUB_PROGRAM) -o ./build/wtb.x $(LIBS) $(OMP) $(COND) $(EXTRA) 
	cp ./build/wtb.x $(DIR)wtb.x

subprograms: lib
	[ -d ./build/subprograms ] || mkdir ./build/subprograms
	$(FOR) -c ./src/subprograms/bands-kpath-tool.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/berry_curvature_bz-tool.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/berry_curvature_kpath-tool.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/bse_diel-tool.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/bse_diel-tool-pol.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/bse_kpath-tool.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/bse_kpath-tool-temp.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/bse_solver-tool-diel.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/bse_solver-tool-diel-temp.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/diel-pp-bse.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/diel-pp-bse-pol.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/tdos-tool-stxt.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/diel-pp.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/diel-pp-pol.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/efmass.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/exciton_lifetime.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/sp_diel-tool.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/sp_diel-tool-pol.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/sp_opt_bz-tool.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/sp_solver-tool-diel.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/boltzmann_transport.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a
	$(FOR) -c ./src/subprograms/emission_PL.F90 $(OMP) $(LIBS) $(COND) $(EXTRA) -L./build/libwtb.a	
	mv *.o ./build/subprograms

lib: subroutines
	ar rcs ./build/libwtb.a $(SUB_OBJ)
	 		
subroutines :
	[ -d ./build/subroutines ] || mkdir ./build/subroutines
	$(FOR) -c ./src/subroutines/berry_curvature_subs.F90 
	$(FOR) -c ./src/subroutines/boltzmann_subs.F90 
	$(FOR) -c ./src/subroutines/bse_subs.F90 
	$(FOR) -c ./src/subroutines/bse_subs_kpath.F90 
	$(FOR) -c ./src/subroutines/bse_subs_temp.F90 
	$(FOR) -c ./src/subroutines/special_funct.F90
	$(FOR) -c ./src/subroutines/ei_spec_funct.F90
	$(FOR) -c ./src/subroutines/coulomb_pot.F90 
	$(FOR) -c ./src/subroutines/diel-pp-subs.F90 
	$(FOR) -c ./src/subroutines/dos_subs.F90 
	$(FOR) -c ./src/subroutines/efmass-subs.F90 
	$(FOR) -c ./src/subroutines/general_subs.F90 $(LIBS)
	$(FOR) -c ./src/subroutines/hamiltonians.F90 
	$(FOR) -c ./src/subroutines/hamiltonian_tb.F90 $(LIBS) $(OMP)
	$(FOR) -c ./src/subroutines/mhkpack_subs.F90 
	$(FOR) -c ./src/subroutines/module_input_read.F90 
	$(FOR) -c ./src/subroutines/optics.F90  
	$(FOR) -c ./src/subroutines/spin_txt_subs.F90 
	$(FOR) -c ./src/subroutines/emission_subs.F90	
	mv *.o ./build/subroutines

clean:
	rm -r ./build/subroutines
	rm -r ./build/subprograms
	rm -r ./build/wtb.x
	rm -r ./build/libwtb.a
	rm ./*.mod
	rm ./bin/*.py
	rm ./bin/*.x

.PHONY : all clean 
