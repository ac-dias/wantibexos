#!/usr/bin/env bash

nthreads=10
wtbexec="wantibexos_folder/bin/wtb.x"

cat > input_mos2_opt.dat << EOF
NTHREADS= $nthreads
SYSDIM= "2D"
DFT= "V"

OUTPUT= "./out/"
CALC_DATA= "./out/"
PARAMS_FILE= "tb_mos2.dat"                                                             
                                                           
                                                       
MESH_TYPE= "RK2D"
RK= 120

BSE= T
SPEC= T
DTDIAG= T

COULOMB_POT= V2DT2
NBANDSC= 2
NBANDSV= 2
LC= 8.00

CSHIFT= 0.08
ENSPECI= 0.0
ENSPECF= 4.0
	
EOF


$wtbexec < input_mos2_opt.dat





