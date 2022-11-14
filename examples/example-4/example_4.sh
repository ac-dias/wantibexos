#!/usr/bin/env bash

temperature="50 100 200 300"
temperature=($temperature)

nthreads=10
wtbexec="wantibexos_folder/bin/wtb.x"

mkdir ./out

for (( i = 0; i<${#temperature[@]}; i++ )); do	

temp="${temperature[$i]}.dat"

mkdir ./out/$temp

cat > input_mos2_$temp.dat << EOF
NTHREADS= $nthreads
SYSDIM= "2D"
DFT= "V"

OUTPUT= "./out/$temp/"
CALC_DATA= "./out/$temp/"
PARAMS_FILE= "tb_mos2.dat"                                                             
                                                      
 
MESH_TYPE= "RK2D"
RK= 120

BSET= T
SPEC= T
DTDIAG= T

COULOMB_POT= V2DT2
NBANDSC= 2
NBANDSV= 2
LC= 8.00

CSHIFT= 0.08
ENSPECI= 0.0
ENSPECF= 4.0

TA= DE
ST= 0.036
PHAVG= 225
TEMP= $temp
	
EOF



$wtbexec < input_mos2_$temp.dat


done

