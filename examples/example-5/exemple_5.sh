#!/usr/bin/env bash

nthreads=10
wtbexec="wantibexos_folder/bin/wtb.x"

mkdir ./out/

cat tmd-kpoints-bse.dat << EOF
6
5
0.0 0.0 0.0               !ponto Gamma
0.6667  -0.3333  0.0      !ponto K
0.6667  -0.3333  0.0      !ponto K
0.3333   0.3333 0.0       !ponto K'
0.3333   0.3333 0.0       !ponto K'
0.0 0.0 0.0               !ponto Gamma
EOF

cat > input_mos2_bseb.dat << EOF
NTHREADS= $nthreads
SYSDIM= "2D"
DFT= "V"

OUTPUT= "./out/"
CALC_DATA= "./out/"
PARAMS_FILE= "tb_mos2.dat"                                                                                                                       
KPATH_BSE= "tmd-kpoints-bse.dat"                                                       
 
MESH_TYPE= "RK2D"
RK= 120

BSE_BND= T


COULOMB_POT= V2DT2
NBANDSC= 2
NBANDSV= 2
LC= 8.00

	
EOF



$wtbexec < input_mos2_bseb.dat





