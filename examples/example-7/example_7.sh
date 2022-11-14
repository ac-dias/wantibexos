#!/usr/bin/env bash

nthreads=10
wtbexec="wantibexos_folder/bin/wtb.x"

mkdir ./out/

cat > tmd-kpoints.dat << EOF
6
50
0.0 0.0 0.0               !Gamma
0.6667  -0.3333  0.0      !K
0.6667  -0.3333  0.0      !K
0.3333   0.3333 0.0       !K'
0.3333   0.3333 0.0       !K'
0.0 0.0 0.0               !Gamma
EOF

cat > input.dat << EOF
NTHREADS= $nthreads
SYSDIM= "2D"
DFT= "V"

OUTPUT= "./out/"
CALC_DATA= "./out/"
PARAMS_FILE= "tb_mos2.dat"                                                             
KPATH_FILE= "tmd-kpoints.dat"                                                           
                                                      
BERRY= T
	
EOF

$wtbexec < input.dat





