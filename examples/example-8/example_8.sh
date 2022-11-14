#!/usr/bin/env bash

nthreads=10
wtbexec="wantibexos_folder/bin/wtb.x"



#to calculate PCE it's necessary to get absorption spectrum, from this data you get EGD and EBGS, 
# EGS you get from exciton band structure and EG from band structure

#first calculation - exciton band structure, this calculation is very long

mkdir ./out-exc-bs/

cat > gaas-kpoints-bse.dat << EOF
12
5
0.50000     0.50000     0.50000 !L  
0.00000     0.00000     0.00000 !G
0.00000     0.00000     0.00000 !G 
0.50000     0.00000     0.50000 !X
0.50000     0.00000     0.50000 !X  
0.50000     0.25000     0.75000 !W
0.50000     0.25000     0.75000 !W
0.50000     0.50000     0.50000 !L
0.50000     0.50000     0.50000 !L  
0.75000     0.37500     0.37500 !K
0.75000     0.37500     0.37500 !K  
0.00000     0.00000     0.00000 !G
EOF

cat > bse_bands_input.dat << EOF

NTHREADS= 1
SYSDIM= "3D"
DFT= "V"

OUTPUT= "./out-exc-bs/"
CALC_DATA= "./out-exc-bs/"
PARAMS_FILE= "tb_gaas.dat"
KPATH_BSE= "gaas-kpoints-bse.dat"                                                       
 

MESH_TYPE= "RK3D"
RK= 80

BSE_BND= T

COULOMB_POT= V3D
NBANDSC= 1
NBANDSV= 3


EOF

$wtbexec < bse_bands_input.dat

#second calculation - band structure and absorption spectrum

mkdir ./out/

cat > gaas-kpoints.dat << EOF
12
50
0.50000     0.50000     0.50000 !L  
0.00000     0.00000     0.00000 !G
0.00000     0.00000     0.00000 !G 
0.50000     0.00000     0.50000 !X
0.50000     0.00000     0.50000 !X  
0.50000     0.25000     0.75000 !W
0.50000     0.25000     0.75000 !W
0.50000     0.50000     0.50000 !L
0.50000     0.50000     0.50000 !L  
0.75000     0.37500     0.37500 !K
0.75000     0.37500     0.37500 !K  
0.00000     0.00000     0.00000 !G
EOF

cat > bse_opt_input.dat << EOF

NTHREADS= 1
SYSDIM= "3D"
DFT= "V"

OUTPUT= "./out/"
CALC_DATA= "./out/"
PARAMS_FILE= "tb_gaas.dat"                                                             
KPATH_FILE= "gaas-kpoints.dat"                                                           
                                                       
 

MESH_TYPE= "RK3D"
RK= 80

BANDS= T
BSE= T
SPEC= T
DTDIAG= T

COULOMB_POT= V3D
NBANDSC= 1
NBANDSV= 3

CSHIFT= 0.08
ENSPECI= 0.0
ENSPECF= 4.0

EOF

$wtbexec < bse_opt_input.dat

#third calculation - PCE

cat > pce_input.dat << EOF

NTHREADS= 1
SYSDIM= "3D"
DFT= "V"

OUTPUT= "./out/"
CALC_DATA= "./out/"
PARAMS_FILE= "tb_gaas.dat"
                                                             
                                                   
MESH_TYPE= "RK3D"
RK= 80

PCE= T
SES= AM15G
CTEMP= 298.15
THMAX= 1E-6
EG= 0.89
EGD= 0.89
EGS= 0.89
EBGS= 0.89

COULOMB_POT= V3D
NBANDSC= 1
NBANDSV= 3

CSHIFT= 0.08
ENSPECI= 0.0
ENSPECF= 4.0

EOF

$wtbexec < pce_input.dat




