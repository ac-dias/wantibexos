#!/usr/bin/env python3
import sisl
import numpy as np
#import matplotlib.pyplot as plt
import os
import sys

############################################################################
def calctype(data):
    
    if data == 'unpolarized':
     outdata = 'NP'
        
    elif data == 'polarized':
     outdata = 'SP'
        
    elif data == 'non-colinear':
     outdata = 'SOC'
        
    elif data == 'spin-orbit':
     outdata = 'SOC'
            
    return outdata   


def ncases(var):
    return "%8f"%var

###############################################################################

#inputfdf=  './teste-honpas/mos2.fdf'

inputfdf=  sys.argv[1]  
#fermi= sys.argv[2]
fermi= 0.00

#os.system("cp ./teste-honpas/mos2.out ./teste-honpas/run.out")
geom = sisl.get_sile(inputfdf).read_geometry()
tshs = sisl.get_sile(inputfdf).read_hamiltonian(geometry=geom) #pegar hamiltoniano do siesta

#fermi = sisl.get_sile(folder).read_fermi_level()

scs = 0.0000 #scissors operator

nbasis= tshs.no
ncell=  tshs.nsc[0]*tshs.nsc[1]*tshs.nsc[2]
sptype=str(tshs.spin)

sptype2=sptype[5:-1]

if sptype2 == 'unpolarized' :

 f = open("system-info-NP.txt", "a")
 print(tshs, file=f)
 f.close()

 nbasis= tshs.no
 ncell=  tshs.nsc[0]*tshs.nsc[1]*tshs.nsc[2]
 sptype=str(tshs.spin)


 f = open("tb-NP.dat", "a")

 print(calctype(sptype2),file=f,flush=True)
 print(ncases(scs),file=f,flush=True)
 #os.system("grep 'siesta:         Fermi =' run.out | awk '{print $4;}' >>  honpas_tb-NP.dat")
 print(fermi,file=f)
 print(ncases(tshs.cell[0,0]),'',ncases(tshs.cell[0,1]),'',ncases(tshs.cell[0,2]),file=f)
 print(ncases(tshs.cell[1,0]),'',ncases(tshs.cell[1,1]),'',ncases(tshs.cell[1,2]),file=f)
 print(ncases(tshs.cell[2,0]),'',ncases(tshs.cell[2,1]),'',ncases(tshs.cell[2,2]),file=f)
 print(nbasis,file=f)
 print(ncell,file=f)
 print('#rcell x',' ','rcell y',' ','rcell z',' ','i',' ','j',' ','ReH',' ','ImH',' ','S',file=f)

 for i in range(0,ncell):
  for j in range(0,nbasis): 
   for k in range(0,nbasis):
  	
    a = tshs.geometry.o2sc(k+i*nbasis)[0]
    b = tshs.geometry.o2sc(k+i*nbasis)[1]
    c = tshs.geometry.o2sc(k+i*nbasis)[2]
		
    d = tshs.H[j,k+i*nbasis]
    e = tshs.S[j,k+i*nbasis]
	
    print(ncases(a),' ',ncases(b), ' ',ncases(c), ' ',j+1,' ',k+1,' ',ncases(d),' ',ncases(0.00),' ',ncases(e),file=f)


 f.close()

 f = open("basis_set-NP", "a")

 print('bindex aspecie ax ay az l m spin',file=f)

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
    xyz = tshs.geometry.xyz[ia, :]
    for orbital in atom:
        print(io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',0,file=f)
        io += 1

 f.close()


####################################################################

if sptype2 == 'polarized' :

 f = open("system-info-sp.txt", "a")
 print(tshs, file=f)
 f.close()

 nbasis= tshs.no
 ncell=  tshs.nsc[0]*tshs.nsc[1]*tshs.nsc[2]
 sptype=str(tshs.spin)


 f = open("basis_set-sp", "a")

 print('#bindex aspecie ax ay az l m spin',file=f)

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
    xyz = tshs.geometry.xyz[ia, :]
    for orbital in atom:
        print(io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',1,file=f)
        io += 1

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
    xyz = tshs.geometry.xyz[ia, :]
    for orbital in atom:
        print(nbasis+io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',-1,file=f)
        io += 1

 f.close()


 f = open("tb-sp.dat", "a")

 print(calctype(sptype2),file=f,flush=True)
 print(ncases(scs),file=f,flush=True)
 #os.system("grep 'siesta:         Fermi =' run.out | awk '{print $4;}' >>  honpas_tb-sp.dat")
 print(fermi,file=f)
 print(ncases(tshs.cell[0,0]),'',ncases(tshs.cell[0,1]),'',ncases(tshs.cell[0,2]),file=f)
 print(ncases(tshs.cell[1,0]),'',ncases(tshs.cell[1,1]),'',ncases(tshs.cell[1,2]),file=f)
 print(ncases(tshs.cell[2,0]),'',ncases(tshs.cell[2,1]),'',ncases(tshs.cell[2,2]),file=f)
 print(2*nbasis,file=f)
 print(ncell,file=f)
 print('#rcell x',' ','rcell y',' ','rcell z',' ','i',' ','j',' ','ReH',' ','ImH',' ','S',file=f)
#np.empty aloca todos os arrays vazios
#alocando todos os arrays com 0

 S = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))
 reH = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))
 imH = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))

 a = np.zeros((ncell+1,(2*nbasis)+2))
 b = np.zeros((ncell+1,(2*nbasis)+2))
 c = np.zeros((ncell+1,(2*nbasis)+2))


 for i in range(0,ncell):
  for j in range(0,nbasis): 
   for k in range(0,nbasis):

#parte up-up
  
    reH[i+1,j+1,k+1] = tshs[j,k+i*nbasis][0]
    #imH[i+1,j+1,k+1] = tshs[j,k+i*nbasis][4]
    S[i+1,j+1,k+1] = tshs.S[j,k+i*nbasis]

    

#parte up-dn

    #reH[i+1,(j+1),(k+1)+nbasis] = tshs[j,k+i*nbasis][2]
    #imH[i+1,(j+1),(k+1)+nbasis] = tshs[j,k+i*nbasis][3]
    #S[i+1,(j+1),(k+1)+nbasis] = 0.0

#parte dn-up

    #reH[i+1,(j+1)+nbasis,(k+1)] = tshs[j,k+i*nbasis][6]
    #imH[i+1,(j+1)+nbasis,(k+1)] = tshs[j,k+i*nbasis][7]
    #S[i+1,(j+1)+nbasis,(k+1)] = 0.0

#parte dn-dn

    reH[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs[j,k+i*nbasis][1]
    #imH[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs[j,k+i*nbasis][5]
    S[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs.S[j,k+i*nbasis]
    
#escrevendo a localização dos atomos da base

 for i in range(0,ncell):
  for k in range(0,nbasis): 

   a[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[0]
   b[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[1]
   c[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[2]
	
   a[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[0]
   b[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[1]
   c[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[2]
	

 for i in range(1,ncell+1):
  for j in range(1,(2*nbasis)+1): 
   for k in range(1,(2*nbasis)+1):
  	
#   a = tshs.geometry.o2sc(k+i*nbasis)[0]
#   b = tshs.geometry.o2sc(k+i*nbasis)[1]
#   c = tshs.geometry.o2sc(k+i*nbasis)[2]
		
#   d = tshs.H[j,k+i*nbasis]
#   e = tshs.S[j,k+i*nbasis]
	
    print(ncases(a[i,k]),' ',ncases(b[i,k]), ' ',ncases(c[i,k]), ' ',j,' ',k,' ',ncases(reH[i,j,k]),' ',ncases(imH[i,j,k]),' ',ncases(S[i,j,k]),file=f)


 f.close()

####################################################################

if sptype2 == 'non-colinear' :

 f = open("system-info-nc.txt", "a")
 print(tshs, file=f)
 f.close()

 nbasis= tshs.no
 ncell=  tshs.nsc[0]*tshs.nsc[1]*tshs.nsc[2]
 sptype=str(tshs.spin)


 f = open("basis_set-nc", "a")

 print('#bindex aspecie ax ay az l m spin',file=f)

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
    xyz = tshs.geometry.xyz[ia, :]
    for orbital in atom:
        print(io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',1,file=f)
        io += 1

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
    xyz = tshs.geometry.xyz[ia, :]
    for orbital in atom:
        print(nbasis+io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',-1,file=f)
        io += 1

 f.close()


 f = open("tb-nc.dat", "a")

 print(calctype(sptype2),file=f,flush=True)
 print(ncases(scs),file=f,flush=True)
 #os.system("grep 'siesta:         Fermi =' run.out | awk '{print $4;}' >>  honpas_tb-nc.dat")
 print(fermi,file=f)
 print(ncases(tshs.cell[0,0]),'',ncases(tshs.cell[0,1]),'',ncases(tshs.cell[0,2]),file=f)
 print(ncases(tshs.cell[1,0]),'',ncases(tshs.cell[1,1]),'',ncases(tshs.cell[1,2]),file=f)
 print(ncases(tshs.cell[2,0]),'',ncases(tshs.cell[2,1]),'',ncases(tshs.cell[2,2]),file=f)
 print(2*nbasis,file=f)
 print(ncell,file=f)
 print('#rcell x',' ','rcell y',' ','rcell z',' ','i',' ','j',' ','ReH',' ','ImH',' ','S',file=f)

#np.empty aloca todos os arrays vazios
#alocando todos os arrays com 0

 S = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))
 reH = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))
 imH = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))

 a = np.zeros((ncell+1,(2*nbasis)+2))
 b = np.zeros((ncell+1,(2*nbasis)+2))
 c = np.zeros((ncell+1,(2*nbasis)+2))


 for i in range(0,ncell):
  for j in range(0,nbasis): 
   for k in range(0,nbasis):

#parte up-up
  
     reH[i+1,j+1,k+1] = tshs[j,k+i*nbasis][0]
    #imH[i+1,j+1,k+1] = tshs[j,k+i*nbasis][4]
     S[i+1,j+1,k+1] = tshs.S[j,k+i*nbasis]

    

#parte up-dn

     reH[i+1,(j+1),(k+1)+nbasis] = tshs[j,k+i*nbasis][2]
    #imH[i+1,(j+1),(k+1)+nbasis] = tshs[j,k+i*nbasis][3]
     S[i+1,(j+1),(k+1)+nbasis] = 0.0

#parte dn-up

     reH[i+1,(j+1)+nbasis,(k+1)] = tshs[j,k+i*nbasis][3]
    #imH[i+1,(j+1)+nbasis,(k+1)] = tshs[j,k+i*nbasis][7]
     S[i+1,(j+1)+nbasis,(k+1)] = 0.0

#parte dn-dn

     reH[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs[j,k+i*nbasis][1]
    #imH[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs[j,k+i*nbasis][5]
     S[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs.S[j,k+i*nbasis]
    
#escrevendo a localização dos atomos da base

 for i in range(0,ncell):
  for k in range(0,nbasis): 

   a[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[0]
   b[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[1]
   c[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[2]
	
   a[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[0]
   b[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[1]
   c[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[2]
	

 for i in range(1,ncell+1):
  for j in range(1,(2*nbasis)+1): 
   for k in range(1,(2*nbasis)+1):
  	
#   a = tshs.geometry.o2sc(k+i*nbasis)[0]
#   b = tshs.geometry.o2sc(k+i*nbasis)[1]
#   c = tshs.geometry.o2sc(k+i*nbasis)[2]
		
#   d = tshs.H[j,k+i*nbasis]
#   e = tshs.S[j,k+i*nbasis]
	
    print(ncases(a[i,k]),' ',ncases(b[i,k]), ' ',ncases(c[i,k]), ' ',j,' ',k,' ',ncases(reH[i,j,k]),' ',ncases(imH[i,j,k]),' ',ncases(S[i,j,k]),file=f)


 f.close()

####################################################################

if sptype2 == 'spin-orbit' :

 f = open("system-info-soc.txt", "a")
 print(tshs, file=f)
 f.close()
 
 nbasis= tshs.no
 ncell=  tshs.nsc[0]*tshs.nsc[1]*tshs.nsc[2]
 sptype=str(tshs.spin)

 f = open("basis_set-soc", "a")

 print('#bindex aspecie ax ay az l m spin',file=f)

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
     xyz = tshs.geometry.xyz[ia, :]
     for orbital in atom:
         print(io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',1,file=f)
         io += 1

 io = 0
 for ia, atom in enumerate(tshs.geometry.atoms):
     xyz = tshs.geometry.xyz[ia, :]
     for orbital in atom:
         print(nbasis+io+1,' ', atom.tag,' ', ncases(xyz[0]),' ', ncases(xyz[1]),' ', ncases(xyz[2]),' ', orbital.l,' ', orbital.m,' ',-1,file=f)
         io += 1

 f.close()


 f = open("tb-soc.dat", "a")

 print(calctype(sptype2),file=f,flush=True)
 print(ncases(scs),file=f,flush=True)
 #os.system("grep 'siesta:         Fermi =' run.out | awk '{print $4;}' >>  honpas_tb-soc.dat")
 print(fermi,file=f)
 print(ncases(tshs.cell[0,0]),'',ncases(tshs.cell[0,1]),'',ncases(tshs.cell[0,2]),file=f)
 print(ncases(tshs.cell[1,0]),'',ncases(tshs.cell[1,1]),'',ncases(tshs.cell[1,2]),file=f)
 print(ncases(tshs.cell[2,0]),'',ncases(tshs.cell[2,1]),'',ncases(tshs.cell[2,2]),file=f)
 print(2*nbasis,file=f)
 print(ncell,file=f)
 print('#rcell x',' ','rcell y',' ','rcell z',' ','i',' ','j',' ','ReH',' ','ImH',' ','S',file=f)

#np.empty aloca todos os arrays vazios
#alocando todos os arrays com 0

 S = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))
 reH = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))
 imH = np.zeros((ncell+1,(2*nbasis)+2,(2*nbasis)+2))

 a = np.zeros((ncell+1,(2*nbasis)+2))
 b = np.zeros((ncell+1,(2*nbasis)+2))
 c = np.zeros((ncell+1,(2*nbasis)+2))


 for i in range(0,ncell):
  for j in range(0,nbasis): 
   for k in range(0,nbasis):

#parte up-up
  
    reH[i+1,j+1,k+1] = tshs[j,k+i*nbasis][0]
    imH[i+1,j+1,k+1] = tshs[j,k+i*nbasis][4]
    S[i+1,j+1,k+1] = tshs.S[j,k+i*nbasis]

    

#parte up-dn

    reH[i+1,(j+1),(k+1)+nbasis] = tshs[j,k+i*nbasis][2]
    imH[i+1,(j+1),(k+1)+nbasis] = tshs[j,k+i*nbasis][3]
    S[i+1,(j+1),(k+1)+nbasis] = 0.0

#parte dn-up

    reH[i+1,(j+1)+nbasis,(k+1)] = tshs[j,k+i*nbasis][6]
    imH[i+1,(j+1)+nbasis,(k+1)] = tshs[j,k+i*nbasis][7]
    S[i+1,(j+1)+nbasis,(k+1)] = 0.0

#parte dn-dn

    reH[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs[j,k+i*nbasis][1]
    imH[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs[j,k+i*nbasis][5]
    S[i+1,(j+1)+nbasis,(k+1)+nbasis] = tshs.S[j,k+i*nbasis]
    
#escrevendo a localização dos atomos da base

 for i in range(0,ncell):
  for k in range(0,nbasis): 

   a[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[0]
   b[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[1]
   c[i+1,k+1] = tshs.geometry.o2sc(k+i*nbasis)[2]
	
   a[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[0]
   b[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[1]
   c[i+1,(k+1)+nbasis] = tshs.geometry.o2sc(k+i*nbasis)[2]
	

 for i in range(1,ncell+1):
  for j in range(1,(2*nbasis)+1): 
   for k in range(1,(2*nbasis)+1):
  	
#   a = tshs.geometry.o2sc(k+i*nbasis)[0]
#   b = tshs.geometry.o2sc(k+i*nbasis)[1]
#   c = tshs.geometry.o2sc(k+i*nbasis)[2]
		
#   d = tshs.H[j,k+i*nbasis]
#   e = tshs.S[j,k+i*nbasis]
	
    print(ncases(a[i,k]),' ',(b[i,k]), ' ',ncases(c[i,k]), ' ',j,' ',k,' ',ncases(reH[i,j,k]),' ',ncases(imH[i,j,k]),' ',ncases(S[i,j,k]),file=f)


 f.close()

#os.system("rm run.out")	
####################################################################	

