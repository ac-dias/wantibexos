!
!   bcast_input_read.F90
!   
!
!   Created by Alexandre Rocha on 01/06/26.
!   Copyright 2026 ___ORGANIZATIONNAME___. All rights reserved.
!

subroutine bcast_input_read()

use input_variables

! Integers
call MPI_BCAST(nthreads,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ngrid,3,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(nc,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(nv,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(excwf0,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(excwff,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(nocpf,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(power,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(nmu,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(nsteps,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(nomega,1,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)

! Reals
call MPI_BCAST(rk,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(numdos,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ebse0,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ebsef,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(numbse,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(sme,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(cshift,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(smeboltz,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(mshift,3,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ediel,3,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ktol,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(ez,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(w,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(lc,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(r0,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(st,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(phavg,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(temp,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(ctemp,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(tmax,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(eg,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(egd,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(egs,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ebgs,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(dk,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(fermishift,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(exc,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(mag,3,MPI_REAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(mu0,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(muf,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(btemp,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(elft,3,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(hlft,3,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(klat,3,MPI_REAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(ni,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ns,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(omegamax,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(ezgw,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(wgw,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(lcgw,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(r0gw,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(edielgw,3,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ktolgw,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(smegw,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ifactor,1,MPI_REAL,root,MPI_COMM_WORLD,ierr)

! Characters
call MPI_BCAST(params,70,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(orbw,70,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(kpaths,70,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(kpathsbse,70,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(outputfolder,70,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(calcparms,70,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(meshtype,70,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(meshgen,70,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(coultype,5,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(sysdim,5,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(dft,1,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ta,2,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(ses,6,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(emfile,70,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(bsealgo,12,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(coultypegw,5,MPI_CHARACTER,root,MPI_COMM_WORLD,ierr)

! Logicals
call MPI_BCAST(bandscalc,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(doscalc,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(bse,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(bsepol,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(bsekpath,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(spdiel,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(spdielpol,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(sppolbz,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(spec,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(berryk,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(berrybz,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(pponly,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(bsewf,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(tmcoef,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(dtfull,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(cpol,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(pce,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(renorm,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(emt,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(spintxt,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(boltz,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(lowdin,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)

call MPI_BCAST(gwmesh,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(gwbnd,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(gwmeshuse,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(gwbnduse,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(gwbsebnd,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
call MPI_BCAST(selfxonly,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)

use bcast_input_read

subroutine bcast_hamil

use hamiltonian_input_variables

! Scalars
call MPI_BCAST(w90basis,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(ntype,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(nocp,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

call MPI_BCAST(efermi,    1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(scs,       1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)

call MPI_BCAST(nvec,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

call MPI_BCAST(rlat,      9, MPI_REAL,    0, MPI_COMM_WORLD, ierr)

call MPI_BCAST(systype, len(4), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

! Spin-polarized parameters
call MPI_BCAST(w90basisu, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(w90basisd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

call MPI_BCAST(nvecu,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
call MPI_BCAST(nvecd,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

! Allocated arrays
call MPI_BCAST(ffactor, size(nvec), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

call MPI_BCAST(rvec, &
     nvec, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

call MPI_BCAST(hopmatrices, &
      nvec*w90basis*w90basis, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

call MPI_BCAST(ihopmatrices, &
      nvec*w90basis*w90basis, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

call MPI_BCAST(ovp, &
      nvec*w90basis*w90basis, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

end subroutine bcasst_hamil
