!############subroutines################
INCLUDE "./subroutines/berry_curvature_subs.f90"
INCLUDE "./subroutines/bse_subs.f90"
INCLUDE "./subroutines/bse_subs_kpath.f90"
INCLUDE "./subroutines/coulomb_pot.f90"
INCLUDE "./subroutines/diel-pp-subs.f90"
INCLUDE "./subroutines/dos_subs.f90"
INCLUDE "./subroutines/general_subs.f90"
INCLUDE "./subroutines/mhkpack_subs.f90"
INCLUDE "./subroutines/module_input_read.f90"
INCLUDE "./subroutines/optics.f90"
INCLUDE "./subroutines/hamiltonian_tb.f90"
INCLUDE "./subroutines/bse_subs_temp.f90"
INCLUDE "./subroutines/pce-subs.f90"

!############subprograms##################

INCLUDE "./subprograms/bands-kpath-tool.f90"

INCLUDE "./subprograms/berry_curvature_bz-tool.f90"
INCLUDE "./subprograms/berry_curvature_kpath-tool.f90"

INCLUDE "./subprograms/bse_diel-tool.f90"
INCLUDE "./subprograms/bse_diel-tool-pol.f90"
INCLUDE "./subprograms/bse_kpath-tool.f90"
INCLUDE "./subprograms/bse_kpath-tool-temp.f90"

#ifdef FAST1
INCLUDE "./subprograms/bse_solver-tool-diel-f.f90"
#else
INCLUDE "./subprograms/bse_solver-tool-diel.f90"
#endif

#ifdef FAST1
INCLUDE "./subprograms/bse_solver-tool-diel-temp-f.f90"
#else
INCLUDE "./subprograms/bse_solver-tool-diel-temp.f90"
#endif

!INCLUDE "./subprograms/bse_solver-tool-pol.f90"

INCLUDE "./subprograms/diel-pp.f90"
INCLUDE "./subprograms/diel-pp-pol.f90"
INCLUDE "./subprograms/diel-pp-bse.f90"
INCLUDE "./subprograms/diel-pp-bse-pol.f90"
INCLUDE "./subprograms/exciton_lifetime.f90"

INCLUDE "./subprograms/sp_diel-tool.f90"
INCLUDE "./subprograms/sp_diel-tool-pol.f90"
INCLUDE "./subprograms/sp_opt_bz-tool.f90"
INCLUDE "./subprograms/sp_solver-tool-diel.f90"
!INCLUDE "./subprograms/sp_solver-tool-pol.f90"


INCLUDE "./subprograms/tdos-tool.f90"
INCLUDE "./subprograms/pce-code.f90"


!ifort tbt_main.f90 -o tbt.x -mkl -qopenmp -shared-intel
program main

	use input_variables
	use hamiltonian_input_variables
	implicit none
	
	integer :: erro
	double precision:: t0,tf
	integer,dimension(8) :: values,values2

	double precision,dimension(3,3) :: rlatv
	double precision :: flag
	character(len=70) ::cflag
	

	call input_read

	OPEN(UNIT=2055, FILE= params,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wannier90 hr input file"	

	read(2055,*) systype
	!read(2055,*) ediel(2)
	read(2055,*) scs
	read(2055,*) efermi
	!read(2055,*) lc


	read(2055,*) rlatv(1,1),rlatv(1,2),rlatv(1,3)
	read(2055,*) rlatv(2,1),rlatv(2,2),rlatv(2,3)
	read(2055,*) rlatv(3,1),rlatv(3,2),rlatv(3,3)


	if (meshtype .eq. "RK3D") then

		call rkmesh(rk,rlatv,ngrid)

	else if (meshtype .eq. "RK2D") then

		call rkmesh2D(rk,rlatv,ngrid)
	else

		continue
	end if

	close(2055)

	OPEN(UNIT=2077, FILE= trim(calcparms)//"log.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log output file "

	call param_out(2077,nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,bandscalc,doscalc,bse,bsepol,bsekpath,spec,&
		     spdiel,spdielpol,sppolbz,berryk,berrybz,pponly,bsewf,excwf0,excwff,&
		     tmcoef,ez,w,lc,r0,sysdim,dtfull,cpol,cshift,dft,bset,bsetbnd,&
		     st,phavg,temp,ta,pce,ses,ctemp,tmax,eg,egd,egs,ebgs,renorm)

	call cpu_time(t0)
	call date_and_time(VALUES=values)

	write(2077,*)
	write(2077,*) 'Begin','  ','day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(2077,*)

		!nthreads,outputfolder,calcparms,ngrid,nc,nv,edos0,edosf,numdos, &
		!ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		!exc,mshift,coultype
		!systype,efermi,edielh,lc,scs,rlat,w90basis
		!hopmatrices,ihopmatrices,rvec,ffactor,nvec

	!calculo estrutura eletronica
	if (pponly) go to 131

	if (bandscalc) then

	 call bandstool(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,rk,dft)
	 write(2077,*) "Band Structure finished"
     call flush(2077)
	end if


	if (doscalc) then

	 call dostool(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,rk,meshtype,dft)

	 write(2077,*) "DOS finished"
     call flush(2077)
	end if

	!calculo optica
	
	if (bse .and. bset) then
	
	write(2077,*) "ERROR: Choose BSE=T or BSET=T not both"
	stop
	
	elseif (bse) then
	
	call bsesolver(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,cshift,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w,r0,lc,rk,meshtype,bsewf,excwf0,excwff,dtfull,cpol,tmcoef)

	 write(2077,*) "BSE and Single particle dielectric calculation finished"
     call flush(2077)
	
	elseif (bset) then
	
	call bsesolvertemp(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,cshift,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w,r0,lc,rk,meshtype,bsewf,excwf0,excwff,dtfull,&
		     cpol,tmcoef,st,phavg,ta,temp)
		     
	 write(2077,*) "BSE and Single particle dielectric calculation, with temperature, finished"	
     call flush(2077)	     	
	
	end if	
	

	if (bsekpath .and. bsetbnd) then
	
	write(2077,*) "ERROR: Choose BSE_BND=T or BSET_BND=T not both"
	stop
	
	elseif (bsekpath) then
	
	 call bsebnds(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w,r0,lc,rk,meshtype,bsewf,excwf0,excwff)	
		     
	 write(2077,*) "BSE exciton band structure finished"
     call flush(2077)		     
	
	elseif (bsetbnd) then
	
	call bsebndstemp(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w,r0,lc,rk,meshtype,bsewf,excwf0,excwff,&
		     st,phavg,ta,temp)	
		     
	 write(2077,*) "BSE exciton band structure, with temperature, finished"
     call flush(2077)
	
	end if		



	if (spdiel) then

	call spoptics(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,cshift,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,rk,meshtype,tmcoef)

	 write(2077,*) "Single particle dielectric calculation finished"
     call flush(2077)

	end if


	if (sppolbz) then

	 call spoptpolbz(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,rk,meshtype)

	 write(2077,*) "Single particle optical activity in BZ finished"
     call flush(2077)

	end if

	!calculo curvatura de berry

	if (berryk) then

	call berrycurv(nthreads,outputfolder,params,kpaths,sme)

	 write(2077,*) "Berry curvature in kpath finished"
     call flush(2077)
	end if


	if (berrybz) then

	call berrycurvbz(nthreads,outputfolder,params,sme,ngrid,mshift)

	 write(2077,*) "Berry curvature in BZ finished"
     call flush(2077)

	end if
	
	
131 continue

	!calculo espectro

	if ((spec) .and. (bse .or. bset)) then

	call bsedielraw(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,cshift) !calculate dielectric constants with BSE

	call bseoptprop(numbse,outputfolder) !calculate abs coeficient and other properties with BSE
	
	if (bset) then
	
		continue
	
	else
	
		call lifetime(sysdim,numbse,ngrid,rlat,nc,nv,outputfolder)
	
	end if
	


	 write(2077,*) "BSE dielectric properties calculated"
	 call flush(2077)

	 if (cpol) then
	
	 call bsedielrawpol(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,cshift) !calculate dielectric constants with BSE for 													different light polarization
	 call bseoptproppol(numbse,outputfolder) !calculate abs coeficient and other properties with BSE
	

	 
	  if (bset) then
	
		continue
	
	 else
	
		call lifetimepol(sysdim,numbse,ngrid,rlat,nc,nv,outputfolder)
	
	 end if

	  write(2077,*) "BSE absorption spectrum with light polarization calculated"
      call flush(2077)
	
	 end if

	end if

	

	if ((spec) .and. ((spdiel) .or. (bse .or. bset) )) then

	call spdielraw(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,cshift) !calculate dielectric constants

	call spoptprop(numbse,outputfolder) !calculate abs coeficient and other properties

	 write(2077,*) "Single particle dielectric properties calculated"
     call flush(2077)

	end if
	
	

	if ((spec) .and. ((cpol) .or. (bse .or. bset) )) then
	
	call spdielrawpol(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,cshift) !calculate dielectric constants different 												       light polarizations
	call spoptproppol(numbse,outputfolder) !calculate abs coeficient and other properties


	 write(2077,*) "Single particle absorption spectrum with light polarization calculated"
     call flush(2077)

	end if


	if (pce) then
	
	call pcecalc(outputfolder,numbse,"IPA",ctemp,ses,tmax,eg,egd)
	
		 write(2077,*) "PCE with single particle absortion spectrum"
         call flush(2077)
		 
	call pcecalc(outputfolder,numbse,"BSE",ctemp,ses,tmax,egs,ebgs)	
		 
		 write(2077,*) "PCE with excitonic effects (BSE)"
         call flush(2077)
	
	end if


	call cpu_time(tf)
	call date_and_time(VALUES=values2)

	write(2077,*)
	write(2077,*) 'End','   ','day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(2077,*)
	!write(2077,*) "Total Time:",(tf-t0)/nthreads,"s"	




	close(2077)

end program main







