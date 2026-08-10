!############subroutines################

!INCLUDE "./subroutines/berry_curvature_subs.F90"
!INCLUDE "./subroutines/bse_subs.F90"
!INCLUDE "./subroutines/bse_subs_kpath.F90"
!INCLUDE "./subroutines/coulomb_pot.F90"
!INCLUDE "./subroutines/diel-pp-subs.F90"
!INCLUDE "./subroutines/dos_subs.F90"
!INCLUDE "./subroutines/general_subs.F90"
!INCLUDE "./subroutines/mhkpack_subs.F90"
!INCLUDE "./subroutines/module_input_read.F90"
!INCLUDE "./subroutines/optics.F90"

!!INCLUDE "./subroutines/lowdin.F90"

!INCLUDE "./subroutines/hamiltonian_tb.F90"
!INCLUDE "./subroutines/bse_subs_temp.F90"
!INCLUDE "./subroutines/pce-subs.F90"
!INCLUDE "./subroutines/efmass-subs.F90"
!INCLUDE "./subroutines/spin_txt_subs.F90"
!INCLUDE "./subroutines/hamiltonians.F90"
!INCLUDE "./subroutines/boltzmann_subs.F90"

!############subprograms##################

!INCLUDE "./subprograms/bands-kpath-tool.F90"

!INCLUDE "./subprograms/berry_curvature_bz-tool.F90"
!INCLUDE "./subprograms/berry_curvature_kpath-tool.F90"

!INCLUDE "./subprograms/bse_diel-tool.F90"
!INCLUDE "./subprograms/bse_diel-tool-pol.F90"
!INCLUDE "./subprograms/bse_kpath-tool.F90"
!INCLUDE "./subprograms/bse_kpath-tool-temp.F90"


!INCLUDE "./subprograms/bse_solver-tool-diel.F90"
!INCLUDE "./subprograms/bse_solver-tool-diel-temp.F90"


!!INCLUDE "./subprograms/bse_solver-tool-pol.F90"

!INCLUDE "./subprograms/diel-pp.F90"
!INCLUDE "./subprograms/diel-pp-pol.F90"
!INCLUDE "./subprograms/diel-pp-bse.F90"
!INCLUDE "./subprograms/diel-pp-bse-pol.F90"
!INCLUDE "./subprograms/exciton_lifetime.F90"

!INCLUDE "./subprograms/sp_diel-tool.F90"
!INCLUDE "./subprograms/sp_diel-tool-pol.F90"
!INCLUDE "./subprograms/sp_opt_bz-tool.F90"
!INCLUDE "./subprograms/sp_solver-tool-diel.F90"

!!INCLUDE "./subprograms/sp_solver-tool-pol.F90"


!INCLUDE "./subprograms/tdos-tool-stxt.F90"
!INCLUDE "./subprograms/pce-code.F90"
!INCLUDE "./subprograms/efmass.F90"
!INCLUDE "./subprograms/boltzmann_transport.F90"

!!INCLUDE "./subprograms/spin_txt.F90"


!ifort tbt_main.F90 -o tbt.x -mkl -qopenmp -shared-intel
program main

#ifdef MPI
!	include 'mpif.h'
        use mpi
#endif
	use input_variables
	use hamiltonian_input_variables
	implicit none
	
	integer :: erro
	real:: t0,tf
	integer,dimension(8) :: values,values2

	real,dimension(3,3) :: rlatv
	real :: flag
	character(len=70) ::cflag

      integer:: MPIError, Node, Nodes

#ifdef MPI
      call MPI_Init(MPIError)
      call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
      call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )
#else
      Node =  0
      Nodes = 1
#endif

    if (Node == 0) then
	  call input_read
    endif

! Broadcast all variables between nodes
#ifdef MPI
    call bcast_input_read
#endif

    if (Node == 0 ) then
	OPEN(UNIT=2055, FILE= params,STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening hamiltonian input file (main)"	

	read(2055,*) systype
	!read(2055,*) ediel(2)
	read(2055,*) scs
	read(2055,*) efermi
	!read(2055,*) lc


	read(2055,*) rlatv(1,1),rlatv(1,2),rlatv(1,3)
	read(2055,*) rlatv(2,1),rlatv(2,2),rlatv(2,3)
	read(2055,*) rlatv(3,1),rlatv(3,2),rlatv(3,3)
    close(2055)
    endif

#ifdef MPI
    call MPI_BCAST(rlatv, 9, MPI_REAL, 0, MPI_COMM_WORLD, MPIError)
#endif

	if (meshtype .eq. "RK3D") then

		call rkmesh(rk,rlatv,ngrid)

	else if (meshtype .eq. "RK2D") then

		call rkmesh2D(rk,rlatv,ngrid)
		
	else if (meshtype .eq. "RK1D") then

		call rkmesh1D(rk,rlatv,ngrid)	
	else

		continue
	end if


    if (Node == 0) then
     OPEN(UNIT=2077, FILE= trim(calcparms)//"log.dat",STATUS='unknown', IOSTAT=erro)
     print*, trim(calcparms)//"log.dat"
    	if (erro/=0) stop "Error opening log output file "

! This subroutine only prints things
    call param_out(2077,nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     mshift,coultype,bandscalc,doscalc,bse,bsepol,bsekpath,spec,&
		     spdiel,spdielpol,sppolbz,berryk,berrybz,pponly,bsewf,excwf0,excwff,&
		     tmcoef,ez,w,lc,r0,sysdim,dtfull,cpol,cshift,dft,&
		     st,phavg,temp,ta,pce,ses,ctemp,tmax,eg,egd,egs,ebgs,renorm,emt,emfile,&
		     dk,nocpf,fermishift,bsealgo,spintxt,meshgen,exc,mag,boltz,nmu,nsteps,mu0,&
		     bsehamwrite,bsehamread,bsehamfile,&
		     muf,btemp,klat,elft,hlft,smeboltz,ni,ns,edielgw,r0gw,lcgw,ezgw,wgw,coultypegw,&
		     nomega,omegamax,gwmesh,gwmeshuse,gwbnd,gwbnduse,gwbsebnd,ktolgw,smegw,selfxonly,ifactor)
		     
   

	call cpu_time(t0)
	call date_and_time(VALUES=values)

	write(2077,*)
	write(2077,*) 'Begin','  ','month',values(2),'day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(2077,*)

    endif

	!calculo estrutura eletronica
	if (pponly) go to 131
	
	if (gwmesh) then
#ifdef MPI
         if (Nodes /= 1) then
          write(2077,*) "Not implemented in parallel yet"
          call mpi_abort(MPI_Comm_World,1, MPIError)
         endif
#endif
	  !call gw_pi0_calc(nthreads,outputfolder,ngrid,smegw,params,&
          !             exc,mshift,nocpf,fermishift,dft,mag,&
          !             nomega,omegamax,sysdim)
	
	 !call gwmeshcalc(nthreads,outputfolder,calcparms,ngrid,smegw,ktolgw,params,edielgw,&
          !             exc,mshift,coultypegw,ezgw,wgw,r0gw,lcgw,rk,meshtype,nocpf,fermishift,dft,mag,&
           !            nomega,omegamax,sysdim,selfxonly)
                       
	 !call sgw(nthreads,outputfolder,ngrid,smegw,ktolgw,params,edielgw,&
          !             exc,mshift,coultypegw,ezgw,wgw,r0gw,lcgw,rk,meshtype,nocpf,fermishift,dft,mag,&
           !            nomega,omegamax,sysdim,selfxonly) 
                       
	 call sgw_mesh(nthreads,outputfolder,ngrid,ifactor,smegw,ktolgw,params,edielgw,&
                       exc,mshift,coultypegw,ezgw,wgw,r0gw,lcgw,rk,meshtype,nocpf,fermishift,dft,mag,&
                       nomega,omegamax,sysdim,selfxonly)                                               
	
	 write(2077,*) "G0W0 k-mesh finished"
         call flush(2077)
	
	end if	
	
	!if (gwbnd) then
	
	  !call gw_pi0_calc(nthreads,outputfolder,ngrid,smegw,params,&
          !             exc,mshift,nocpf,fermishift,dft,mag,&
          !             nomega,omegamax,sysdim)
	
	 !call gwmeshcalc(nthreads,outputfolder,calcparms,ngrid,smegw,ktolgw,params,edielgw,&
          !             exc,mshift,coultypegw,ezgw,wgw,r0gw,lcgw,rk,meshtype,nocpf,fermishift,dft,mag,&
           !            nomega,omegamax,sysdim,selfxonly)
                       
	 !call sgw_bnd(nthreads,outputfolder,ngrid,smegw,ktolgw,params,kpaths,edielgw,&
          !             exc,mshift,coultypegw,ezgw,wgw,r0gw,lcgw,rk,meshtype,nocpf,fermishift,dft,mag,&
           !            nomega,omegamax,sysdim,selfxonly)                       
	
	 !write(2077,*) "G0W0 k-path finished"
         !call flush(2077)
	
	!end if		
	


	if (bandscalc) then
#ifdef MPI
        if (Nodes /= 1) then
         write(2077,*) "Not implemented in parallel yet"
         call mpi_abort(MPI_Comm_World,1, MPIError)
        endif
#endif
	 call bandstool(nthreads,outputfolder,params,kpaths,orbw, &
		     exc,mag,mshift,dft,nocpf,fermishift)
	 write(2077,*) "Band Structure finished"
         call flush(2077)
	end if


	if (doscalc) then
#ifdef MPI
        if (Nodes /= 1) then
         write(2077,*) "Not implemented in parallel yet"
         call mpi_abort(MPI_Comm_World,1, MPIError)
        endif
#endif

	 call dostool(nthreads,outputfolder,ngrid,numdos, &
		      sme,params,orbw,exc,mag,mshift,dft,nocpf,fermishift,spintxt)

	 write(2077,*) "DOS finished"
         call flush(2077)
	end if
	
	if (spintxt) then
#ifdef MPI
        if (Nodes /= 1) then
         write(2077,*) "Not implemented in parallel yet"
         call mpi_abort(MPI_Comm_World,1, MPIError)
        endif
#endif
	write(2077,*) "Spin Texture finished"
	call flush(2077)
	end if
	
		
	if (emt) then
#ifdef MPI
        if (Nodes /= 1) then
         write(2077,*) "Not implemented in parallel yet"
         call mpi_abort(MPI_Comm_World,1, MPIError)
        endif
#endif

	 !call efmass(nthreads,dft,outputfolder,params,emfile,dk,nocpf,fermishift,exc,mag,sysdim)
	 
	 call efmass_num(nthreads,dft,outputfolder,params,emfile,dk,nocpf,fermishift,exc,mag,sysdim)

	 write(2077,*) "Effective Mass Tensor finished"
         call flush(2077)
	end if	
	
	if (boltz) then
#ifdef MPI
        if (Nodes /= 1) then
         write(2077,*) "Not implemented in parallel yet"
         call mpi_abort(MPI_Comm_World,1, MPIError)
        endif
#endif
		call boltztransport(nthreads,outputfolder,ngrid,nsteps,smeboltz,params,exc,mag,mshift,dft,nocpf,fermishift, &
                           mu0,muf,nmu,btemp,klat,elft,hlft)
                           
               write(2077,*) "Boltzmann transport calculations finished"
               call flush(2077)            
	
	end if

	!calculo optica
	
	if (bse) then
         if (Node==0) then
	  write(2077,*) "BSE and Single particle dielectric calculation started"
          call flush(2077)
         endif
		if (temp .eq. 0.0) then
		
		call bsesolver(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,cshift,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w,r0,lc,rk,meshtype,bsewf,excwf0,excwff,dtfull,cpol,tmcoef,&
		     nocpf,fermishift,bsealgo,bsehamwrite,bsehamread,bsehamfile,dft,mag)

                if (Node == 0) then
	 	  write(2077,*) "BSE and Single particle dielectric calculation finished"
     		  call flush(2077)
		endif

		else
#ifdef MPI
            if (Nodes /= 1) then
             write(2077,*) "Not implemented in parallel yet"
             call mpi_abort(MPI_Comm_World,1, MPIError)
            endif
#endif

		call bsesolvertemp(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,cshift,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w,r0,lc,rk,meshtype,bsewf,excwf0,excwff,dtfull,&
		     cpol,tmcoef,st,phavg,ta,temp,nocpf,fermishift,bsealgo,bsehamwrite,bsehamread,bsehamfile,dft,mag)
		
                if (Node==0) then      
	 	 write(2077,*) "BSE and Single particle dielectric calculation, with temperature, finished"	
     		 call flush(2077)	     	
		endif

		end if
	
	end if
	

	if (bsekpath) then
	
		if (temp .eq. 0.0) then
		
		call bsebnds(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w,r0,lc,rk,meshtype,bsewf,excwf0,excwff,&
		     nocpf,fermishift,bsealgo,dft,mag)	
		     
		 write(2077,*) "BSE exciton band structure finished"
     		 call flush(2077)	
		
		else 
		
		call bsebndstemp(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w,r0,lc,rk,meshtype,bsewf,excwf0,excwff,&
		     st,phavg,ta,temp,nocpf,fermishift,bsealgo,dft,mag)	
		     
	 	write(2077,*) "BSE exciton band structure, with temperature, finished"
     		call flush(2077)
		
		end if
	
	end if
	
	
	if (spdiel .and. Node == 0) then

	call spoptics(nthreads,dft,outputfolder,ngrid,nc,nv, &
		     cshift,params,exc,mshift,tmcoef,nocpf,fermishift,mag)

	 write(2077,*) "Single particle dielectric calculation finished"
         call flush(2077)

	end if


	if (sppolbz) then

	 call spoptpolbz(nthreads,dft,outputfolder,ngrid,nc,nv, &
		     cshift,params,exc,mshift,nocpf,fermishift,meshgen,mag)

	 write(2077,*) "Single particle optical activity in BZ finished"
         call flush(2077)

	end if

	!calculo curvatura de berry

	if (berryk) then

	 call berrycurv(nthreads,dft,outputfolder,params,kpaths,sme,nocpf,fermishift,exc,mag)

	 write(2077,*) "Berry curvature in kpath finished"
         call flush(2077)
	end if


	if (berrybz) then

	 call berrycurvbz(nthreads,dft,outputfolder,params,sme,ngrid,mshift,nocpf,fermishift,meshgen,exc,mag)

	 write(2077,*) "Berry curvature in BZ finished"
         call flush(2077)

	end if
	
		
	
131 continue

	! The spectrum post-processing routines below use ordinary formatted files.
	! Wait until all MPI ranks have finished producing their distributed results,
	! then let only rank zero create and consume the shared post-processing files.
#ifdef MPI
	call MPI_BARRIER(MPI_COMM_WORLD,MPIError)
#endif

	if (Node == 0) then

	!calculo espectro

	if ((spec) .and. (bse)) then

	call bsedielraw(nthreads,dft,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,cshift) !calculate dielectric constants with BSE

	call bseoptprop(numbse,outputfolder,tmax,ni,ns) !calculate abs coeficient and other properties with BSE
		
	call lifetime(sysdim,numbse,ngrid,rlat,nc,nv,outputfolder)
	
	call emissionopt(nthreads,"BSE",outputfolder,numbse,ngrid,nc,nv,ctemp,cshift,egs)
	
        if (Node == 0) then
	 write(2077,*) "BSE dielectric properties calculated"
	 call flush(2077)
	endif
 
	 if (cpol) then
	
	 call bsedielrawpol(nthreads,dft,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,cshift) !calculate dielectric constants with BSE for 													different light polarization
	 call bseoptproppol(numbse,outputfolder,tmax,ni,ns) !calculate abs coeficient and other properties with BSE
	
	 call lifetimepol(sysdim,numbse,ngrid,rlat,nc,nv,outputfolder)
	 
	 call emissionopt(nthreads,"BSP",outputfolder,numbse,ngrid,nc,nv,ctemp,cshift,egs)
	
         if (Node == 0) then
	  write(2077,*) "BSE absorption/PL spectrum with light polarization calculated"
          call flush(2077)
	 endif

	end if
	
	end if

	if ((spec) .and. ((spdiel) .or. (bse) )) then

	call spdielraw(nthreads,dft,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,cshift) !calculate dielectric constants

	call spoptprop(numbse,outputfolder,tmax,ni,ns) !calculate abs coeficient and other properties
	
	call emissionopt(nthreads,"IPA",outputfolder,numbse,ngrid,nc,nv,ctemp,cshift,eg)
       
        if (Node == 0) then
	  write(2077,*) "Single particle dielectric properties calculated"
          call flush(2077)
        endif
     		if (cpol) then
     		
     			call spdielrawpol(nthreads,dft,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,cshift) !calculate dielectric constants different 												       light polarizations
			call spoptproppol(numbse,outputfolder,tmax,ni,ns) !calculate abs coeficient and other properties

			call emissionopt(nthreads,"IPP",outputfolder,numbse,ngrid,nc,nv,ctemp,cshift,eg)
                     if (Node == 0) then
	 		write(2077,*) "Single particle absorption/PL spectrum with light polarization calculated"
     			call flush(2077)
                     endif
     		
     		end if

	end if

	end if

#ifdef MPI
	! Do not let the remaining ranks finalize while rank zero still owns the files.
	call MPI_BARRIER(MPI_COMM_WORLD,MPIError)
#endif






	!if (pce) then
	
	 !call pcecalc(outputfolder,numbse,"IPA",ctemp,ses,tmax,eg,egd)
	
	!	 write(2077,*) "PCE with single particle absortion spectrum"
        ! call flush(2077)
		 
	!call pcecalc(outputfolder,numbse,"BSE",ctemp,ses,tmax,egs,ebgs)	
		 
	!	 write(2077,*) "PCE with excitonic effects (BSE)"
        ! call flush(2077)
	
	!end if


	call cpu_time(tf)
	call date_and_time(VALUES=values2)
    if (Node == 0) then
	write(2077,*)
	write(2077,*) 'End','   ','month',values2(2),'day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(2077,*)
	!write(2077,*) "Total Time:",(tf-t0)/nthreads,"s"	




	close(2077)
    endif

#ifdef MPI
          call mpi_barrier(MPI_COMM_WORLD,MPIError)
          call mpi_finalize(MPIError)
#endif

end program main
