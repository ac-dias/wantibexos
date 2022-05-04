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
INCLUDE "./subroutines/sistema_tb.f90"

!############subprograms##################

INCLUDE "./subprograms/bands-kpath-tool.f90"

INCLUDE "./subprograms/berry_curvature_bz-tool.f90"
INCLUDE "./subprograms/berry_curvature_kpath-tool.f90"

INCLUDE "./subprograms/bse_diel-tool.f90"
INCLUDE "./subprograms/bse_diel-tool-pol.f90"
INCLUDE "./subprograms/bse_kpath-tool.f90"
INCLUDE "./subprograms/bse_solver-tool-diel.f90"
INCLUDE "./subprograms/bse_solver-tool-pol.f90"

INCLUDE "./subprograms/diel-pp.f90"
INCLUDE "./subprograms/diel-pp-pol.f90"
INCLUDE "./subprograms/diel-pp-bse.f90"
INCLUDE "./subprograms/diel-pp-bse-pol.f90"
INCLUDE "./subprograms/exciton_lifetime.f90"

INCLUDE "./subprograms/sp_diel-tool.f90"
INCLUDE "./subprograms/sp_diel-tool-pol.f90"
INCLUDE "./subprograms/sp_opt_bz-tool.f90"
INCLUDE "./subprograms/sp_solver-tool-diel.f90"
INCLUDE "./subprograms/sp_solver-tool-pol.f90"


INCLUDE "./subprograms/tdos-tool.f90"

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
	
	!add variables for teste
	integer :: exc0,excf
	

	call input_read

	OPEN(UNIT=2055, FILE= params,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de entrada parametros"	

	read(2055,*) systype
	read(2055,*) ediel(2)
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

	OPEN(UNIT=2077, FILE= trim(calcparms)//"log-test.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida parametros calculo"

	call param_out(2077,nthreads,outputfolder,calcparms,ngrid,nc,nv,edos0,edosf,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,bandscalc,doscalc,bse,bsepol,bsekpath,spec,&
		     spdiel,spdielpol,sppolbz,berryk,berrybz,pponly,bsewf,excwf0,excwff,&
		     tmcoef,ez,w,lc,r0,sysdim)

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


	exc0 = 1
	excf = 3


	 call bsewf(nthreads,outputfolder,calcparms,ngrid,nc,nv,edos0,edosf,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,rk,meshtype,exc0,excf)



	call cpu_time(tf)
	call date_and_time(VALUES=values2)

	write(2077,*)
	write(2077,*) 'End','   ','day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(2077,*)




	close(2077)

end program main







