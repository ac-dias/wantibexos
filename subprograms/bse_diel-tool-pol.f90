

!ifort bse_diel-tool.f90 -o spec_bse_diel.x -qopenmp -mkl

subroutine bsedielrawpol(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,sme)

	use omp_lib
	use hamiltonian_input_variables
	implicit none

	integer :: i,j,erro
	integer :: dimbse     

	double precision,allocatable,dimension(:) :: esp
	double precision,allocatable,dimension(:) :: exxf,eyyf,ezzf 
	double precision,allocatable,dimension(:) :: espf,esmf 

	double precision :: flag
	double precision :: vc
	double precision :: elux

	double precision :: rxx,ryy,rzz
	double precision :: rsp,rsm

	double precision :: ixx,iyy,izz
	double precision :: isp,ism
	
	double precision,allocatable,dimension(:,:) :: dielfxx,dielfyy,dielfzz
	double precision,allocatable,dimension(:,:) :: dielfsp,dielfsm	

	character(len=2) :: aread

	double precision :: spinf


	!modificacoes versao 2.1

	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv
	!integer :: ncrpa,nvrpa
	!integer :: ncbz,nvbz
	double precision :: edos0,edosf,numdos
	double precision :: ebse0,ebsef,numbse
	double precision :: sme,exc,rk
	double precision,dimension(3) :: mshift
	double precision :: ktol
	character(len=70) :: params   !parametros TB
	character(len=70) :: orbw     !peso orbitais lcount
	character(len=70) :: kpaths    !kpath
	character(len=70) :: kpathsbse    !kpath
	!character(len=70) :: diein    !ambiente dieletrico
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: calcparms
	character(len=70) :: meshtype
	character(len=5) :: coultype
	double precision,dimension(3) :: ediel
	logical :: renorm

	call hamiltonian_input_read(200,params)

	!ediel(2) = edielh

	if ( systype .eq. "NP" ) then

	spinf = 2.0

	else

	spinf = 1.0

	end if



	dimbse=ngrid(1)*ngrid(2)*ngrid(3)*nc*nv

	!arquivo de entrada

	OPEN(UNIT=100, FILE=trim(outputfolder)//"bse_opt_diel-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_opt_diel-pol input file"


	call OMP_SET_NUM_THREADS(nthreads)


	!arquivos de saida

	OPEN(UNIT=300, FILE=trim(outputfolder)//"bse_diel-pol_x.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_x output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"bse_diel-pol_y.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_y output file"
	OPEN(UNIT=302, FILE=trim(outputfolder)//"bse_diel-pol_z.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_z output file"
	OPEN(UNIT=303, FILE=trim(outputfolder)//"bse_diel-pol_sp.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_sp output file"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"bse_diel-pol_sm.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening bse_diel-pol_sm output file"



	!parametros do calculo



	!termino parametros calculo


	!call vcell(systype,rlat,vc)
	call vcell3D(rlat,vc)

	allocate(esp(dimbse))
	allocate(exxf(dimbse),eyyf(dimbse),ezzf(dimbse))
	allocate(espf(dimbse),esmf(dimbse))

	read(100,*) aread
		
	do j=1,dimbse

	read(100,*) esp(j),exxf(j),eyyf(j),ezzf(j),espf(j),esmf(j)

	end do

	allocate(dielfxx(int(numbse),3),dielfyy(int(numbse),3),dielfzz(int(numbse),3))
	allocate(dielfsp(int(numbse),3),dielfsm(int(numbse),3))


	write(300,*) "#","  ","energy","  ","real","  ","imag"
	write(301,*) "#","  ","energy","  ","real","  ","imag"
	write(302,*) "#","  ","energy","  ","real","  ","imag"
	write(303,*) "#","  ","energy","  ","real","  ","imag"
	write(304,*) "#","  ","energy","  ","real","  ","imag"


	!$omp do 
	!ordered
	do j=1,int(numbse)

		elux= ebse0 + dble(((ebsef-ebse0)*(j-1))/(numbse-1.))

		call optdiel(1,vc,dimbse,ngrid,elux,esp,exxf,sme,rxx,ixx)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,eyyf,sme,ryy,iyy)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,ezzf,sme,rzz,izz)
		
		call optdiel(1,vc,dimbse,ngrid,elux,esp,espf,sme,rsp,isp)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,esmf,sme,rsm,ism)
		
		dielfxx(j,1) = elux
		dielfxx(j,2) = rxx
		dielfxx(j,3) = ixx
		
		dielfyy(j,1) = elux
		dielfyy(j,2) = ryy
		dielfyy(j,3) = iyy
		
		dielfzz(j,1) = elux
		dielfzz(j,2) = rzz
		dielfzz(j,3) = izz
		
		dielfsp(j,1) = elux
		dielfsp(j,2) = rsp
		dielfsp(j,3) = isp
		
		dielfsm(j,1) = elux
		dielfsm(j,2) = rsm
		dielfsm(j,3) = ism
		
			


		!!$omp ordered
		!write(300,*) real(elux),real(spinf*rxx),real(spinf*ixx)
		!write(301,*) real(elux),real(spinf*ryy),real(spinf*iyy)
		!write(302,*) real(elux),real(spinf*rzz),real(spinf*izz)
		
		!write(303,*) real(elux),real(spinf*rsp),real(spinf*isp)
		!write(304,*) real(elux),real(spinf*rsm),real(spinf*ism)
	
		!call flush(300)
		!!$omp end ordered


	end do
	!$omp end  do
	
	if (renorm) then
	
	call imagrenorm(dimbse,int(numbse),esp,dielfxx)
	call imagrenorm(dimbse,int(numbse),esp,dielfyy)
	call imagrenorm(dimbse,int(numbse),esp,dielfzz)
	call imagrenorm(dimbse,int(numbse),esp,dielfsp)
	call imagrenorm(dimbse,int(numbse),esp,dielfsm)
	
	else
		continue
	end if		
	
	do i=1,int(numbse)
	
		write(300,*) dielfxx(i,1),real(spinf*dielfxx(i,2)),real(spinf*dielfxx(i,3))
		write(301,*) dielfyy(i,1),real(spinf*dielfyy(i,2)),real(spinf*dielfyy(i,3))
		write(302,*) dielfzz(i,1),real(spinf*dielfzz(i,2)),real(spinf*dielfzz(i,3))
		
		write(303,*) dielfsp(i,1),real(spinf*dielfsp(i,2)),real(spinf*dielfsp(i,3))
		write(304,*) dielfsm(i,1),real(spinf*dielfsm(i,2)),real(spinf*dielfsm(i,3))


														
		!write(301,*) elux,real(spinf*rxy),real(spinf*ixy)
		!write(302,*) elux,real(spinf*rxz),real(spinf*ixz)
		!write(303,*) elux,real(spinf*ryy),real(spinf*iyy)
		!write(304,*) elux,real(spinf*ryz),real(spinf*iyz)
		!write(305,*) elux,real(spinf*rzz),real(spinf*izz)
	
	end do	

	deallocate(esp)
	deallocate(exxf,eyyf,ezzf)
	deallocate(espf,esmf)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)
	
	deallocate(dielfxx,dielfyy,dielfzz)
	deallocate(dielfsp,dielfsm)	

	close(100)

	close(200)

	close(300)
	close(301)
	close(302)
	close(303)
	close(304)
	




	

end subroutine bsedielrawpol


