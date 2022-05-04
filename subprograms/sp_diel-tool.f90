
!ifort rpa_diel-tool.f90 -o spec_rpa_diel.x -qopenmp -mkl

subroutine spdielraw(nthreads,outputfolder,renorm,params,ngrid,nc,nv,ebse0,ebsef,numbse,sme)

	use omp_lib
	use hamiltonian_input_variables
	implicit none

	integer :: i,j,erro
	integer :: dimbse     

	double precision,allocatable,dimension(:) :: esp
	double precision,allocatable,dimension(:) :: exxf,exyf,exzf 
	double precision,allocatable,dimension(:) :: eyyf,eyzf,ezzf 

	double precision :: flag
	double precision :: vc
	double precision :: elux

	double precision :: rxx,rxy,rxz
	double precision :: ryy,ryz,rzz

	double precision :: ixx,ixy,ixz
	double precision :: iyy,iyz,izz

	character(len=2) :: aread

	double precision :: spinf
	
	double precision,allocatable,dimension(:,:) :: dielfxx,dielfyy,dielfzz
	double precision,allocatable,dimension(:,:) :: dielfxy,dielfxz,dielfyz	

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


	!call input_read

	dimbse=ngrid(1)*ngrid(2)*ngrid(3)*nc*nv

	!arquivo de entrada

	OPEN(UNIT=100, FILE=trim(outputfolder)//"sp_optics.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening sp_optics input file"

	call OMP_SET_NUM_THREADS(nthreads)

	call hamiltonian_input_read(200,params)
	!ediel(2) = edielh

	if ( systype .eq. "NP" ) then

	spinf = 2.0

	else

	spinf = 1.0

	end if

	!arquivos de saida

	OPEN(UNIT=300, FILE=trim(outputfolder)//"sp_diel_xx.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_xx output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"sp_diel_xy.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_xy output file"
	OPEN(UNIT=302, FILE=trim(outputfolder)//"sp_diel_xz.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_xz output file"
	OPEN(UNIT=303, FILE=trim(outputfolder)//"sp_diel_yy.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_yy output file"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"sp_diel_yz.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_yz output file"
	OPEN(UNIT=305, FILE=trim(outputfolder)//"sp_diel_zz.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Error opening sp_diel_zz output file"


	!parametros do calculo


	!termino parametros calculo


	!call vcell(systype,rlat,vc)
	call vcell3D(rlat,vc)

	allocate(esp(dimbse))
	allocate(exxf(dimbse),exyf(dimbse),exzf(dimbse),eyyf(dimbse))
	allocate(eyzf(dimbse),ezzf(dimbse))

	read(100,*) aread
		
	do j=1,dimbse

	read(100,*) esp(j),exxf(j),eyyf(j),ezzf(j),exyf(j),exzf(j),eyzf(j)

	end do


	allocate(dielfxx(int(numbse),3),dielfyy(int(numbse),3),dielfzz(int(numbse),3))
	allocate(dielfxy(int(numbse),3),dielfxz(int(numbse),3),dielfyz(int(numbse),3))

	write(300,*) "#","  ","energy","  ","real","  ","imag"
	write(301,*) "#","  ","energy","  ","real","  ","imag"
	write(302,*) "#","  ","energy","  ","real","  ","imag"
	write(303,*) "#","  ","energy","  ","real","  ","imag"
	write(304,*) "#","  ","energy","  ","real","  ","imag"
	write(305,*) "#","  ","energy","  ","real","  ","imag"

	!$omp do 
	!ordered
	do j=1,int(numbse)

		elux= ebse0 + dble(((ebsef-ebse0)*(j-1))/(numbse-1.))

		call optdiel(1,vc,dimbse,ngrid,elux,esp,exxf,sme,rxx,ixx)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,eyyf,sme,ryy,iyy)
		call optdiel(1,vc,dimbse,ngrid,elux,esp,ezzf,sme,rzz,izz)


		call optdiel(0,vc,dimbse,ngrid,elux,esp,exyf,sme,rxy,ixy)
		call optdiel(0,vc,dimbse,ngrid,elux,esp,exzf,sme,rxz,ixz)
		call optdiel(0,vc,dimbse,ngrid,elux,esp,eyzf,sme,ryz,iyz)
		
		dielfxx(j,1) = elux
		dielfxx(j,2) = rxx
		dielfxx(j,3) = ixx
		
		dielfyy(j,1) = elux
		dielfyy(j,2) = ryy
		dielfyy(j,3) = iyy
		
		dielfzz(j,1) = elux
		dielfzz(j,2) = rzz
		dielfzz(j,3) = izz
		
		dielfxy(j,1) = elux
		dielfxy(j,2) = rxy
		dielfxy(j,3) = ixy
		
		dielfxz(j,1) = elux
		dielfxz(j,2) = rxz
		dielfxz(j,3) = ixz
		
		dielfyz(j,1) = elux
		dielfyz(j,2) = ryz
		dielfyz(j,3) = iyz			

		!!$omp ordered
		!write(300,*) real(elux),real(spinf*rxx),real(spinf*ixx)
		!write(301,*) real(elux),real(spinf*rxy),real(spinf*ixy)
		!write(302,*) real(elux),real(spinf*rxz),real(spinf*ixz)
		!write(303,*) real(elux),real(spinf*ryy),real(spinf*iyy)
		!write(304,*) real(elux),real(spinf*ryz),real(spinf*iyz)
		!write(305,*) real(elux),real(spinf*rzz),real(spinf*izz)
	
		!call flush(300)
		!!$omp end ordered


	end do
	!$omp end  do
	
	if (renorm) then
	
	call imagrenorm(dimbse,int(numbse),esp,dielfxx)
	call imagrenorm(dimbse,int(numbse),esp,dielfyy)
	call imagrenorm(dimbse,int(numbse),esp,dielfzz)
	call imagrenorm(dimbse,int(numbse),esp,dielfxy)
	call imagrenorm(dimbse,int(numbse),esp,dielfxz)
	call imagrenorm(dimbse,int(numbse),esp,dielfyz)
	
	else
		continue
	end if		
	
	do i=1,int(numbse)
	
		write(300,*) dielfxx(i,1),real(spinf*dielfxx(i,2)),real(spinf*dielfxx(i,3))
		write(301,*) dielfxy(i,1),real(spinf*dielfxy(i,2)),real(spinf*dielfxy(i,3))
		write(302,*) dielfxz(i,1),real(spinf*dielfxz(i,2)),real(spinf*dielfxz(i,3))
		write(303,*) dielfyy(i,1),real(spinf*dielfyy(i,2)),real(spinf*dielfyy(i,3))
		write(304,*) dielfyz(i,1),real(spinf*dielfyz(i,2)),real(spinf*dielfyz(i,3))
		write(305,*) dielfzz(i,1),real(spinf*dielfzz(i,2)),real(spinf*dielfzz(i,3))
														
		!write(301,*) elux,real(spinf*rxy),real(spinf*ixy)
		!write(302,*) elux,real(spinf*rxz),real(spinf*ixz)
		!write(303,*) elux,real(spinf*ryy),real(spinf*iyy)
		!write(304,*) elux,real(spinf*ryz),real(spinf*iyz)
		!write(305,*) elux,real(spinf*rzz),real(spinf*izz)
	
	end do	

	deallocate(esp)
	deallocate(exxf,exyf,exzf,eyyf)
	deallocate(eyzf,ezzf)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)
	
	deallocate(dielfxx,dielfyy,dielfzz)
	deallocate(dielfxy,dielfxz,dielfyz)	

	close(100)

	close(200)

	close(300)
	close(301)
	close(303)
	close(304)
	close(305)	




	

end subroutine spdielraw


