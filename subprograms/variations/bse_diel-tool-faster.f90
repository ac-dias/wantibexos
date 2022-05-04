INCLUDE "./subroutines/module_input_read.f90"
INCLUDE "./subroutines/general_subs.f90"
INCLUDE "./subroutines/mhkpack_subs.f90"


!ifort bse_diel-tool-faster.f90 -o spec_bse_diel-f.x -qopenmp -mkl

program main
	use omp_lib
	use input_variables
	use hamiltonian_input_variables
	implicit none

	integer :: i,j,erro
	integer :: dimbse     

	double precision,allocatable,dimension(:) :: erpa
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


	call input_read
	call hamiltonian_input_read(200,params)

	if (meshtype .eq. "RK3D") then

		call rkmesh(rk,rlat,ngrid)

	else if (meshtype .eq. "RK2D") then

		call rkmesh2D(rk,rlat,ngrid)
	else

		continue
	end if

	dimbse=ngrid(1)*ngrid(2)*ngrid(3)*nc*nv

	!arquivo de entrada

	OPEN(UNIT=100, FILE=trim(outputfolder)//"bse_opt_diel-faster.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de entrada bse xx yy zz"


	call OMP_SET_NUM_THREADS(nthreads)


	!arquivos de saida

	OPEN(UNIT=300, FILE=trim(outputfolder)//"bse_diel_xx-f.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Erro na abertura do arquivo de saida rpa diel xx"
	OPEN(UNIT=303, FILE=trim(outputfolder)//"bse_diel_yy-f.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Erro na abertura do arquivo de saida rpa diel yy"
	OPEN(UNIT=305, FILE=trim(outputfolder)//"bse_diel_zz-f.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Erro na abertura do arquivo de saida rpa diel zz"


	!parametros do calculo

	OPEN(UNIT=2077, FILE= trim(calcparms)//"spec_bse_calc-diel.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida parametros calculo"

	call param_out(2077,nthreads,outputfolder,calcparms,ngrid,nc,nv,edos0,edosf,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype)



	close(2077)

	!termino parametros calculo


	call vcell(systype,rlat,vc)

	allocate(erpa(dimbse))
	allocate(exxf(dimbse),exyf(dimbse),exzf(dimbse),eyyf(dimbse))
	allocate(eyzf(dimbse),ezzf(dimbse))

	read(100,*) aread
	!read(101,*) aread
		
	do j=1,dimbse

	read(100,*) erpa(j),exxf(j),eyyf(j),ezzf(j),exyf(j),exzf(j),eyzf(j)
	!read(101,*) erpa(j),exyf(j),exzf(j),eyzf(j)

	end do



	write(300,*) "#","  ","energy","  ","real","  ","imag"
	write(303,*) "#","  ","energy","  ","real","  ","imag"
	write(305,*) "#","  ","energy","  ","real","  ","imag"

	!$omp do ordered
	do j=1,int(numbse)

		elux= ebse0 + dble(((ebsef-ebse0)*(j-1))/(numbse-1.))

		call rpadiel2(1,vc,dimbse,ngrid,elux,erpa,exxf,sme,rxx,ixx)
		call rpadiel2(1,vc,dimbse,ngrid,elux,erpa,eyyf,sme,ryy,iyy)
		call rpadiel2(1,vc,dimbse,ngrid,elux,erpa,ezzf,sme,rzz,izz)


		!call rpadiel2(0,vc,dimbse,ngrid,elux,erpa,exyf,sme,rxy,ixy)
		!call rpadiel2(0,vc,dimbse,ngrid,elux,erpa,exzf,sme,rxz,ixz)
		!call rpadiel2(0,vc,dimbse,ngrid,elux,erpa,eyzf,sme,ryz,iyz)

		!$omp ordered
		write(300,*) elux,real(rxx),real(ixx)
		!write(301,*) elux,rxy,ixy
		!write(302,*) elux,rxz,ixz
		write(303,*) elux,real(ryy),real(iyy)
		!write(304,*) elux,ryz,iyz
		write(305,*) elux,real(rzz),real(izz)
	
		!$omp end ordered


	end do
	!$omp end  do

	deallocate(erpa)
	deallocate(exxf,exyf,exzf,eyyf)
	deallocate(eyzf,ezzf)

	close(100)
	close(101)

	close(200)

	close(300)
	!close(301)
	!close(302)
	close(303)
	!close(304)
	close(305)	




	

end program main

subroutine rpadiel2(rtype,vc,dimse,ngrid,elux,exciton,fosc,sme,rpart,ipart)


	implicit none

	double precision :: elux
	integer :: dimse
	double precision :: sme
	double precision :: rpart,ipart
	integer :: i,j

	integer :: rtype
	double precision :: delta

	integer,dimension(3) :: ngrid

	double precision :: vc

	double precision :: raux,iaux

	double precision,dimension(dimse) :: exciton
	double precision,dimension(dimse) :: fosc
	
	double precision :: caux,eaux

	double precision,parameter :: pi=acos(-1.)

	double precision,parameter :: gama= 180.90!180.90

	double complex,parameter :: imag= cmplx(0.0,1.0)

	double precision :: baux

	baux = real(ngrid(1)*ngrid(2)*ngrid(3))*vc

	caux=gama/(baux*4.0*pi)


	if ( rtype .eq. 1) then

		delta = 1.0
	else 

		delta = 0.0
	end if 


	rpart= 0.0
	ipart= 0.0


	do i=1,dimse

		eaux = ((elux-exciton(i))**2)+(sme**2)



		raux = fosc(i)*((exciton(i)-elux)/(eaux))

		iaux = fosc(i)*((sme)/(eaux))

		!raux = raux/eaux
		
		!iaux = iaux/eaux
	

		rpart=rpart+raux

		ipart=ipart+iaux


	end do

	!write(*,*) caux

	ipart= ipart*caux

	rpart= delta+(rpart*caux)


end subroutine rpadiel2
