
!gfortran spec_rpa-tool.f90 -o spec_rpa-pp.x

subroutine sprawpol(nthreads,outputfolder,params,ngrid,nc,nv,ebse0,ebsef,numbse,sme)

	use omp_lib
	use hamiltonian_input_variables
	implicit none

	integer :: i,j,erro
	integer :: dimbse     

	double precision,allocatable,dimension(:) :: exlinear,fosclinearx 
	double precision,allocatable,dimension(:) :: fosclineary 
	double precision,allocatable,dimension(:) :: fosclinearz 
	double precision,allocatable,dimension(:) :: foscsigmap 
	double precision,allocatable,dimension(:) :: foscsigmam 


	double precision :: flag
	double precision :: vc
	double precision :: elux

	double precision :: intensidadelinearx
	double precision :: intensidadelineary
	double precision :: intensidadelinearz
	double precision :: intensidadesigmap
	double precision :: intensidadesigmam

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


	dimbse=ngrid(1)*ngrid(2)*ngrid(3)*nc*nv

	!arquivo de entrada
	OPEN(UNIT=100, FILE=trim(outputfolder)//"sp_optics-pol.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Erro na abertura do arquivo de entrada sp optics"



	call OMP_SET_NUM_THREADS(nthreads)

	call hamiltonian_input_read(200,params)
	ediel(2) = edielh

	if ( systype .eq. "NP" ) then

	spinf = 2.0

	else

	spinf = 1.0

	end if



	!arquivos de saida


	OPEN(UNIT=300, FILE=trim(outputfolder)//"sp_spectrum.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Erro na abertura do arquivo de saida spec sp "

	!parametros do calculo



	!termino parametros calculo


	call vcell(systype,rlat,vc)

	allocate(exlinear(dimbse))
	allocate(fosclinearx(dimbse),fosclineary(dimbse),foscsigmap(dimbse),foscsigmam(dimbse))
	allocate(fosclinearz(dimbse))



	read(100,*) aread
		
	do j=1,dimbse

	read(100,*) exlinear(j),fosclinearx(j),fosclineary(j),fosclinearz(j),foscsigmap(j),foscsigmam(j)

	end do



	write(300,*) "#","  ","x","  ","y","  ","z","  ","sigma_plus","  ","sigma_minus"

	!$omp do ordered
	do j=1,int(numbse)

		elux= ebse0 + dble(((ebsef-ebse0)*(j-1))/(numbse-1.))

		call optintensitypol(vc,dimbse,ngrid,elux,exlinear,fosclinearx,sme,intensidadelinearx)
		call optintensitypol(vc,dimbse,ngrid,elux,exlinear,fosclineary,sme,intensidadelineary)
		call optintensitypol(vc,dimbse,ngrid,elux,exlinear,fosclinearz,sme,intensidadelinearz)
		call optintensitypol(vc,dimbse,ngrid,elux,exlinear,foscsigmap,sme,intensidadesigmap)
		call optintensitypol(vc,dimbse,ngrid,elux,exlinear,foscsigmam,sme,intensidadesigmam)
		!$omp ordered
	write(300,"(6F15.4)") elux,spinf*intensidadelinearx,spinf*intensidadelineary&
		     ,spinf*intensidadelinearz,spinf*intensidadesigmap,spinf*intensidadesigmam
	
		call flush(300)
		!$omp end ordered


	end do
	!$omp end  do

	deallocate(exlinear)
	deallocate(fosclinearx,fosclineary,fosclinearz,foscsigmap,foscsigmam)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)

	close(100)
	close(200)

	close(300)



	

end subroutine sprawpol


