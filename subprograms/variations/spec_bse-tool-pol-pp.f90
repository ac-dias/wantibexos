!gfortran spec_bse-tool.f90 -o spec_bse-pp.x

subroutine bserawpol(nthreads,outputfolder,params,ngrid,nc,nv,ebse0,ebsef,numbse,sme)

	use omp_lib
	use hamiltonian_input_variables
	implicit none

	integer :: i,j,erro
	integer :: dimbse     

	double precision,allocatable,dimension(:) :: exlinearx,fosclinearx 
	double precision,allocatable,dimension(:) :: exlineary,fosclineary 
	double precision,allocatable,dimension(:) :: exlinearz,fosclinearz 
	double precision,allocatable,dimension(:) :: exsigmap,foscsigmap 
	double precision,allocatable,dimension(:) :: exsigmam,foscsigmam 


	double precision :: flag
	double precision :: vc
	double precision :: elux

	double precision :: intensidadelinearx
	double precision :: intensidadelineary
	double precision :: intensidadelinearz
	double precision :: intensidadesigmap
	double precision :: intensidadesigmam

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
	double precision :: nflag


	!call input_read

	dimbse=ngrid(1)*ngrid(2)*ngrid(3)*nc*nv

	!arquivo de entrada
	OPEN(UNIT=100, FILE=trim(outputfolder)//"bse_opt_pol_x.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Erro na abertura do arquivo de entrada bse linear x"

	OPEN(UNIT=101, FILE=trim(outputfolder)//"bse_opt_pol_y.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Erro na abertura do arquivo de entrada bse linear y"

	OPEN(UNIT=102, FILE=trim(outputfolder)//"bse_opt_pol_z.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Erro na abertura do arquivo de entrada bse linear z"

	OPEN(UNIT=103, FILE=trim(outputfolder)//"bse_opt_pol_sp.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Erro na abertura do arquivo de entrada bse sigma +"

	OPEN(UNIT=104, FILE=trim(outputfolder)//"bse_opt_pol_sm.dat",STATUS='old', IOSTAT=erro)
	if (erro/=0) stop "Erro na abertura do arquivo de entrada bse sigma -"

	call OMP_SET_NUM_THREADS(nthreads)

	call hamiltonian_input_read(200,params)
	ediel(2) = edielh

	if ( systype .eq. "NP" ) then

	spinf = 2.0

	else

	spinf = 1.0

	end if

	!arquivos de saida


	OPEN(UNIT=300, FILE= trim(outputfolder)//"bse_spectrum.dat",STATUS='unknown', IOSTAT=erro)
	if (erro/=0) stop "Erro na abertura do arquivo de saida spec bse "


	!parametros do calculo


	!termino parametros calculo


	call vcell3D(rlat,vc)

	allocate(exlinearx(dimbse),exlineary(dimbse),exsigmap(dimbse),exsigmam(dimbse))
	allocate(fosclinearx(dimbse),fosclineary(dimbse),foscsigmap(dimbse),foscsigmam(dimbse))
	allocate(exlinearz(dimbse),fosclinearz(dimbse))




		
	do j=1,dimbse

		read(100,*) exlinearx(j),fosclinearx(j),nflag

		read(101,*) exlineary(j),fosclineary(j),nflag

		read(102,*) exlinearz(j),fosclinearz(j),nflag

		read(103,*) exsigmap(j),foscsigmap(j),nflag

		read(104,*) exsigmam(j),foscsigmam(j),nflag

	end do



	write(300,*) "#","  ","x","  ","y","  ","z","  ","sigma_plus","  ","sigma_minus"

	!$omp do ordered
	do j=1,int(numbse)

		elux= ebse0 + dble(((ebsef-ebse0)*(j-1))/(numbse-1.))

		call optintensitypol(vc,dimbse,ngrid,elux,exlinearx,fosclinearx,sme,intensidadelinearx)
		call optintensitypol(vc,dimbse,ngrid,elux,exlineary,fosclineary,sme,intensidadelineary)
		call optintensitypol(vc,dimbse,ngrid,elux,exlinearz,fosclinearz,sme,intensidadelinearz)
		call optintensitypol(vc,dimbse,ngrid,elux,exsigmap,foscsigmap,sme,intensidadesigmap)
		call optintensitypol(vc,dimbse,ngrid,elux,exsigmam,foscsigmam,sme,intensidadesigmam)

		!$omp ordered
	write(300,"(6F15.4)") elux,spinf*intensidadelinearx,spinf*intensidadelineary,spinf*intensidadelinearz&
                    ,spinf*intensidadesigmap,spinf*intensidadesigmam

	!write(300,*) elux,intensidadelinearx,intensidadelineary
	
		call flush(300)
		!$omp end ordered


	end do
	!$omp end do

	deallocate(exlinearx,exlineary,exlinearz,exsigmap,exsigmam)
	deallocate(fosclinearx,fosclineary,fosclinearz,foscsigmap,foscsigmam)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ffactor)

	close(100)
	close(101)
	close(102)
	close(103)

	close(200)


	close(300)




	

end subroutine bserawpol


