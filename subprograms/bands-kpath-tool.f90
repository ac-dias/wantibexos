

!gfortran bands-kpath-tool.f90 -o bands_tool.x -llapack95 -lopenblas -fopenmp

subroutine bandstool(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,rk,dft)

	use omp_lib
	use hamiltonian_input_variables
	implicit none

	integer :: i,j,k,erro,j2
	double precision,parameter :: pi=acos(-1.)

	!double complex,allocatable,dimension(:,:) :: htb
	double precision,allocatable,dimension(:) :: eigv

	double precision,allocatable,dimension(:,:,:) :: ebands

	double complex,allocatable,dimension(:,:) :: autovetores

	double precision,allocatable,dimension(:,:) :: kpts

	double precision,allocatable,dimension(:) :: orbweight

	double precision :: kx,ky,kz,kp

	double precision :: lco

	double complex :: spz

	!variaveis kpoints

	integer :: nks
	integer :: nkpts
	double precision,allocatable,dimension(:,:) :: ks

	!call input_read

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
	character(len=1) :: dft	

	!fim modificacoes versao 2.1

	! INPUT : lendo os parametros do modelo de tight-binding
	
	OPEN(UNIT=202, FILE= kpaths,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening kpath input file"


	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"bands.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bands output file"

	!OPEN(UNIT=500, FILE=trim(outputfolder)//"nocp.dat",STATUS='unknown', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de saida nocp"

	call OMP_SET_NUM_THREADS(nthreads)

	call hamiltonian_input_read(200,params)

	!ediel(2) = edielh


	

	read(202,*) nks
	read(202,*) nkpts

	allocate(ks(nks,3))

	do i=1,nks

		read(202,*) ks(i,1),ks(i,2),ks(i,3)
	
	end do

	allocate(orbweight(w90basis))

	OPEN(UNIT=203, FILE= orbw,STATUS='old', IOSTAT=erro)
    	  	if (erro/=0) then
	
	orbweight= 1.0

	else

	do i=1,w90basis
		read(203,*) orbweight(i)
	end do

	end if



	!termino leitura parametros

	!parametros do calculo

	!OPEN(UNIT=2077, FILE= trim(calcparms)//"bands_calc.dat",STATUS='unknown', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de saida parametros calculo"

	!call param_out(2077,nthreads,outputfolder,calcparms,ngrid,nc,nv,edos0,edosf,numdos, &
	!	     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
	!	     exc,mshift,coultype)

	!close(2077)

	!termino parametros calculo
 

	allocate(autovetores(w90basis,w90basis),eigv(w90basis))
	allocate(ebands(w90basis,nkpts*(nks-1),4))
	!allocate(kpts(nkpts*(nks-1),4))
	allocate(kpts((nks/2)*nkpts,4))

	!definindo kpath

	call kpath2(outputfolder,rlat(1,:),rlat(2,:),rlat(3,:),nks,ks,nkpts,kpts)

	!termino definicao kpath

	ebands=0.0

	!do j=1,nkpts*(nks-1)
	do j=1,(nks/2)*nkpts

	        kp= kpts(j,1)
		kx= kpts(j,2)
		ky= kpts(j,3)
		kz= kpts(j,4)

		call eigsys(nthreads,scs,exc,nocp,ffactor,kx,ky,kz,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,efermi,eigv,autovetores)

		!write(500,*) nocp,kx,ky,kz 


		do i=1,w90basis

			call layercont(w90basis,autovetores(i,:),orbweight,lco)
			call spinvl(w90basis,autovetores(i,:),dft,systype,spz)

			ebands(i,j,1) = kp
			ebands(i,j,2) = eigv(i)
			ebands(i,j,3) = real(spz)
			ebands(i,j,4) = lco
			


		end do





	end do



	do i=1,w90basis

			write(300,*) "#band",i

		!do j=1,nkpts*(nks-1)
		do j=1,(nks/2)*nkpts

			write(300,*) real(ebands(i,j,1)),real(ebands(i,j,2)),real(ebands(i,j,3)),real(ebands(i,j,4))

		end do

			write(300,*)

	end do

	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)

	deallocate(eigv)
	deallocate(kpts)
	deallocate(orbweight)
	deallocate(autovetores)
	deallocate(ebands)
	deallocate(ks)


	close(200)
	close(300)
	close(202)
	close(203)

	!close(500)


end subroutine bandstool
