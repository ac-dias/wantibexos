
!gfortran berry_curvature_kpath-tool.f90 -o berry_ct-tool.x -llapack95 -lopenblas -fopenmp

subroutine berrycurv(nthreads,outputfolder,params,kpaths,sme)

	use omp_lib
	use hamiltonian_input_variables

	implicit none

	integer :: i,j,k,erro,j2,m
	double precision,parameter :: pi=acos(-1.)


	double complex,allocatable,dimension(:,:) :: htb
	double precision,allocatable,dimension(:) :: eigv

	double complex,allocatable,dimension(:,:) :: autovetores

	double precision,allocatable,dimension(:,:) :: kpts

	double precision :: kx,ky,kz,kp

	double precision,allocatable,dimension(:,:) :: eigvf

	double complex,allocatable,dimension(:,:,:) :: vector

	double complex :: bxx,bxy,bxz,byy,byz,bzz
	double precision,parameter :: gammas= 0.001

	integer,allocatable,dimension(:) :: nocpk 

	!variaveis kpoints

	integer :: nks,nocp1
	integer :: nkpts
	double precision,allocatable,dimension(:,:) :: ks

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

	!fim modificacoes versao 2.1

	!call input_read


	OPEN(UNIT=202, FILE= kpaths,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening kpath input file"

	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"berry_curv_kpath.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening berry_curv_kpath output file"

	call hamiltonian_input_read(200,params)

	!ediel(2) = edielh


	read(202,*) nks
	read(202,*) nkpts

	allocate(ks(nks,3))

	do i=1,nks

		read(202,*) ks(i,1),ks(i,2),ks(i,3)
	
	end do



	!termino leitura parametros

	!parametros do calculo


	!termino parametros calculo 

	allocate(autovetores(w90basis,w90basis),eigv(w90basis))
	
	!allocate(kpts(nkpts*(nks-1),4))
	!allocate(nocpk(nkpts*(nks-1)))
	
	allocate(kpts((nks/2)*nkpts,4))
	allocate(nocpk((nks/2)*nkpts))

	!definindo kpath


	call kpath2(outputfolder,rlat(1,:),rlat(2,:),rlat(3,:),nks,ks,nkpts,kpts)

	allocate(eigvf(nkpts*(nks-1),w90basis),vector(nkpts*(nks-1),w90basis,w90basis))


	! $omp parallel do
	!do i=1,nkpts*(nks-1)
	do i=1,(nks/2)*nkpts

		call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpts(i,2),kpts(i,3),kpts(i,4),w90basis,nvec,rlat,rvec,hopmatrices,&
		   ihopmatrices,efermi,eigv,autovetores)

			do k=1,w90basis


				eigvf(i,k)=eigv(k)

				do m=1,w90basis

				vector(i,k,m)=autovetores(k,m)

				end do

				
			end do

			

	end do
	! $omp end parallel do
	!termino definicao kpath
	
	write(300,*) "kp yz xz xy"

	!$omp do ordered
	!do j=1,nkpts*(nks-1)
	do j=1,(nks/2)*nkpts

	        kp= kpts(j,1)
		kx= kpts(j,2)
		ky= kpts(j,3)
		kz= kpts(j,4)


		call berryct2(nocpk(j),kx,ky,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		  	ihopmatrices,efermi,eigvf(j,:),vector(j,:,:),nthreads,gammas,bxx,bxy,bxz,byy,byz,bzz)

		!$omp ordered
			 write(300,"(4F15.4)") kp,aimag(byz),aimag(bxz),aimag(bxy)
		!$omp end ordered


	end do
	!$omp end do




	deallocate(rvec,hopmatrices,ihopmatrices,eigv)
	deallocate(ffactor)
	deallocate(kpts)
	deallocate(autovetores)
	deallocate(eigvf,vector)
	deallocate(nocpk)



	close(200)
	close(202)
	close(300)


end subroutine berrycurv
