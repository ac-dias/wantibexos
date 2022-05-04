
!gfortran berry_curvature_kpath-tool.f90 -o berry_ct-tool.x -llapack95 -lopenblas -fopenmp

subroutine berrycurvbz(nthreads,outputfolder,params,sme,ngrid,mshift)

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

	double precision,allocatable,dimension(:,:) :: output

	!double precision,dimension(3) :: shift !shift no mhkpack

	!variaveis kpoints
	integer :: nk
	integer :: nks
	integer :: nkpts
	double precision,allocatable,dimension(:,:) :: ks

	integer,allocatable,dimension(:) :: nocpk 

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


	!OPEN(UNIT=202, FILE= kpaths,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada kpath"

	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"berry_curv_bz.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening berry_curv_bz output file"

	call OMP_SET_NUM_THREADS(nthreads)
	call hamiltonian_input_read(200,params)

	!ediel(2) = edielh

	!read(202,*) nks
	!read(202,*) nkpts

	!allocate(ks(nks,3))

	!do i=1,nks

	!	read(202,*) ks(i,1),ks(i,2),ks(i,3)
	
	!end do


	!termino leitura parametros

	!parametros do calculo


	!termino parametros calculo

 
	allocate(autovetores(w90basis,w90basis),eigv(w90basis))


	allocate(kpts(ngrid(1)*ngrid(2)*ngrid(3),3))

	!shift = 0.0
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpts)

	!definindo kpath


	allocate(eigvf(ngrid(1)*ngrid(2)*ngrid(3),w90basis),vector(ngrid(1)*ngrid(2)*ngrid(3),w90basis,w90basis))

	allocate(nocpk(ngrid(1)*ngrid(2)*ngrid(3)))

	! $omp parallel do
	do i=1,ngrid(1)*ngrid(2)*ngrid(3)

		call eigsys(nthreads,scs,exc,nocpk(i),ffactor,kpts(i,1),kpts(i,2),kpts(i,3),w90basis,nvec,rlat,rvec,hopmatrices,&
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

	allocate(output(ngrid(1)*ngrid(2)*ngrid(3),9))
	
	write(300,*) "kx ky kz yz xz xy"	

	!$omp parallel do private(j,bxx,bxy,bxz,byy,byz,bzz)    

	do j=1,ngrid(1)*ngrid(2)*ngrid(3)




		call berryct2(nocpk(j),kpts(j,1),kpts(j,2),ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		  	ihopmatrices,efermi,eigvf(j,:),vector(j,:,:),nthreads,gammas,bxx,bxy,bxz,byy,byz,bzz)
		! $omp ordered
			 !write(300,*) kpts(j,1),kpts(j,2),aimag(berry)
		! $omp end ordered
			output(j,1) = kpts(j,1)
			output(j,2) = kpts(j,2)
			output(j,3) = kpts(j,3)
			output(j,4) = aimag(bxx)
			output(j,5) = aimag(bxy)
			output(j,6) = aimag(bxz)
			output(j,7) = aimag(byy)
			output(j,8) = aimag(byz)
			output(j,9) = aimag(bzz)
	

	end do
	!$omp end parallel do


	do j=1,ngrid(1)*ngrid(2)*ngrid(3)
		write(300,"(6F15.4)") real(output(j,1)),real(output(j,2)),real(output(j,3)),output(j,8),output(j,6),output(j,5)
	end do


	deallocate(rvec,hopmatrices,ihopmatrices,eigv)
	deallocate(ffactor)
	deallocate(kpts)
	deallocate(autovetores)
	deallocate(eigvf,vector)
	deallocate(output)


	close(200)
	close(201)
	!close(202)

	close(300)


end subroutine berrycurvbz
