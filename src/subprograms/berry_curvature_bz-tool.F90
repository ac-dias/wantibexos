
!gfortran berry_curvature_kpath-tool.f90 -o berry_ct-tool.x -llapack95 -lopenblas -fopenmp

subroutine berrycurvbz(nthreads,dft,outputfolder,params,sme,ngrid,mshift,nocpf,fermishift,meshgen,exc,mag)

	use omp_lib
	use hamiltonian_input_variables

	implicit none

	integer :: i,j,k,erro,j2,m
	real,parameter :: pi=acos(-1.)


	complex,allocatable,dimension(:,:) :: htb
	real,allocatable,dimension(:) :: eigv

	complex,allocatable,dimension(:,:) :: autovetores

	real,allocatable,dimension(:,:) :: kpts

	real :: kx,ky,kz,kp

	real,allocatable,dimension(:,:) :: eigvf

	complex,allocatable,dimension(:,:,:) :: vector

	complex ::bxy,bxz,byz
	real,parameter :: gammas= 1.0E-37

	!real,allocatable,dimension(:,:) :: output

	!real,dimension(3) :: shift !shift no mhkpack

	!variaveis kpoints
	integer :: nk
	integer :: nks
	integer :: nkpts
	real,allocatable,dimension(:,:) :: ks

	integer,allocatable,dimension(:) :: nocpk 
	
	real,allocatable,dimension(:,:) :: output

	!modificacoes versao 2.1

	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv
	!integer :: ncrpa,nvrpa
	!integer :: ncbz,nvbz
	real :: edos0,edosf,numdos
	real :: ebse0,ebsef,numbse
	real :: sme,exc,rk
	real,dimension(3) :: mshift
	real :: ktol
	character(len=70) :: params   !parametros TB
	character(len=70) :: orbw     !peso orbitais lcount
	character(len=70) :: kpaths    !kpath
	character(len=70) :: kpathsbse    !kpath
	!character(len=70) :: diein    !ambiente dieletrico
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: calcparms
	character(len=70) :: meshtype
	character(len=70) :: meshgen	
	character(len=5) :: coultype
	real,dimension(3) :: ediel,mag
	character(len=1) :: dft
	
	
	integer :: nocpf,ngkpt
	real :: fermishift

	!fim modificacoes versao 2.1

	!call input_read


	!OPEN(UNIT=202, FILE= kpaths,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada kpath"

	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"berry_curv_bz.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening berry_curv_bz output file"

	call OMP_SET_NUM_THREADS(nthreads)
	
	
	
	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	allocate(ovp(nvec,w90basis,w90basis))

	call hamiltonian_input_read(200,params)
	

	end select



 
	allocate(autovetores(w90basis,w90basis),eigv(w90basis))
	
	
	select case (meshgen)
	
	case ("2DHEX")
	
	call gridgen2dhexcount(ngrid(1),ngrid(2),rlat(1,1),ngkpt)
	
	allocate(kpts(ngkpt,3))
	
	call gridgen2dhex(ngrid(1),ngrid(2),rlat(1,1),ngkpt,kpts)
	
	
	case ("2DRET")
	
	ngkpt = ngrid(1)*ngrid(2)
	
	allocate(kpts(ngkpt,3))
	
	call gridgen2dret(ngrid(1),ngrid(2),rlat(1,1),rlat(2,2),kpts)
	
	case default

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	
	allocate(kpts(ngkpt,3))

	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpts)
	
	end select






	allocate(eigvf(ngkpt,w90basis),vector(ngkpt,w90basis,w90basis))

	allocate(nocpk(ngkpt))

	!$omp parallel do default(shared) private(i,k,autovetores,eigv)
	do i=1,ngkpt
	
#ifdef MKL
		call MKL_SET_NUM_THREADS(1)
#endif		

#ifdef AOCL		
		call bli_thread_set_num_threads(1)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(1)
		
#endif
	
		call eigsys(nthreads,dft,systype,scs,exc,nocpk(i),ffactor,kpts(i,1),kpts(i,2),kpts(i,3),w90basis,nvec,&
		           rlat,rvec,hopmatrices,ihopmatrices,ovp,efermi,eigv,autovetores,nocpf,fermishift,mag)
#ifdef MKL
		call MKL_SET_NUM_THREADS(nthreads)
#endif		

#ifdef AOCL		
		call bli_thread_set_num_threads(nthreads)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(nthreads)
		
#endif	

			do k=1,w90basis


				eigvf(i,k)=eigv(k)

				do m=1,w90basis

				vector(i,k,m)=autovetores(k,m)

				end do

				
			end do

			

	end do
	!$omp end parallel do

	!termino definicao kpath

	allocate(output(ngkpt,6))
	
	output= 0.0	

	!$omp parallel do default(shared) private(j,bxy,bxz,byz)    
	do j=1,ngkpt




		call berryct2(nocpk(j),kpts(j,1),kpts(j,2),kpts(j,3),ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		  	ihopmatrices,efermi,eigvf(j,:),vector(j,:,:),nthreads,bxy,bxz,byz)
		  	
		! $omp ordered
			 !write(300,*) kpts(j,1),kpts(j,2),aimag(berry)
		! $omp end ordered
			output(j,1) = kpts(j,1)
			output(j,2) = kpts(j,2)
			output(j,3) = kpts(j,3)
			!output(j,4) = aimag(bxx)
			output(j,4) = aimag(bxy)
			output(j,5) = aimag(bxz)
			!output(j,7) = aimag(byy)
			output(j,6) = aimag(byz)
			!output(j,9) = aimag(bzz)
			
			!write(300,"(6F15.4)") kpts(j,1),kpts(j,2),kpts(j,3),aimag(byz),aimag(bxz),aimag(bxy)
			!call flush(300)

	end do
	!$omp end parallel do


	write(300,*) "#kx ky kz yz xz xy"

	do j=1,ngkpt
		write(300,"(6F15.4)") output(j,1),output(j,2),output(j,3),output(j,6),output(j,5),output(j,4)
		!write(300,*) output(j,1),output(j,2),output(j,3),output(j,8),output(j,6),output(j,5)
	
	end do



	deallocate(rvec,hopmatrices,ihopmatrices,eigv)
	deallocate(ovp)
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
