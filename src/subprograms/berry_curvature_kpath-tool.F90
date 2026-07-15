
!gfortran berry_curvature_kpath-tool.f90 -o berry_ct-tool.x -llapack95 -lopenblas -fopenmp

subroutine berrycurv(nthreads,dft,outputfolder,params,kpaths,sme,nocpf,fermishift,exc,mag)

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

	complex :: bxy,bxz,byz
	real,parameter :: gammas= 1.0E-37

	integer,allocatable,dimension(:) :: nocpk 

	!variaveis kpoints

	integer :: nks,nocp1
	integer :: nkpts
	real,allocatable,dimension(:,:) :: ks

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
	character(len=5) :: coultype
	real,dimension(3) :: ediel,mag
	character(len=1) :: dft
	
	integer :: nocpf
	real :: fermishift
	
	real,allocatable,dimension(:) :: gwcor
	integer :: gwin
	real :: gwaux
	character(len=1) :: gwchar	

	!fim modificacoes versao 2.1

	!call input_read


	OPEN(UNIT=202, FILE= kpaths,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening kpath input file"

	OPEN(UNIT=304, FILE=trim(outputfolder)//"gw_qp_energy_cor_avg.dat",STATUS='old', IOSTAT=gwin)
    	

	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"berry_curv_kpath.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening berry_curv_kpath output file"

	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	allocate(ovp(nvec,w90basis,w90basis))

	call hamiltonian_input_read(200,params)
	

	end select

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


	call kpath(outputfolder,rlat(1,:),rlat(2,:),rlat(3,:),nks,ks,nkpts,kpts)

	allocate(eigvf(nkpts*(nks-1),w90basis),vector(nkpts*(nks-1),w90basis,w90basis))


    	if (gwin .eq. 0) then
    	
    		allocate(gwcor(w90basis))
    		
    			read(304,*) gwchar
    		   		
    		do i=1,w90basis
    		
    			read(304,*) gwcor(i),gwaux,gwaux,gwaux,gwaux
    			
    		
    		end do
    		
    		write(2077,*) "G0W0 correction applied in kpath Berry curvature"
    	
    	else
    	
    		continue
    	
    	end if  

	!$omp parallel do default(shared) private(i,k,eigv,autovetores)
	!do i=1,nkpts*(nks-1)
	do i=1,(nks/2)*nkpts

#ifdef MKL
		call MKL_SET_NUM_THREADS(1)
#endif		

#ifdef AOCL		
		call bli_thread_set_num_threads(1)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(1)
		
#endif
	
		call eigsys(nthreads,dft,systype,scs,exc,nocpk(i),ffactor,kpts(i,2),kpts(i,3),kpts(i,4),w90basis,nvec,&
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

				if (gwin .eq. 0) then

				  eigvf(i,k)=eigv(k)+gwcor(k)
				
				else
				
				  eigvf(i,k)=eigv(k)

				end if

				do m=1,w90basis

				vector(i,k,m)=autovetores(k,m)

				end do

				
			end do

			

	end do
	!$omp end parallel do
	!termino definicao kpath
	
	write(300,*) "kp yz xz xy"

	if (gwin .eq. 0) then
    	
    		deallocate(gwcor)
        end if	

	
	!do j=1,nkpts*(nks-1)
	do j=1,(nks/2)*nkpts

	        kp= kpts(j,1)
		kx= kpts(j,2)
		ky= kpts(j,3)
		kz= kpts(j,4)


		call berryct2(nocpk(j),kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		  	ihopmatrices,efermi,eigvf(j,:),vector(j,:,:),nthreads,bxy,bxz,byz)

		
			 write(300,"(4F15.4)") kp,aimag(byz),aimag(bxz),aimag(bxy)
		


	end do
	




	deallocate(rvec,hopmatrices,ihopmatrices,eigv)
	deallocate(ovp)
	deallocate(ffactor)
	deallocate(kpts)
	deallocate(autovetores)
	deallocate(eigvf,vector)
	deallocate(nocpk)



	close(200)
	close(202)
	close(300)
	close(304)


end subroutine berrycurv
