
#ifdef ELPA
#ifndef ELPA_API_VERSION
#define ELPA_API_VERSION 20211125
#endif
#endif

subroutine bsesolver(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w1,r0,lc,rk,meshtype,bsewf,excwf0,excwff,&
		     dtfull,cpol,tmcoef,nocpf,fermishift,bsealgo,bsehamwrite,bsehamread,bsehamfile,dft,mag)

#ifdef MPI
	use mpi
!include 'mpif.h'
#endif
#ifdef ELPA
	use elpa
#endif
	use omp_lib
	use hamiltonian_input_variables
	use input_variables, only: bsermatfile
	use bse_q_optics, only: rmn_data, rmn_read, rmn_destroy, rmn_apply_q0_optical_correction

	implicit none

	real,parameter:: pi=acos(-1.)

	integer :: dimbse !=ngrid*ngrid*nc*nv ! dimensão da matriz bse

	real,allocatable,dimension(:,:) :: eigv
	complex,allocatable,dimension(:,:,:) :: vector

	real,allocatable,dimension(:,:) :: kpt !pontos k do grid
	real,allocatable,dimension(:,:) :: kpt_bse

	integer,allocatable,dimension(:,:) :: stt
	integer,allocatable,dimension(:,:) :: stt_bse

	real,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores

	!real,allocatable,dimension(:,:) :: lcount !pertencimento a camada de um dado estado

	!real,allocatable,dimension(:) :: orbweight

	integer:: counter,c,v,i,j,f,l,h,k,kl,erro

	integer :: ncaux,nvaux

	complex:: matrizelbse


	complex,allocatable,dimension(:) :: hrx,hry,hrz,hrsp,hrsm 

	real,allocatable,dimension(:) :: actxx,actyy,actzz,actxy,actxz,actyz,actsp,actsm

	complex,parameter :: imag=cmplx(0.0,1.0)

	real :: a,r0,ed,ec,ev,egap

	integer :: ngkpt
	real,allocatable,dimension(:,:) :: vecres
		
	!real,dimension(3) :: shift !shift no mhkpack

	!variaveis relacionadas a marcacao do tempo

	real:: t0,tf
	double precision :: task_start,task_elapsed,task_elapsed_max
	integer,dimension(8) :: values,values2

	integer,allocatable,dimension(:) :: nocpk
	
	integer :: nocpf
	real :: fermishift
	character(len=12) :: bsealgo
	character(len=70) :: bsehamfile
	logical :: bsehamwrite,bsehamread,bseham_ok
	character(len=1) :: dft
	real,dimension(3) :: mag	
	
	complex,allocatable,dimension(:,:,:) :: sk
	type(rmn_data) :: rmn
	logical :: use_rmn,rmn_ok
	character(len=256) :: rmn_message

	!definicoes diagonalizacao 
	INTEGER   ::       ifail
	integer,allocatable,dimension(:) :: ifailx
	real,parameter :: ABSTOL=1.0e-6
	INTEGER  ::        INFO
	real,allocatable,dimension(:) :: W,RWORK
		COMPLEX,allocatable,dimension(:,:) :: hbse,hbse_dist,hbse_eigenvectors
		COMPLEX,allocatable,dimension(:) :: eigvec_dist

        INTEGER :: LWMAX
   	INTEGER :: LWORK
	INTEGER :: LIWORK, LRWORK
	INTEGER,allocatable,dimension(:) :: IWORK
        complex,allocatable,dimension (:) :: WORK

#ifdef ELPA
	class(elpa_t), pointer :: elpa_instance
	complex,allocatable,dimension(:,:) :: elpa_eigenvectors
	integer :: elpa_status
#endif

        !zheevr definitions
        real :: VL,VU
        integer :: IL,IU,M
        COMPLEX,allocatable,dimension(:,:) :: Z
        integer :: LDZ
        integer,allocatable,dimension(:) :: ISUPPZ

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
	character(len=7) :: coultype
	real,dimension(3) :: ediel
	logical :: bsewf
	integer :: excwf0,excwff
	real :: ez,w1,lc
	
	logical :: cpol,dtfull
	logical :: tmcoef	

	integer :: MPIError, Node, Nodes
	integer :: blacs_ctxt, nprow, npcol, myrow, mycol
	integer :: mb, nb, locr, locc, lld, ig, jg, li, lj
	integer :: desca(9), descz(9)
	integer :: numroc, indxl2g
	integer,dimension(3) :: bseham_metadata
	character(len=160) :: bseham_path
#ifdef MPI
	integer,parameter :: bse_mpi_block_elements=16777216
	integer(kind=8) :: mpi_transfer_elements
#endif

	!fim modificacoes versao 2.1

	!call input_read

	! INPUT : lendo os parametros do modelo de tight-binding
	!OPEN(UNIT=203, FILE= orbw,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada orb weight"

#ifdef MPI
	call MPI_COMM_RANK(MPI_COMM_WORLD, Node, MPIError)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, Nodes, MPIError)
#else
	Node = 0
	Nodes = 1
	MPIError = 0
#endif

	!OUTPUT
    
    if (Node == 0) then
	OPEN(UNIT=300, FILE=trim(outputfolder)//"log_bse_optics.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_bse-diel output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"bse_oscf.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_oscf.dat output file"
    	OPEN(UNIT=302, FILE=trim(outputfolder)//"bse_oscf-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse_oscf-pol.dat output file"	


	OPEN(UNIT=401, FILE=trim(outputfolder)//"ipa_oscf.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening ipa_oscf.dat output file"
    	
    	
    	if (tmcoef) then
    	OPEN(UNIT=402, FILE=trim(outputfolder)//"tm_coef.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef output file"
    	OPEN(UNIT=404, FILE=trim(outputfolder)//"tm_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef-pol output file"    	
    	else
    	 continue
    	end if    
    	
    	OPEN(UNIT=403, FILE=trim(outputfolder)//"ipa_oscf-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening ipa_oscf-pol.dat output file"	

    endif

	!OPEN(UNIT=500, FILE=trim(outputfolder)//"log_bse-matrix.dat",STATUS='unknown', IOSTAT=erro)
    	!if (erro/=0) stop "Error opening bse hamiltonian matrix progress output file"

	call cpu_time(t0)
	call date_and_time(VALUES=values)
	call OMP_SET_NUM_THREADS(nthreads)
!#ifdef gnu
!	call OPENBLAS_SET_NUM_THREADS(nthreads)
!#endif


    if (Node == 0) then
	!inicio leitura parametros

	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	allocate(ovp(nvec,w90basis,w90basis))

	call hamiltonian_input_read(200,params)
	

	end select

	!ediel(2) = edielh


	!termino leitura parametros

	!parametros do calculo


	!dtfull = .false.
	!cpol = .true.


	write(301,*) "#","  ", "exciton energy","  ","xx","  ","yy","  ","zz"," ","xy","  ","xz","  ","yz"
	write(302,*) "#","  ", "exciton energy","  ","xx","  ","yy","  ","zz"," ","sp","  ","sm"	

    endif

#ifdef MPI
      call mpi_barrier(MPI_COMM_WORLD,MPIError)
      call bcast_hamil()
#endif

		use_rmn=len_trim(bsermatfile) > 0
		if (use_rmn) then
			if (dft == 'S') stop 'BSE_RMAT_FILE currently requires the orthonormal Wannier representation'
			call rmn_read(trim(bsermatfile),rmn,rmn_ok,rmn_message)
			if (.not. rmn_ok) then
				write(*,*) trim(rmn_message)
				stop 'Unable to initialize the Q=0 Wannier position matrix'
			end if
			if (rmn%num_wann /= w90basis) then
				write(*,*) 'Wannier position matrix basis:',rmn%num_wann,' Hamiltonian basis:',w90basis
				stop 'Incompatible Wannier position and Hamiltonian bases'
			end if
			if (Node == 0) then
				write(300,*) 'Q=0 optical position matrix:',trim(bsermatfile)
				write(300,*) 'Q=0 position-matrix treatment: dH/dk - i[A,H]'
			end if
		end if

	!termino parametros calculo 

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimbse = ngkpt*nc*nv

	!call alat(systype,rlat,a)

	!if (systype .eq. "2D") then
	!	ed = (ediel(1)+ediel(3))/2.
	!	r0= ((ediel(2)-1.0)*lc)/(ediel(1)+ediel(3))
	!else
	!	r0 = 1.0
	!	ed = ediel(2)
		
	!end if

    if (Node == 0) then
	!Informações para o arquivo de log do calculo
	write(300,*)
	write(300,*)
	write(300,*)
	write(300,*) 'threads:', nthreads
	write(300,*)
	write(300,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
	write(300,*)
	write(300,"(A14,1E15.4)") 'ktol-coulomb:', ktol
	write(300,*)
	write(300,"(A13,3F15.4)") 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	!write(300,*) 'angulo theta:',beta,'rad'

	write(300,*)
	write(300,*) 'conduction bands, nc:',nc,'  ','valence bands, nv:',nv
	write(300,*)
	write(300,*) 'Coulomb Potential:',coultype
	write(300,*)
	write(300,*) 'BSE_ALGO:',bsealgo

	!write(300,*)
	!write(300,*) 'intensidade termo exchange:',exc
	!write(300,*)
	!write(300,*) 'intensidade termo exchange l:',excl
	!write(300,*)
	!write(300,*) 'intensidade campo eletrico:',ez
	!write(300,*)
	!write(300,*) 'direção da magnetização do termo de exchange:',mag(1),'x',mag(2),'y',mag(3),'z'
	!write(300,*)

	!write(300,*) 'arquivo de saida 1:',logbse
	!write(300,*) 'arquivo de saida 2:',bseoptx
	!write(300,*) 'arquivo de saida 3:',bseopty
	!write(300,*) 'arquivo de saida 4:',bseoptsp
	!write(300,*) 'arquivo de saida 5:',bseoptsm

	
	write(300,*)
	write(300,*) 'begin','  ','month',values(2),'day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(300,*) 

	call flush(300)
    endif

    allocate(kpt(ngkpt,3))

	!shift = 0.0
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)

	! Keep the orbital index first so every eigenvector passed to the BSE
	! kernel is a contiguous column.
	allocate(eigv(ngkpt,nc+nv),vector(w90basis,nc+nv,ngkpt))
	allocate(nocpk(ngkpt))

	allocate (stt(dimbse,4))
	allocate(hrx(dimbse),hry(dimbse),hrz(dimbse),hrsp(dimbse),hrsm(dimbse))

	allocate(actxx(dimbse),actyy(dimbse),actzz(dimbse),actxy(dimbse))
	allocate(actxz(dimbse),actyz(dimbse))
	allocate(actsp(dimbse),actsm(dimbse))	
	
	
	if (Nodes == 1) then
	 select case (bsealgo)
		
		case ("cheev")
		
		LWMAX = 2*dimbse-1
				
		LWORK = -1

		allocate (RWORK(3*dimbse-2))
		allocate(W(dimbse),WORK(LWMAX))

		case ("cheevr")
		
		!LWMAX = (1 + 5*dimbse + 2*dimbse**2)
				
		LWORK = -1
      		LIWORK = -1
      		LRWORK = -1
		
		allocate(W(dimbse))
		allocate(WORK(2*dimbse))
		allocate(ISUPPZ(2*dimbse))
		allocate(IWORK(10*dimbse))
		allocate(RWORK(24*dimbse))
		allocate(Z(dimbse,dimbse)) 
		
		case ("cheevx")
		
		!LWMAX = (1 + 5*dimbse + 2*dimbse**2)
				
		LWORK = -1
 
		
		allocate(W(dimbse))
		allocate(WORK(2*dimbse))
		allocate(IWORK(5*dimbse))
		allocate(RWORK(7*dimbse))
		allocate(Z(dimbse,dimbse)) 
		allocate(ifailx(dimbse))
		
		case ("cheevd")
		
		LWMAX = (2*dimbse-1)**2
		
		LWORK = -1
      		LIWORK = -1
      		LRWORK = -1
      		!LWORK = 2*dimbse + dimbse**2
      		!LIWORK = 3 + 5*dimbse
      		!LRWORK = 1 + 5*dimbse + 2*dimbse**2
		allocate(W(dimbse),WORK(2*dimbse + dimbse**2))
		allocate (IWORK(3 + 5*dimbse))
		allocate (RWORK(1 + 5*dimbse + 2*dimbse**2))

		case ("elpa")

		write(*,*) "BSE_ALGO=elpa requires an MPI run with at least two ranks"
		stop
		
		case default
		
		 write(*,*) "please select a subroutine for diagonalization"
		 stop
		
	 end select

	    endif ! end the check if we are going to perform a parallel calculation

	    if (Nodes > 1) then
	    	allocate(W(dimbse))
	    end if

	!allocate(orbweight(w90basis))
	
	!do i=1,w90basis
	!	read(203,*) orbweight(i)
	!end do

	!allocate(lcount(ngkpt,w90basis))

	egap = 50.0

    if (Node == 0) then ! I am going to run this in serial for now
     allocate(eaux(w90basis),vaux(w90basis,w90basis))
     write(300,*) "We are going to perform independent particle calculations in serial for now"

	!$omp parallel do default(shared) private(i,j,l,h,eaux,vaux)
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
	
		call eigsys(nthreads,dft,systype,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,ovp,efermi,eaux,vaux,nocpf,fermishift,mag)
#ifdef MKL
		call MKL_SET_NUM_THREADS(nthreads)
#endif		

#ifdef AOCL		
		call bli_thread_set_num_threads(nthreads)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(nthreads)
		
#endif	

			do j=1,nc+nv
	
				eigv(i,j)= eaux(nocpk(i)-nv+j)

			end do
			
			!if (ntype .eq. 2) then

			!do kl = 1,w90basis

			!	call layercont(w90basis,vaux(kl,:),orbweight,lcount(i,kl))


			!end do

			!else 

			!	continue

			!end if

			if ((eigv(i,nv+1)-eigv(i,nv)) .le. egap) then

				egap = eigv(i,nv+1)-eigv(i,nv)

			else

				continue

			end if

			


			do l=1,nc+nv


				do h=1,w90basis

					vector(h,l,i)=vaux(nocpk(i)-nv+l,h)


				end do
			

			end do

	
	end do
	!$omp end parallel do



	deallocate(eaux,vaux)
	write(300,*) 'direct gap:', egap
	write(300,*) 'eigenvalues and eigenvectors calculated'
	call flush(300)
    endif

#ifdef MPI
	call MPI_BCAST(nocpk, ngkpt, MPI_INTEGER, 0, MPI_COMM_WORLD, MPIError)
	if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('single-particle band indices',MPIError)
	mpi_transfer_elements=int(ngkpt,kind=8)*int(nc+nv,kind=8)
	call bse_mpi_bcast_real_blocks(eigv,mpi_transfer_elements,0,MPI_COMM_WORLD,MPIError,bse_mpi_block_elements)
	if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('single-particle eigenvalues',MPIError)
	mpi_transfer_elements=mpi_transfer_elements*int(w90basis,kind=8)
	call bse_mpi_bcast_complex_blocks(vector,mpi_transfer_elements,0,MPI_COMM_WORLD,MPIError,bse_mpi_block_elements)
	if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('single-particle eigenvectors',MPIError)
	call MPI_BCAST(egap, 1, MPI_REAL, 0, MPI_COMM_WORLD, MPIError)
	if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('direct gap',MPIError)
#endif

	!write(*,*) "autovetores e autovalores"


	! Store each k-point overlap matrix contiguously: sk(:,:,ik).
	allocate(sk(w90basis,w90basis,ngkpt))
	if (dft .eq. "S") then
		
		do i=1,ngkpt
		
			call overlap(w90basis,nvec,rvec,ovp,kpt(i,1),kpt(i,2),kpt(i,3),sk(:,:,i))
		
		end do
		
		if (Node == 0) then
			write(300,*) 'overlap matrices calculated'
			call flush(300)
		end if
	else
		sk = 0.0
		do i=1,ngkpt
			do j=1,w90basis
				sk(j,j,i) = 1.0
			end do
		end do
	end if


	!definindo os numeros quanticos dos estados


	!allocate (stt(ngkpt*nc*nv,4))
    
	call quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpk,stt)
	!allocate(hrx(dimbse),hry(dimbse),hrz(dimbse))

    allocate(vecres(dimbse,17))

    if (Node == 0) then
	write(300,*) 'quantum numbers for exciton basis set finished'
	call flush(300)

	write(401,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(403,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","sp","  ","sm"
	
	if (tmcoef) then	
	write(402,*) "number of kpoints:",ngkpt
	write(402,*) "number of conduction states",nc
	write(402,*) "number of valence states",nv
	write(402,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ","ec "," ","ev"," ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","xy","  ","xz","  ","yz"
		      
	write(404,*) "number of kpoints:",ngkpt
	write(404,*) "number of conduction states",nc
	write(404,*) "number of valence states",nv
	write(404,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ","ec "," ","ev","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","sp","  ","sm"		   
		      
	else
	
	 continue
	end if	
    endif

		if (Node == 0) then
			write(300,*) 'IPA transitions: begin'
			call flush(300)
		end if
#ifdef MPI
		call MPI_BARRIER(MPI_COMM_WORLD,MPIError)
		task_start = MPI_WTIME()
#else
		call cpu_time(task_start)
#endif

		 ! $omp parallel default(shared) private(i,w90basis,rlat,rvec,hopmatrices,ihopmatrices)

			 !$omp parallel do private(rmn_ok,rmn_message)
	do i=1,dimbse

		!ec = eigv(stt(i,4),stt(i,3))
		!ev = eigv(stt(i,4),stt(i,2))

			call optsp(eigv(stt(i,4),stt(i,2)),vector(:,stt(i,2),stt(i,4)),&
			     eigv(stt(i,4),stt(i,3)),vector(:,stt(i,3),stt(i,4)),&
			     kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),ffactor,sme,&
			     w90basis,nvec,rlat,rvec,hopmatrices,&
			     ihopmatrices,hrx(i),hry(i),hrz(i))
			if (use_rmn) then
				call rmn_apply_q0_optical_correction(rmn,kpt(stt(i,4),:),rlat,&
					eigv(stt(i,4),stt(i,2)),vector(:,stt(i,2),stt(i,4)),&
					eigv(stt(i,4),stt(i,3)),vector(:,stt(i,3),stt(i,4)),sme,&
					hrx(i),hry(i),hrz(i),rmn_ok,rmn_message)
				if (.not. rmn_ok) stop 'Unable to evaluate the Q=0 Wannier position matrix'
			end if
		     
		     hrsp(i) = (hrx(i)+cmplx(0.,1.)*hry(i))*(1.0/sqrt(2.))
		     hrsm(i) = (hrx(i)-cmplx(0.,1.)*hry(i))*(1.0/sqrt(2.))
		     
		vecres(i,1) = real(kpt(stt(i,4),1))
		vecres(i,2) = real(kpt(stt(i,4),2))
		vecres(i,3) = real(kpt(stt(i,4),3))
		vecres(i,4) = real(nocpk(stt(i,4)))
		vecres(i,5) = real(nocpk(stt(i,4))-nv+stt(i,3))
		vecres(i,6) = real(nocpk(stt(i,4))-nv+stt(i,2))		
		vecres(i,7) = real(eigv(stt(i,4),stt(i,3))-eigv(stt(i,4),stt(i,2)))		
		vecres(i,8) = real(hrx(i)*conjg(hrx(i)))		
		vecres(i,9) = real(hry(i)*conjg(hry(i)))		
		vecres(i,10) = real(hrz(i)*conjg(hrz(i)))
		vecres(i,11) = real(hrx(i)*conjg(hry(i)))
		vecres(i,12) = real(hrx(i)*conjg(hrz(i)))
		vecres(i,13) = real(hry(i)*conjg(hrz(i)))		 
		vecres(i,14) = real(hrsp(i)*conjg(hrsp(i)))
		vecres(i,15) = real(hrsm(i)*conjg(hrsm(i)))				     
		vecres(i,16) = real(eigv(stt(i,4),stt(i,3)))
		vecres(i,17) =	real(eigv(stt(i,4),stt(i,2)))

	end do
	 !$omp end parallel do


	 ! $omp end parallel
	
		call Bubblem(7,17,vecres, dimbse)

#ifdef MPI
		task_elapsed = MPI_WTIME() - task_start
		call MPI_REDUCE(task_elapsed,task_elapsed_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,MPIError)
#else
		call cpu_time(task_elapsed)
		task_elapsed_max = task_elapsed - task_start
#endif

	    if (Node == 0) then
	do i=1,dimbse
	
		write(401,"(7F15.6)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,11),vecres(i,12),vecres(i,13)
		write(403,"(6F15.6)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,14),vecres(i,15)
		
		if (tmcoef) then
	        write(402,"(3F10.6,3I10.0,9F10.6)")  vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,16),vecres(i,17),vecres(i,7),&
						  vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,11),vecres(i,12),vecres(i,13)
	        
	        write(404,"(3F10.6,3I10.0,8F10.6)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,16),vecres(i,17),vecres(i,7),vecres(i,8),&
						  vecres(i,9),vecres(i,10),&
						  vecres(i,14),vecres(i,15)						  
						  
		end if
	
	end do	

		write(300,*) 'IPA transitions: finished'
		write(300,"(A,F12.3,A)") 'IPA transitions wall time: ',task_elapsed_max,' s'
		write(300,*) 'IPA single particle optics finished'
	call flush(300)
    endif

	deallocate(vecres)
	! These compact views make the remaining BSE-kernel arguments contiguous
	! without duplicating the eigenvector storage.
	allocate(stt_bse(4,dimbse),kpt_bse(3,ngkpt))
	do i=1,dimbse
		do j=1,4
			stt_bse(j,i) = stt(i,j)
		end do
	end do
	do i=1,ngkpt
		do j=1,3
			kpt_bse(j,i) = kpt(i,j)
		end do
	end do
	
	!go to 789


    if (Node == 0) then
	 if (bsehamread) then
		write(300,*) 'BSE Hamiltonian restart read: begin'
	 else
		write(300,*) 'BSE Hamiltonian construction: begin'
	 end if
	 call flush(300)
    endif
#ifdef MPI
	call MPI_BARRIER(MPI_COMM_WORLD,MPIError)
	task_start = MPI_WTIME()
#else
	call cpu_time(task_start)
#endif

	if (Node == 0) then
		call bse_hamiltonian_memory_report(300,'BSE Hamiltonian, global dense equivalent',dimbse,dimbse)
	end if

    if (Nodes == 1) then
		bseham_metadata = (/ dimbse,1,1 /)
		bseham_path = trim(outputfolder)//trim(bsehamfile)
		allocate(hbse(dimbse,dimbse))
		if (bsehamread) then
			call bse_hamiltonian_read(bseham_path,bseham_metadata,dimbse,dimbse,hbse,bseham_ok)
			if (.not. bseham_ok) then
				write(*,*) 'Unable to read a compatible BSE Hamiltonian: ',trim(bseham_path)
				stop
			end if
		else
			hbse=0.0
			!$omp parallel do private(i,j) schedule(dynamic)
			do j=1,dimbse
				do i=1,j
					hbse(i,j)= matrizelbse(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,rlat,stt_bse(:,i),eigv(stt_bse(4,i)&
					    ,stt_bse(3,i)),eigv(stt_bse(4,i),stt_bse(2,i)),vector(:,stt_bse(3,i),stt_bse(4,i)),&
					    vector(:,stt_bse(2,i),stt_bse(4,i)),kpt_bse(:,stt_bse(4,i)),stt_bse(:,j),eigv(stt_bse(4,j),stt_bse(3,j))&
					    ,eigv(stt_bse(4,j),stt_bse(2,j)) &
					    ,vector(:,stt_bse(3,j),stt_bse(4,j)),vector(:,stt_bse(2,j),stt_bse(4,j)),kpt_bse(:,stt_bse(4,j)),dft,nvec,rvec,&
					    sk(:,:,stt_bse(4,i)),sk(:,:,stt_bse(4,j)))
				end do
			end do
			!$omp end parallel do
			if (bsehamwrite) then
				call bse_hamiltonian_write(bseham_path,bseham_metadata,dimbse,dimbse,hbse,bseham_ok)
				if (.not. bseham_ok) then
					write(*,*) 'Unable to save BSE Hamiltonian: ',trim(bseham_path)
					stop
				end if
			end if
		end if
    else
#ifdef MPI
		! hbse_dist is assembled and diagonalized locally on the BLACS grid.
		! Its N**2 size is never passed as one MPI count; transfers use blocks.
		mb = 64
		nb = 64
		nprow = int(sqrt(real(Nodes)))
		do while (mod(Nodes,nprow) /= 0)
			nprow = nprow - 1
		end do
		npcol = Nodes/nprow

		call BLACS_GET(-1, 0, blacs_ctxt)
		call BLACS_GRIDINIT(blacs_ctxt, 'R', nprow, npcol)
		call BLACS_GRIDINFO(blacs_ctxt, nprow, npcol, myrow, mycol)

		locr = numroc(dimbse, mb, myrow, 0, nprow)
		locc = numroc(dimbse, nb, mycol, 0, npcol)
		lld = max(1,locr)
		if (Node == 0) then
			call bse_hamiltonian_memory_report(300,'BSE Hamiltonian on MPI rank 0',lld,max(1,locc))
		end if
		allocate(hbse_dist(lld,max(1,locc)))
		call DESCINIT(desca, dimbse, dimbse, mb, nb, 0, 0, blacs_ctxt, lld, INFO)
		call DESCINIT(descz, dimbse, dimbse, mb, nb, 0, 0, blacs_ctxt, lld, INFO)
		bseham_metadata = (/ dimbse,1,1 /)
		if (trim(bsealgo) == 'elpa') bseham_metadata(3) = 2
		bseham_path = trim(outputfolder)//trim(bsehamfile)

		if (bsehamread) then
			call bse_hamiltonian_read_parallel(bseham_path,bseham_metadata,dimbse,hbse_dist,locr,locc,lld, &
										  mb,nb,nprow,npcol,myrow,mycol,MPI_COMM_WORLD,bseham_ok)
			if (.not. bseham_ok) then
				write(*,*) 'Unable to read a compatible BSE Hamiltonian: ',trim(bseham_path)
				call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)
			end if
		else
			hbse_dist = 0.0
			!$omp parallel do private(li,ig,lj,jg) schedule(dynamic)
			do lj=1,locc
				jg = indxl2g(lj, nb, mycol, 0, npcol)
				do li=1,locr
					ig = indxl2g(li, mb, myrow, 0, nprow)

#ifdef ELPA
					! ELPA requires the complete Hermitian matrix.  PCHEEV only
					! reads its upper triangle, so preserve the cheaper construction
					! for every non-ELPA distributed run.
					if (ig <= jg .or. trim(bsealgo) == "elpa") then
#else
					if (ig <= jg) then
#endif
						hbse_dist(li,lj)= matrizelbse(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,rlat,stt_bse(:,ig),&
						    eigv(stt_bse(4,ig),stt_bse(3,ig)),eigv(stt_bse(4,ig),stt_bse(2,ig)),vector(:,stt_bse(3,ig),stt_bse(4,ig)),&
						    vector(:,stt_bse(2,ig),stt_bse(4,ig)),kpt_bse(:,stt_bse(4,ig)),stt_bse(:,jg),eigv(stt_bse(4,jg),stt_bse(3,jg)),&
						    eigv(stt_bse(4,jg),stt_bse(2,jg)),vector(:,stt_bse(3,jg),stt_bse(4,jg)),vector(:,stt_bse(2,jg),stt_bse(4,jg)),&
						    kpt_bse(:,stt_bse(4,jg)),dft,nvec,rvec,sk(:,:,stt_bse(4,ig)),sk(:,:,stt_bse(4,jg)))
					end if
				end do
			end do
			!$omp end parallel do
			if (bsehamwrite) then
				call bse_hamiltonian_write_parallel(bseham_path,bseham_metadata,dimbse,hbse_dist,locr,locc,lld, &
										   mb,nb,nprow,npcol,myrow,mycol,MPI_COMM_WORLD,bseham_ok)
				if (.not. bseham_ok) then
					write(*,*) 'Unable to save BSE Hamiltonian: ',trim(bseham_path)
					call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)
				end if
			end if
		end if
#else
		stop "MPI/ScaLAPACK path requested without MPI support"
#endif
    end if

#ifdef MPI
	task_elapsed = MPI_WTIME() - task_start
	call MPI_REDUCE(task_elapsed,task_elapsed_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,MPIError)
#else
	call cpu_time(task_elapsed)
	task_elapsed_max = task_elapsed - task_start
#endif


    if (Node == 0) then
	 if (bsehamread) then
		write(300,*) 'BSE Hamiltonian restart read: finished'
		write(300,"(A,F12.3,A)") 'BSE Hamiltonian restart read wall time: ',task_elapsed_max,' s'
	 else
		write(300,*) 'BSE Hamiltonian construction: finished'
		write(300,"(A,F12.3,A)") 'BSE Hamiltonian construction wall time: ',task_elapsed_max,' s'
	 end if
	 call flush(300)
    endif

	!call OMP_SET_NUM_THREADS(nthreads)

	if (Node == 0) then
		write(300,*) 'BSE diagonalization: begin'
		call flush(300)
	end if
#ifdef MPI
	call MPI_BARRIER(MPI_COMM_WORLD,MPIError)
	task_start = MPI_WTIME()
#else
	call cpu_time(task_start)
#endif

    if (Nodes == 1 ) then
		select case (bsealgo)
		
		case ("cheev")

     	
      		CALL CHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      		LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      		!LWORK = MIN( 2*dimbse + dimbse**2, INT( WORK( 1 ) ) )
      		CALL CHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      		IF(INFO .GT. 0 ) THEN
        	WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         	STOP
      		END IF   

		case ("cheevr")
		
		CALL CHEEVR( 'Vectors','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
		              M, W, Z, dimbse, ISUPPZ, WORK, LWORK,RWORK, LRWORK, IWORK, LIWORK,INFO )
		LWORK = MIN( 2*dimbse, INT( WORK( 1 ) ) )
      		LRWORK = MIN( 24*dimbse, INT( RWORK( 1 ) ) )
      		LIWORK = MIN( 10*dimbse, IWORK( 1 ) )	
		CALL CHEEVR( 'Vectors','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
		              M, W, Z, dimbse, ISUPPZ, WORK, LWORK,RWORK, LRWORK, IWORK, LIWORK,INFO )
       	        IF( INFO.GT. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues.'
                STOP
         	END IF  		
		
		IF( W(dimbse) .eq. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues. Change BSE_ALGO to cheev'
                STOP
         	END IF  
         	
         	hbse = Z
         	
         	deallocate(Z)
         	
         	case ("cheevx")
		
		CALL CHEEVX( 'Vectors','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
		              M, W, Z, dimbse, WORK, LWORK,RWORK, IWORK, ifailx,INFO )
		LWORK = MIN( 2*dimbse, INT( WORK( 1 ) ) )

		CALL CHEEVX( 'Vectors','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
		              M, W, Z, dimbse, WORK, LWORK,RWORK, IWORK, ifailx,INFO )
       	        IF( INFO.GT. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues.'
                STOP
         	END IF  		
		
		IF( W(dimbse) .eq. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues. Change BSE_ALGO to cheev'
                STOP
         	END IF  
         	
         	hbse = Z
         	
         	deallocate(Z)
		
		case ("cheevd")
		
		CALL CHEEVD( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,& 	      
		              RWORK, LRWORK, IWORK, LIWORK, INFO )
      		LWORK = MIN( 2*dimbse + dimbse**2, INT( WORK( 1 ) ) )
      		LRWORK = MIN( 1 + 5*dimbse + 2*dimbse**2, INT( RWORK( 1 ) ) )
      		LIWORK = MIN( 3 + 5*dimbse, IWORK( 1 ) )
      		!LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      		!LRWORK = MIN( LWMAX, INT( RWORK( 1 ) ) )
      		!LIWORK = MIN( LWMAX, IWORK( 1 ) )
      		
		CALL CHEEVD( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,& 	      
		              RWORK, LRWORK, IWORK, LIWORK, INFO )
       	        IF(INFO .GT. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues.'
                STOP
         	END IF  
		
		IF( W(dimbse) .eq. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues. Change BSE_ALGO to cheev'
                STOP
         	END IF  		

		case ("elpa")

		write(*,*) "BSE_ALGO=elpa requires an MPI run with at least two ranks"
		stop
		
		case default
		
		 write(*,*) "please select a subroutine for diagonalization"
		 stop
		
		end select
	else
#ifdef MPI

#ifdef ELPA
		if (trim(bsealgo) == "elpa") then
			allocate(elpa_eigenvectors(lld,max(1,locc)),stat=erro)
			if (erro /= 0) then
				write(*,*) 'Could not allocate ELPA eigenvector storage'
				call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)
			end if

			elpa_status = elpa_init(ELPA_API_VERSION)
			call elpa_check(elpa_status, 'elpa_init', MPI_COMM_WORLD)
			elpa_instance => elpa_allocate(elpa_status)
			call elpa_check(elpa_status, 'elpa_allocate', MPI_COMM_WORLD)

			call elpa_instance%set('na', dimbse, elpa_status)
			call elpa_check(elpa_status, 'set na', MPI_COMM_WORLD)
			call elpa_instance%set('nev', dimbse, elpa_status)
			call elpa_check(elpa_status, 'set nev', MPI_COMM_WORLD)
			call elpa_instance%set('local_nrows', locr, elpa_status)
			call elpa_check(elpa_status, 'set local_nrows', MPI_COMM_WORLD)
			call elpa_instance%set('local_ncols', locc, elpa_status)
			call elpa_check(elpa_status, 'set local_ncols', MPI_COMM_WORLD)
			call elpa_instance%set('nblk', mb, elpa_status)
			call elpa_check(elpa_status, 'set nblk', MPI_COMM_WORLD)
			call elpa_instance%set('mpi_comm_parent', MPI_COMM_WORLD, elpa_status)
			call elpa_check(elpa_status, 'set mpi_comm_parent', MPI_COMM_WORLD)
			call elpa_instance%set('process_row', myrow, elpa_status)
			call elpa_check(elpa_status, 'set process_row', MPI_COMM_WORLD)
			call elpa_instance%set('process_col', mycol, elpa_status)
			call elpa_check(elpa_status, 'set process_col', MPI_COMM_WORLD)
			! Do not force ELPA's optional omp_threads setting here.  ELPA builds
			! without OpenMP expose no such configurable option (status -3), while
			! OpenMP-enabled ELPA obtains the caller's OpenMP thread count during
			! setup when this setting is left unset.
			elpa_status = elpa_instance%setup()
			call elpa_check(elpa_status, 'elpa setup', MPI_COMM_WORLD)
			! The one-stage complex solver is the robust ELPA path on the local
			! ARM/NEON build; the two-stage complex kernel returned invalid BSE
			! eigenvalues for a valid 3600-state Hermitian matrix.
			call elpa_instance%set('solver', ELPA_SOLVER_1STAGE, elpa_status)
			call elpa_check(elpa_status, 'set ELPA 1-stage solver', MPI_COMM_WORLD)
			call elpa_instance%eigenvectors(hbse_dist, W, elpa_eigenvectors, elpa_status)
			call elpa_check(elpa_status, 'ELPA eigenvectors', MPI_COMM_WORLD)
			call elpa_deallocate(elpa_instance, elpa_status)
			call elpa_check(elpa_status, 'elpa_deallocate', MPI_COMM_WORLD)
			call elpa_uninit(elpa_status)
			call elpa_check(elpa_status, 'elpa_uninit', MPI_COMM_WORLD)

			deallocate(hbse_dist)
			call move_alloc(elpa_eigenvectors, hbse_dist)
		else
#else
		if (trim(bsealgo) == "elpa") then
			write(*,*) 'BSE_ALGO=elpa was requested, but this executable was built without ELPA support'
			call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)
		end if
#endif
		LWORK = -1
		LRWORK = -1
		! PCHEEV overwrites Z before the back transformation is complete, so A and Z
		! must be distinct distributed arrays (unlike the one-process CHEEV path).
		allocate(hbse_eigenvectors(lld,max(1,locc)),stat=erro)
		if (erro /= 0) then
			write(*,*) 'Could not allocate distributed BSE eigenvector storage'
			call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)
		end if
		allocate(WORK(1),RWORK(1))
		call PCHEEV('V','U',dimbse,hbse_dist,1,1,desca,W,hbse_eigenvectors,1,1,descz,WORK,LWORK,RWORK,LRWORK,INFO)
		LWORK = max(1,int(real(WORK(1))))
		LRWORK = max(1,int(RWORK(1)))
		deallocate(WORK,RWORK)
		allocate(WORK(LWORK),RWORK(LRWORK))
		call PCHEEV('V','U',dimbse,hbse_dist,1,1,desca,W,hbse_eigenvectors,1,1,descz,WORK,LWORK,RWORK,LRWORK,INFO)
		if (INFO .ne. 0) then
			write(*,*) 'ScaLAPACK PCHEEV failed with INFO = ', INFO
			call MPI_ABORT(MPI_COMM_WORLD, INFO, MPIError)
		end if
		deallocate(hbse_dist)
		call move_alloc(hbse_eigenvectors,hbse_dist)
#ifdef ELPA
		end if
#endif
#endif
	endif

#ifdef MPI
	task_elapsed = MPI_WTIME() - task_start
	call MPI_REDUCE(task_elapsed,task_elapsed_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,MPIError)
#else
	call cpu_time(task_elapsed)
	task_elapsed_max = task_elapsed - task_start
#endif
	if (Node == 0) then
		write(300,*) 'BSE diagonalization: finished'
		write(300,"(A,F12.3,A)") 'BSE diagonalization wall time: ',task_elapsed_max,' s'
		call flush(300)
	end if

    !allocate(pinter(dimbse),pintra(dimbse))

	!if (ntype .eq. 2) then

	!do i=1,dimbse

	!	call excitonil(i,w90basis,nc,nv,ngrid*ngrid,stt,dimbse,hbse(:,i),lcount,pinter(i),pintra(i))

	!	write(306,*) W(i),pinter(i),pintra(i)
	!	call flush(306)

	!end do 

	!else 

	!	continue

	!end if 
	
		if (Nodes == 1) then
		select case (bsealgo)
	
		case ("cheev")
		
		deallocate(RWORK,WORK)
		
		case ("cheevr")
		
		deallocate(WORK,ISUPPZ,IWORK,RWORK)
		
		case ("cheevx")
		
		deallocate(WORK,IWORK,RWORK,ifailx)
		
		case ("cheevd")		
		
		deallocate(WORK,IWORK,RWORK)

		case default 
		
		 write(*,*) "please select a subroutine for diagonalization"
		 go to 789	
		 	
		end select
		else
			if (allocated(WORK)) deallocate(WORK)
			if (allocated(RWORK)) deallocate(RWORK)
		end if

		if (Nodes == 1) then
			if (Node == 0) then
				if (bsewf) then
					do i=excwf0,excwff
						call excwf(outputfolder,ngkpt,kpt,nc,nv,nocpk,stt,W(i),i,hbse(:,i))
					end do
				end if

				call dielbsev(nthreads,dimbse,hbse,hrx,actxx)
				call dielbsev(nthreads,dimbse,hbse,hry,actyy)
				call dielbsev(nthreads,dimbse,hbse,hrz,actzz)

				if (dtfull) then
					call dielbsep(nthreads,dimbse,hbse,hrx,hry,actxy)
					call dielbsep(nthreads,dimbse,hbse,hrx,hrz,actxz)
					call dielbsep(nthreads,dimbse,hbse,hry,hrz,actyz)
				else
					actxy = 0.0
					actxz = 0.0
					actyz = 0.0
				end if

				if (cpol) then
					call dielbsev(nthreads,dimbse,hbse,hrsp,actsp)
					call dielbsev(nthreads,dimbse,hbse,hrsm,actsm)
				else
					actsp = 0.0
					actsm = 0.0
				end if
			end if
		else
#ifdef MPI
			if (bsewf) then
				allocate(eigvec_dist(dimbse))
				do i=excwf0,excwff
					call bse_eigenvector_column_dist(dimbse,hbse_dist,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,i,eigvec_dist,MPIError)
					if (Node == 0) call excwf(outputfolder,ngkpt,kpt,nc,nv,nocpk,stt,W(i),i,eigvec_dist)
				end do
				deallocate(eigvec_dist)
			end if

			call dielbsev_dist(dimbse,hbse_dist,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hrx,actxx,MPIError)
			call dielbsev_dist(dimbse,hbse_dist,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hry,actyy,MPIError)
			call dielbsev_dist(dimbse,hbse_dist,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hrz,actzz,MPIError)

			if (dtfull) then
				call dielbsep_dist(dimbse,hbse_dist,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hrx,hry,actxy,MPIError)
				call dielbsep_dist(dimbse,hbse_dist,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hrx,hrz,actxz,MPIError)
				call dielbsep_dist(dimbse,hbse_dist,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hry,hrz,actyz,MPIError)
			else
				actxy = 0.0
				actxz = 0.0
				actyz = 0.0
			end if

			if (cpol) then
				call dielbsev_dist(dimbse,hbse_dist,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hrsp,actsp,MPIError)
				call dielbsev_dist(dimbse,hbse_dist,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hrsm,actsm,MPIError)
			else
				actsp = 0.0
				actsm = 0.0
			end if
#endif
		end if

		if (Node == 0) then
			write(300,*) 'exciton Hamiltonian diagonalized'
			call flush(300)
			write(300,*) "exciton ground state",W(1)
			write(300,*) 'xx tensor component'
			write(300,*) 'yy tensor component'
			write(300,*) 'zz tensor component'
			if (dtfull) then
				write(300,*) 'xy tensor component'
				write(300,*) 'xz tensor component'
				write(300,*) 'yz tensor component'
			else
				write(300,*) 'xy tensor component set to 0'
				write(300,*) 'xz tensor component set to 0'
				write(300,*) 'yz tensor component set to 0'
			end if
			if (cpol) then
				write(300,*) 'sp polarization'
				write(300,*) 'sm polarization'
			else
				write(300,*) 'sp polarization set to 0'
				write(300,*) 'sm polarization set to 0'
			end if
			write(300,*) 'optics finished'
			call flush(300)

			do i=1,dimbse
				write(301,"(7F15.6)") W(i),actxx(i),actyy(i),actzz(i),actxy(i),actxz(i),actyz(i)
				call flush(301)
				if (cpol) then
					write(302,"(6F15.6)") W(i),actxx(i),actyy(i),actzz(i),actsp(i),actsm(i)
					call flush(302)
				end if
			end do
		end if
			
		deallocate(eigv,vector)
		deallocate(rvec,hopmatrices)
		deallocate(ihopmatrices,ffactor)
		deallocate(hrx,hry,hrz,hrsp,hrsm)
	deallocate(actxx,actxy,actxz,actyy,actyz,actzz)
	deallocate(actsp,actsm)
		if (Nodes == 1) then
			deallocate(hbse)
		else
#ifdef MPI
			call BLACS_GRIDEXIT(blacs_ctxt)
			deallocate(hbse_dist)
#endif
		end if
		deallocate(W,stt,stt_bse,nocpk)
	deallocate(ovp)

	deallocate(kpt,kpt_bse)
	
		deallocate(sk)
		call rmn_destroy(rmn)


789     continue

	call cpu_time(tf)
	call date_and_time(VALUES=values2)

		if (Node == 0) then
			write(300,*)
			write(300,*) 'end','   ','month',values2(2),'day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
			write(300,*)
		end if




		if (Node == 0) then
		close(200)
	!close(203)



	close(300)
	close(301)
	close(302)
	
	close(401)
	close(402)
	close(403)
		close(404)
		end if
	
	!close(500)			


	


end subroutine bsesolver

#ifdef MPI
#ifdef ELPA
subroutine elpa_check(status, operation, comm)
	use mpi
	use elpa

	implicit none

	integer, intent(in) :: status, comm
	character(len=*), intent(in) :: operation
	integer :: mpi_ierr

	if (status /= ELPA_OK) then
		write(*,*) 'ELPA failed during ', trim(operation), '; status = ', status
		call MPI_ABORT(comm, 1, mpi_ierr)
		stop
	end if
end subroutine elpa_check
#endif
#endif
