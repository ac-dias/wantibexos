!gfortran -mcmodel=large bse_kpath-tool.f90 -o bse_kpath-tool.x -llapack95 -lopenblas -fopenmp 

subroutine bsebndstemp(nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     exc,mshift,coultype,ez,w1,r0,lc,rk,meshtype,bsewf,excwf0,excwff,&
		     st,phavg,ta,temp,nocpf,fermishift,bsealgo,bsekpathmpi,bsekpathcheckpoint,bsekpathcheckpointfile,dft,mag)

#ifdef MPI
	use mpi
#endif
	use omp_lib
	use hamiltonian_input_variables

	implicit none

	real,parameter:: pi=acos(-1.)
	integer :: dimbse
	real,allocatable,dimension(:,:) :: exk !(dimbse,nqpts*npath)


	integer :: ncaux,nvaux

	integer,allocatable,dimension(:,:) :: stt,stt_bse !(ngrid*ngrid*nc*nv,4)
	real,allocatable,dimension(:,:) :: kpt,qpt !pontos k do grid (ngrid*ngrid,2)
	real,allocatable,dimension(:,:) :: kpt_bse

	real,dimension(4) :: q

	real,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores

	real,allocatable,dimension(:,:) :: energy !variavel que guarda os autovalores para cada ponto k
	complex,allocatable,dimension(:,:,:) :: vector !variavel que guarda os autovetores para cada ponto k

	real,allocatable,dimension(:,:) :: energyq !variavel que guarda os autovalores para cada ponto k
	complex,allocatable,dimension(:,:,:) :: vectorq !variavel que guarda os autovetores para cada ponto k

	real,allocatable,dimension(:,:) :: qauxv !(nqpts*npath,4) 

	integer:: counter,c,v,i,j,f,l,h,erro,i2,k

	complex:: matrizelbsekqtemp

	real :: a,r0,ed

	integer :: ngkpt,nqpath,checkpoints_loaded,checkpoints_restored
	integer,dimension(3) :: nqgrid
	real,dimension(3) :: qshift

	!variaveis relacionadas a marcacao do tempo

	real:: t0,tf
	integer,dimension(8) :: values,values2

	!variaveis kpoints

	integer :: nks
	integer :: nkpts
	real,allocatable,dimension(:,:) :: ks

	integer,allocatable,dimension(:) :: nocpk,nocpq
	
	integer :: nocpf
	real :: fermishift
	character(len=12) :: bsealgo,bsekpathmpi
	character(len=70) :: bsekpathcheckpointfile
	character(len=1) :: dft
	integer :: Node,Nodes,MPIError
	logical :: bsekpathcheckpoint,checkpoint_ok,checkpoint_loaded,qgrid
	character(len=240) :: checkpoint_path
	character(len=8) :: checkpoint_label
	character(len=256) :: qinput
	character(len=16) :: qsampling
	
	real,dimension(3) :: mag
	complex,allocatable,dimension(:,:,:) :: sk,skq		


	!definicoes diagonalizacao 
	INTEGER   ::       ifail
	integer,allocatable,dimension(:) :: ifailx
	real,parameter :: ABSTOL=1.0e-6
	INTEGER          INFO
	real,allocatable,dimension(:) :: W,RWORK
	COMPLEX,allocatable,dimension(:,:) :: hbse

        INTEGER ::          LWMAX
   	INTEGER ::         LWORK
	INTEGER ::         LIWORK, LRWORK
	INTEGER,allocatable,dimension(:) :: IWORK
        complex,allocatable,dimension (:) :: WORK
        
        !zheevr definitions
        real :: VL,VU
        integer :: IL,IU,M
        complex,allocatable,dimension(:,:) :: Z
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
	character(len=2) :: ta
	real,dimension(3) :: ediel
	real :: ez,w1,lc
	logical :: bsewf
	integer :: excwf0,excwff
	
	real :: st,phavg,temp
	real :: gapcortemp,gapcortemp2,tcor
		
	

	!call input_read

	Node=0
	Nodes=1
	MPIError=0
#ifdef MPI
	call MPI_COMM_RANK(MPI_COMM_WORLD,Node,MPIError)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,Nodes,MPIError)
#endif

	if (trim(bsekpathmpi) .ne. "Q") then
		if (Node .eq. 0) write(*,*) "Unsupported BSE_KPATH_MPI; use Q (HYBRID is reserved)"
#ifdef MPI
		call MPI_ABORT(MPI_COMM_WORLD,1,MPIError)
#endif
		stop
	end if

	! INPUT 
	OPEN(UNIT=500, FILE= kpathsbse,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bse-kpath input file"
	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"

	! Only rank zero owns the shared formatted output files.
	if (Node .eq. 0) then
		OPEN(UNIT=300, FILE=trim(outputfolder)//"log_bse_kpath.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_bse_kpath output file"
		OPEN(UNIT=400, FILE=trim(outputfolder)//"bands_bse.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bands_bse output file"
	end if


	call cpu_time(t0)
	call date_and_time(VALUES=values)
	call OMP_SET_NUM_THREADS(nthreads)

	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	allocate(ovp(nvec,w90basis,w90basis))

	call hamiltonian_input_read(200,params)
	

	end select
	!ediel(2) = edielh

	! KPATH_BSE accepts its original path format, or the explicit GRID format:
	!   GRID
	!   nq1 nq2 nq3
	!   shift1 shift2 shift3     ! optional, defaults to 0 0 0
	! A grid shift is expressed in units of a grid spacing along each reciprocal
	! lattice vector, so shift=0 includes Gamma and shift=0.5 is half-shifted.
	qgrid=.false.
	qshift=0.0
	read(500,'(A)',IOSTAT=erro) qinput
	if (erro .ne. 0) stop "Error reading BSE q-point sampling"
	read(qinput,*,IOSTAT=erro) qsampling
	if (erro .ne. 0) stop "Error reading BSE q-point sampling mode"

	select case (trim(qsampling))
	case ("GRID","grid")
		qgrid=.true.
		read(500,*,IOSTAT=erro) nqgrid
		if (erro .ne. 0 .or. any(nqgrid .le. 0)) stop "Invalid BSE q-grid dimensions"
		read(500,*,IOSTAT=erro) qshift
		if (erro .gt. 0) stop "Invalid BSE q-grid shift"
		if (erro .lt. 0) qshift=0.0
	case ("PATH","path")
		read(500,*,IOSTAT=erro) nks
		if (erro .ne. 0) stop "Error reading BSE q-path endpoints"
	case default
		! Backward-compatible input: the first record is the endpoint count.
		read(qinput,*,IOSTAT=erro) nks
		if (erro .ne. 0) stop "Unknown BSE q-point sampling mode"
	end select

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	if (qgrid) then
		nqpath=nqgrid(1)*nqgrid(2)*nqgrid(3)
	else
		read(500,*,IOSTAT=erro) nkpts
		if (erro .ne. 0 .or. nks .le. 0 .or. mod(nks,2) .ne. 0 .or. nkpts .lt. 2) then
			stop "Invalid BSE q-path; use an even endpoint count and at least two points"
		end if
		allocate(ks(nks,3))
		do i=1,nks
			read(500,*,IOSTAT=erro) ks(i,1),ks(i,2),ks(i,3)
			if (erro .ne. 0) stop "Error reading BSE q-path endpoint"
		end do
		nqpath=(nks/2)*nkpts
	end if
	dimbse = ngkpt*nc*nv

	!call alat(systype,rlat,a)
	!lc = rlat3(3)

	!if (systype .eq. "2D") then
	!	ed = (ediel(1)+ediel(3))/2.
	!	r0= ((ediel(2)-1.0)*lc)/(ediel(1)+ediel(3))
	!else
	!	r0 = 1.0
	!	ed = ediel(2)
		
	!end if

	allocate(qauxv(nqpath,4))

	if (Node .eq. 0) then
		if (qgrid) then
			call qgridbse(rlat(1,:),rlat(2,:),rlat(3,:),nqgrid,qshift,qauxv)
		else
			call kpathbse(outputfolder,rlat(1,:),rlat(2,:),rlat(3,:),nks,ks,nkpts,qauxv)
		end if
	end if
#ifdef MPI
	call MPI_BCAST(qauxv,4*nqpath,MPI_REAL,0,MPI_COMM_WORLD,MPIError)
#endif

	!Informações para o arquivo de log do calculo

	if (Node .eq. 0) then
	write(300,*) 'threads:', nthreads
	write(300,*)
	write(300,*) 'temperature', temp
	write(300,*)
	write(300,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
	write(300,*)
	write(300,*) 'ktol-coulomb:', ktol
	write(300,*)
	write(300,*) 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	write(300,*)
	write(300,*) 'conduction bands, nc:',nc,'  ','valence bands, nv:',nv
	write(300,*)
	write(300,*) 'Coulomb Potential:',coultype
	write(300,*)
	write(300,*) 'BSE_ALGO:',bsealgo	
	write(300,*) 'BSE_KPATH_MPI:',bsekpathmpi
	if (qgrid) then
		write(300,*) 'BSE q-point sampling: GRID',nqgrid(1),nqgrid(2),nqgrid(3)
		write(300,*) 'BSE q-grid shift:',qshift(1),qshift(2),qshift(3)
	else
		write(300,*) 'BSE q-point sampling: PATH'
	end if
	write(300,*)
	write(300,*) 'number of q points:','   ',nqpath
	write(300,*)
	write(300,*)
	write(300,*) 'begin','  ','month',values(2),'day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(300,*) 

	call flush(300)
	end if

	allocate(kpt(ngkpt,3))
	
	!shift= 0.0
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)
	!call gridgenmhp(ngrid,rlat,kpt)
	allocate(kpt_bse(3,ngkpt))
	do i=1,ngkpt
		do j=1,3
		kpt_bse(j,i) = kpt(i,j)
		end do
	end do

	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	! Store the orbital index first; BSE eigenvector sections are then contiguous.
	allocate(energy(ngkpt,nc+nv),vector(w90basis,nc+nv,ngkpt))

	allocate(nocpk(ngkpt))
	allocate(nocpq(ngkpt))

	select case (ta)
	
	case("FA")
	
		tcor = 0.0
	
	case("VE")
	
		tcor = gapcortemp(st,phavg,temp)
	
	case("BE")
	
		tcor = gapcortemp2(st,phavg,temp)
	
	case default
	
		tcor= 0.00
	
	end select

	!$omp parallel do default(shared) private(i,eaux,vaux,j,l,h)
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

		call eigsys(nthreads,dft,systype,scs+tcor,exc,nocpk(i),&
			    ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
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
	
				energy(i,j)= eaux(nocpk(i)-nv+j)

			end do


			do l=1,nc+nv


				do h=1,w90basis

					vector(h,l,i)=vaux(nocpk(i)-nv+l,h)


				end do
			

			end do

	
	end do
	!$omp end parallel do

	deallocate(eaux,vaux)

	if (dft .eq. "S") then
	
		! Keep each overlap matrix contiguous when passed into the BSE kernel.
		allocate(sk(w90basis,w90basis,ngkpt))
		
		do i=1,ngkpt
		
			call overlap(w90basis,nvec,rvec,ovp,kpt(i,1),kpt(i,2),kpt(i,3),sk(:,:,i))
		
		end do
		!write(300,*) 'overlap matrices calculated'
		!call flush(300)	
	
	end if	

	!definindo os numeros quanticos dos estados




	if (Node .eq. 0) then
		write(300,*) "quantum numbers for exciton basis set finished"
		call flush(300)
	end if

	counter=counter-1 !numero total de estados para equação bse

	allocate(qpt(ngkpt,3))
	allocate(exk(dimbse,nqpath))
	exk=0.0
	checkpoints_loaded=0
	if (Node .eq. 0) then
		call bse_hamiltonian_memory_report(300,'BSE Q-point Hamiltonian per active MPI rank',dimbse,dimbse)
	end if



	q_loop: do i=Node+1,nqpath,Nodes

		
		!definindo os pontos q

		q(1)= qauxv(i,1)
		q(2)= qauxv(i,2) 
		q(3)= qauxv(i,3)
		q(4)= qauxv(i,4)

		if (bsekpathcheckpoint) then
			write(checkpoint_label,"(I8.8)") i
			checkpoint_path=trim(outputfolder)//trim(bsekpathcheckpointfile)//"_q"//checkpoint_label//".bin"
			call bse_kpath_checkpoint_read(checkpoint_path,i,dimbse,q,exk(:,i),checkpoint_loaded)
			if (checkpoint_loaded) then
				checkpoints_loaded=checkpoints_loaded+1
				cycle q_loop
			end if
		end if

		!gerando o grid k+q
		
		call monhkhorst_packq(q,ngrid(1),ngrid(2),ngrid(3),mshift,&
				      rlat(1,:),rlat(2,:),rlat(3,:),qpt)


		allocate(eaux(w90basis),vaux(w90basis,w90basis))

		allocate(energyq(ngkpt,nc+nv),vectorq(w90basis,nc+nv,ngkpt))

	allocate (stt(ngkpt*nc*nv,4))


	!$omp parallel do default(shared) private(i2,j,l,h,eaux,vaux)
		do i2=1,ngkpt


#ifdef MKL
		call MKL_SET_NUM_THREADS(1)
#endif		

#ifdef AOCL		
		call bli_thread_set_num_threads(1)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(1)
		
#endif	




			call eigsys(nthreads,dft,systype,scs+tcor,exc,nocpq(i2),&
			            ffactor,qpt(i2,1),qpt(i2,2),qpt(i2,3),w90basis,nvec,&
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
	
					energyq(i2,j)= eaux(nocpq(i2)-nv+j) 

				end do

	
				do l=1,nc+nv


					do h=1,w90basis

						vectorq(h,l,i2)=vaux(nocpq(i2)-nv+l,h)

	
					end do
			

				end do



	
		end do
	!$omp end parallel do
	
	if (dft .eq. "S") then
	
		allocate(skq(w90basis,w90basis,ngkpt))
		
		 
		do i2=1,ngkpt
		
			call overlap(w90basis,nvec,rvec,ovp,qpt(i2,1),qpt(i2,2),qpt(i2,3),skq(:,:,i2))
		
		end do
		if (Node .eq. 0) then
			write(300,*) 'overlap matrices calculated'
			call flush(300)
		end if
	
	
	end if	

		call quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpq,stt)
		allocate(stt_bse(4,dimbse))
		do i2=1,dimbse
			do j=1,4
				stt_bse(j,i2) = stt(i2,j)
			end do
		end do

		deallocate(eaux,vaux)


		allocate(hbse(dimbse,dimbse),W(dimbse))

		hbse=0.

	!$omp parallel do 
        !collapse(2)

		do i2=1,dimbse



			do j=i2,dimbse




hbse(i2,j)= matrizelbsekqtemp(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,q,rlat,stt_bse(:,i2),&
          energyq(stt_bse(4,i2),stt_bse(3,i2))&
          ,energy(stt_bse(4,i2),stt_bse(2,i2)),vectorq(:,stt_bse(3,i2),stt_bse(4,i2)),vector(:,stt_bse(2,i2),stt_bse(4,i2)),&
          kpt_bse(:,stt_bse(4,i2)),stt_bse(:,j),energyq(stt_bse(4,j),stt_bse(3,j)),energy(stt_bse(4,j),stt_bse(2,j))&
          ,vectorq(:,stt_bse(3,j),stt_bse(4,j)),vector(:,stt_bse(2,j),stt_bse(4,j)),kpt_bse(:,stt_bse(4,j)),temp,dft,nvec,rvec,&
          sk(:,:,stt_bse(4,i2)),sk(:,:,stt_bse(4,j)),skq(:,:,stt_bse(4,i2)),skq(:,:,stt_bse(4,j)))

			end do



		end do

	!$omp end parallel do

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
		
		case default
		
		 write(*,*) "please select a subroutine for diagonalization"
		 stop
		
	end select
	
	if (bsewf) then
	
	select case (bsealgo)
		
		case ("cheev")

     	
      		CALL CHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      		LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      		CALL CHEEV( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      		IF( INFO.GT. 0 ) THEN
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
		              RWORK, LRWORK, IWORK, LIWORK,INFO )
      		LWORK = MIN( 2*dimbse + dimbse**2, INT( WORK( 1 ) ) )
      		LRWORK = MIN( 1 + 5*dimbse + 2*dimbse**2, INT( RWORK( 1 ) ) )
      		LIWORK = MIN( 3 + 5*dimbse, IWORK( 1 ) )
		CALL CHEEVD( 'Vectors', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,& 	      
		              RWORK, LRWORK, IWORK, LIWORK,INFO )
       	        IF( INFO.GT. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues.'
                STOP
         	END IF 
         	
         	IF( W(dimbse) .eq. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues. Change BSE_ALGO to cheev'
                STOP
         	END IF           	 
		
		case default
		
		 write(*,*) "please select a subroutine for diagonalization"
		 stop
		
	end select
	
	
	else

 	select case (bsealgo)
		
		case ("cheev")

     	
      		CALL CHEEV( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK, INFO )
      		LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      		CALL CHEEV( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK, RWORK,INFO )
      		IF( INFO.GT. 0 ) THEN
        	WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         	STOP
      		END IF   

		case ("cheevr")
		
		CALL CHEEVR( 'N','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
		              M, W, Z, dimbse, ISUPPZ, WORK, LWORK,RWORK, LRWORK, IWORK, LIWORK,INFO )
		LWORK = MIN( 2*dimbse, INT( WORK( 1 ) ) )
      		LRWORK = MIN( 24*dimbse, INT( RWORK( 1 ) ) )
      		LIWORK = MIN( 10*dimbse, IWORK( 1 ) )	
		CALL CHEEVR( 'N','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
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
		
		CALL CHEEVX( 'N','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
		              M, W, Z, dimbse, WORK, LWORK,RWORK, IWORK, ifailx,INFO )
		LWORK = MIN( 2*dimbse, INT( WORK( 1 ) ) )

		CALL CHEEVX( 'N','A', 'U', dimbse, hbse, dimbse,VL,VU,IL,IU,ABSTOL,&
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
		
		CALL CHEEVD( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,& 	      
		              RWORK, LRWORK, IWORK, LIWORK,INFO )
      		LWORK = MIN( 2*dimbse + dimbse**2, INT( WORK( 1 ) ) )
      		LRWORK = MIN( 1 + 5*dimbse + 2*dimbse**2, INT( RWORK( 1 ) ) )
      		LIWORK = MIN( 3 + 5*dimbse, IWORK( 1 ) )
		CALL CHEEVD( 'N', 'U', dimbse, hbse, dimbse, W, WORK, LWORK,& 	      
		              RWORK, LRWORK, IWORK, LIWORK,INFO )
       	        IF( INFO.GT. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues.'
                STOP
         	END IF 
         	
         	IF( W(dimbse) .eq. 0 ) THEN
                WRITE(*,*)'The algorithm failed to compute eigenvalues. Change BSE_ALGO to cheev'
                STOP
         	END IF           	 
		
		case default
		
		 write(*,*) "please select a subroutine for diagonalization"
		 stop
		
	end select
      	
      	end if  

		!write(*,*) W(1)


		do i2=1,dimbse

			exk(i2,i)=W(i2)

		!	write(400,*) qauxv(i,1),W(i2)
		!	call flush(400)

		end do


	if (bsewf) then
	
	      	do i2=excwf0,excwff
      	
      			call excwfi(outputfolder,ngkpt,kpt,q,nc,nv,nocpk,stt,W(i2),i2,i,hbse(:,i2))
      	
      		end do
	
	else
	
	 continue
	
	end if

		if (bsekpathcheckpoint) then
			call bse_kpath_checkpoint_write(checkpoint_path,i,dimbse,q,exk(:,i),checkpoint_ok)
			if (.not. checkpoint_ok) then
				write(*,*) 'Unable to write BSE q-path checkpoint: ',trim(checkpoint_path)
#ifdef MPI
				call MPI_ABORT(MPI_COMM_WORLD,1,MPIError)
#endif
				stop
			end if
		end if

		deallocate(hbse,W)
		deallocate(energyq,vectorq)
		deallocate(stt,stt_bse)

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
		 go to 580	
		 	
	end select
		!deallocate (IWORK)


		qpt=0.0
		nocpq=0

		if (Node .eq. 0) then
			write(300,*) 'progress:',i,'/',nqpath
			call flush(300)
		end if
		
		if (dft .eq. "S") then
		
			deallocate(skq)
		
		end if		

	end do q_loop

	checkpoints_restored=checkpoints_loaded
#ifdef MPI
	if (Nodes .gt. 1) then
		call MPI_REDUCE(checkpoints_loaded,checkpoints_restored,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,MPIError)
	end if
#endif
	if (Node .eq. 0 .and. bsekpathcheckpoint) then
		write(300,*) 'BSE q-path checkpoints restored:',checkpoints_restored
		call flush(300)
	end if

#ifdef MPI
	if (Nodes .gt. 1) then
		if (Node .eq. 0) then
			call MPI_REDUCE(MPI_IN_PLACE,exk,dimbse*nqpath,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,MPIError)
		else
			call MPI_REDUCE(exk,exk,dimbse*nqpath,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,MPIError)
		end if
	end if
#endif

	if (Node .eq. 0) then
	do i2=1,dimbse

		do i=1,nqpath

			if (qgrid) then
				write(400,*) qauxv(i,2),qauxv(i,3),qauxv(i,4),exk(i2,i)
			else
				write(400,*) real(qauxv(i,1)),exk(i2,i)
			end if
			call flush(400)

		end do

		write(400,*)

	end do
	end if

580 continue

	deallocate(energy,vector)
	deallocate(kpt,kpt_bse)
	deallocate(qpt)
	deallocate(exk)
	deallocate(nocpq,nocpk)
	deallocate(rvec,hopmatrices)
	deallocate(ihopmatrices,ovp,ffactor)
	
		if (dft .eq. "S") then
		
			deallocate(sk)
		
		end if		

	call cpu_time(tf)
	call date_and_time(VALUES=values2)

	if (Node .eq. 0) then
		write(300,*)
		write(300,*) 'end','   ','month',values2(2),'day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
		write(300,*)
	end if

	close(200)



	if (Node .eq. 0) then
		close(300)
		close(400)
	end if
	close(500)



end subroutine bsebndstemp
