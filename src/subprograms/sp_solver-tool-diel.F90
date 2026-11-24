!gfortran -mcmodel=large rpa_optics-tool-diel.f90 -o rpa_opt-tool-diel.x -llapack95 -lopenblas -fopenmp

subroutine spoptics(nthreads,dft,outputfolder,ngrid,nc,nv, &
		     sme,params,exc,mshift,tmcoef,nocpf,fermishift,mag)

#ifdef MPI
	use mpi
#endif
	use omp_lib
	use hamiltonian_input_variables
	use input_variables, only: bsecenterfile
	use bse_q_optics, only: rmn_data, rmn_center_setup, rmn_destroy, rmn_dfts_dipoles, &
		rmn_apply_q0_optical_correction

	implicit none

	real,parameter:: pi=acos(-1.)

	integer :: dimsp !=ngrid*ngrid*nc*nv ! dimensão da matriz bse


	real,allocatable,dimension(:,:) :: eigt !band energies, eigt(band,k)
	real,allocatable,dimension(:,:) :: tau !DFT=S orbital centres
	real,allocatable,dimension(:,:,:,:) :: dmat !DFT=S orbital dipoles (BSE_CENTER_FILE)
	type(rmn_data) :: rmn
	logical :: use_dmat,use_rmn,rmn_ok
	character(len=256) :: rmn_message
	complex,allocatable,dimension(:,:,:) :: vector !eigenvectors, vector(orbital,band,k)
	complex,allocatable,dimension(:) :: hrx,hry,hrz !optical vertices of the transitions
	complex,allocatable,dimension(:,:,:) :: opt_hx,opt_hy,opt_hz !DFT=S vertices <c|dH/dk-ebar dS/dk|v>/(Ec-Ev), (nc,nv,k)

	real,allocatable,dimension(:,:) :: kpt !pontos k do grid

	integer,allocatable,dimension(:,:) :: stt

	real,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores

	!real,allocatable,dimension(:) :: esp

	integer:: counter,c,v,i,j,f,l,h,k,kl,erro,iko

	! The k points are shared between the MPI ranks (all of them in a serial run): block kfirst..klast of this rank, its
	! transitions ifirst..ilast, and the wall times of the stage on this rank and on the slowest one (1, 2 eigensystems:
	! computing, sharing; 3, 4 vertices: computing, sharing)
	integer :: Node,Nodes,MPIError
	integer :: kfirst,klast,nt_k,ifirst,ilast
	integer,parameter :: ipa_mpi_block_elements=16777216
	double precision :: ipa_t(4),ipa_tmax(4),ipa_tw

	integer :: ncaux,nvaux

	complex:: matrizelbse
	
	complex :: hxsp,hysp,hzsp,hrsp,hrsm

	!complex :: exx,exy,exz,eyy,eyz,ezz


	real :: ec,ev

	complex,parameter :: imag=cmplx(0.0,1.0)


	integer :: ngkpt
	!real,allocatable,dimension(:) :: auxx,auxy,auxz,auyy,auyz,auzz
	real,allocatable,dimension(:,:) :: vecres

	integer,allocatable,dimension(:) :: nocpk
	
	integer :: nocpf
	real :: fermishift
	
	real:: t0,tf
	integer,dimension(8) :: values,values2	



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
	logical :: tmcoef

	!fim modificacoes versao 2.1
	!call input_read

	! INPUT : lendo os parametros do modelo de tight-binding
	!OPEN(UNIT=201, FILE= diein,STATUS='old', IOSTAT=erro)
    	!if (erro/=0) stop "Erro na abertura do arquivo de entrada ambiente dieletrico"


	Node = 0
	Nodes = 1
#ifdef MPI
	call MPI_COMM_RANK(MPI_COMM_WORLD, Node, MPIError)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, Nodes, MPIError)
#endif

	!OUTPUT, written by rank 0 alone
	if (Node == 0) then
	OPEN(UNIT=301, FILE=trim(outputfolder)//"ipa_oscf.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening ipa_oscf output file"
    	
    	
    	if (tmcoef) then
    	OPEN(UNIT=302, FILE=trim(outputfolder)//"tm_coef.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef output file"
    	OPEN(UNIT=304, FILE=trim(outputfolder)//"tm_coef-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening tm_coef-pol output file"    	
    	else
    	 continue
    	end if    
    	
    	OPEN(UNIT=303, FILE=trim(outputfolder)//"ipa_oscf-pol.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening ipa_oscf-pol output file"	
    	
	OPEN(UNIT=400, FILE=trim(outputfolder)//"log_ipa-diel.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_ipa-diel output file"  	


	end if


	call cpu_time(t0)
	call date_and_time(VALUES=values)
	call OMP_SET_NUM_THREADS(nthreads)



	!inicio leitura parametros

	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	call hamiltonian_input_read(200,params)

	allocate(ovp(nvec,w90basis,w90basis))
	ovp = 0.0

	end select
	!ediel(2) = edielh


	!termino leitura parametros

	!parametros do calculo


	!termino parametros calculo

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimsp = ngkpt*nc*nv
	
	if (Node == 0) then
	write(400,*)
	write(400,*)
	write(400,*)
	write(400,*) 'threads:', nthreads
	write(400,*)
	write(400,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
	write(400,*)
	write(400,"(A13,3F15.4)") 'kmesh shift:',mshift(1),mshift(2),mshift(3)


	write(400,*)
	write(400,*) 'conduction bands, nc:',nc,'  ','valence bands, nv:',nv
	write(400,*)
	write(400,*) 'begin','  ','month',values(2),'day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(400,*) 

	call flush(400)	
	end if


	allocate(kpt(ngkpt,3))
	allocate(nocpk(ngkpt))

	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)

	! Every rank takes a block of consecutive k points (ipa_kblock) and diagonalizes H(k) c = E S(k) c for them with OpenMP,
	! then receives the other blocks; a serial run is one rank with every k point.
	call ipa_kblock(Nodes,Node,ngkpt,kfirst,klast)
	nt_k = nc*nv
	ifirst = (kfirst-1)*nt_k+1
	ilast = klast*nt_k
	ipa_t = 0.0d0
	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(eigt(nc+nv,ngkpt),vector(w90basis,nc+nv,ngkpt))
	eigt = 0.0
	if (Node == 0) then
		write(400,"(A,I0,A,I0,A,I0,A)") ' IPA: ',ngkpt,' k points, shared between ',Nodes,' MPI ranks of ',nthreads,' threads'
		call flush(400)
	end if
	ipa_tw = omp_get_wtime()

	!$omp parallel do default(shared) private(i,j,l,h,eaux,vaux)
	do i=kfirst,klast

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
	
				eigt(j,i)= eaux(nocpk(i)-nv+j)

			end do
			
			do l=1,nc+nv

				do h=1,w90basis

					vector(h,l,i)=vaux(nocpk(i)-nv+l,h)

				end do

			end do

	
	end do
	!$omp end parallel do
	deallocate(eaux,vaux)
	ipa_t(1) = omp_get_wtime()-ipa_tw

	! every rank receives the band indices, energies and eigenvectors of the k points of the other ranks
#ifdef MPI
	if (Nodes > 1) then
		ipa_tw = omp_get_wtime()
		call ipa_share_kblocks_integer(nocpk,ngkpt,MPI_COMM_WORLD,MPIError)
		if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('single-particle band indices',MPIError)
		call ipa_share_kblocks_real(eigt,nc+nv,ngkpt,MPI_COMM_WORLD,MPIError,ipa_mpi_block_elements)
		if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('single-particle eigenvalues',MPIError)
		call ipa_share_kblocks_complex(vector,(nc+nv)*w90basis,ngkpt,MPI_COMM_WORLD,MPIError,ipa_mpi_block_elements)
		if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('single-particle eigenvectors',MPIError)
		ipa_t(2) = omp_get_wtime()-ipa_tw
	end if
	call MPI_REDUCE(ipa_t(1:2),ipa_tmax(1:2),2,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,MPIError)
#else
	ipa_tmax(1:2) = ipa_t(1:2)
#endif
	if (Node == 0) then
		write(400,*) 'eigenvalues and eigenvectors calculated'
		write(400,"(A,F12.3,A,F12.3,A)") ' IPA eigensystems wall time (slowest rank): ',ipa_tmax(1),' s, sharing: ',ipa_tmax(2),' s'
		call flush(400)
	end if


	!definindo os numeros quanticos dos estados



	allocate (stt(ngkpt*nc*nv,4))

	call quantumnumbers2(w90basis,ngkpt,nc,nv,nocpk,nocpk,stt)

	!allocate(esp(dimsp))
	

	!allocate(auxx(dimsp),auxy(dimsp),auxz(dimsp))
	!allocate(auyy(dimsp),auyz(dimsp),auzz(dimsp))
	
	allocate(hrx(dimsp),hry(dimsp),hrz(dimsp))
	if (Node == 0) then
	allocate(vecres(dimsp,17))

	write(301,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","xy","  ","xz","  ","yz"
	write(303,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz","  ","sp","  ","sm"
	
	if (tmcoef) then	
	write(302,*) "number of kpoints:",ngkpt
	write(302,*) "number of conduction states",nc
	write(302,*) "number of valence states",nv
	write(302,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ","ec "," ","ev"," ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","xy","  ","xz","  ","yz"
		      
	write(304,*) "number of kpoints:",ngkpt
	write(304,*) "number of conduction states",nc
	write(304,*) "number of valence states",nv
	write(304,*) "#"," ", "kx", " ", "ky"," ", "kz"," ","nocp"," ", "nc"," ", "nv","  ","ec "," ","ev","  ", "energy","  ","xx","  ",&
		      "yy","  ","zz","  ","sp","  ","sm"		   
		      
	else
	
	 continue
	end if
	end if
	
	! DFT=S: orbital centres for the atomic-gauge optical vertex
	allocate(tau(w90basis,3))
	tau = 0.0
	if (dft .eq. "S") call orbital_centres(w90basis,tau)
	! BSE_CENTER_FILE, the position matrix of the basis, as in bsesolver:
	! DFT=S: the dipoles between basis orbitals (siesta2wtb.py --rmatrix), for
	! optdip_s; DFT=W: the -i[A,H] term of the Wannier position matrix
	use_rmn = len_trim(bsecenterfile) > 0
	use_dmat = dft .eq. "S" .and. use_rmn
	if (use_rmn) then
		call rmn_center_setup(trim(bsecenterfile),dft,w90basis,rmn,rmn_ok,rmn_message)
		if (rmn_ok .and. use_dmat) then
			allocate(dmat(nvec,w90basis,w90basis,3))
			call rmn_dfts_dipoles(rmn,w90basis,nvec,rvec,rlat,ovp,tau,dmat,rmn_ok,rmn_message)
			call rmn_destroy(rmn)
		end if
		if (.not. rmn_ok) then
			write(*,*) trim(rmn_message)
			if (use_dmat) stop 'Unable to build the DFT=S orbital dipoles'
			stop 'Unable to initialize the Wannier position data'
		end if
	end if

	! DFT=S vertices also restore efermi+fermishift: H and dH/dk retain
	! the input energy reference, while eigsys returns Fermi-shifted energies.
	! The optical vertex takes the DFT energies: the scissor scs, which eigsys
	! adds to the conduction bands, shifts the transition energies (BSE
	! diagonal, IPA spectra) but not the eigenvectors, so it must not change
	! r_cv. With the shifted energies <c|dH/dk|v>/(Ec-Ev) would shrink by
	! (Ec-Ev)/(Ec-Ev+scs), and for DFT=S the -(Ec+Ev)/2 dS/dk term would change.

	! DFT=S, hoisted: dH/dk, dS/dk and D(k) depend on k only. optsp_s_kall and optdip_s_kall build those N x N matrices once per
	! k point and contract them with all (c,v) pairs of that k by matrix products (one k per thread, serial BLAS inside), instead
	! of once per transition. Each rank does the k points of its block (those of the eigensystems) and the transitions of
	! those k points (consecutive: the states are ordered by k); the vertices are shared afterwards.
	ipa_tw = omp_get_wtime()
	if (dft .eq. "S") then
		allocate(opt_hx(nc,nv,kfirst:klast),opt_hy(nc,nv,kfirst:klast),opt_hz(nc,nv,kfirst:klast))
		!$omp parallel do default(shared) private(iko) schedule(dynamic,1)
		do iko=kfirst,klast
#ifdef MKL
			call MKL_SET_NUM_THREADS(1)
#endif
#ifdef AOCL
			call bli_thread_set_num_threads(1)
#endif
#ifdef OPENBLAS
			call OPENBLAS_SET_NUM_THREADS(1)
#endif
			call optsp_s_kall(nv,nc,eigt(1:nv,iko),eigt(nv+1:nv+nc,iko),scs,efermi,fermishift,&
			     vector(:,1:nv,iko),vector(:,nv+1:nv+nc,iko),&
			     kpt(iko,1),kpt(iko,2),kpt(iko,3),sme,&
			     w90basis,nvec,rvec,hopmatrices,ihopmatrices,ovp,tau,&
			     opt_hx(:,:,iko),opt_hy(:,:,iko),opt_hz(:,:,iko))
			if (use_dmat) call optdip_s_kall(nv,nc,eigt(1:nv,iko),eigt(nv+1:nv+nc,iko),scs,sme,&
			     vector(:,1:nv,iko),vector(:,nv+1:nv+nc,iko),&
			     kpt(iko,1),kpt(iko,2),kpt(iko,3),w90basis,nvec,rvec,dmat,&
			     opt_hx(:,:,iko),opt_hy(:,:,iko),opt_hz(:,:,iko))
		end do
		!$omp end parallel do
#ifdef MKL
		call MKL_SET_NUM_THREADS(nthreads)
#endif
#ifdef AOCL
		call bli_thread_set_num_threads(nthreads)
#endif
#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(nthreads)
#endif
	end if

	! the transitions of the k points of this rank
	if (ilast >= ifirst) then
		if (stt(ifirst,4) /= kfirst .or. stt(ilast,4) /= klast) then
			write(*,*) 'IPA: the transitions of the k points ',kfirst,' to ',klast,' are not the rows ',ifirst,' to ',ilast,' of stt'
			stop 'Unexpected order of the exciton basis states'
		end if
	end if

	!$omp parallel do default(shared) private(i,rmn_ok,rmn_message)

	do i=ifirst,ilast

		if (dft .eq. "S") then
		hrx(i) = opt_hx(stt(i,3)-nv,stt(i,2),stt(i,4))
		hry(i) = opt_hy(stt(i,3)-nv,stt(i,2),stt(i,4))
		hrz(i) = opt_hz(stt(i,3)-nv,stt(i,2),stt(i,4))
		else
	       call optsp(eigt(stt(i,2),stt(i,4)),vector(:,stt(i,2),stt(i,4)),&
		     eigt(stt(i,3),stt(i,4))-scs,vector(:,stt(i,3),stt(i,4)),&
		     kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),ffactor,sme,&
		     w90basis,nvec,rlat,rvec,hopmatrices,&
		     ihopmatrices,hrx(i),hry(i),hrz(i))
		end if
		if (use_rmn .and. .not. use_dmat) then
			call rmn_apply_q0_optical_correction(rmn,kpt(stt(i,4),:),rlat,&
				eigt(stt(i,2),stt(i,4)),vector(:,stt(i,2),stt(i,4)),&
				eigt(stt(i,3),stt(i,4))-scs,vector(:,stt(i,3),stt(i,4)),sme,&
				hrx(i),hry(i),hrz(i),rmn_ok,rmn_message)
			if (.not. rmn_ok) stop 'Unable to evaluate the Q=0 Wannier position matrix'
		end if

	end do

	!$omp end parallel do

	if (dft .eq. "S") deallocate(opt_hx,opt_hy,opt_hz)
	ipa_t(3) = omp_get_wtime()-ipa_tw

	! every rank receives the vertices of the other ranks (the transitions of a k point are consecutive); only rank 0 writes
#ifdef MPI
	if (Nodes > 1) then
		ipa_tw = omp_get_wtime()
		call ipa_share_kblocks_complex(hrx,nt_k,ngkpt,MPI_COMM_WORLD,MPIError,ipa_mpi_block_elements)
		if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('optical vertices (x)',MPIError)
		call ipa_share_kblocks_complex(hry,nt_k,ngkpt,MPI_COMM_WORLD,MPIError,ipa_mpi_block_elements)
		if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('optical vertices (y)',MPIError)
		call ipa_share_kblocks_complex(hrz,nt_k,ngkpt,MPI_COMM_WORLD,MPIError,ipa_mpi_block_elements)
		if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('optical vertices (z)',MPIError)
		ipa_t(4) = omp_get_wtime()-ipa_tw
	end if
	call MPI_REDUCE(ipa_t(3:4),ipa_tmax(3:4),2,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,MPIError)
#else
	ipa_tmax(3:4) = ipa_t(3:4)
#endif

	if (Node == 0) then

	!$omp parallel do default(shared) private(i,hxsp,hysp,hzsp,hrsp,hrsm)
	do i=1,dimsp

		hxsp = hrx(i)
		hysp = hry(i)
		hzsp = hrz(i)
		hrsp = (hxsp+cmplx(0.,1.)*hysp)*(1.0/sqrt(2.))
		hrsm = (hxsp-cmplx(0.,1.)*hysp)*(1.0/sqrt(2.))

		vecres(i,1) = real(kpt(stt(i,4),1))
		vecres(i,2) = real(kpt(stt(i,4),2))
		vecres(i,3) = real(kpt(stt(i,4),3))
		vecres(i,4) = real(nocpk(stt(i,4)))
		vecres(i,5) = real(nocpk(stt(i,4))-nv+stt(i,3))
		vecres(i,6) = real(nocpk(stt(i,4))-nv+stt(i,2))		
		vecres(i,7) = real(eigt(stt(i,3),stt(i,4))-eigt(stt(i,2),stt(i,4)))		
		vecres(i,8) = real(hxsp*conjg(hxsp))		
		vecres(i,9) = real(hysp*conjg(hysp))		
		vecres(i,10) = real(hzsp*conjg(hzsp))
		vecres(i,11) = real(hxsp*conjg(hysp))
		vecres(i,12) = real(hxsp*conjg(hzsp))
		vecres(i,13) = real(hysp*conjg(hzsp))		 
		vecres(i,14) = real(hrsp*conjg(hrsp))
		vecres(i,15) = real(hrsm*conjg(hrsm))
		vecres(i,16) = real(eigt(stt(i,3),stt(i,4)))
		vecres(i,17) =	real(eigt(stt(i,2),stt(i,4)))			

	end do
	!$omp end parallel do

	call Sortrows(7,17,vecres, dimsp)
	
	do i=1,dimsp
	
		write(301,"(7F15.6)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,11),vecres(i,12),vecres(i,13)
		write(303,"(6F15.6)") vecres(i,7),vecres(i,8),vecres(i,9),vecres(i,10),vecres(i,14),vecres(i,15)
		
		if (tmcoef) then
	        write(302,"(3F10.6,3I10.0,9F10.6)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,16),vecres(i,17),vecres(i,7),&
						  vecres(i,8),vecres(i,9),vecres(i,10),&
						  vecres(i,11),vecres(i,12),vecres(i,13)
	        
	        write(304,"(3F10.6,3I10.0,8F10.6)") vecres(i,1),vecres(i,2),vecres(i,3),int(vecres(i,4)),int(vecres(i,5)),&
						  int(vecres(i,6)),vecres(i,16),vecres(i,17),vecres(i,7),vecres(i,8),&
						  vecres(i,9),vecres(i,10),&
						  vecres(i,14),vecres(i,15)						  
						  
		end if
	
	end do

	write(400,"(A,F12.3,A,F12.3,A)") ' IPA vertices wall time (slowest rank): ',ipa_tmax(3),' s, sharing: ',ipa_tmax(4),' s'
	end if

	if (Node == 0) then
	write(400,*) 'IPA single particle optics finished'
	call flush(400)	
	end if

	deallocate(eigt,vector,hrx,hry,hrz)
	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)
	deallocate(tau)
	if (use_dmat) deallocate(dmat)
	if (use_rmn .and. .not. use_dmat) call rmn_destroy(rmn)
	deallocate(ovp)

	!deallocate(auxx,auxy,auxz,auyy,auyz,auzz)

	deallocate(nocpk)
	deallocate(stt)
	!deallocate(esp)
	deallocate(kpt)
	if (allocated(vecres)) deallocate(vecres)


	call cpu_time(tf)
	call date_and_time(VALUES=values2)

	close(200)

	if (Node == 0) then

	write(400,*)
	write(400,*) 'end','   ','month',values2(2),'day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(400,*)

	close(400)



	close(300)
	close(301)
	close(303)	
	
	if (tmcoef) then
	close(302)
	close(304)
	else
	 continue
	end if

	end if
	

	


end subroutine spoptics
