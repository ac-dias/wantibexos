
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
	use input_variables, only: bsecenterfile, bsecentgrad, bseff, bsefffile, bseffexch, bseffecut, sysdim, bse_qfold, bse_nstates, bse_real, &
		bse_eph_file, bse_eph_eta, bse_eph_temp
	use bse_realbasis
	use bse_symmetry
	use bse_q_optics, only: rmn_data, rmn_read, rmn_destroy, rmn_centers, &
		center_phase_build, rmn_apply_q0_optical_correction, rmn_dfts_dipoles
	use bse_formfactor, only: ff_data, ff_read, ff_destroy, ff_setup_check, ff_mesh_check, &
		ff_overlap_check, ff_exchange_setup
	use bse_fastkernel, only: fk_active, fk_vertex, fk_setup, fk_free
	use bse_eph, only: eph_data, eph_read, eph_destroy, eph_build_kplusq, &
	                   exciton_phonon_linewidth, exciton_phonon_phassist

	implicit none

	real,parameter:: pi=acos(-1.)

	integer :: dimbse !=ngrid*ngrid*nc*nv ! dimensão da matriz bse

	real,allocatable,dimension(:,:) :: eigv
	real,allocatable,dimension(:,:) :: tau !DFT=S orbital centres
	complex,allocatable,dimension(:,:,:) :: vector
	complex,allocatable,dimension(:,:) :: center_phase

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
	complex,allocatable,dimension(:,:,:) :: opt_hx,opt_hy,opt_hz !DFT=S optical vertices <c|dH/dk-ebar dS/dk|v>/(Ec-Ev), (nc,nv,k)
	integer :: iko

	real,allocatable,dimension(:) :: actxx,actyy,actzz,actxy,actxz,actyz,actsp,actsm

	complex,parameter :: imag=cmplx(0.0,1.0)

	real :: a,r0,ed,ec,ev,egap

	integer :: ngkpt
	real,allocatable,dimension(:,:) :: vecres
		
	!real,dimension(3) :: shift !shift no mhkpack

	!variaveis relacionadas a marcacao do tempo

	real:: t0,tf
	double precision :: task_start,task_elapsed,task_elapsed_max
	! independent-particle stage shared between the MPI ranks: the block kfirst..klast of k points of this rank, the
	! transitions ifirst..ilast of those k points, band energies eigt(band,k) (contiguous per k) and the wall times of the
	! stage on this rank and on the slowest one: 1, 2 eigensystems (computing, sharing them with the other ranks), 3, 4 vertices
	! (computing, sharing)
	integer :: kfirst,klast,nt_k,ifirst,ilast
	real,allocatable,dimension(:,:) :: eigt
	real :: egap_loc
	double precision :: ipa_t(4),ipa_tmax(4),ipa_tw
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
	logical :: use_rmn,use_center_phase,use_center_grad,rmn_ok,use_dmat
	real,allocatable,dimension(:,:,:,:) :: dmat !DFT=S orbital dipoles (BSE_CENTER_FILE)
	character(len=256) :: rmn_message
	type(ff_data) :: ffdat
	logical :: use_ff,ff_ok
	character(len=256) :: ff_message
	! WTB_KERNEL= fast (default) | old | check: the DFT=S direct kernel from the vertex table of bse_fastkernel,
	! the per-element sandwich of sandwich_phase, or both on a sample of elements (then fast)
	character(len=16) :: kmode
	logical :: fastk,kcheck
	integer :: nsk,fk_nimgtot,ipass,ichk,jchk,nchk,stepchk
	! k points of the rows and columns of the Hamiltonian elements this rank evaluates (the vertex table holds those pairs only)
	logical,allocatable,dimension(:) :: fk_needr,fk_needc
	logical :: fk_upper
	integer :: fk_ilo,fk_jhi
	double precision :: fk_mib,fk_mibmax
	! sk(:,:,isk) and sk(:,:,jsk) of the two states: the index min(k,nsk) is set before
	! the call, since gfortran 15 bounds-checks (-fcheck=bounds) an array section passed
	! as an argument with a MIN subscript against a temporary that is never set
	integer :: isk,jsk
	double precision :: fk_t1,fk_t2
	complex,dimension(2) :: ekern
	real :: dchk,echk
	real :: ffspin,ffdev
	real,allocatable,dimension(:,:) :: wannier_centers

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
	real(kind=8) :: elpa_emin,elpa_emax,elpa_trace
	logical :: elpa_good
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
	character(len=10) :: coultype
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
	! ELPA solver (BSE_ALGO= elpa, elpa1, elpa2) and the number of lowest eigenpairs it computes (BSE_NSTATES=)
	logical :: use_elpa
	integer :: elpa_stage, nstates, elpa_hold
	integer,dimension(3) :: bseham_metadata
	character(len=160) :: bseham_path
#ifdef MPI
	integer,parameter :: bse_mpi_block_elements=16777216
	integer(kind=8) :: mpi_transfer_elements
#endif

	!fim modificacoes versao 2.1

	! EPH variables
	type(eph_data) :: eph
	logical :: eph_ok
	character(len=256) :: eph_message
	integer, allocatable :: kplusq_idx(:,:)
	real,    allocatable :: Gamma_eph(:)
	real,    allocatable :: chi2_xx_eph(:), chi2_yy_eph(:), chi2_zz_eph(:)
	real,    allocatable :: omega_grid_eph(:)
	complex, allocatable :: mu_x_eph(:), mu_y_eph(:), mu_z_eph(:)
	complex, allocatable :: hbse_work(:,:)
	integer :: io_eph, iomega

	! BSE_REAL= T: the Hamiltonian in the time-reversal real basis (bse_realbasis), a real symmetric matrix
	logical :: use_real,rb_ok,rb_full
	! symmetry of the system, written to the log
	type(sym_group) :: sym_g
	real(kind=8),allocatable,dimension(:,:) :: sym_tau
	character(len=8),allocatable,dimension(:) :: sym_spec
	integer,allocatable,dimension(:) :: sym_lq
	logical :: sym_have_pos,sym_have_lab
	character(len=256) :: symflag
	real(kind=8) :: symtol,symclu
	integer :: symerr
	real,allocatable,dimension(:,:) :: hbse_r,hbse_dist_r,hbse_eigenvectors_r,Zr
	real,allocatable,dimension(:) :: WORKr
	character(len=8) :: bse_eigname
	real,allocatable,dimension(:,:) :: hv_re,hv_im
	real,allocatable,dimension(:) :: rb_rvec,rb_w
	complex,allocatable,dimension(:) :: rb_cvec
	complex,allocatable,dimension(:,:) :: rb_vsc
	integer,allocatable,dimension(:) :: rb_kb
	integer :: rb_nself,rb_nmis,rbu,rbv,rbk,rbnr,rbnc
	real :: rb_devmax,rb_tgdev,rb_eref
	! time-reversal representation (rb_tr_rep): mode 1 U = 1, mode 2 spin flip; Kramers partners of the transitions
	integer :: rb_mode
	integer,allocatable,dimension(:) :: rb_pmap,rb_tbar,rb_sig
	real,allocatable,dimension(:) :: rb_psgn
	real :: rb_dev1,rb_dev2,rb_devo,rb_deve_k,rb_devo_k,rb_kdev
	integer :: rbj,rbic,rbiv
	real :: blk(2,2),blkt(2,2)
#ifdef ELPA
	real,allocatable,dimension(:,:) :: elpa_eigenvectors_r
#endif

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
	
	call hamiltonian_input_read(200,params)

	allocate(ovp(nvec,w90basis,w90basis))
	ovp = 0.0
	

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

	! BSE_CENTER_FILE, the position matrix <m,0|r|n,R> of the basis (DFT=W:
	! Wannier90's seedname_r.dat; DFT=S: siesta2wtb.py --rmatrix): the optical
	! vertex, and for DFT=W the Wannier centres of the kernel (DFT=S takes the
	! orbital centres of basis_set-*, with or without the file). BSE_FF puts
	! the form factors of the orbitals in the kernel instead.
	use_rmn=len_trim(bsecenterfile) > 0
	use_dmat=use_rmn .and. dft == 'S'
	use_ff=bseff
	use_center_phase=use_rmn .and. dft /= 'S' .and. .not. use_ff
	use_center_grad=bsecentgrad .and. use_center_phase
	if (bsecentgrad .and. .not. use_rmn) stop 'BSE_CENTER_GRAD requires BSE_CENTER_FILE'
	if (bsecentgrad .and. dft == 'S') stop 'BSE_CENTER_GRAD is implemented for DFT=W; for DFT=S use BSE_FF'
	if (bsecentgrad .and. use_ff) stop 'BSE_CENTER_GRAD and BSE_FF both correct the kernel vertices: use one of them'
	if (use_ff .and. len_trim(bsefffile) == 0) stop 'BSE_FF= T requires BSE_FF_FILE'
	if (bseffexch .and. .not. use_ff) stop 'BSE_FF_EXCHANGE= T requires BSE_FF= T'
	call bse_qfold_check(bse_qfold,rlat,use_ff)
	allocate(wannier_centers(3,w90basis))
	wannier_centers=0.0
	if (use_rmn) then
		call rmn_read(trim(bsecenterfile),rmn,rmn_ok,rmn_message)
		if (.not. rmn_ok) then
			write(*,*) trim(rmn_message)
			stop 'Unable to initialize the Wannier position data'
		end if
		if (rmn%num_wann /= w90basis) then
			write(*,*) 'Wannier position matrix basis:',rmn%num_wann,' Hamiltonian basis:',w90basis
			stop 'Incompatible Wannier position and Hamiltonian bases'
		end if
		if (use_center_phase) then
			call rmn_centers(rmn,wannier_centers,rmn_ok,rmn_message)
			if (.not. rmn_ok) then
				write(*,*) trim(rmn_message)
				stop 'Unable to initialize the Wannier centres'
			end if
		end if
		if (Node == 0 .and. use_dmat) then
			write(300,*) 'DFT=S position matrix:',trim(bsecenterfile)
			write(300,*) 'Q=0 optical vertex: orbital dipoles <c|D|v> added'
		else if (Node == 0) then
			write(300,*) 'Wannier position data:',trim(bsecenterfile)
			write(300,*) 'Q=0 position-matrix treatment: dH/dk - i[A,H]'
			if (use_center_grad) then
				write(300,*) 'G=0 direct Coulomb embedding: Wannier-centre phases + pair-midpoint dipoles'
			else if (use_center_phase) then
				write(300,*) 'G=0 direct Coulomb embedding: Wannier-centre phases enabled'
			end if
		end if
	end if
	! BSE_FF: the form factors F_mn(R; Q) of the orbitals of the basis
	! (utils/wtb_formfactor.py) in the kernel instead of point centres, for
	! DFT=W and DFT=S; BSE_FF_EXCHANGE adds the exchange term
	if (use_ff) then
		call ff_read(trim(bsefffile),ffdat,ff_ok,ff_message)
		if (ff_ok) call ff_setup_check(ffdat,w90basis,rlat,ngrid,ff_ok,ff_message)
		if (.not. ff_ok) then
			write(*,*) trim(ff_message)
			stop 'Unable to initialize the form factors (BSE_FF_FILE)'
		end if
		ffdat%direct=.true.
		if (Node == 0) then
			write(300,'(3A,I0,A,I0,A)') ' Form factors: ',trim(bsefffile),' (',ffdat%norb,' orbitals, ', &
				ffdat%nr,' lattice vectors)'
			write(300,'(A,I0,A,3(1X,I0),A,I0,A,F0.1,A)') ' Form factors: ',ffdat%ndir,' direct Q (mesh', &
				ffdat%mesh,'), ',ffdat%nexc,' exchange G up to ',ffdat%ecut,' eV'
			write(300,*) 'G=0 direct Coulomb embedding: orbital form factors (BSE_FF)'
		end if
	else if (Node == 0 .and. len_trim(bsefffile) > 0) then
		write(300,*) 'BSE_FF_FILE not read: BSE_FF= F'
	end if

	!termino parametros calculo

	ngkpt = ngrid(1)*ngrid(2)*ngrid(3)
	dimbse = ngkpt*nc*nv

	! BSE_ALGO= elpa1 (elpa is the same), elpa2: ELPA 1-stage or 2-stage solver of the distributed Hamiltonian.
	! BSE_NSTATES= n makes the ELPA solver compute only the n lowest eigenpairs (0 = all of them).
	use_elpa = .false.
	elpa_stage = 0
	select case (trim(bsealgo))
	case ("elpa","elpa1")
		use_elpa = .true.
		elpa_stage = 1
	case ("elpa2")
		use_elpa = .true.
		elpa_stage = 2
	end select
	nstates = dimbse
	if (bse_nstates > 0) then
		if (.not. use_elpa) then
			if (Node == 0) write(*,*) 'BSE_NSTATES= needs BSE_ALGO= elpa1 or elpa2; BSE_ALGO= ',trim(bsealgo),' computes all the states'
#ifdef MPI
			call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)
#endif
			stop
		end if
		nstates = min(bse_nstates,dimbse)
	end if
	if (bsewf .and. nstates < dimbse .and. excwff > nstates) then
		if (Node == 0) write(*,*) 'EXC_WF_F= ',excwff,' is larger than the ',nstates,' states of BSE_NSTATES='
#ifdef MPI
		call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)
#endif
		stop
	end if

	! BSE_REAL= T: real symmetric Hamiltonian of the time-reversal real basis (bse_realbasis).  It needs real H(R) and S(R)
	! (time-reversal invariant, no spin exchange term), and it leaves out the options that need the complex eigenvectors
	use_real = bse_real
	rb_full = .false.
	if (use_real) then
		if (bsehamwrite .or. bsehamread) call rb_refuse('BSE_REAL= T does not support BSE_HAM_SAVE / BSE_HAM_READ yet')
		if (len_trim(bse_eph_file) > 0) call rb_refuse('BSE_REAL= T does not support BSE_EPH_FILE')
		if (abs(exc) > 0.0) call rb_refuse('BSE_REAL= T needs a time-reversal invariant Hamiltonian: SPIN_EXC= must be 0')
		! time reversal on the orbitals: T = K (real H(R)) or T = i sigma_y K (spin-orbit coupling)
		allocate(rb_pmap(w90basis),rb_psgn(w90basis))
		call rb_tr_rep(dft,systype,w90basis,nvec,hopmatrices,ihopmatrices,ovp,rb_mode,rb_pmap,rb_psgn,rb_dev1,rb_dev2)
		if (rb_mode == 0) then
			if (Node == 0) then
				write(*,*) 'BSE_REAL= T: the imaginary part of H(R) is up to',rb_dev1,'eV'
				if (rb_dev2 >= 0.0) write(*,*) 'BSE_REAL= T: H(R) and S(R) differ from the spin-flip form of time reversal by',rb_dev2,'eV'
			end if
			if (rb_dev2 >= 0.0) call rb_refuse('BSE_REAL= T needs a time-reversal invariant tight-binding model: a real H(R), or the &
				&spin-orbit form H(R) = U H(R)* U^+ with the spin flip U of the spinor layout of the file')
			call rb_refuse('BSE_REAL= T needs a real H(R) (time-reversal invariant tight-binding model)')
		end if
		rb_full = use_elpa
	end if

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
	if (use_elpa) then
		write(300,*) 'ELPA solver: ',elpa_stage,'-stage'
		if (nstates < dimbse) then
			write(300,*) 'BSE_NSTATES: the',nstates,'lowest of the',dimbse,'exciton states are computed'
		else
			write(300,*) 'BSE_NSTATES: all',dimbse,'exciton states are computed'
		end if
	else
		! BSE_ALGO= cheev, cheevd, cheevr, cheevx name the serial LAPACK driver of the complex Hermitian problem.  With
		! BSE_REAL= T the real symmetric driver of the same family is called (the matrix is real), and with MPI ranks
		! ScaLAPACK is used whatever the keyword.
		select case (trim(bsealgo))
		case ("cheev")
			bse_eigname = merge('SSYEV ','CHEEV ',use_real)
		case ("cheevd")
			bse_eigname = merge('SSYEVD','CHEEVD',use_real)
		case ("cheevr")
			bse_eigname = merge('SSYEVR','CHEEVR',use_real)
		case ("cheevx")
			bse_eigname = merge('SSYEVX','CHEEVX',use_real)
		case default
			bse_eigname = '?'
		end select
		if (Nodes > 1) then
			bse_eigname = merge('PSSYEV','PCHEEV',use_real)
		end if
		if (use_real) then
			write(300,*) 'Eigensolver: ',trim(bse_eigname),' (real symmetric, BSE_REAL= T)'
		else
			write(300,*) 'Eigensolver: ',trim(bse_eigname),' (complex Hermitian)'
		end if
	end if

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
	allocate(rb_kb(ngkpt))
	rb_kb = 0
	if (use_real) then
		call rb_partners(ngkpt,kpt,rlat,1.0e-4,rb_kb,rb_nself,rb_ok)
		if (.not. rb_ok) call rb_refuse('BSE_REAL= T needs a k mesh that is symmetric under k -> -k (an unshifted mesh or a shift of 1/2)')
		! partner of the transitions at a self-conjugate k point: itself (T^2 = +1), or its Kramers partner (T^2 = -1): the
		! band j goes to j+1 (j odd) or j-1 (j even), sigma = +1 (odd), -1 (even), sigma_a = sigma_c sigma_v
		allocate(rb_tbar(nc*nv),rb_sig(nc*nv))
		do rbic=1,nc
			do rbiv=1,nv
				rbj = (rbic-1)*nv+rbiv
				rb_tbar(rbj) = rbj
				rb_sig(rbj) = 1
			end do
		end do
		if (rb_mode == 2 .and. rb_nself > 0) then
			if (mod(nv,2) /= 0 .or. mod(nc,2) /= 0) call rb_refuse('BSE_REAL= T with spin-orbit coupling needs NBANDSV and NBANDSC &
				&even (whole Kramers doublets) on a k mesh with self-conjugate k points (an unshifted mesh); a mesh shifted by half &
				&a step has none')
			do rbic=1,nc
				do rbiv=1,nv
					rbj = (rbic-1)*nv+rbiv
					rb_tbar(rbj) = (rbic-1+2*mod(rbic,2)-1)*nv+(rbiv-1+2*mod(rbiv,2))
					rb_sig(rbj) = (2*mod(rbic,2)-1)*(2*mod(rbiv,2)-1)
				end do
			end do
		end if
		call rb_setup(ngkpt,nc*nv,rb_kb,rb_tbar,rb_sig)
		if (Node == 0) then
			write(300,"(A,I0,A,I0,A)") ' BSE_REAL: time-reversal real basis, ',(ngkpt-rb_nself)/2,' pairs of k points (k, -k) and ', &
				rb_nself,' self-conjugate k points'
			if (rb_mode == 2) write(300,"(A,ES9.2,A)") ' BSE_REAL: spin-orbit coupling, time reversal T = i sigma_y K (T^2 = -1), deviation of H(R), S(R) from this form ',rb_dev2,' eV'
			write(300,*) 'BSE_REAL: the Hamiltonian is a real symmetric matrix of dimension',dimbse
			call flush(300)
		end if
	end if
	allocate(center_phase(w90basis,ngkpt))
	center_phase=cmplx(1.0,0.0)
	if (use_center_phase) then
		call center_phase_build(wannier_centers,kpt,center_phase,rmn_ok,rmn_message)
		if (.not. rmn_ok) then
			write(*,*) trim(rmn_message)
			stop 'Unable to build Wannier-centre phase table'
		end if
	end if

	! Keep the orbital index first so every eigenvector passed to the BSE
	! kernel is a contiguous column.
	allocate(eigv(ngkpt,nc+nv),vector(w90basis,nc+nv,ngkpt))
	allocate(nocpk(ngkpt))

	allocate (stt(dimbse,4))
	allocate(hrx(dimbse),hry(dimbse),hrz(dimbse),hrsp(dimbse),hrsm(dimbse))

	allocate(actxx(dimbse),actyy(dimbse),actzz(dimbse),actxy(dimbse))
	allocate(actxz(dimbse),actyz(dimbse))
	allocate(actsp(dimbse),actsm(dimbse))	
	
	
	if (Nodes == 1 .and. use_real) then
		! real symmetric solvers: the workspace is allocated at the diagonalization
		select case (bsealgo)
		case ("cheev","cheevd","cheevr","cheevx")
			allocate(W(dimbse))
		case ("elpa","elpa1","elpa2")
			write(*,*) "BSE_ALGO=",trim(bsealgo)," requires an MPI run with at least two ranks"
			stop 1
		case default
			write(*,*) "please select a subroutine for diagonalization"
			stop
		end select
	else if (Nodes == 1) then
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

		case ("elpa","elpa1","elpa2")

		write(*,*) "BSE_ALGO=",trim(bsealgo)," requires an MPI run with at least two ranks"
		stop 1
		
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

	! Independent-particle stage, shared between the MPI ranks (a serial run is one rank with every k point): H(k) c = E S(k) c
	! of different k points is independent, so each rank takes a block of consecutive k points (ipa_kblock), diagonalizes
	! them with OpenMP over the block, and then receives the other blocks. The vertices below are done the same way.
	call ipa_kblock(Nodes,Node,ngkpt,kfirst,klast)
	nt_k = nc*nv
	ifirst = (kfirst-1)*nt_k+1
	ilast = klast*nt_k
	ipa_t = 0.0d0
	allocate(eigt(nc+nv,ngkpt))
	eigt = 0.0
	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	if (Node == 0) then
		write(300,"(A,I0,A,I0,A,I0,A)") ' IPA: ',ngkpt,' k points, shared between ',Nodes,' MPI ranks of ',nthreads,' threads'
		call flush(300)
	end if
	ipa_tw = omp_get_wtime()

	!$omp parallel do default(shared) private(i,j,l,h,eaux,vaux) reduction(min:egap)
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
			
			!if (ntype .eq. 2) then

			!do kl = 1,w90basis

			!	call layercont(w90basis,vaux(kl,:),orbweight,lcount(i,kl))


			!end do

			!else 

			!	continue

			!end if

			if ((eigt(nv+1,i)-eigt(nv,i)) .le. egap) then

				egap = eigt(nv+1,i)-eigt(nv,i)

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
	ipa_t(1) = omp_get_wtime()-ipa_tw

	! BSE_REAL: the vectors of the self-conjugate k points (k = -k + G): real vectors (T^2 = +1, bse_realbasis), or Kramers
	! pairs T psi_{2j-1} = psi_{2j}, T psi_{2j} = -psi_{2j-1} (T^2 = -1, spin-orbit coupling)
	rb_devmax = 0.0
	rb_devo = 0.0
	if (use_real) then
		do i=kfirst,klast
			if (rb_kb(i) /= i) cycle
			if (rb_mode == 2) then
				if (mod(nocpk(i),2) /= 0) call rb_refuse('BSE_REAL= T with spin-orbit coupling: an odd number of occupied bands at a &
					&self-conjugate k point (the Fermi level cuts a Kramers doublet)')
				call rb_kramers_gauge(dft,w90basis,nvec,rvec,ovp,kpt(i,1),kpt(i,2),kpt(i,3),rb_pmap,rb_psgn,nc+nv,eigt(:,i), &
					vector(:,:,i),rb_deve_k,rb_devo_k)
				rb_devmax = max(rb_devmax,rb_deve_k)
				rb_devo = max(rb_devo,rb_devo_k)
				cycle
			end if
			allocate(rb_vsc(w90basis,nc+nv),rb_w(nc+nv))
			call rb_selfconj_vectors(dft,w90basis,nvec,rvec,hopmatrices,ihopmatrices,ffactor,ovp,kpt(i,1),kpt(i,2),kpt(i,3), &
				nc+nv,nocpk(i)-nv+1,rb_w,rb_vsc,erro)
			if (erro /= 0) call rb_refuse('BSE_REAL: the real eigensolver failed at a self-conjugate k point')
			do j=1,nc+nv
				rb_eref = rb_w(j)-(efermi+fermishift)
				if (nocpk(i)-nv+j > nocpk(i)) rb_eref = rb_eref+scs
				rb_devmax = max(rb_devmax,abs(rb_eref-eigt(j,i)))
				vector(:,j,i) = rb_vsc(:,j)
			end do
			deallocate(rb_vsc,rb_w)
		end do
	end if

	! every rank receives the band indices, energies and eigenvectors of the k points of the other ranks, and the smallest direct gap
#ifdef MPI
	if (Nodes > 1) then
		ipa_tw = omp_get_wtime()
		call ipa_share_kblocks_integer(nocpk,ngkpt,MPI_COMM_WORLD,MPIError)
		if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('single-particle band indices',MPIError)
		call ipa_share_kblocks_real(eigt,nc+nv,ngkpt,MPI_COMM_WORLD,MPIError,bse_mpi_block_elements)
		if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('single-particle eigenvalues',MPIError)
		call ipa_share_kblocks_complex(vector,(nc+nv)*w90basis,ngkpt,MPI_COMM_WORLD,MPIError,bse_mpi_block_elements)
		if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('single-particle eigenvectors',MPIError)
		egap_loc = egap
		call MPI_ALLREDUCE(egap_loc, egap, 1, MPI_REAL, MPI_MIN, MPI_COMM_WORLD, MPIError)
		if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('direct gap',MPIError)
		ipa_t(2) = omp_get_wtime()-ipa_tw
	end if
#endif
	if (use_real) then
#ifdef MPI
		call MPI_ALLREDUCE(MPI_IN_PLACE,rb_devmax,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,MPIError)
		call MPI_ALLREDUCE(MPI_IN_PLACE,rb_devo,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,MPIError)
#endif
		call rb_tr_gauge(ngkpt,nc+nv,w90basis,rb_kb,rb_pmap,rb_psgn,vector,eigt,nocpk,rb_tgdev,rb_nmis)
		if (Node == 0) then
			write(300,"(A,ES9.2,A,ES9.2,A)") ' BSE_REAL: max |E(k)-E(-k)| of the eigensolver ',rb_tgdev, &
				' eV, real-solver check at the self-conjugate k points ',rb_devmax,' eV'
			call flush(300)
		end if
		if (rb_mode == 2 .and. rb_devo > 1.0e-3) then
			if (Node == 0) write(*,*) 'BSE_REAL: overlap of the Kramers-paired vectors at a self-conjugate k point deviates by',rb_devo
			call rb_refuse('BSE_REAL= T with spin-orbit coupling: the Kramers partners are not orthogonal to the other bands of the &
				&window at a self-conjugate k point (an accidental degeneracy of two doublets)')
		end if
		if (rb_tgdev > 1.0e-3 .or. rb_devmax > 1.0e-3 .or. rb_nmis > 0) then
			if (Node == 0) write(*,*) 'BSE_REAL: bands of k and -k differ by',rb_tgdev,'eV (real solver:',rb_devmax,'eV, Fermi-level band mismatches:',rb_nmis,')'
			call rb_refuse('BSE_REAL= T: the band structure is not time-reversal symmetric, E(k) /= E(-k) (with spin-orbit coupling: &
				&a Kramers doublet at a self-conjugate k point is not degenerate)')
		end if
	end if
	do i=1,ngkpt
		do j=1,nc+nv
			eigv(i,j) = eigt(j,i)
		end do
	end do

#ifdef MPI
	call MPI_REDUCE(ipa_t(1:2),ipa_tmax(1:2),2,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,MPIError)
#else
	ipa_tmax(1:2) = ipa_t(1:2)
#endif
	if (Node == 0) then
		write(300,*) 'direct gap:', egap
		write(300,*) 'eigenvalues and eigenvectors calculated'
		write(300,"(A,F12.3,A,F12.3,A)") ' IPA eigensystems wall time (slowest rank): ',ipa_tmax(1),' s, sharing: ',ipa_tmax(2),' s'
		call flush(300)
	end if

	!write(*,*) "autovetores e autovalores"


	! Store each k-point overlap matrix contiguously: sk(:,:,ik).
	kmode = 'fast'
	call get_environment_variable('WTB_KERNEL',kmode)
	if (len_trim(kmode) == 0) kmode = 'fast'
	fastk = (dft .eq. "S") .and. (.not. use_ff) .and. (trim(kmode) .ne. 'old')
	kcheck = fastk .and. (trim(kmode) .eq. 'check')
	if (fastk .and. .not. kcheck) then
		! the vertex table of bse_fastkernel needs one S(k) at a time (scratch), not the table of all k
		nsk = 1
	else
		nsk = ngkpt
	end if
	allocate(sk(w90basis,w90basis,nsk))
	if (dft .eq. "S") then

		if (nsk .eq. ngkpt) then
		do i=1,ngkpt

			call overlap(w90basis,nvec,rvec,ovp,kpt(i,1),kpt(i,2),kpt(i,3),sk(:,:,i))

		end do
		end if

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

    if (Node == 0) allocate(vecres(dimbse,17))

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

	! DFT=S: orbital centres for the atomic-gauge optical vertex
	allocate(tau(w90basis,3))
	tau = 0.0
	if (dft .eq. "S") call orbital_centres(w90basis,tau)
	! the symmetry of the system, to the log: point group, operations, the part that keeps the k mesh and the Hamiltonian
	! (WTB_SYMMETRY=off in the environment skips it, WTB_SYM_TOL= sets the tolerance on the positions, Angstrom)
	if (Node == 0) then
		symflag = ''
		call get_environment_variable('WTB_SYMMETRY',symflag)
		if (trim(symflag) /= 'off') then
			symtol = 3.0d-3
			symflag = ''
			call get_environment_variable('WTB_SYM_TOL',symflag)
			if (len_trim(symflag) > 0) read(symflag,*,iostat=symerr) symtol
			symclu = symtol
			symflag = ''
			call get_environment_variable('WTB_SYM_CLUSTER',symflag)
			if (len_trim(symflag) > 0) read(symflag,*,iostat=symerr) symclu
			allocate(sym_tau(w90basis,3),sym_spec(w90basis),sym_lq(w90basis))
			sym_tau = 0.0d0
			sym_spec = ''
			sym_lq = -1
			sym_have_pos = .false.
			sym_have_lab = .false.
			if (dft .eq. "S") then
				sym_tau = dble(tau)
				sym_have_pos = any(abs(tau) > 0.0)
				call orbital_labels(w90basis,sym_spec,sym_lq,sym_have_lab)
			else if (use_center_phase) then
				sym_tau = transpose(dble(wannier_centers))
				sym_have_pos = .true.
			end if
			write(300,*)
			call sym_analyze(300,dble(rlat),w90basis,sym_tau,sym_have_pos,sym_spec,sym_lq,sym_have_lab,nvec,rvec, &
				hopmatrices,ihopmatrices,ovp,.true.,ngrid,dble(mshift),symtol,symclu,sym_g)
			write(300,*)
			deallocate(sym_tau,sym_spec,sym_lq)
		end if
	end if
	! DFT=S with BSE_CENTER_FILE: the dipoles between basis orbitals, from
	! the position matrix of siesta2wtb.py --rmatrix, for optdip_s
	if (use_dmat) then
		allocate(dmat(nvec,w90basis,w90basis,3))
		call rmn_dfts_dipoles(rmn,w90basis,nvec,rvec,rlat,ovp,tau,dmat,rmn_ok,rmn_message)
		if (.not. rmn_ok) then
			write(*,*) trim(rmn_message)
			stop 'Unable to build the DFT=S orbital dipoles'
		end if
	end if

		if (Node == 0) then
			write(300,*) 'IPA transitions: begin'
			call flush(300)
		end if
#ifdef MPI
		call MPI_BARRIER(MPI_COMM_WORLD,MPIError)
		task_start = MPI_WTIME()
#else
		task_start = omp_get_wtime()
#endif

		 ! $omp parallel default(shared) private(i,w90basis,rlat,rvec,hopmatrices,ihopmatrices)

	! DFT=S vertices also restore efermi+fermishift: H and dH/dk retain
	! the input energy reference, while eigsys returns Fermi-shifted energies.
	! The optical vertex takes the DFT energies: the scissor scs, which eigsys
	! adds to the conduction bands, shifts the transition energies (BSE
	! diagonal, IPA spectra) but not the eigenvectors, so it must not change
	! r_cv. With the shifted energies <c|dH/dk|v>/(Ec-Ev) would shrink by
	! (Ec-Ev)/(Ec-Ev+scs), and for DFT=S the -(Ec+Ev)/2 dS/dk term would change.

	! DFT=S, hoisted: dH/dk, dS/dk and D(k) depend on k only. optsp_s_kall and optdip_s_kall build those N x N matrices once per
	! k point and contract them with all (c,v) pairs of that k by matrix products (one k per thread, serial BLAS inside), instead
	! of once per transition. Each rank does the k points of its block (the same k points as the eigensystems above) and the
	! transitions of those k points (they are consecutive: the states are ordered by k); the vertices are shared afterwards.
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
			! the residual position term of BSE_CENTER_FILE, added to the vertices of this k point
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
			 !$omp parallel do private(rmn_ok,rmn_message)
	do i=ifirst,ilast

		!ec = eigv(stt(i,4),stt(i,3))
		!ev = eigv(stt(i,4),stt(i,2))

		if (dft .eq. "S") then
		hrx(i) = opt_hx(stt(i,3)-nv,stt(i,2),stt(i,4))
		hry(i) = opt_hy(stt(i,3)-nv,stt(i,2),stt(i,4))
		hrz(i) = opt_hz(stt(i,3)-nv,stt(i,2),stt(i,4))
		else
		call optsp(eigv(stt(i,4),stt(i,2)),vector(:,stt(i,2),stt(i,4)),&
		     eigv(stt(i,4),stt(i,3))-scs,vector(:,stt(i,3),stt(i,4)),&
		     kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),ffactor,sme,&
		     w90basis,nvec,rlat,rvec,hopmatrices,&
		     ihopmatrices,hrx(i),hry(i),hrz(i))
		end if
		if (use_rmn .and. .not. use_dmat) then
			call rmn_apply_q0_optical_correction(rmn,kpt(stt(i,4),:),rlat,&
				eigv(stt(i,4),stt(i,2)),vector(:,stt(i,2),stt(i,4)),&
				eigv(stt(i,4),stt(i,3))-scs,vector(:,stt(i,3),stt(i,4)),sme,&
				hrx(i),hry(i),hrz(i),rmn_ok,rmn_message)
			if (.not. rmn_ok) stop 'Unable to evaluate the Q=0 Wannier position matrix'
		end if

	end do
	 !$omp end parallel do

	if (dft .eq. "S") deallocate(opt_hx,opt_hy,opt_hz)
	deallocate(eigt)
	ipa_t(3) = omp_get_wtime()-ipa_tw

	! every rank receives the vertices of the transitions of the other ranks (the rows of a k point are consecutive in hrx)
#ifdef MPI
	if (Nodes > 1) then
		ipa_tw = omp_get_wtime()
		call ipa_share_kblocks_complex(hrx,nt_k,ngkpt,MPI_COMM_WORLD,MPIError,bse_mpi_block_elements)
		if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('optical vertices (x)',MPIError)
		call ipa_share_kblocks_complex(hry,nt_k,ngkpt,MPI_COMM_WORLD,MPIError,bse_mpi_block_elements)
		if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('optical vertices (y)',MPIError)
		call ipa_share_kblocks_complex(hrz,nt_k,ngkpt,MPI_COMM_WORLD,MPIError,bse_mpi_block_elements)
		if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('optical vertices (z)',MPIError)
		ipa_t(4) = omp_get_wtime()-ipa_tw
	end if
#endif

	! circular components, for the BSE spectra of every rank
	do i=1,dimbse
		hrsp(i) = (hrx(i)+cmplx(0.,1.)*hry(i))*(1.0/sqrt(2.))
		hrsm(i) = (hrx(i)-cmplx(0.,1.)*hry(i))*(1.0/sqrt(2.))
	end do

	! the sorted IPA files are written by rank 0 alone: it makes the table of the transitions
	if (Node == 0) then
	 !$omp parallel do
	 do i=1,dimbse
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
	 call Sortrows(7,17,vecres, dimbse)
	end if

#ifdef MPI
		task_elapsed = MPI_WTIME() - task_start
		call MPI_REDUCE(task_elapsed,task_elapsed_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,MPIError)
		call MPI_REDUCE(ipa_t(3:4),ipa_tmax(3:4),2,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,MPIError)
#else
		task_elapsed = omp_get_wtime()
		task_elapsed_max = task_elapsed - task_start
		ipa_tmax(3:4) = ipa_t(3:4)
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
		write(300,"(A,F12.3,A,F12.3,A)") ' IPA vertices wall time (slowest rank): ',ipa_tmax(3),' s, sharing: ',ipa_tmax(4),' s'
		write(300,"(A,F12.3,A)") 'IPA transitions wall time: ',task_elapsed_max,' s'
		write(300,*) 'IPA single particle optics finished'
	call flush(300)
    endif

	if (allocated(vecres)) deallocate(vecres)
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
	
	fk_active = .false.
	if (Node == 0) then
		if (bse_qfold > 1) then
			write(300,"(A,I0,A,I0,A)") ' BSE_QFOLD= ',bse_qfold,': direct kernel = sum over the ',bse_qfold**2, &
				' classes q+g+Gamma_UC of the unit-cell zone, shortest image of each (weight 1 per class)'
		end if
		if (fastk) then
			write(300,*) 'BSE kernel: vertex table from per-state vectors (bse_fastkernel), WTB_KERNEL= ',trim(kmode)
		else
			write(300,*) 'BSE kernel: sandwich_phase per matrix element, WTB_KERNEL= ',trim(kmode)
		end if
		call flush(300)
	end if
#ifdef MPI
	! BLACS grid of the distributed Hamiltonian, made here because the vertex table needs to know which rows and
	! columns of the matrix this rank holds
	if (Nodes > 1) then
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
#ifdef ELPA
		! ELPA stops with status -13 when a rank holds no row or no column of the matrix, which happens when the matrix has
		! fewer blocks than the grid has process rows or columns (N <= mb*(npcol-1)), for BSE_REAL= F as well
		if (use_elpa) then
			elpa_hold = min(locr,locc)
			call MPI_ALLREDUCE(MPI_IN_PLACE, elpa_hold, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, MPIError)
			if (elpa_hold == 0) then
				if (Node == 0) then
					write(*,*) 'BSE_ALGO= ',trim(bsealgo),': the BSE matrix is too small for ELPA on ',Nodes,' MPI ranks: ',dimbse,' x ',dimbse, &
						' elements on a ',nprow,' x ',npcol,' process grid with blocks of ',mb,' leave at least one rank without any part of it.'
					write(*,*) 'ELPA needs every rank to hold rows and columns of the matrix: use fewer MPI ranks, or BSE_ALGO= cheev (ScaLAPACK).'
					call flush(6)
				end if
				call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)
			end if
		end if
#endif
	end if
#endif

	! BSE_REAL: time-reversal symmetry of the kernel and of the orbital centres, tested on a spread of states before the
	! build (the fast kernel is not active yet): H(atilde,btilde) = conjg(H(a,b)), atilde = sigma_a abar
	if (use_real) then
		if (rb_mode == 2 .and. dft .eq. "S") then
			rb_kdev = 0.0
			do rbj=1,w90basis
				rb_kdev = max(rb_kdev,maxval(abs(tau(rbj,:)-tau(rb_pmap(rbj),:))))
			end do
			if (rb_kdev > 1.0e-3) call rb_refuse('BSE_REAL= T with spin-orbit coupling: the two spin orbitals of an orbital have different &
				&centres in the basis_set-* file')
		end if
		call rb_kernel_check(rb_kdev)
		if (Node == 0) then
			write(300,"(A,ES9.2)") ' BSE_REAL: time-reversal check of the kernel, max |H(a~,b~) - H(a,b)*| / (1 + |H|) ',rb_kdev
			call flush(300)
		end if
		if (rb_kdev > 1.0e-4) then
			if (Node == 0) write(*,*) 'BSE_REAL: the kernel is not time-reversal symmetric, deviation',rb_kdev
			call rb_refuse('BSE_REAL= T: the kernel H(a~,b~) differs from H(a,b)*: the Coulomb potential or the form factors are not &
				&time-reversal symmetric (not even in q), use BSE_REAL= F')
		end if
	end if

	! a Hamiltonian read from the file needs no kernel, except for the check
	if (fastk .and. (kcheck .or. .not. bsehamread)) then
		allocate(fk_needr(ngkpt),fk_needc(ngkpt))
		fk_needr = .true.
		fk_needc = .true.
		! only the upper triangle i <= j is built (not the full matrix of ELPA); with the states ordered by k
		! that is k1 <= k2
		fk_upper = all(stt_bse(4,2:dimbse) >= stt_bse(4,1:dimbse-1))
#ifdef ELPA
		if (Nodes > 1 .and. use_elpa) fk_upper = .false.
#endif
#ifdef MPI
		! the check samples elements of all k on rank 0, which then builds its own block as the others do
		if (Nodes > 1 .and. .not. kcheck) then
			fk_needr = .false.
			fk_needc = .false.
			if (locr > 0 .and. locc > 0) then
				fk_ilo = indxl2g(1, mb, myrow, 0, nprow)
				fk_jhi = indxl2g(locc, nb, mycol, 0, npcol)
				do li=1,locr
					ig = indxl2g(li, mb, myrow, 0, nprow)
					if (ig <= fk_jhi .or. .not. fk_upper) fk_needr(stt_bse(4,ig)) = .true.
				end do
				do lj=1,locc
					jg = indxl2g(lj, nb, mycol, 0, npcol)
					if (jg >= fk_ilo .or. .not. fk_upper) fk_needc(stt_bse(4,jg)) = .true.
				end do
			end if
		end if
#endif
		if (use_real) then
			! BSE_REAL: the rows of the real basis are the representative states a (k <= kbar) and the columns the states
			! b and bbar, so the (k1, k2) pairs are not the upper triangle of the k indices
			fk_upper = .false.
			if (.not. kcheck) then
				fk_needr = .false.
				fk_needc = .false.
				if (Nodes == 1) then
					do rbk=1,ngkpt
						fk_needr(rbk) = (rb_kb(rbk) >= rbk)
					end do
					fk_needc = .true.
				else
#ifdef MPI
					do li=1,locr
						ig = indxl2g(li, mb, myrow, 0, nprow)
						rbk = stt_bse(4,rb_row(ig))
						fk_needr(rbk) = .true.
						fk_needc(rbk) = .true.
						fk_needc(rb_kb(rbk)) = .true.
					end do
					do lj=1,locc
						jg = indxl2g(lj, nb, mycol, 0, npcol)
						rbk = stt_bse(4,rb_row(jg))
						fk_needr(rbk) = .true.
						fk_needc(rbk) = .true.
						fk_needc(rb_kb(rbk)) = .true.
					end do
#endif
				end if
			end if
		end if
		call fk_setup(w90basis,nvec,nc+nv,ngkpt,rvec,ovp,kpt_bse,rlat,vector,tau,max(1,bse_qfold),fk_needr,fk_needc, &
			fk_upper,fk_nimgtot,fk_t1,fk_t2)
		fk_active = .true.
		deallocate(fk_needr,fk_needc)
		fk_mib = dble(size(fk_vertex,kind=8))*8.d0/1048576.d0
		fk_mibmax = fk_mib
#ifdef MPI
		call MPI_REDUCE(fk_mib,fk_mibmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,MPIError)
#endif
		if (Node == 0) then
			write(300,"(A,F10.3,A,F10.3,A)") ' fast kernel: per-state vectors ',fk_t1,' s, vertex table ',fk_t2,' s'
			write(300,"(A,I0,A,F12.1,A,F12.1,A)") ' fast kernel: ',fk_nimgtot,' (k1,k2,image) on rank 0, ',fk_mib, &
				' MiB on rank 0, ',fk_mibmax,' MiB on the largest rank'
			call flush(300)
		end if
	end if

	! WTB_KERNEL=check: the same elements from the vertex table and from sandwich_phase (all of the
	! elements when dimbse <= 60, otherwise a regular sample), reported on rank 0
	if (kcheck .and. Node == 0) then
		stepchk = max(1,dimbse/60)
		nchk = 0
		dchk = 0.0
		echk = 0.0
		do ichk=1,dimbse,stepchk
			do jchk=ichk,dimbse,stepchk
				isk = min(stt_bse(4,ichk),nsk)
				jsk = min(stt_bse(4,jchk),nsk)
				do ipass=1,2
					fk_active = (ipass == 2)
					ekern(ipass)= matrizelbse(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,rlat,stt_bse(:,ichk),&
					    eigv(stt_bse(4,ichk),stt_bse(3,ichk)),eigv(stt_bse(4,ichk),stt_bse(2,ichk)),&
					    vector(:,stt_bse(3,ichk),stt_bse(4,ichk)),vector(:,stt_bse(2,ichk),stt_bse(4,ichk)),&
					    kpt_bse(:,stt_bse(4,ichk)),stt_bse(:,jchk),eigv(stt_bse(4,jchk),stt_bse(3,jchk)),&
					    eigv(stt_bse(4,jchk),stt_bse(2,jchk)),vector(:,stt_bse(3,jchk),stt_bse(4,jchk)),&
					    vector(:,stt_bse(2,jchk),stt_bse(4,jchk)),kpt_bse(:,stt_bse(4,jchk)),dft,nvec,rvec,&
					    sk(:,:,isk),sk(:,:,jsk),use_center_phase, &
					    center_phase(:,stt_bse(4,ichk)),center_phase(:,stt_bse(4,jchk)), &
					    use_center_grad,rmn,wannier_centers,ffdat)
				end do
				fk_active = .true.
				nchk = nchk + 1
				dchk = max(dchk,abs(ekern(1)-ekern(2)))
				echk = max(echk,abs(ekern(1)))
			end do
		end do
		write(300,"(A,I0,A,ES10.3,A,ES10.3,A,ES10.3)") ' WTB_KERNEL=check: ',nchk,' elements, max |old - fast| = ',dchk, &
			', max |old| = ',echk,', relative ',dchk/max(echk,tiny(1.0))
		write(*,"(A,I0,A,ES10.3,A,ES10.3)") ' WTB_KERNEL=check: ',nchk,' elements, max |old - fast| = ',dchk, &
			', max |old| = ',echk
		call flush(300)
	end if

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
	task_start = omp_get_wtime()
#endif

	if (Node == 0) then
		if (use_real) then
			call bse_hamiltonian_memory_report_real(300,'BSE Hamiltonian, global dense equivalent (real symmetric)',dimbse,dimbse)
		else
			call bse_hamiltonian_memory_report(300,'BSE Hamiltonian, global dense equivalent',dimbse,dimbse)
		end if
	end if

	! BSE_FF: every q of the mesh in the file, the overlaps of the BSE bands
	! from F(R; 0), and the vectors of the exchange term
	if (use_ff) then
		call ff_mesh_check(ffdat,kpt_bse,rlat,ff_ok,ff_message)
		if (.not. ff_ok) then
			write(*,*) trim(ff_message)
			stop 'The form factors do not hold the q of the BSE mesh (BSE_FF_FILE)'
		end if
		call ff_overlap_check(ffdat,stt_bse,vector,kpt_bse,rlat,ffdev)
		if (Node == 0) write(300,'(A,ES9.2)') ' Form factors: max |<n k|m k> - delta_nm| of the BSE bands '// &
			'from F(R;0): ',ffdev
		if (bseffexch .and. .not. bsehamread) then
			ffspin=1.0
			if (systype == "NP") ffspin=2.0
			call ff_exchange_setup(ffdat,stt_bse,vector,kpt_bse,rlat,ngrid,sysdim,ffspin,bseffecut,ff_ok,ff_message)
			if (.not. ff_ok) then
				write(*,*) trim(ff_message)
				stop 'Unable to set up the exchange term (BSE_FF_EXCHANGE)'
			end if
			ffdat%exchange=.true.
			if (Node == 0) write(300,'(A,I0,A,F0.1,A,F0.1)') ' Exchange term (BSE_FF_EXCHANGE): ',ffdat%ng, &
				' G up to ',ffdat%gmax,' eV, spin factor ',ffspin
		end if
	end if

    if (Nodes == 1 .and. use_real) then
		! BSE_REAL: real symmetric matrix of the real basis, by units (an s,t couple or a single function), upper triangle
		allocate(hbse_r(dimbse,dimbse))
		hbse_r = 0.0
		!$omp parallel do private(rbu,blk) schedule(dynamic)
		do rbv=1,rb_nunit
			do rbu=1,rbv
				call rb_block(rb_ustart(rbu),rb_ustart(rbv),blk)
				hbse_r(rb_ustart(rbu):rb_ustart(rbu)+rb_ulen(rbu)-1,rb_ustart(rbv):rb_ustart(rbv)+rb_ulen(rbv)-1) = &
					blk(1:rb_ulen(rbu),1:rb_ulen(rbv))
			end do
		end do
		!$omp end parallel do
    else if (Nodes == 1) then
		bseham_metadata = (/ dimbse,1,1 /)
		if (use_center_phase) bseham_metadata(3) = 3
		if (use_center_grad) bseham_metadata(3) = 7   ! 5: the kernel before the pair-midpoint term
		if (use_ff) bseham_metadata(3) = bseham_metadata(3)+16   ! BSE_FF
		if (use_ff .and. bseffexch) bseham_metadata(3) = bseham_metadata(3)+32   ! BSE_FF_EXCHANGE
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
			!$omp parallel do private(i,j,isk,jsk) schedule(dynamic)
			do j=1,dimbse
				jsk = min(stt_bse(4,j),nsk)
				do i=1,j
					isk = min(stt_bse(4,i),nsk)
					hbse(i,j)= matrizelbse(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,rlat,stt_bse(:,i),eigv(stt_bse(4,i)&
					    ,stt_bse(3,i)),eigv(stt_bse(4,i),stt_bse(2,i)),vector(:,stt_bse(3,i),stt_bse(4,i)),&
					    vector(:,stt_bse(2,i),stt_bse(4,i)),kpt_bse(:,stt_bse(4,i)),stt_bse(:,j),eigv(stt_bse(4,j),stt_bse(3,j))&
					    ,eigv(stt_bse(4,j),stt_bse(2,j)) &
					    ,vector(:,stt_bse(3,j),stt_bse(4,j)),vector(:,stt_bse(2,j),stt_bse(4,j)),kpt_bse(:,stt_bse(4,j)),dft,nvec,rvec,&
					    sk(:,:,isk),sk(:,:,jsk),use_center_phase, &
					    center_phase(:,stt_bse(4,i)),center_phase(:,stt_bse(4,j)), &
					    use_center_grad,rmn,wannier_centers,ffdat)
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
		if (use_real) then
		! BSE_REAL: the same blocks, on the BLACS grid; an s,t couple never splits across blocks (even block size)
		if (Node == 0) then
			call bse_hamiltonian_memory_report_real(300,'BSE Hamiltonian on MPI rank 0 (real symmetric)',lld,max(1,locc))
		end if
		allocate(hbse_dist_r(lld,max(1,locc)))
		call DESCINIT(desca, dimbse, dimbse, mb, nb, 0, 0, blacs_ctxt, lld, INFO)
		call DESCINIT(descz, dimbse, dimbse, mb, nb, 0, 0, blacs_ctxt, lld, INFO)
		hbse_dist_r = 0.0
		!$omp parallel do private(li,ig,lj,jg,rbnr,rbnc,blk,blkt) schedule(dynamic)
		do lj=1,locc
			jg = indxl2g(lj, nb, mycol, 0, npcol)
			if (rb_kind(jg) == 2) cycle
			rbnc = 1
			if (rb_kind(jg) == 1) rbnc = 2
			do li=1,locr
				ig = indxl2g(li, mb, myrow, 0, nprow)
				if (rb_kind(ig) == 2) cycle
				rbnr = 1
				if (rb_kind(ig) == 1) rbnr = 2
				if (ig <= jg) then
					call rb_block(ig,jg,blk)
					hbse_dist_r(li:li+rbnr-1,lj:lj+rbnc-1) = blk(1:rbnr,1:rbnc)
				else if (rb_full) then
					call rb_block(jg,ig,blkt)
					hbse_dist_r(li:li+rbnr-1,lj:lj+rbnc-1) = transpose(blkt(1:rbnc,1:rbnr))
				end if
			end do
		end do
		!$omp end parallel do
		else
		! the BLACS grid and the local sizes (mb, nb, locr, locc, lld) are set before the vertex table
		if (Node == 0) then
			call bse_hamiltonian_memory_report(300,'BSE Hamiltonian on MPI rank 0',lld,max(1,locc))
		end if
		allocate(hbse_dist(lld,max(1,locc)))
		call DESCINIT(desca, dimbse, dimbse, mb, nb, 0, 0, blacs_ctxt, lld, INFO)
		call DESCINIT(descz, dimbse, dimbse, mb, nb, 0, 0, blacs_ctxt, lld, INFO)
		bseham_metadata = (/ dimbse,1,1 /)
		if (use_elpa) bseham_metadata(3) = 2
		if (use_center_phase) bseham_metadata(3) = bseham_metadata(3)+2
		if (use_center_grad) bseham_metadata(3) = bseham_metadata(3)+4   ! +2: before the pair-midpoint term
		if (use_ff) bseham_metadata(3) = bseham_metadata(3)+16   ! BSE_FF
		if (use_ff .and. bseffexch) bseham_metadata(3) = bseham_metadata(3)+32   ! BSE_FF_EXCHANGE
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
			!$omp parallel do private(li,ig,lj,jg,isk,jsk) schedule(dynamic)
			do lj=1,locc
				jg = indxl2g(lj, nb, mycol, 0, npcol)
				jsk = min(stt_bse(4,jg),nsk)
				do li=1,locr
					ig = indxl2g(li, mb, myrow, 0, nprow)
					isk = min(stt_bse(4,ig),nsk)

#ifdef ELPA
					! ELPA requires the complete Hermitian matrix.  PCHEEV only
					! reads its upper triangle, so preserve the cheaper construction
					! for every non-ELPA distributed run.
					if (ig <= jg .or. use_elpa) then
#else
					if (ig <= jg) then
#endif
						hbse_dist(li,lj)= matrizelbse(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,rlat,stt_bse(:,ig),&
						    eigv(stt_bse(4,ig),stt_bse(3,ig)),eigv(stt_bse(4,ig),stt_bse(2,ig)),vector(:,stt_bse(3,ig),stt_bse(4,ig)),&
						    vector(:,stt_bse(2,ig),stt_bse(4,ig)),kpt_bse(:,stt_bse(4,ig)),stt_bse(:,jg),eigv(stt_bse(4,jg),stt_bse(3,jg)),&
						    eigv(stt_bse(4,jg),stt_bse(2,jg)),vector(:,stt_bse(3,jg),stt_bse(4,jg)),vector(:,stt_bse(2,jg),stt_bse(4,jg)),&
						    kpt_bse(:,stt_bse(4,jg)),dft,nvec,rvec,sk(:,:,isk),sk(:,:,jsk), &
						    use_center_phase,center_phase(:,stt_bse(4,ig)),center_phase(:,stt_bse(4,jg)), &
						    use_center_grad,rmn,wannier_centers,ffdat)
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
		end if
#else
		stop "MPI/ScaLAPACK path requested without MPI support"
#endif
    end if
	if (use_ff) call ff_destroy(ffdat)
	! the vertex table is not needed by the diagonalization, which takes the memory
	call fk_free()

#ifdef MPI
	task_elapsed = MPI_WTIME() - task_start
	call MPI_REDUCE(task_elapsed,task_elapsed_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,MPIError)
#else
	task_elapsed = omp_get_wtime()
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
	task_start = omp_get_wtime()
#endif

    if (use_real) then
		if (Nodes == 1) then
			select case (bsealgo)
			case ("cheev")
				LWORK = -1
				allocate(WORKr(1))
				CALL SSYEV('Vectors','U',dimbse,hbse_r,dimbse,W,WORKr,LWORK,INFO)
				LWORK = rb_lwork(WORKr(1),max(1_8,3_8*dimbse-1_8))
				if (LWORK < 0) then
					write(*,*) 'SSYEV needs a workspace that does not fit an INTEGER(4): use BSE_ALGO= elpa1 or run with MPI ranks'
					stop
				end if
				deallocate(WORKr)
				allocate(WORKr(LWORK))
				CALL SSYEV('Vectors','U',dimbse,hbse_r,dimbse,W,WORKr,LWORK,INFO)
				if (INFO .ne. 0) then
					write(*,*) 'The algorithm failed to compute eigenvalues (SSYEV, INFO =',INFO,')'
					stop
				end if
				deallocate(WORKr)
			case ("cheevd")
				LWORK = -1
				LIWORK = -1
				allocate(WORKr(1),IWORK(1))
				CALL SSYEVD('Vectors','U',dimbse,hbse_r,dimbse,W,WORKr,LWORK,IWORK,LIWORK,INFO)
				! documented minimum for 'Vectors': 1+6N+2N**2 and 3+5N; the query returns the first in a REAL (see rb_lwork)
				LWORK = rb_lwork(WORKr(1),1_8+6_8*dimbse+2_8*dimbse*dimbse)
				LIWORK = max(1,IWORK(1),3+5*dimbse)
				if (LWORK < 0) then
					write(*,*) 'SSYEVD needs a workspace that does not fit an INTEGER(4) (N > 32766): use BSE_ALGO= elpa1 or run with MPI ranks'
					stop
				end if
				deallocate(WORKr,IWORK)
				allocate(WORKr(LWORK),IWORK(LIWORK))
				CALL SSYEVD('Vectors','U',dimbse,hbse_r,dimbse,W,WORKr,LWORK,IWORK,LIWORK,INFO)
				if (INFO .ne. 0) then
					write(*,*) 'The algorithm failed to compute eigenvalues (SSYEVD, INFO =',INFO,')'
					stop
				end if
				if (W(dimbse) .eq. 0) then
					write(*,*) 'The algorithm failed to compute eigenvalues. Change BSE_ALGO to cheev'
					stop
				end if
				deallocate(WORKr,IWORK)
			case ("cheevr")
				LWORK = -1
				LIWORK = -1
				allocate(WORKr(1),IWORK(1),ISUPPZ(2*dimbse),Zr(dimbse,dimbse))
				CALL SSYEVR('Vectors','A','U',dimbse,hbse_r,dimbse,VL,VU,IL,IU,ABSTOL,M,W,Zr,dimbse,ISUPPZ,WORKr,LWORK,IWORK,LIWORK,INFO)
				LWORK = rb_lwork(WORKr(1),26_8*dimbse)
				LIWORK = max(1,IWORK(1),10*dimbse)
				if (LWORK < 0) then
					write(*,*) 'SSYEVR needs a workspace that does not fit an INTEGER(4): use BSE_ALGO= elpa1 or run with MPI ranks'
					stop
				end if
				deallocate(WORKr,IWORK)
				allocate(WORKr(LWORK),IWORK(LIWORK))
				CALL SSYEVR('Vectors','A','U',dimbse,hbse_r,dimbse,VL,VU,IL,IU,ABSTOL,M,W,Zr,dimbse,ISUPPZ,WORKr,LWORK,IWORK,LIWORK,INFO)
				if (INFO .ne. 0) then
					write(*,*) 'The algorithm failed to compute eigenvalues (SSYEVR, INFO =',INFO,')'
					stop
				end if
				if (W(dimbse) .eq. 0) then
					write(*,*) 'The algorithm failed to compute eigenvalues. Change BSE_ALGO to cheev'
					stop
				end if
				hbse_r = Zr
				deallocate(WORKr,IWORK,ISUPPZ,Zr)
			case ("cheevx")
				LWORK = -1
				allocate(WORKr(1),IWORK(5*dimbse),ifailx(dimbse),Zr(dimbse,dimbse))
				CALL SSYEVX('Vectors','A','U',dimbse,hbse_r,dimbse,VL,VU,IL,IU,ABSTOL,M,W,Zr,dimbse,WORKr,LWORK,IWORK,ifailx,INFO)
				LWORK = rb_lwork(WORKr(1),8_8*dimbse)
				if (LWORK < 0) then
					write(*,*) 'SSYEVX needs a workspace that does not fit an INTEGER(4): use BSE_ALGO= elpa1 or run with MPI ranks'
					stop
				end if
				deallocate(WORKr)
				allocate(WORKr(LWORK))
				CALL SSYEVX('Vectors','A','U',dimbse,hbse_r,dimbse,VL,VU,IL,IU,ABSTOL,M,W,Zr,dimbse,WORKr,LWORK,IWORK,ifailx,INFO)
				if (INFO .ne. 0) then
					write(*,*) 'The algorithm failed to compute eigenvalues (SSYEVX, INFO =',INFO,')'
					stop
				end if
				if (W(dimbse) .eq. 0) then
					write(*,*) 'The algorithm failed to compute eigenvalues. Change BSE_ALGO to cheev'
					stop
				end if
				hbse_r = Zr
				deallocate(WORKr,IWORK,ifailx,Zr)
			end select
		else
#ifdef MPI
#ifdef ELPA
			if (use_elpa) then
allocate(elpa_eigenvectors_r(lld,max(1,locc)),stat=erro)
			if (erro /= 0) then
				write(*,*) 'Could not allocate ELPA eigenvector storage'
				call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)
			end if
			! with BSE_NSTATES= below the dimension ELPA fills only the first nstates columns and eigenvalues
			elpa_eigenvectors_r = 0.0
			W = 0.0
			! spectrum bounds and trace of the Hamiltonian, to check what ELPA returns
			call bse_elpa_bounds_real(dimbse,hbse_dist_r,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,elpa_emin,elpa_emax,elpa_trace,MPIError)

			elpa_status = elpa_init(ELPA_API_VERSION)
			call elpa_check(elpa_status, 'elpa_init', MPI_COMM_WORLD)
			elpa_instance => elpa_allocate(elpa_status)
			call elpa_check(elpa_status, 'elpa_allocate', MPI_COMM_WORLD)

			call elpa_instance%set('na', dimbse, elpa_status)
			call elpa_check(elpa_status, 'set na', MPI_COMM_WORLD)
			call elpa_instance%set('nev', nstates, elpa_status)
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
			call elpa_instance%set('omp_threads', nthreads, elpa_status)
			call elpa_check(elpa_status, 'set omp_threads', MPI_COMM_WORLD)
			elpa_status = elpa_instance%setup()
			call elpa_check(elpa_status, 'elpa setup', MPI_COMM_WORLD)
			if (elpa_stage == 2) then
				call elpa_instance%set('solver', ELPA_SOLVER_2STAGE, elpa_status)
				call elpa_check(elpa_status, 'set ELPA 2-stage solver', MPI_COMM_WORLD)
			else
				call elpa_instance%set('solver', ELPA_SOLVER_1STAGE, elpa_status)
				call elpa_check(elpa_status, 'set ELPA 1-stage solver', MPI_COMM_WORLD)
			end if
			call elpa_instance%eigenvectors(hbse_dist_r, W, elpa_eigenvectors_r, elpa_status)
			call elpa_check(elpa_status, 'ELPA eigenvectors', MPI_COMM_WORLD)
			call elpa_deallocate(elpa_instance, elpa_status)
			call elpa_check(elpa_status, 'elpa_deallocate', MPI_COMM_WORLD)
			call elpa_uninit(elpa_status)
			call elpa_check(elpa_status, 'elpa_uninit', MPI_COMM_WORLD)
			call bse_elpa_verify(dimbse,nstates,W,elpa_emin,elpa_emax,elpa_trace,Node,elpa_good)
			! every rank decides alike, and none aborts before rank 0 has written (and flushed) the reason: an abort of
			! another rank kills it with the message still in its buffer
			call MPI_ALLREDUCE(MPI_IN_PLACE, elpa_good, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, MPIError)
			if (.not. elpa_good) call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)

			! only the nstates lowest eigenpairs are defined: clear the other columns so that the
			! optical amplitudes below get nothing from them
			if (nstates < dimbse) then
				W(nstates+1:dimbse) = 0.0
				do lj=1,locc
					jg = indxl2g(lj, nb, mycol, 0, npcol)
					if (jg > nstates) elpa_eigenvectors_r(:,lj) = 0.0
				end do
			end if

			deallocate(hbse_dist_r)
			call move_alloc(elpa_eigenvectors_r, hbse_dist_r)
			else
#else
			if (use_elpa) then
				write(*,*) 'BSE_ALGO=',trim(bsealgo),' was requested, but this executable was built without ELPA support'
				call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)
			end if
#endif
			! PSSYEV overwrites Z before the back transformation is complete: A and Z are distinct distributed arrays
			allocate(hbse_eigenvectors_r(lld,max(1,locc)),stat=erro)
			if (erro /= 0) then
				write(*,*) 'Could not allocate distributed BSE eigenvector storage'
				call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)
			end if
			LWORK = -1
			allocate(WORKr(1))
			call PSSYEV('V','U',dimbse,hbse_dist_r,1,1,desca,W,hbse_eigenvectors_r,1,1,descz,WORKr,LWORK,INFO)
			LWORK = rb_lwork(WORKr(1),1_8)
			if (LWORK < 0) then
				write(*,*) 'ScaLAPACK PSSYEV needs a workspace that does not fit an INTEGER(4)'
				call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)
			end if
			deallocate(WORKr)
			allocate(WORKr(LWORK))
			call PSSYEV('V','U',dimbse,hbse_dist_r,1,1,desca,W,hbse_eigenvectors_r,1,1,descz,WORKr,LWORK,INFO)
			if (INFO .ne. 0) then
				write(*,*) 'ScaLAPACK PSSYEV failed with INFO = ', INFO
				call MPI_ABORT(MPI_COMM_WORLD, INFO, MPIError)
			end if
			deallocate(WORKr,hbse_dist_r)
			call move_alloc(hbse_eigenvectors_r,hbse_dist_r)
#ifdef ELPA
			end if
#endif
#endif
		end if
    else if (Nodes == 1 ) then
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
		LWORK = 2*dimbse + dimbse**2
		LRWORK = 1 + 5*dimbse + 2*dimbse**2
		LIWORK = 3 + 5*dimbse
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

		case ("elpa","elpa1","elpa2")

		write(*,*) "BSE_ALGO=",trim(bsealgo)," requires an MPI run with at least two ranks"
		stop 1
		
		case default
		
		 write(*,*) "please select a subroutine for diagonalization"
		 stop
		
		end select
	else
#ifdef MPI

#ifdef ELPA
		if (use_elpa) then
			allocate(elpa_eigenvectors(lld,max(1,locc)),stat=erro)
			if (erro /= 0) then
				write(*,*) 'Could not allocate ELPA eigenvector storage'
				call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)
			end if
			! with BSE_NSTATES= below the dimension ELPA fills only the first nstates columns and eigenvalues
			elpa_eigenvectors = cmplx(0.0,0.0)
			W = 0.0
			! spectrum bounds and trace of the Hamiltonian, to check what ELPA returns
			call bse_elpa_bounds(dimbse,hbse_dist,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,elpa_emin,elpa_emax,elpa_trace,MPIError)

			elpa_status = elpa_init(ELPA_API_VERSION)
			call elpa_check(elpa_status, 'elpa_init', MPI_COMM_WORLD)
			elpa_instance => elpa_allocate(elpa_status)
			call elpa_check(elpa_status, 'elpa_allocate', MPI_COMM_WORLD)

			call elpa_instance%set('na', dimbse, elpa_status)
			call elpa_check(elpa_status, 'set na', MPI_COMM_WORLD)
			call elpa_instance%set('nev', nstates, elpa_status)
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
			call elpa_instance%set('omp_threads', nthreads, elpa_status)
			call elpa_check(elpa_status, 'set omp_threads', MPI_COMM_WORLD)
			elpa_status = elpa_instance%setup()
			call elpa_check(elpa_status, 'elpa setup', MPI_COMM_WORLD)
			if (elpa_stage == 2) then
				call elpa_instance%set('solver', ELPA_SOLVER_2STAGE, elpa_status)
				call elpa_check(elpa_status, 'set ELPA 2-stage solver', MPI_COMM_WORLD)
			else
				call elpa_instance%set('solver', ELPA_SOLVER_1STAGE, elpa_status)
				call elpa_check(elpa_status, 'set ELPA 1-stage solver', MPI_COMM_WORLD)
			end if
			call elpa_instance%eigenvectors(hbse_dist, W, elpa_eigenvectors, elpa_status)
			call elpa_check(elpa_status, 'ELPA eigenvectors', MPI_COMM_WORLD)
			call elpa_deallocate(elpa_instance, elpa_status)
			call elpa_check(elpa_status, 'elpa_deallocate', MPI_COMM_WORLD)
			call elpa_uninit(elpa_status)
			call elpa_check(elpa_status, 'elpa_uninit', MPI_COMM_WORLD)
			call bse_elpa_verify(dimbse,nstates,W,elpa_emin,elpa_emax,elpa_trace,Node,elpa_good)
			! every rank decides alike, and none aborts before rank 0 has written (and flushed) the reason: an abort of
			! another rank kills it with the message still in its buffer
			call MPI_ALLREDUCE(MPI_IN_PLACE, elpa_good, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, MPIError)
			if (.not. elpa_good) call MPI_ABORT(MPI_COMM_WORLD, 1, MPIError)

			! only the nstates lowest eigenpairs are defined: clear the other columns so that the
			! optical amplitudes below get nothing from them
			if (nstates < dimbse) then
				W(nstates+1:dimbse) = 0.0
				do lj=1,locc
					jg = indxl2g(lj, nb, mycol, 0, npcol)
					if (jg > nstates) elpa_eigenvectors(:,lj) = cmplx(0.0,0.0)
				end do
			end if

			deallocate(hbse_dist)
			call move_alloc(elpa_eigenvectors, hbse_dist)
		else
#else
		if (use_elpa) then
			write(*,*) 'BSE_ALGO=',trim(bsealgo),' was requested, but this executable was built without ELPA support'
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
	task_elapsed = omp_get_wtime()
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
	
		if (use_real) then
		continue
		else if (Nodes == 1) then
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

		if (use_real) then
			! optical amplitudes in the real basis: sum_r a_r h'_r with h' = U^+ h (bse_realbasis)
			allocate(hv_re(dimbse,5),hv_im(dimbse,5))
			call rb_vector(hrx,hv_re(:,1),hv_im(:,1))
			call rb_vector(hry,hv_re(:,2),hv_im(:,2))
			call rb_vector(hrz,hv_re(:,3),hv_im(:,3))
			call rb_vector(hrsp,hv_re(:,4),hv_im(:,4))
			call rb_vector(hrsm,hv_re(:,5),hv_im(:,5))
			if (Nodes == 1) then
				if (Node == 0) then
					if (bsewf) then
						allocate(rb_cvec(dimbse))
						do i=excwf0,excwff
							call rb_expand(hbse_r(:,i),rb_cvec)
							call excwf(outputfolder,ngkpt,kpt,nc,nv,nocpk,stt,W(i),i,rb_cvec)
						end do
						deallocate(rb_cvec)
					end if

					call dielbsev_real(nthreads,dimbse,hbse_r,hv_re(:,1),hv_im(:,1),actxx)
					call dielbsev_real(nthreads,dimbse,hbse_r,hv_re(:,2),hv_im(:,2),actyy)
					call dielbsev_real(nthreads,dimbse,hbse_r,hv_re(:,3),hv_im(:,3),actzz)

					if (dtfull) then
						call dielbsep_real(nthreads,dimbse,hbse_r,hv_re(:,1),hv_im(:,1),hv_re(:,2),hv_im(:,2),actxy)
						call dielbsep_real(nthreads,dimbse,hbse_r,hv_re(:,1),hv_im(:,1),hv_re(:,3),hv_im(:,3),actxz)
						call dielbsep_real(nthreads,dimbse,hbse_r,hv_re(:,2),hv_im(:,2),hv_re(:,3),hv_im(:,3),actyz)
					else
						actxy = 0.0
						actxz = 0.0
						actyz = 0.0
					end if

					if (cpol) then
						call dielbsev_real(nthreads,dimbse,hbse_r,hv_re(:,4),hv_im(:,4),actsp)
						call dielbsev_real(nthreads,dimbse,hbse_r,hv_re(:,5),hv_im(:,5),actsm)
					else
						actsp = 0.0
						actsm = 0.0
					end if
				end if
			else
#ifdef MPI
				if (bsewf) then
					allocate(rb_rvec(dimbse),rb_cvec(dimbse))
					do i=excwf0,excwff
						call bse_eigenvector_column_dist_real(dimbse,hbse_dist_r,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,i,rb_rvec,MPIError)
						if (Node == 0) then
							call rb_expand(rb_rvec,rb_cvec)
							call excwf(outputfolder,ngkpt,kpt,nc,nv,nocpk,stt,W(i),i,rb_cvec)
						end if
					end do
					deallocate(rb_rvec,rb_cvec)
				end if

				call dielbsev_dist_real(dimbse,hbse_dist_r,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hv_re(:,1),hv_im(:,1),actxx,MPIError)
				call dielbsev_dist_real(dimbse,hbse_dist_r,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hv_re(:,2),hv_im(:,2),actyy,MPIError)
				call dielbsev_dist_real(dimbse,hbse_dist_r,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hv_re(:,3),hv_im(:,3),actzz,MPIError)

				if (dtfull) then
					call dielbsep_dist_real(dimbse,hbse_dist_r,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hv_re(:,1),hv_im(:,1),hv_re(:,2),hv_im(:,2),actxy,MPIError)
					call dielbsep_dist_real(dimbse,hbse_dist_r,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hv_re(:,1),hv_im(:,1),hv_re(:,3),hv_im(:,3),actxz,MPIError)
					call dielbsep_dist_real(dimbse,hbse_dist_r,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hv_re(:,2),hv_im(:,2),hv_re(:,3),hv_im(:,3),actyz,MPIError)
				else
					actxy = 0.0
					actxz = 0.0
					actyz = 0.0
				end if

				if (cpol) then
					call dielbsev_dist_real(dimbse,hbse_dist_r,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hv_re(:,4),hv_im(:,4),actsp,MPIError)
					call dielbsev_dist_real(dimbse,hbse_dist_r,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hv_re(:,5),hv_im(:,5),actsm,MPIError)
				else
					actsp = 0.0
					actsm = 0.0
				end if
#endif
			end if
		else if (Nodes == 1) then
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
			if (nstates < dimbse) then
				write(300,*) 'BSE_NSTATES: bse_oscf.dat has the',nstates,'lowest excitons, up to',W(nstates),'eV;'
				write(300,*) 'the spectrum is complete only below that energy'
				write(300,*) 'the spectra add 10 x the mean level spacing of these states, ', &
					10.0*(W(nstates)-W(1))/real(nstates),' eV, to CSHIFT (avgsme in optics.F90)'
			end if
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

			do i=1,nstates
				write(301,"(7F15.6)") W(i),actxx(i),actyy(i),actzz(i),actxy(i),actxz(i),actyz(i)
				call flush(301)
				if (cpol) then
					write(302,"(6F15.6)") W(i),actxx(i),actyy(i),actzz(i),actsp(i),actsm(i)
					call flush(302)
				end if
			end do
		end if

		! --- Exciton-phonon coupling post-processing ---
		! MPI: q-points are distributed round-robin inside the EPH subroutines.
		! OMP: inner lam loops are threaded inside the EPH subroutines.
		! The eigenvector matrix is first replicated on all ranks (one MPI_ALLREDUCE)
		! so the existing EPH routines (which expect a full hbse) work unchanged.
		if (len_trim(bse_eph_file) > 0) then
			call eph_read(trim(bse_eph_file), eph, eph_ok, eph_message)
			if (.not. eph_ok) then
				if (Node == 0) write(*,*) trim(eph_message)
				stop
			end if
			allocate(kplusq_idx(ngkpt, eph%nq))
			call eph_build_kplusq(kpt_bse, ngkpt, eph%qpts, eph%nq, rlat, ktol, kplusq_idx)

			! Gather distributed eigenvectors to all ranks
			allocate(hbse_work(dimbse, dimbse))
			if (Nodes == 1) then
				hbse_work = hbse
			else
#ifdef MPI
				hbse_work = cmplx(0.0, 0.0)
				do lj = 1, locc
					jg = indxl2g(lj, nb, mycol, 0, npcol)
					do li = 1, locr
						ig = indxl2g(li, mb, myrow, 0, nprow)
						if (ig <= dimbse .and. jg <= dimbse) hbse_work(ig, jg) = hbse_dist(li, lj)
					end do
				end do
				call bse_mpi_allreduce_complex_sum_blocks(hbse_work(1,1), &
					int(dimbse,kind=8)*int(dimbse,kind=8), MPI_COMM_WORLD, MPIError, 16777216)
#endif
			end if

			! Exciton optical matrix elements mu(lam) = sum_j A^lam_j * hopt_j
			allocate(mu_x_eph(dimbse), mu_y_eph(dimbse), mu_z_eph(dimbse))
			do i = 1, dimbse
				mu_x_eph(i) = sum(hbse_work(:,i) * hrx)
				mu_y_eph(i) = sum(hbse_work(:,i) * hry)
				mu_z_eph(i) = sum(hbse_work(:,i) * hrz)
			end do

			! Linewidth: Gamma(lam) = pi * sum_{q,nu,lp} |gex(lam,lp)|^2 * Bose-Einstein weights
			allocate(Gamma_eph(dimbse))
			Gamma_eph = 0.0
			call exciton_phonon_linewidth(eph, vector, w90basis, nc, nv, ngkpt, &
			                              hbse_work, stt_bse, dimbse, W, kplusq_idx, &
			                              bse_eph_temp, bse_eph_eta, Gamma_eph, Node, Nodes)
			if (Node == 0) then
				open(newunit=io_eph, file=trim(outputfolder)//'exc_eph_linewidth.dat', status='replace', action='write')
				write(io_eph,'(A)') '# state   E_lambda(eV)   Gamma(meV)   lifetime(ps)'
				do i = 1, dimbse
					write(io_eph,'(I8,3F18.8)') i, W(i), Gamma_eph(i)*1000.0, &
					                             0.6582119569 / (Gamma_eph(i) + 1.0e-30)
				end do
				close(io_eph)
				write(300,*) 'EPH linewidth written to exc_eph_linewidth.dat'
				call flush(300)
			end if
			deallocate(Gamma_eph)

			! Phonon-assisted absorption spectrum
			allocate(omega_grid_eph(int(numbse)), chi2_xx_eph(int(numbse)), &
			         chi2_yy_eph(int(numbse)), chi2_zz_eph(int(numbse)))
			do iomega = 1, int(numbse)
				omega_grid_eph(iomega) = ebse0 + (ebsef - ebse0) * (iomega - 1) / max(int(numbse) - 1, 1)
			end do
			chi2_xx_eph = 0.0;  chi2_yy_eph = 0.0;  chi2_zz_eph = 0.0
			call exciton_phonon_phassist(eph, vector, w90basis, nc, nv, ngkpt, &
			                             hbse_work, stt_bse, dimbse, W, &
			                             mu_x_eph, mu_y_eph, mu_z_eph, kplusq_idx, &
			                             bse_eph_temp, bse_eph_eta, &
			                             omega_grid_eph, int(numbse), &
			                             chi2_xx_eph, chi2_yy_eph, chi2_zz_eph, Node, Nodes)
			if (Node == 0) then
				open(newunit=io_eph, file=trim(outputfolder)//'exc_eph_phassist.dat', status='replace', action='write')
				write(io_eph,'(A)') '# omega(eV)   Im_chi_xx   Im_chi_yy   Im_chi_zz'
				do iomega = 1, int(numbse)
					write(io_eph,'(4F18.8)') omega_grid_eph(iomega), &
					    chi2_xx_eph(iomega), chi2_yy_eph(iomega), chi2_zz_eph(iomega)
				end do
				close(io_eph)
				write(300,*) 'Phonon-assisted spectrum written to exc_eph_phassist.dat'
				call flush(300)
			end if
			deallocate(omega_grid_eph, chi2_xx_eph, chi2_yy_eph, chi2_zz_eph)

			deallocate(mu_x_eph, mu_y_eph, mu_z_eph)
			deallocate(hbse_work)
			deallocate(kplusq_idx)
			call eph_destroy(eph)
		end if

		deallocate(eigv,vector)
		deallocate(rvec,hopmatrices)
		deallocate(ihopmatrices,ffactor)
		deallocate(hrx,hry,hrz,hrsp,hrsm)
	deallocate(actxx,actxy,actxz,actyy,actyz,actzz)
	deallocate(actsp,actsm)
		if (use_real) then
			deallocate(hv_re,hv_im)
			call rb_free()
			if (allocated(rb_pmap)) deallocate(rb_pmap,rb_psgn)
			if (allocated(rb_tbar)) deallocate(rb_tbar,rb_sig)
			if (Nodes == 1) then
				deallocate(hbse_r)
			else
#ifdef MPI
				call BLACS_GRIDEXIT(blacs_ctxt)
				deallocate(hbse_dist_r)
#endif
			end if
		else if (Nodes == 1) then
			deallocate(hbse)
		else
#ifdef MPI
			call BLACS_GRIDEXIT(blacs_ctxt)
			deallocate(hbse_dist)
#endif
		end if
		deallocate(W,stt,stt_bse,nocpk,rb_kb)
	deallocate(tau)
	if (use_dmat) deallocate(dmat)
	deallocate(ovp)

	deallocate(kpt,kpt_bse,center_phase,wannier_centers)
	
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


	


contains

	subroutine rb_refuse(msg)
		! BSE_REAL= T cannot be used with this input: every rank stops
		character(len=*),intent(in) :: msg
		if (Node == 0) write(*,*) trim(msg)
#ifdef MPI
		call MPI_ABORT(MPI_COMM_WORLD,1,MPIError)
#endif
		stop 1
	end subroutine rb_refuse

	function rb_melem(i,j) result(h)
		! kernel element H(i,j) of the states i and j of the complex basis of the mesh
		integer,intent(in) :: i,j
		complex :: h
		integer :: isk,jsk
		isk = min(stt_bse(4,i),nsk)
		jsk = min(stt_bse(4,j),nsk)
		h = rb_melem_sk(i,j,sk(:,:,isk),sk(:,:,jsk))
	end function rb_melem

	function rb_melem_sk(i,j,sk1,sk2) result(h)
		! the same, with the overlap matrices S(k) of the two states given (the table sk holds one scratch matrix in the
		! fast kernel mode)
		integer,intent(in) :: i,j
		complex,dimension(w90basis,w90basis),intent(in) :: sk1,sk2
		complex :: h
		h = matrizelbse(coultype,ktol,w90basis,ediel,lc,ez,w1,r0,ngrid,rlat,stt_bse(:,i),&
		    eigv(stt_bse(4,i),stt_bse(3,i)),eigv(stt_bse(4,i),stt_bse(2,i)),vector(:,stt_bse(3,i),stt_bse(4,i)),&
		    vector(:,stt_bse(2,i),stt_bse(4,i)),kpt_bse(:,stt_bse(4,i)),stt_bse(:,j),eigv(stt_bse(4,j),stt_bse(3,j)),&
		    eigv(stt_bse(4,j),stt_bse(2,j)),vector(:,stt_bse(3,j),stt_bse(4,j)),vector(:,stt_bse(2,j),stt_bse(4,j)),&
		    kpt_bse(:,stt_bse(4,j)),dft,nvec,rvec,sk1,sk2,use_center_phase, &
		    center_phase(:,stt_bse(4,i)),center_phase(:,stt_bse(4,j)), &
		    use_center_grad,rmn,wannier_centers,ffdat)
	end function rb_melem_sk

	subroutine rb_block(ru,rc,b)
		! the real matrix elements between the real basis functions of the unit that starts at ru (rows) and the one
		! that starts at rc (columns): at most 2 x 2, from H(a,b) and, for two pairs, H(a,bbar) (bse_realbasis)
		integer,intent(in) :: ru,rc
		real,intent(out) :: b(2,2)
		complex :: h1,h2
		integer :: ku,kc,ii,jj,kki,kkj,nrow,ncol
		ku = rb_kind(ru)
		kc = rb_kind(rc)
		nrow = 1
		if (ku == 1) nrow = 2
		ncol = 1
		if (kc == 1) ncol = 2
		h1 = rb_melem(rb_row(ru),rb_row(rc))
		h2 = cmplx(0.0,0.0)
		if (ku /= 0 .and. kc /= 0) h2 = real(rb_sgn(rc))*rb_melem(rb_row(ru),rb_rowbar(rc))
		b = 0.0
		do jj=1,ncol
			kkj = 0
			if (kc == 1) kkj = jj
			do ii=1,nrow
				kki = 0
				if (ku == 1) kki = ii
				b(ii,jj) = rb_element(kki,kkj,h1,h2)
			end do
		end do
	end subroutine rb_block

	subroutine rb_kernel_check(dev)
		! time-reversal symmetry of the kernel, H(atilde,btilde) = conjg(H(a,b)) with atilde = sigma_a abar, for the
		! states of three functions of the real basis table (the first, a middle one, the last: a function of a
		! self-conjugate k point when there is one) in all combinations; dev is the largest
		! |H(atilde,btilde) - conjg(H(a,b))| / (1 + |H(a,b)|).  The overlap matrices S(k) of the six states are built here
		! (the fast kernel keeps one scratch matrix), so the check does not depend on the kernel mode.
		real,intent(out) :: dev
		integer :: fu(3),rows(6),fi,fj,i1,i2,m
		complex,allocatable :: skc(:,:,:)
		complex :: h1,h2
		fu = (/1,max(1,rb_n/2),rb_n/)
		do fi=1,3
			rows(2*fi-1) = rb_row(fu(fi))
			rows(2*fi) = rb_rowbar(fu(fi))
		end do
		allocate(skc(w90basis,w90basis,6))
		do m=1,6
			if (dft .eq. "S") then
				call overlap(w90basis,nvec,rvec,ovp,kpt(stt_bse(4,rows(m)),1),kpt(stt_bse(4,rows(m)),2),kpt(stt_bse(4,rows(m)),3),skc(:,:,m))
			else
				skc(:,:,m) = cmplx(0.0,0.0)
				do i1=1,w90basis
					skc(i1,i1,m) = cmplx(1.0,0.0)
				end do
			end if
		end do
		dev = 0.0
		do fi=1,3
			do fj=1,3
				i1 = 2*fi-1
				i2 = 2*fj-1
				h1 = rb_melem_sk(rows(i1),rows(i2),skc(:,:,i1),skc(:,:,i2))
				h2 = real(rb_sgn(fu(fi))*rb_sgn(fu(fj)))*rb_melem_sk(rows(i1+1),rows(i2+1),skc(:,:,i1+1),skc(:,:,i2+1))
				dev = max(dev,abs(h2-conjg(h1))/(1.0+abs(h1)))
			end do
		end do
		deallocate(skc)
	end subroutine rb_kernel_check


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
