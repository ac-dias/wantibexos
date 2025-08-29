module prec

	implicit none
	integer,parameter :: qprec= 8
	integer,parameter :: cqprec=4

end module prec

module hamiltonian_input_variables

	implicit none
	
	integer :: w90basis,ntype,nocp
	real :: efermi,scs
	!real:: edielh
	integer :: nvec
	real,dimension(3,3) :: rlat

	integer,allocatable,dimension(:,:) :: rveci
	real,allocatable,dimension(:,:) :: rvec
	real,allocatable,dimension(:,:,:) :: hopmatrices,ihopmatrices
	real,allocatable,dimension(:,:,:) :: ovp
	integer,allocatable,dimension(:) :: ffactor

	character(len=4) :: systype

	!real :: lc
	
	!variables for spin polarized hamiltonian w90
	real,allocatable,dimension(:,:,:) :: hopmatricesu,ihopmatricesu
	real,allocatable,dimension(:,:,:) :: hopmatricesd,ihopmatricesd
	integer,allocatable,dimension(:) :: ffactoru,ffactord
	integer,allocatable,dimension(:,:) :: rvecu,rvecd
	integer :: w90basisu,w90basisd
	integer :: nvecu,nvecd	
					

	

end module hamiltonian_input_variables

subroutine hamiltonian_nort_input_read(unidade,inputfile)

	use hamiltonian_input_variables
	implicit none
	integer :: i,j,k
	integer :: unidade,erro
	integer :: flagaux
	character(len=70) :: inputfile,charflag
	character(len=1) :: cflag

	OPEN(UNIT=unidade, FILE= inputfile,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening non-ortogonal hr input file"
    	
    	read(unidade,*) systype
	read(unidade,*) scs
	read(unidade,*) efermi 
	
	read(unidade,*) rlat(1,1),rlat(1,2),rlat(1,3)
	read(unidade,*) rlat(2,1),rlat(2,2),rlat(2,3)
	read(unidade,*) rlat(3,1),rlat(3,2),rlat(3,3)	
	
	read(unidade,*) w90basis
	read(unidade,*) nvec	
	
	read(unidade,*) cflag,cflag,cflag,cflag,cflag,cflag,cflag,cflag,cflag,cflag,cflag
	
	allocate(ffactor(nvec))	
	
	allocate(rvec(nvec,3),hopmatrices(nvec,w90basis,w90basis),ihopmatrices(nvec,w90basis,w90basis))
	allocate(ovp(nvec,w90basis,w90basis))
	
	do i=1,nvec

		do j=1,w90basis

			do k=1,w90basis

			read(unidade,*) rvec(i,1),rvec(i,2),rvec(i,3),flagaux,flagaux,hopmatrices(i,k,j),ihopmatrices(i,k,j),ovp(i,k,j)
			
			end do
		end do

	end do	
	
	ffactor = 1  	
    	
end subroutine hamiltonian_nort_input_read

subroutine hamiltonian_input_read(unidade,inputfile)
	
	use hamiltonian_input_variables
	implicit none
	integer :: i,j,k
	integer :: unidade,erro
	integer :: flagaux
	character(len=70) :: inputfile,charflag

	OPEN(UNIT=unidade, FILE= inputfile,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wannier90 hr input file"
	
	read(unidade,*) systype
	!read(unidade,*) edielh
	read(unidade,*) scs
	read(unidade,*) efermi


	read(unidade,*) rlat(1,1),rlat(1,2),rlat(1,3)
	read(unidade,*) rlat(2,1),rlat(2,2),rlat(2,3)
	read(unidade,*) rlat(3,1),rlat(3,2),rlat(3,3)
	
	select case (systype)
	
	case("SP")
	
	!lendo parte spin up
	
	read(unidade,*) charflag

	read(unidade,*) w90basisu
	read(unidade,*) nvecu

	allocate(ffactoru(nvecu))

	read(unidade,*)(ffactoru(i), i=1, nvecu)


	allocate(rvecu(nvecu,3),hopmatricesu(nvecu,w90basisu,w90basisu),ihopmatricesu(nvecu,w90basisu,w90basisu))
	

	do i=1,nvecu

		do j=1,w90basisu

			do k=1,w90basisu

			read(unidade,*) rvecu(i,1),rvecu(i,2),rvecu(i,3),flagaux,flagaux,hopmatricesu(i,k,j),ihopmatricesu(i,k,j)
			
			end do
		end do

	end do	
	
	!lendo parte spin down	
	
	read(unidade,*) charflag

	read(unidade,*) w90basisd
	read(unidade,*) nvecd

	allocate(ffactord(nvecd))

	read(unidade,*)(ffactord(i), i=1, nvecd)



	allocate(rvecd(nvecd,3),hopmatricesd(nvecd,w90basisd,w90basisd),ihopmatricesd(nvecd,w90basisd,w90basisd))
	

	do i=1,nvecd

		do j=1,w90basisd

			do k=1,w90basisd

			read(unidade,*) rvecd(i,1),rvecd(i,2),rvecd(i,3),flagaux,flagaux,hopmatricesd(i,k,j),ihopmatricesd(i,k,j)
			
			end do
		end do

	end do	
	
	!juntando as partes up e down
	
	w90basis = w90basisd+w90basisu
	nvec = nvecd+nvecu
	
	allocate(ffactor(nvec))	
	allocate(rveci(nvec,3))
	allocate(rvec(nvec,3),hopmatrices(nvec,w90basis,w90basis),ihopmatrices(nvec,w90basis,w90basis))	
	
	ffactor(1:nvecu) = ffactoru
	ffactor(nvecu+1:nvec) = ffactord
	
	rveci(1:nvecu,:) = rvecu
	rveci(nvecu+1:nvec,:) = rvecd	
	
	hopmatrices = 0.0
	ihopmatrices = 0.0 
	
	hopmatrices(1:nvecu,1:w90basisu,1:w90basisu) = hopmatricesu
	ihopmatrices(1:nvecu,1:w90basisu,1:w90basisu) = ihopmatricesu
	
	hopmatrices(nvecu+1:nvec,w90basisu+1:w90basis,w90basisu+1:w90basis) = hopmatricesd	
	ihopmatrices(nvecu+1:nvec,w90basisu+1:w90basis,w90basisu+1:w90basis) = ihopmatricesd				
	
	deallocate(ffactord)
	deallocate(rvecd,hopmatricesd,ihopmatricesd)	
	deallocate(ffactoru)
	deallocate(rvecu,hopmatricesu,ihopmatricesu)		
	
	case default


	read(unidade,*) charflag

	read(unidade,*) w90basis
	read(unidade,*) nvec

	allocate(ffactor(nvec))

	read(unidade,*)(ffactor(i), i=1, nvec)


	allocate(rveci(nvec,3))
	allocate(rvec(nvec,3),hopmatrices(nvec,w90basis,w90basis),ihopmatrices(nvec,w90basis,w90basis))
	

	do i=1,nvec

		do j=1,w90basis

			do k=1,w90basis

			read(unidade,*) rveci(i,1),rveci(i,2),rveci(i,3),flagaux,flagaux,hopmatrices(i,k,j),ihopmatrices(i,k,j)
			
			end do
		end do

	end do	

	end select
	
	
	do i=1,nvec
	
		rvec(i,1) = rveci(i,1)*rlat(1,1)+rveci(i,2)*rlat(2,1)+rveci(i,3)*rlat(3,1)
		rvec(i,2) = rveci(i,1)*rlat(1,2)+rveci(i,2)*rlat(2,2)+rveci(i,3)*rlat(3,2)
		rvec(i,3) = rveci(i,1)*rlat(1,3)+rveci(i,2)*rlat(2,3)+rveci(i,3)*rlat(3,3)
	
	end do
	
	deallocate(rveci)

end subroutine hamiltonian_input_read


module input_variables

	implicit none 

	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv
	!integer :: ncrpa,nvrpa
	!integer :: ncbz,nvbz
	real :: rk
	real :: numdos
	real :: ebse0,ebsef,numbse
	real :: sme,cshift,smeboltz
	real,dimension(3) :: mshift, ediel
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
	character(len=5) :: sysdim
	character(len=1) :: dft
	character(len=2) :: ta
	character(len=6) :: ses
	character(len=70) :: emfile
		

	logical :: bandscalc,doscalc
	logical :: bse,bsepol,bsekpath
	logical :: spdiel,spdielpol,sppolbz
	logical :: spec
	logical :: berryk,berrybz
	logical :: pponly
	logical :: bsewf
	logical :: tmcoef
	logical :: dtfull,cpol
	!logical :: bset,bsetbnd
	logical :: pce 
	logical :: renorm
	logical :: emt
	logical :: spintxt
	logical :: boltz
	
	integer :: excwf0,excwff
	
	real :: ez,w,lc,r0
	real :: st,phavg,temp
	
	real :: ctemp,tmax,eg,egd
	real :: egs,ebgs
	real :: dk
	
	integer :: nocpf
	real :: fermishift
	character(len=12) :: bsealgo
	
	real :: exc
	real,dimension(3) :: mag
	
	integer :: power
	logical :: lowdin
	
	integer :: nmu
	integer :: nsteps
	real :: mu0,muf
	real :: btemp,klat
	real,dimension(3) :: elft,hlft


	

end module input_variables

subroutine input_read
	
	use input_variables
	use hamiltonian_input_variables
	implicit none

	character(len=70) :: a,b
	character(len=70) :: add
	integer :: erro


	!default values

	nthreads = 1
	
	nocpf = 0
	fermishift = 0.0

	ngrid(1)= 1
	ngrid(2)= 1
	ngrid(3)= 1

	rk= 0.0

	nc = 1
	nv = 1

	!edos0 = -5.0
	!edosf = 5.0
	numdos = 6001

	ebse0 = 0.0
	ebsef = 3.0
	numbse = 6001

	sme = 0.08
	cshift = 0.01
	smeboltz = 0.01

	exc = 0.000
	mag(1) = 0.000
	mag(2) = 0.000
	mag(3) = 1.000

	mshift(1) = 0.0
	mshift(2) = 0.0
	mshift(3) = 0.0

	ktol = 0.001 


	outputfolder = "./"

	calcparms = "./"

	meshtype = "MKH"
	
	meshgen = "MKH"

	coultype = "V3D"

	orbw = "unknown"

	kpaths = "unknown"

	kpathsbse = "unknown"  
	
	emfile = "unknown"

	bandscalc = .false.

	doscalc = .false.

	bse = .false.

	bsepol = .false.

	bsekpath = .false.

	spdiel = .false.

	spdielpol = .false.

	sppolbz = .false.

	spec = .false.

	berryk = .false.

	berrybz = .false.
	
	pponly = .false.
	
	bsewf = .false.
	
	tmcoef = .false.
	
	dtfull = .false.
	
	cpol = .false.
	

	
	excwf0 = 1
	excwff = 2
	


	sysdim = "3D"	
	dft = "W"
	ta = "FA"
	
	bsealgo = "cheev"
	!bset = .false.
	!bsetbnd = .false.
	renorm = .true.
	
	ediel(1) = 1.0
	ediel(2) = 1.0
	ediel(3) = 1.0
	ez = 1.0
	w = 0.0
	lc = 1.0
	r0 = 1.0	
	
	st = 0.0 
	phavg = 0.0 
	temp = 0.0
	
	
	pce = .false.
	ses = "AM15G"
	ctemp = 298.15 
	tmax = 5E-06
	eg = 0.00
	egd = 0.00
	egs = 0.00
	ebgs = 0.00
	
	emt = .false.
	dk = 0.14
	
	spintxt = .false.
	
	boltz = .false.
	nmu = 601
	nsteps = 6001
	mu0 = -2.0
	muf = 2.0
	btemp = 300
	klat = 0.0
	
	elft = 1.0
	hlft = 1.0
	
	!lowdin = .false.
	!power = 5

	!end default values

	do 

	read(*,*,iostat=erro) a,b
	if (erro/=0) exit

	select case (a)

	case ("MAG_X=")

		read(b,*) mag(1)
		
	case ("MAG_Y=")

		read(b,*) mag(2)
		
	case ("MAG_Z=")

		read(b,*) mag(3)					
	
	case ("SPIN_EXC=")

		read(b,*) exc
	
	case ("KMESH_GEN=")

		read(b,*) meshgen	
	
	case ("SPIN_TXT=")

		read(b,*) spintxt
	
	case ("FERMI_SHIFT=")

		read(b,*) fermishift	
	
	case ("NOCPF=")

		read(b,*) nocpf	
	
	case ("BSE_ALGO=")

		read(b,*) bsealgo
	
	case ("dK=")

		read(b,*) dk
	
	case ("EM_TENSOR=")

		read(b,*) emt 	
	
	case ("EM_FILE=")

		read(b,*) emfile 	
	
	case ("RNMD=")

		read(b,*) renorm 	
	
	case ("PCE=")

		read(b,*) pce 	
		
	case ("SES=")

		read(b,*) ses
			
	case ("CTEMP=")

		read(b,*) ctemp
		
	case ("THMAX=")

		read(b,*) tmax	
		
	case ("EG=")

		read(b,*) eg	
		
	case ("EGD=")

		read(b,*) egd	
		
	case ("EGS=")

		read(b,*) egs	
		
	case ("EBGS=")

		read(b,*) ebgs														
	
	case ("TA=")

		read(b,*) ta 	
	
	!case ("BSET=")

	!	read(b,*) bset 	
	
	!case ("BSET_BND=")

	!	read(b,*) bsetbnd 	
	
	case ("TEMP=")

		read(b,*) temp 	
	
	case ("PHAVG=")

		read(b,*) phavg 	
	
	case ("ST=")

		read(b,*) st 	
	
	case ("DFT=")

		read(b,*) dft 	
	
	case ("CSHIFT=")

		read(b,*) cshift 	
	
	case ("DTDIAG=")

		read(b,*) dtfull 
		
	case ("CPOL=")

		read(b,*) cpol 			
	
	case ("SYSDIM=")

		read(b,*) sysdim 	
	
	case ("R_0=")
	
		read(b,*) r0	
	
	case ("LC=")
	
		read(b,*) lc	
	
	case ("EDIEL_Z=")
	
		read(b,*) ez
	
	case ("W_COUL=")
	
		read(b,*) w
	
	case ("TMCOEF=")
	
		read(b,*) tmcoef
	
	case ("EXC_WF_F=")
	
		read(b,*) excwff

	case ("EXC_WF_I=")
	
		read(b,*) excwf0
		
	case ("BSE_WF=")
	
		read(b,*) bsewf					
	
	case ("PP_ONLY=")
	
		read(b,*) pponly

	case ("BERRY_BZ=")

		read(b,*) berrybz

	case ("BERRY=")

		read(b,*) berryk

	case ("SPEC=")

		read(b,*) spec


	case ("OPT_BZ=")

		read(b,*) sppolbz

	!case ("DIEL_POL=")

	!	read(b,*) spdielpol

	case ("IPA=")

		read(b,*) spdiel

	case ("BSE_BND=")

		read(b,*) bsekpath

	!case ("BSE_POL=")

	!	read(b,*) bsepol

	case ("BSE=")

		read(b,*) bse

	case ("DOS=")

		read(b,*) doscalc

	case ("BANDS=")

		read(b,*) bandscalc

	case ("RK=")

		read(b,*) rk

	case ("COULOMB_POT=")

		coultype = b

	case ("MESH_TYPE=")

		meshtype = b

	case ("OUTPUT=")
		
		outputfolder = b

	case ("CALC_DATA=")
		
		calcparms = b

	case ("NGX=")

		read(b,*) ngrid(1) 

	case ("NGY=")

		read(b,*) ngrid(2) 

	case ("NGZ=")

		read(b,*) ngrid(3) 

	case ("NBANDSC=")

		read(b,*) nc 

	case ("NBANDSV=")

		read(b,*) nv 


	case ('NEDOS=')

		read(b,*) numdos 

	case ("ENSPECI=")

		read(b,*) ebse0 

	case ("ENSPECF=")

		read(b,*) ebsef 

	case ("NESPEC=")

		read(b,*) numbse 

	case ("SIGMA=")

		read(b,*) sme 		

	case ("KTOL=")

		read(b,*) ktol 

	case ("PARAMS_FILE=")

		params = b

	case ("KPATH_FILE=")

		kpaths = b

	case ("KPATH_BSE=")

		kpathsbse = b

	case ("ORB_W=")

		orbw = b

	case ("EDIEL_T=")

		read(b,*) ediel(1) 

	case ("EDIEL=")

		read(b,*) ediel(2) 		

	case ("EDIEL_B=")

		read(b,*) ediel(3) 


	case ("NTHREADS=")

		read(b,*) nthreads 

	case ("SHIFT_1=")

		read(b,*) mshift(1) 

	case ("SHIFT_2=")

		read(b,*) mshift(2) 

	case ("SHIFT_3=")

		read(b,*) mshift(3)
		
	case ("BOLTZ=")

		read(b,*) boltz 
		
	case ("BTEMP=")

		read(b,*) btemp
		
	case ("NEBOLTZ=")

		read(b,*) nsteps
		
	case ("NEMU=")

		read(b,*) nmu
		
	case ("MU_I=")

		read(b,*) mu0	
		
	case ("MU_F=")

		read(b,*) muf
	
	case ("KappaL=")
		
		read(b,*) klat
															 
	case ("ELFT_X=")
	
		read(b,*) elft(1)

	case ("ELFT_Y=")
	
		read(b,*) elft(2)

	case ("ELFT_Z=")
	
		read(b,*) elft(3)				
	
	case("HLFT_X=")
	
		read(b,*) hlft(1)
		
	case("HLFT_Y=")
	
		read(b,*) hlft(2)
		
	case("HLFT_Z=")
	
		read(b,*) hlft(3)
		
	case("SIGMA_BOLTZ=")
	
		read(b,*) smeboltz						

	case default
	 continue

	end select


	end do




end subroutine input_read

subroutine param_out(unitout,nthreads,outputfolder,calcparms,ngrid,nc,nv,numdos, &
		     ebse0,ebsef,numbse,sme,ktol,params,kpaths,kpathsbse,orbw,ediel, &
		     mshift,coultype,bandscalc,doscalc,bse,bsepol,bsekpath,spec,&
		     spdiel,spdielpol,sppolbz,berryk,berrybz,pponly,bsewf,excwf0,excwff,&
		     tmcoef,ez,w,lc,r0,sysdim,dtfull,cpol,cshift,dft,&
		     st,phavg,temp,ta,pce,ses,ctemp,tmax,eg,egd,egs,ebgs,renorm,emt,emfile,&
		     dk,nocpf,fermishift,bsealgo,spintxt,meshgen,exc,mag,boltz,nmu,nsteps,mu0,muf,btemp,klat,elft,hlft,smeboltz)

	implicit none
	integer :: unitout

	real,dimension(3) :: ediel

	integer :: nthreads
	integer,dimension(3) :: ngrid
	integer :: nc,nv,excwf0,excwff
	integer :: nocpf
	!integer :: power
	!integer :: ncrpa,nvrpa
	!integer :: ncbz,nvbz
	real :: numdos
	real :: ebse0,ebsef,numbse
	real :: sme,cshift,smeboltz
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
	character(len=5) :: sysdim
	character(len=1) :: dft
	character(len=2) :: ta	
	character(len=6) :: ses
	character(len=70) :: emfile
	character(len=12) :: bsealgo		

	logical :: bandscalc,doscalc
	logical :: bse,bsepol,bsekpath
	logical :: spdiel,spdielpol,sppolbz
	logical :: spec,tmcoef
	logical :: berryk,berrybz,pponly,bsewf
	logical :: dtfull,cpol	
	!logical :: bset,bsetbnd
	logical :: pce,renorm
	logical :: emt,spintxt	
	!logical :: lowdin
	logical :: boltz
	
	real :: fermishift
	real :: ez,w,lc,r0
	real :: st,phavg,temp
	
	real :: ctemp,tmax,eg,egd
	real :: egs,ebgs
	real :: dk	
	
	real :: exc
	real,dimension(3) :: mag,elft,hlft
	
	integer :: nmu
	integer :: nsteps
	real :: mu0,muf
	real :: btemp,klat			

	write(unitout,"(A10,I0)") "NTHREADS= ",nthreads
	write(unitout,"(A8,A2)") "SYSDIM= ",sysdim
	write(unitout,"(A5,A1)") "DFT= ",dft		
	write(unitout,*)
	write(unitout,*) "INPUT/OUTPUT FILES"
	write(unitout,*)
	write(unitout,"(A8,A70)") "OUTPUT= ", outputfolder
	write(unitout,"(A11,A70)") "CALC_DATA= ", calcparms
	write(unitout,"(A13,A70)") "PARAMS_FILE= ", params
	write(unitout,"(A12,A70)") "KPATH_FILE= ", kpaths
	write(unitout,"(A11,A70)") "KPATH_BSE= ", kpathsbse
	write(unitout,"(A7,A70)") "ORB_W= ", orbw
	write(unitout,"(A9,A70)") "EM_FILE= ", emfile	
	write(unitout,*)	
	write(unitout,*) "CONTROL"
	write(unitout,*)
	write(unitout,"(A7,L1)")  "BANDS= ",bandscalc
	write(unitout,"(A5,L1)")  "DOS= ",doscalc	
	write(unitout,"(A5,L1)")  "BSE= ",bse
	!write(unitout,"(A6,L1)")  "BSET= ",bset	
	write(unitout,"(A9,L1)")  "BSE_BND= ",bsekpath
	!write(unitout,"(A10,L1)")  "BSET_BND= ",bsetbnd	
	write(unitout,"(A5,L1)")  "IPA= ",spdiel
	write(unitout,"(A8,L1)")  "OPT_BZ= ",sppolbz
	write(unitout,"(A6,L1)")  "SPEC= ",spec
	write(unitout,"(A10,L1)") "BERRY_BZ= ",berrybz
	write(unitout,"(A7,L1)")  "BERRY= ",berryk
	write(unitout,"(A9,L1)")  "PP_ONLY= ",pponly
	write(unitout,"(A8,L1)") "BSE_WF= ",bsewf	
	write(unitout,"(A8,L1)") "TMCOEF= ",tmcoef
	write(unitout,"(A8,L1)") "DTDIAG= ",dtfull
	write(unitout,"(A6,L1)") "CPOL= ",cpol
	!write(unitout,"(A5,L1)") "PCE= ",pce
	write(unitout,"(A11,L1)") "EM_TENSOR= ",emt
	write(unitout,"(A10,L1)") "SPIN_TXT= ",spintxt
	write(unitout,"(A7,L1)") "BOLTZ= ",boltz					
	write(unitout,*)
	write(unitout,*) "K-MESH"
	write(unitout,*)
	write(unitout,"(A11,A15)") "KMESH_GEN= ", meshgen
	write(unitout,"(A5,I0)") "NGX= ", ngrid(1)
	write(unitout,"(A5,I0)") "NGY= ", ngrid(2)
	write(unitout,"(A5,I0)") "NGZ= ", ngrid(3)
	write(unitout,"(A9,F8.4)") "SHIFT_1= ", mshift(1)
	write(unitout,"(A9,F8.4)") "SHIFT_2= ", mshift(2)
	write(unitout,"(A9,F8.4)") "SHIFT_3= ", mshift(3)
	write(unitout,*)
	write(unitout,*) "DOS"
	write(unitout,*)
	!write(unitout,"(A8,F8.4)") "ENDOSI= ", edos0
	!write(unitout,"(A8,F8.4)") "ENDOSF= ", edosf
	write(unitout,"(A7,I0)") "NEDOS= ", int(numdos)
	write(unitout,"(A7,F8.4)") "SIGMA= ", sme
	write(unitout,*)
	write(unitout,*) "BSE/OPTICAL PROPERTIES"
	write(unitout,*)
	write(unitout,"(A10,A12)") "BSE_ALGO= ", bsealgo
	write(unitout,"(A9,I0)") "NBANDSC= ", nc
	write(unitout,"(A9,I0)") "NBANDSV= ", nv
	write(unitout,"(A9,F8.4)") "ENSPECI= ", ebse0
	write(unitout,"(A9,F8.4)") "ENSPECF= ", ebsef
	write(unitout,"(A8,I0)") "NESPEC= ", int(numbse)
	write(unitout,"(A6,F8.4)") "KTOL= ", ktol
	write(unitout,"(A10,I0)") "EXC_WF_I= ",excwf0
	write(unitout,"(A10,I0)") "EXC_WF_F= ",excwff
	write(unitout,"(A13,A5)") "COULOMB_POT= ",coultype
	write(unitout,"(A8,F8.4)") "CSHIFT= ", cshift	
	write(unitout,"(A6,L1)") "RNMD= ",renorm		
	write(unitout,*)
	write(unitout,*) "PARAMETERS FOR COULOMB POTENTIALS"
	write(unitout,*)
	write(unitout,"(A9,F8.4)") "EDIEL_Z= ",ez
	write(unitout,"(A8,F8.4)") "W_COUL= ",w
	write(unitout,"(A4,F8.4)") "LC= ",lc
	write(unitout,"(A5,F8.4)") "R_0= ",r0
	write(unitout,"(A9,F8.4)") "EDIEL_T= ", ediel(1)
	write(unitout,"(A9,F8.4)") "EDIEL_B= ", ediel(3)
	write(unitout,"(A7,F8.4)") "EDIEL= ", ediel(2)	
	write(unitout,*)	
	write(unitout,*) "PARAMETERS FOR BSE with Temperature Effects"
	write(unitout,*)
	write(unitout,"(A4,A2)") "TA= ",ta			
	write(unitout,"(A4,F8.4)") "ST= ",st
	write(unitout,"(A7,F8.4)") "PHAVG= ",phavg
	write(unitout,"(A6,F8.4)") "TEMP= ",temp
	write(unitout,*)	
	write(unitout,*) "PARAMETERS FOR PCE"
	write(unitout,*)
	write(unitout,"(A5,A6)") "SES= ",ses
	write(unitout,"(A7,F8.4)") "CTEMP= ",ctemp
	write(unitout,"(A7,E15.4)") "THMAX= ",tmax				
	write(unitout,"(A4,F8.4)") "EG= ",eg
	write(unitout,"(A5,F8.4)") "EGD= ",egd
	write(unitout,"(A5,F8.4)") "EGS= ",egs
	write(unitout,"(A6,F8.4)") "EBGS= ",ebgs
	write(unitout,*)
	write(unitout,*) "Effective Mass"
	write(unitout,*)
	write(unitout,"(A4,F8.4)") "dK= ",dk					
	write(unitout,*)
	write(unitout,*) "OTHERS"
	write(unitout,*)
	write(unitout,"(A13,F8.4)") "FERMI_SHIFT= ",fermishift
	write(unitout,"(A7,I0)") "NOCPF= ", nocpf	
	write(unitout,*)
	write(unitout,*) "PARAMETERS FOR BOLTZMANN TRANSPORT"
	write(unitout,*)
	write(unitout,"(A7,F8.4)") "BTEMP= ",btemp
	write(unitout,"(A9,I0)") "NEBOLTZ= ",nsteps
	write(unitout,"(A6,I0)") "NEMU= ",nmu
	write(unitout,"(A6,F8.4)") "MU_I= ",mu0
	write(unitout,"(A6,F8.4)") "MU_F= ",muf
	write(unitout,"(A8,F8.4)") "KappaL= ",klat
	write(unitout,"(A8,E15.4)") "ELFT_X= ",elft(1)
	write(unitout,"(A8,E15.4)") "ELFT_Y= ",elft(2)
	write(unitout,"(A8,E15.4)") "ELFT_Z= ",elft(3)		
	write(unitout,"(A8,E15.4)") "HLFT_X= ",hlft(1)
	write(unitout,"(A8,E15.4)") "HLFT_Y= ",hlft(2)
	write(unitout,"(A8,E15.4)") "HLFT_Z= ",hlft(3)
	write(unitout,"(A13,E15.4)") "SIGMA_BOLTZ= ",smeboltz				
	write(unitout,*)								
	write(unitout,*) "EXTERNAL HAMILTONIANS"
	write(unitout,*)
	write(unitout,"(A10,F8.4)") "SPIN_EXC= ", exc
	write(unitout,"(A7,F8.4)") "MAG_X= ", mag(1)
	write(unitout,"(A7,F8.4)") "MAG_Y= ", mag(2)
	write(unitout,"(A7,F8.4)") "MAG_Z= ", mag(3)	
	write(unitout,*)	
	!write(unitout,*) "LOWDIN ORBITALS"
	!write(unitout,*)
	!write(unitout,"(A8,L1)") "LOWDIN= ", lowdin
	!write(unitout,"(A7,I0)") "POWER= ", power
			
	!write(unitout,*)

	write(unitout,*)
    call flush(unitout)

end subroutine param_out


