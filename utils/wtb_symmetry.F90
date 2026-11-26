program wtb_symmetry

	! Symmetry group and operations of the system of a BeXos input file.
	!
	!    bin/wtb_symmetry.x < input.dat
	!
	! It reads the input file as wtb.x does (the tight-binding file PARAMS_FILE=, DFT=, the k mesh NGX= NGY= NGZ= and
	! SHIFT_1..3=, BSE_CENTER_FILE= for DFT=W) and prints the point group, the operations, which of them map the k mesh
	! onto itself, and how well the Hamiltonian obeys them (src/subroutines/bse_symmetry.F90).  For DFT=S the orbital
	! centres, species and shells come from the basis_set-* file of the run directory; for DFT=W the Wannier centres of
	! BSE_CENTER_FILE are the orbital centres (without them only the symmetry of the lattice is given).
	! The tolerance on the positions is 3e-3 Angstrom, WTB_SYM_TOL= (environment variable) changes it; orbitals whose
	! centres are closer than WTB_SYM_CLUSTER (Angstrom, default the tolerance) form one site, for Wannier functions of an atom.

	use input_variables
	use hamiltonian_input_variables
	use bse_symmetry
	use bse_q_optics, only: rmn_data, rmn_read, rmn_destroy, rmn_centers

	implicit none

	type(rmn_data) :: rmn
	type(sym_group) :: g
	real,allocatable,dimension(:,:) :: tau,centers
	real(kind=8),allocatable,dimension(:,:) :: tau8
	character(len=8),allocatable,dimension(:) :: orbspec
	integer,allocatable,dimension(:) :: orblq
	logical :: have_pos,have_lab,rmn_ok
	character(len=256) :: rmn_message,tolstring
	real(kind=8) :: tol,rclu,lat(3,3)
	integer :: i,ierr

	call input_read

	select case (dft)
	case ("S")
		call hamiltonian_nort_input_read(200,params)
	case default
		call hamiltonian_input_read(200,params)
		allocate(ovp(nvec,w90basis,w90basis))
		ovp = 0.0
	end select

	tol = 3.0d-3
	call get_environment_variable('WTB_SYM_TOL',tolstring)
	if (len_trim(tolstring) > 0) read(tolstring,*,iostat=ierr) tol
	rclu = tol
	tolstring = ''
	call get_environment_variable('WTB_SYM_CLUSTER',tolstring)
	if (len_trim(tolstring) > 0) read(tolstring,*,iostat=ierr) rclu

	allocate(tau(w90basis,3),tau8(w90basis,3),orbspec(w90basis),orblq(w90basis))
	tau = 0.0
	orbspec = ''
	orblq = -1
	have_pos = .false.
	have_lab = .false.
	if (dft == "S") then
		call orbital_centres(w90basis,tau)
		have_pos = any(abs(tau) > 0.0)
		call orbital_labels(w90basis,orbspec,orblq,have_lab)
	else if (len_trim(bsecenterfile) > 0) then
		call rmn_read(trim(bsecenterfile),rmn,rmn_ok,rmn_message)
		if (.not. rmn_ok) then
			write(*,*) trim(rmn_message)
			stop 1
		end if
		allocate(centers(3,w90basis))
		call rmn_centers(rmn,centers,rmn_ok,rmn_message)
		if (.not. rmn_ok) then
			write(*,*) trim(rmn_message)
			stop 1
		end if
		do i=1,w90basis
			tau(i,:) = centers(:,i)
		end do
		have_pos = .true.
		call rmn_destroy(rmn)
	end if
	tau8 = dble(tau)
	lat = dble(rlat)

	call sym_analyze(6,lat,w90basis,tau8,have_pos,orbspec,orblq,have_lab,nvec,rvec,hopmatrices,ihopmatrices,ovp,.true., &
		ngrid,dble(mshift),tol,rclu,g)

end program wtb_symmetry
