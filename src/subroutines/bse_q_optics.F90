module bse_q_optics

	implicit none
	private

	type,public :: rmn_data
		integer :: num_wann=0
		integer :: nrpts=0
		integer,allocatable :: rvec(:,:)
		complex,allocatable :: rmn(:,:,:,:)
	end type rmn_data

	public :: rmn_read
	public :: rmn_destroy
	public :: rmn_bloch
	public :: rmn_centers
	public :: center_phase_build
	public :: center_phase_direct_vertices
	public :: rmn_apply_q0_optical_correction

contains

subroutine rmn_destroy(data)

	type(rmn_data),intent(inout) :: data

	if (allocated(data%rvec)) deallocate(data%rvec)
	if (allocated(data%rmn)) deallocate(data%rmn)
	data%num_wann=0
	data%nrpts=0

end subroutine rmn_destroy


subroutine rmn_read(filename,data,ok,message)

	! Read Wannier90 seedname_r.dat. The file stores
	! r_mn(R)=<m,0|r|n,R>, with the x, y, and z complex components on each row.
	! We intentionally retain the complete Wigner-Seitz table exactly as written:
	! unlike seedname_hr.dat, seedname_r.dat has no degeneracy block to divide by.
	! W90's use_ws_distance representation needs seedname_wsvec.dat-dependent
	! corrections. Those are deliberately not approximated by this first path.

	character(len=*),intent(in) :: filename
	type(rmn_data),intent(inout) :: data
	logical,intent(out) :: ok
	character(len=*),intent(out) :: message

	integer :: unit,ios,ir,m,n,mfile,nfile
	integer :: rfile(3)
	logical,allocatable :: seen(:,:)
	logical :: wsvec_exists
	character(len=256) :: header
	character(len=512) :: wsvecfile
	real :: xre,xim,yre,yim,zre,zim

	ok=.false.
	message=''
	call rmn_destroy(data)

	wsvec_exists=.false.
	if (len_trim(filename) > 6 .and. filename(len_trim(filename)-5:len_trim(filename)) == '_r.dat') then
		wsvecfile=filename(:len_trim(filename)-6)//'_wsvec.dat'
		inquire(file=trim(wsvecfile),exist=wsvec_exists)
	end if
	if (wsvec_exists) then
		message='Wannier90 use_ws_distance data require wsvec support, which is not implemented'
		return
	end if

	open(newunit=unit,file=trim(filename),status='old',action='read',iostat=ios)
	if (ios /= 0) then
		message='Unable to open Wannier90 r-matrix file: '//trim(filename)
		return
	end if

	read(unit,'(A)',iostat=ios) header
	if (ios == 0) read(unit,*,iostat=ios) data%num_wann
	if (ios == 0) read(unit,*,iostat=ios) data%nrpts
	if (ios /= 0 .or. data%num_wann <= 0 .or. data%nrpts <= 0) then
		message='Invalid Wannier90 r-matrix header'
		close(unit)
		call rmn_destroy(data)
		return
	end if

	allocate(data%rvec(3,data%nrpts),data%rmn(3,data%num_wann,data%num_wann,data%nrpts), &
		seen(data%num_wann,data%num_wann),stat=ios)
	if (ios /= 0) then
		message='Unable to allocate Wannier90 r-matrix table'
		close(unit)
		call rmn_destroy(data)
		return
	end if
	data%rmn=cmplx(0.0,0.0)

	do ir=1,data%nrpts
		seen=.false.
		do m=1,data%num_wann
			do n=1,data%num_wann
				read(unit,*,iostat=ios) rfile(1),rfile(2),rfile(3),mfile,nfile, &
					xre,xim,yre,yim,zre,zim
				if (ios /= 0) then
					message='Unexpected or malformed record in Wannier90 r-matrix file'
					close(unit)
					deallocate(seen)
					call rmn_destroy(data)
					return
				end if
				if (mfile < 1 .or. mfile > data%num_wann .or. &
					nfile < 1 .or. nfile > data%num_wann) then
					message='Wannier90 r-matrix record has an out-of-range Wannier index'
					close(unit)
					deallocate(seen)
					call rmn_destroy(data)
					return
				end if
				if (m == 1 .and. n == 1) then
					data%rvec(:,ir)=rfile
				else if (any(data%rvec(:,ir) /= rfile)) then
					message='Wannier90 r-matrix records are not grouped by lattice vector'
					close(unit)
					deallocate(seen)
					call rmn_destroy(data)
					return
				end if
				if (seen(mfile,nfile)) then
					message='Duplicate Wannier-index pair in a Wannier90 r-matrix block'
					close(unit)
					deallocate(seen)
					call rmn_destroy(data)
					return
				end if
				seen(mfile,nfile)=.true.
				data%rmn(1,mfile,nfile,ir)=cmplx(xre,xim)
				data%rmn(2,mfile,nfile,ir)=cmplx(yre,yim)
				data%rmn(3,mfile,nfile,ir)=cmplx(zre,zim)
			end do
		end do
		if (.not. all(seen)) then
			message='Incomplete Wannier-index block in Wannier90 r-matrix file'
			close(unit)
			deallocate(seen)
			call rmn_destroy(data)
			return
		end if
	end do

	close(unit)
	deallocate(seen)
	ok=.true.
	message=''

end subroutine rmn_read


subroutine rmn_centers(data,centers,ok,message)

	! Extract the Wannier centres t_m=<0,m|r|0,m> from the R=0 block of
	! seedname_r.dat.  The full r-matrix remains available for the optical
	! position correction; this routine only provides its diagonal centres to
	! the G=0 Coulomb embedding.

	type(rmn_data),intent(in) :: data
	real,dimension(:,:),intent(out) :: centers
	logical,intent(out) :: ok
	character(len=*),intent(out) :: message

	integer :: ir,irzero,m
	real,parameter :: imaginary_tolerance=1.0e-5

	ok=.false.
	message=''
	centers=0.0

	if (.not. allocated(data%rvec) .or. .not. allocated(data%rmn) .or. &
		data%num_wann <= 0 .or. data%nrpts <= 0) then
		message='Wannier90 r-matrix table has not been initialized'
		return
	end if
	if (size(centers,1) /= 3 .or. size(centers,2) /= data%num_wann) then
		message='Incompatible output dimensions for Wannier centres'
		return
	end if

	irzero=0
	do ir=1,data%nrpts
		if (all(data%rvec(:,ir) == 0)) then
			if (irzero /= 0) then
				message='Wannier90 r-matrix table contains duplicate R=0 blocks'
				return
			end if
			irzero=ir
		end if
	end do
	if (irzero == 0) then
		message='Wannier90 r-matrix table has no R=0 block for Wannier centres'
		return
	end if

	do m=1,data%num_wann
		if (maxval(abs(aimag(data%rmn(:,m,m,irzero)))) > imaginary_tolerance) then
			message='Wannier90 R=0 diagonal position matrix is not real within tolerance'
			return
		end if
		centers(:,m)=real(data%rmn(:,m,m,irzero))
	end do

	ok=.true.

end subroutine rmn_centers


subroutine center_phase_build(centers,kpoints,phases,ok,message)

	! P_m(k)=exp(+i k.t_m), using the same positive Fourier convention as
	! hamiltonian() and rmn_bloch().  Precomputing this table avoids evaluating
	! exponentials inside every BSE Hamiltonian matrix element.

	real,dimension(:,:),intent(in) :: centers,kpoints
	complex,dimension(:,:),intent(out) :: phases
	logical,intent(out) :: ok
	character(len=*),intent(out) :: message

	integer :: ik,m
	real :: angle

	ok=.false.
	message=''
	phases=cmplx(1.0,0.0)

	if (size(centers,1) /= 3 .or. size(kpoints,2) /= 3 .or. &
		size(phases,1) /= size(centers,2) .or. size(phases,2) /= size(kpoints,1)) then
		message='Incompatible dimensions for Wannier-centre phase table'
		return
	end if

	do ik=1,size(kpoints,1)
		do m=1,size(centers,2)
			angle=dot_product(kpoints(ik,:),centers(:,m))
			phases(m,ik)=cmplx(cos(angle),sin(angle))
		end do
	end do

	ok=.true.

end subroutine center_phase_build


subroutine center_phase_direct_vertices(c1,c2,v1,v2,phase1,phase2,vc,vv)

	! For q=k1-k2, construct the two vertices whose product gives
	!
	!   sum_ij c1_i^* c2_i v1_j^* v2_j
	!          exp(+i q.t_i) exp(-i q.t_j).
	!
	! This is the G=0 centre-resolved direct kernel.  The screened Coulomb
	! scalar, normalization, and sign remain in the existing Coulomb routines.

	complex,dimension(:),intent(in) :: c1,c2,v1,v2,phase1,phase2
	complex,intent(out) :: vc,vv

	integer :: m
	complex :: relative_phase

	vc=cmplx(0.0,0.0)
	vv=cmplx(0.0,0.0)
	do m=1,size(c1)
		relative_phase=phase1(m)*conjg(phase2(m))
		vc=vc+conjg(c1(m))*c2(m)*relative_phase
		vv=vv+conjg(v1(m))*v2(m)*conjg(relative_phase)
	end do

end subroutine center_phase_direct_vertices


subroutine rmn_bloch(data,kpoint,rlat,connection,ok,message)

	! Evaluate the finite Fourier representation at one Cartesian k point:
	! A^W_alpha(k)=sum_R exp(+i k.R) <0,m|r_alpha|R,n>.
	! This is an exact reconstruction of the supplied finite Wannier90 table,
	! not a fitted or truncated interpolation. The positive phase is the same
	! convention used for H(k) by hamiltonian() in this code.

	type(rmn_data),intent(in) :: data
	real,dimension(3),intent(in) :: kpoint
	real,dimension(3,3),intent(in) :: rlat
	complex,dimension(:,:,:),intent(out) :: connection
	logical,intent(out) :: ok
	character(len=*),intent(out) :: message

	integer :: ir
	real,dimension(3) :: rcart
	real :: phase
	complex :: factor

	ok=.false.
	message=''
	connection=cmplx(0.0,0.0)

	if (.not. allocated(data%rvec) .or. .not. allocated(data%rmn) .or. &
		data%num_wann <= 0 .or. data%nrpts <= 0) then
		message='Wannier90 r-matrix table has not been initialized'
		return
	end if
	if (size(connection,1) /= 3 .or. size(connection,2) /= data%num_wann .or. &
		size(connection,3) /= data%num_wann) then
		message='Incompatible output dimensions for the Wannier connection matrix'
		return
	end if

	do ir=1,data%nrpts
		rcart=real(data%rvec(1,ir))*rlat(1,:)+real(data%rvec(2,ir))*rlat(2,:) + &
			real(data%rvec(3,ir))*rlat(3,:)
		phase=dot_product(kpoint,rcart)
		factor=cmplx(cos(phase),sin(phase))
		connection=connection+factor*data%rmn(:,:,:,ir)
	end do

	ok=.true.

end subroutine rmn_bloch


subroutine rmn_apply_q0_optical_correction(data,kpoint,rlat,ev,vv,ec,vc,sme,hrx,hry,hrz,ok,message)

	! Add the <0|r|R> contribution to one vertical (Q=0) optical matrix element.
	! The legacy optical path returns <c|dH/dk|v>/(Ec-Ev+i*sme). In the
	! Wannier gauge the corresponding velocity numerator is
	!
	!   <c| dH/dk - i [A,H] |v> = <c|dH/dk|v> + i (Ec-Ev) <c|A|v>,
	!
	! where A(k) is reconstructed by rmn_bloch. We apply the added numerator
	! with the same finite-sme convention as the existing path. At sme=0 this
	! is simply +i<c|A|v>, so no new approximation is introduced.

	type(rmn_data),intent(in) :: data
	real,dimension(3),intent(in) :: kpoint
	real,dimension(3,3),intent(in) :: rlat
	real,intent(in) :: ev,ec,sme
	complex,dimension(:),intent(in) :: vv,vc
	complex,intent(inout) :: hrx,hry,hrz
	logical,intent(out) :: ok
	character(len=*),intent(out) :: message

	complex :: connection(3,data%num_wann,data%num_wann)
	complex :: acv(3),factor
	integer :: alpha,m,n

	ok=.false.
	message=''
	if (size(vv) /= data%num_wann .or. size(vc) /= data%num_wann) then
		message='Incompatible eigenvector dimensions for the Wannier connection matrix'
		return
	end if

	call rmn_bloch(data,kpoint,rlat,connection,ok,message)
	if (.not. ok) return

	acv=cmplx(0.0,0.0)
	do alpha=1,3
		do m=1,data%num_wann
			do n=1,data%num_wann
				acv(alpha)=acv(alpha)+conjg(vc(m))*connection(alpha,m,n)*vv(n)
			end do
		end do
	end do

	factor=cmplx(0.0,ec-ev)/cmplx(ec-ev,sme)
	hrx=hrx+factor*acv(1)
	hry=hry+factor*acv(2)
	hrz=hrz+factor*acv(3)

	! factor is i*(Ec-Ev)/(Ec-Ev+i*sme), matching the legacy convention.
	ok=.true.

end subroutine rmn_apply_q0_optical_correction

end module bse_q_optics
