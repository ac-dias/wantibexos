module bse_realbasis

	! Time-reversal real basis of the Q = 0 BSE Hamiltonian (BSE_REAL= T).
	!
	! Basis.  The exciton basis is |a> = |v c k>, a = (k, transition), the states of every k point of the mesh.  For a
	! time-reversal invariant crystal (real H(R) and S(R)) conjg(psi_nk) is an eigenvector of H(-k), so the vectors of the
	! k point kbar = -k (mod G) are taken as the complex conjugates of those of k (rb_tr_gauge), and the state |abar> of
	! the same bands at kbar is the time-reversed partner of |a>.  The kernel and the exchange term then obey
	!
	!    H(abar, bbar) = conjg( H(a, b) )                                                                       (1)
	!
	! The states of a pair k < kbar are replaced by the two real combinations
	!
	!    |s> = ( |a> + |abar> ) / sqrt(2) ,       |t> = ( |a> - |abar> ) / ( i sqrt(2) ) ,
	!
	! and a state of a self-conjugate k point (kbar = k, real vectors, rb_selfconj_vectors) stays as it is, |e> = |a>.
	! With (1) the matrix elements in this basis are all real (h1 = H(a,b), h2 = H(a,bbar), a and b of pairs):
	!
	!    <s_a|H|s_b> = Re h1 + Re h2       <s_a|H|t_b> = Im h1 - Im h2
	!    <t_a|H|s_b> = -(Im h1 + Im h2)    <t_a|H|t_b> = Re h1 - Re h2
	!    <e_a|H|s_b> = sqrt2 Re h1         <e_a|H|t_b> = sqrt2 Im h1        <e_a|H|e_b> = Re h1
	!    <s_a|H|e_b> = sqrt2 Re h1         <t_a|H|e_b> = -sqrt2 Im h1
	!
	! so the Hamiltonian is a real symmetric matrix of the same dimension N, U^+ H U with U unitary (U has 2 x 2 blocks
	! [1, -i; 1, i]/sqrt2 for the pairs).  It takes half of the memory, is diagonalized by the real symmetric solvers, and
	! two complex kernel elements give four real ones.  The eigenvectors a of the real matrix are real; the optical
	! amplitude of an exciton is sum_j conjg(A_j) h_j = sum_r a_r h'_r with h' = U^+ h (rb_vector), and the complex
	! amplitudes A = U a, in the basis of the states of the mesh, are rebuilt with rb_expand when they are needed.
	!
	! Real basis function r of the table: rb_kind(r) = 1 (s) and 2 (t), adjacent, for each transition of each pair
	! k < kbar (the pairs first, an even number of functions, so that a block of an even block size of the distributed
	! matrix never splits an s,t couple), then rb_kind(r) = 0 (e) for the transitions of the self-conjugate k points.
	! rb_row(r) is the row of the representative state a in the complex basis and rb_rowbar(r) that of abar.
	!
	! Spin-orbit models (mode 2, T^2 = -1).  Time reversal is T = U K with the spin flip U = i sigma_y x 1 (a signed
	! permutation of the spinor components, rb_tr_rep), U U* = -1, and H(R) = U H(R)* U^+ is complex.  The pair operator of the
	! excitons, Theta = T_e x T_h, has Theta^2 = (-1)(-1) = +1, so a real basis exists all the same.  The gauge is
	! psi_{n kbar} = U psi*_{n k} (rb_tr_gauge) for the pairs k < kbar.  At a self-conjugate k point the partner of a band
	! is its Kramers partner, T psi_{2j-1} = psi_{2j}, T psi_{2j} = -psi_{2j-1} (rb_kramers_gauge), so
	!
	!    Theta |a> = sigma_a |abar> ,     abar = (vbar, cbar, k) ,     sigma_a = sigma_c sigma_v ,    sigma_n = +1 (odd), -1 (even)
	!
	! and Theta^2 |a> = |a>.  The state paired with |a> is |atilde> = Theta |a> = sigma_a |abar> (sigma = +1 for the pairs
	! k < kbar), and everything above holds with abar -> atilde:  H(atilde, btilde) = conjg(H(a,b)), h2 = H(a, btilde) =
	! sigma_b H(a, bbar), h'_s = ( h_a + sigma_a h_abar )/sqrt2, h'_t = i ( h_a - sigma_a h_abar )/sqrt2.  A self-conjugate k
	! point then has no unpaired function e (the transition partner tbar is never the transition itself).
	! rb_sgn(r) = sigma_a of the function r (1 for the functions of the pairs k < kbar and for e).

	implicit none

	logical,save :: rb_active = .false.
	integer,save :: rb_nk = 0, rb_nt = 0, rb_n = 0, rb_nunit = 0
	integer,allocatable,save,dimension(:) :: rb_kbar
	integer,allocatable,save,dimension(:) :: rb_kind, rb_row, rb_rowbar, rb_sgn
	! units of the real basis: one s,t couple or one single function; first function and length of each
	integer,allocatable,save,dimension(:) :: rb_ustart, rb_ulen
	real,parameter :: rb_r2 = 0.70710678118654752440
	real,parameter :: rb_s2 = 1.41421356237309504880

contains

	integer function rb_lwork(query,lwmin)

		! Workspace size LWORK from the result of a LAPACK workspace query, which comes back in a REAL.  Single precision holds
		! every integer only up to 2**24 = 16777216: above it the spacing is 2, and a size such as 1+6N+2N**2 of SSYEVD, which
		! is always odd, is rounded DOWN by one for every N >= 2895.  LAPACK 3.10 and OpenBLAS round the value up before they
		! return it; MKL does not, and the routine then refuses LWORK as one too small (INFO < 0).  nearest() is the next
		! single precision number above the query, so the result is never below the size the routine asked for.
		! lwmin: the minimum documented for the routine (1 when the query is all there is).  The result is -1 when the size
		! does not fit an INTEGER(4), which LAPACK takes as LWORK (SSYEVD above N = 32766: use ScaLAPACK or ELPA).

		implicit none

		real,intent(in) :: query
		integer(8),intent(in) :: lwmin

		integer(8) :: l

		l = max(int(nearest(query,1.0),8),lwmin,1_8)
		if (l > int(huge(rb_lwork),8)) then
			rb_lwork = -1
		else
			rb_lwork = int(l)
		end if

	end function rb_lwork

	subroutine rb_partners(nk,kpt,rlat,tol,kbar,nself,ok)

		! kbar(i): the k point of the mesh equal to -kpt(i) plus a reciprocal lattice vector, within tol in the
		! fractional coordinates; nself: how many k points are their own partner.  ok = .false. when a k point has none
		! (a mesh that is not symmetric under k -> -k, for instance a shifted one).

		implicit none

		integer,intent(in) :: nk
		real,dimension(nk,3),intent(in) :: kpt
		real,dimension(3,3),intent(in) :: rlat
		real,intent(in) :: tol
		integer,dimension(nk),intent(out) :: kbar
		integer,intent(out) :: nself
		logical,intent(out) :: ok

		real(kind=8),allocatable,dimension(:,:) :: f
		real(kind=8) :: twopi,d
		real,dimension(3) :: b1,b2,b3
		integer :: i,j,m
		logical :: hit

		twopi = 8.0d0*atan(1.0d0)
		call recvec(rlat(1,:),rlat(2,:),rlat(3,:),b1,b2,b3)

		! fractional coordinates of the k points, k = f1 b1 + f2 b2 + f3 b3
		allocate(f(nk,3))
		do i=1,nk
			do m=1,3
				f(i,m) = (dble(kpt(i,1))*dble(rlat(m,1))+dble(kpt(i,2))*dble(rlat(m,2))+dble(kpt(i,3))*dble(rlat(m,3)))/twopi
			end do
		end do

		kbar = 0
		nself = 0
		ok = .true.
		do i=1,nk
			if (kbar(i) /= 0) cycle
			hit = .false.
			do j=i,nk
				if (kbar(j) /= 0) cycle
				hit = .true.
				do m=1,3
					d = f(i,m)+f(j,m)
					if (abs(d-anint(d)) > dble(tol)) hit = .false.
				end do
				if (hit) then
					kbar(i) = j
					kbar(j) = i
					if (i == j) nself = nself+1
					exit
				end if
			end do
			if (.not. hit) then
				ok = .false.
				deallocate(f)
				return
			end if
		end do
		deallocate(f)

	end subroutine rb_partners


	subroutine rb_input_check(ngrid,mshift,temp,exc,ham_io,eph,ok,message)

		! What BSE_REAL= T needs and can take, from the input alone, so that the program can stop before it reads the
		! tight-binding file or computes anything.  The first of them that is not met is returned in message.
		! The mesh k = b_i (i/n_i + s_i), i = 0..n_i-1, contains -k for every k exactly when 2 s_i n_i is an integer for
		! each direction: s_i = 0 (Gamma-centred) or a half step, s_i = 1/(2 n_i), and their multiples by 1/(2 n_i).

		implicit none

		integer,dimension(3),intent(in) :: ngrid
		real,dimension(3),intent(in) :: mshift
		real,intent(in) :: temp,exc
		logical,intent(in) :: ham_io,eph
		logical,intent(out) :: ok
		character(len=*),intent(out) :: message

		integer :: i
		real(kind=8) :: x
		character(len=1),parameter :: axis(3) = (/ 'X','Y','Z' /)
		character(len=12) :: s1,s2

		ok = .true.
		message = ''
		do i=1,3
			x = 2.0d0*dble(mshift(i))*dble(ngrid(i))
			if (abs(x-anint(x)) > 1.0d-3) then
				ok = .false.
				write(s1,'(F12.6)') mshift(i)
				write(s2,'(F12.6)') 0.5d0/dble(ngrid(i))
				write(message,'(A,I0,A,I0,A)') 'BSE_REAL= T needs a k mesh that is symmetric under k -> -k, and SHIFT_'//achar(48+i)// &
					'= '//trim(adjustl(s1))//' with NG'//axis(i)//'= ',ngrid(i),' is not: the shift must be a multiple of 1/(2 NG'// &
					axis(i)//') = '//trim(adjustl(s2))//' (SHIFT_',i,'= 0 or the half step '//trim(adjustl(s2))//' keep the mesh symmetric); set BSE_REAL= F to keep this mesh'
				return
			end if
		end do
		if (temp /= 0.0) then
			ok = .false.
			message = 'BSE_REAL= T is not implemented for TEMP > 0 (the finite-temperature BSE solver); set BSE_REAL= F'
		else if (ham_io) then
			ok = .false.
			message = 'BSE_REAL= T does not support BSE_HAM_SAVE / BSE_HAM_READ yet'
		else if (eph) then
			ok = .false.
			message = 'BSE_REAL= T does not support BSE_EPH_FILE'
		else if (abs(exc) > 0.0) then
			ok = .false.
			message = 'BSE_REAL= T needs a time-reversal invariant Hamiltonian: SPIN_EXC= must be 0'
		end if

	end subroutine rb_input_check


	subroutine rb_setup(nk,nt,kbar,tbar,sig)

		! tables of the real basis for nk k points of nt transitions each (row of transition it of k: (k-1) nt + it).
		! At a self-conjugate k point the partner of transition it is transition tbar(it) of the same k point, with the
		! sign sig(it) = sigma_a (tbar = identity and sig = 1 without spin-orbit coupling: the single functions e;
		! tbar without fixed points with spin-orbit coupling: only pairs).  The pairs k < kbar have the same transition
		! and sigma = 1.

		implicit none

		integer,intent(in) :: nk,nt
		integer,dimension(nk),intent(in) :: kbar
		integer,dimension(nt),intent(in) :: tbar,sig

		integer :: ik,it,r,u,a,ab

		call rb_free()
		rb_nk = nk
		rb_nt = nt
		rb_n = nk*nt
		allocate(rb_kbar(nk),rb_kind(rb_n),rb_row(rb_n),rb_rowbar(rb_n),rb_sgn(rb_n))
		rb_kbar = kbar
		rb_sgn = 1

		r = 0
		do ik=1,nk
			if (kbar(ik) > ik) then
				do it=1,nt
					a = (ik-1)*nt+it
					ab = (kbar(ik)-1)*nt+it
					r = r+1
					rb_kind(r) = 1
					rb_row(r) = a
					rb_rowbar(r) = ab
					r = r+1
					rb_kind(r) = 2
					rb_row(r) = a
					rb_rowbar(r) = ab
				end do
			end if
		end do
		! pairs of transitions inside a self-conjugate k point (Kramers partners)
		do ik=1,nk
			if (kbar(ik) == ik) then
				do it=1,nt
					if (tbar(it) > it) then
						a = (ik-1)*nt+it
						ab = (ik-1)*nt+tbar(it)
						r = r+1
						rb_kind(r) = 1
						rb_row(r) = a
						rb_rowbar(r) = ab
						rb_sgn(r) = sig(it)
						r = r+1
						rb_kind(r) = 2
						rb_row(r) = a
						rb_rowbar(r) = ab
						rb_sgn(r) = sig(it)
					end if
				end do
			end if
		end do
		do ik=1,nk
			if (kbar(ik) == ik) then
				do it=1,nt
					if (tbar(it) == it) then
						a = (ik-1)*nt+it
						r = r+1
						rb_kind(r) = 0
						rb_row(r) = a
						rb_rowbar(r) = a
					end if
				end do
			end if
		end do
		if (r /= rb_n) stop 'rb_setup: the real basis does not have the dimension of the exciton basis'

		! units
		rb_nunit = 0
		r = 1
		do while (r <= rb_n)
			rb_nunit = rb_nunit+1
			if (rb_kind(r) == 1) then
				r = r+2
			else
				r = r+1
			end if
		end do
		allocate(rb_ustart(rb_nunit),rb_ulen(rb_nunit))
		u = 0
		r = 1
		do while (r <= rb_n)
			u = u+1
			rb_ustart(u) = r
			if (rb_kind(r) == 1) then
				rb_ulen(u) = 2
				r = r+2
			else
				rb_ulen(u) = 1
				r = r+1
			end if
		end do
		rb_active = .true.

	end subroutine rb_setup


	subroutine rb_free()

		implicit none

		if (allocated(rb_kbar)) deallocate(rb_kbar)
		if (allocated(rb_kind)) deallocate(rb_kind)
		if (allocated(rb_row)) deallocate(rb_row)
		if (allocated(rb_rowbar)) deallocate(rb_rowbar)
		if (allocated(rb_sgn)) deallocate(rb_sgn)
		if (allocated(rb_ustart)) deallocate(rb_ustart)
		if (allocated(rb_ulen)) deallocate(rb_ulen)
		rb_active = .false.
		rb_n = 0
		rb_nunit = 0

	end subroutine rb_free


	subroutine rb_tr_gauge(nk,nb,norb,kbar,pmap,psgn,vector,eigt,nocpk,devmax,nmis)

		! time-reversal gauge: the band energies, the band index of the Fermi level and the eigenvectors of the k point
		! kbar > k are those of k, time reversed: psi_kbar(i) = psgn(i) conjg(psi_k(pmap(i))), i.e. U psi*_k with the
		! signed permutation U of rb_tr_rep (the identity without spin-orbit coupling: the complex conjugate).  The
		! vectors of the self-conjugate k points were made real by rb_selfconj_vectors, or Kramers-paired by
		! rb_kramers_gauge.  Done after the k points are shared, the same way on every rank.
		! devmax: the largest difference between the band energies of k and of kbar, as the eigensolver gave them (they
		! are equal for a time-reversal invariant Hamiltonian); nmis: the pairs with different Fermi level band index.

		implicit none

		integer,intent(in) :: nk,nb,norb
		integer,dimension(nk),intent(in) :: kbar
		integer,dimension(norb),intent(in) :: pmap
		real,dimension(norb),intent(in) :: psgn
		complex,dimension(norb,nb,nk),intent(inout) :: vector
		real,dimension(nb,nk),intent(inout) :: eigt
		integer,dimension(nk),intent(inout) :: nocpk
		real,intent(out) :: devmax
		integer,intent(out) :: nmis

		integer :: ik,jk,i

		devmax = 0.0
		nmis = 0
		do ik=1,nk
			jk = kbar(ik)
			if (jk > ik) then
				devmax = max(devmax,maxval(abs(eigt(:,jk)-eigt(:,ik))))
				if (nocpk(jk) /= nocpk(ik)) nmis = nmis+1
				do i=1,norb
					vector(i,:,jk) = psgn(i)*conjg(vector(pmap(i),:,ik))
				end do
				eigt(:,jk) = eigt(:,ik)
				nocpk(jk) = nocpk(ik)
			end if
		end do

	end subroutine rb_tr_gauge


	subroutine rb_tr_rep(dft,systype,nbasis,nvec,hop,ihop,ovp,mode,pmap,psgn,dev1,dev2)

		! The representation of time reversal on the orbitals, T = U K with U a signed permutation: mode 1, U = 1
		! (T^2 = +1, real H(R)), or mode 2, the spin flip of a spinor basis, U = i sigma_y x 1 (T^2 = -1): orbital i goes
		! to pmap(i) with the sign psgn(i) (+1 for a spin up orbital, -1 for spin down), U psi*  has components
		! psgn(i) conjg(psi(pmap(i))).  The layout is that of spin_z_matrix: DFT=S, systype SOC or NC: the first half of the
		! orbitals is spin up and the second half spin down; DFT=W, systype SOC or NC: orbitals 2m-1 (up) and 2m (down).
		! Time-reversal invariance of H and S is H(R) = U H(R)* U^+ for every R,
		!
		!    H_ij(R) = psgn(i) psgn(j) conjg( H_{pmap(i) pmap(j)}(R) ) ,
		!
		! tested here for the elements of the file: dev1 is the largest imaginary part of H(R) (mode 1) and dev2 the
		! largest deviation of H and S from the spin-flip form (mode 2, -1 when it was not tried); mode = 0 when neither
		! holds within 1e-4.  Mode 1 is tried first: for a spinor model without spin-orbit coupling, or any other real
		! H(R), it is the one in use.

		implicit none

		character(len=1),intent(in) :: dft
		character(len=4),intent(in) :: systype
		integer,intent(in) :: nbasis,nvec
		real,dimension(nvec,nbasis,nbasis),intent(in) :: hop,ihop,ovp
		integer,intent(out) :: mode
		integer,dimension(nbasis),intent(out) :: pmap
		real,dimension(nbasis),intent(out) :: psgn
		real,intent(out) :: dev1,dev2

		integer :: i,j,w2
		real,parameter :: tol = 1.0e-4
		real :: d,devh,devs

		! mode 1
		do i=1,nbasis
			pmap(i) = i
			psgn(i) = 1.0
		end do
		dev1 = maxval(abs(ihop))
		dev2 = -1.0
		mode = 1
		if (dev1 <= tol) return

		! mode 2
		mode = 0
		if (mod(nbasis,2) /= 0) return
		if (trim(systype) /= 'SOC' .and. trim(systype) /= 'NC') return
		w2 = nbasis/2
		do i=1,nbasis
			if (dft == 'W') then
				if (mod(i,2) == 1) then
					pmap(i) = i+1
					psgn(i) = 1.0
				else
					pmap(i) = i-1
					psgn(i) = -1.0
				end if
			else
				if (i <= w2) then
					pmap(i) = i+w2
					psgn(i) = 1.0
				else
					pmap(i) = i-w2
					psgn(i) = -1.0
				end if
			end if
		end do
		devh = 0.0
		devs = 0.0
		!$omp parallel do default(shared) private(i,j,d) reduction(max:devh,devs)
		do j=1,nbasis
			do i=1,nbasis
				d = psgn(i)*psgn(j)
				devh = max(devh,maxval(abs(hop(:,i,j)-d*hop(:,pmap(i),pmap(j)))), &
					maxval(abs(ihop(:,i,j)+d*ihop(:,pmap(i),pmap(j)))))
				if (dft == 'S') devs = max(devs,maxval(abs(ovp(:,i,j)-d*ovp(:,pmap(i),pmap(j)))))
			end do
		end do
		!$omp end parallel do
		dev2 = max(devh,devs)
		if (dev2 <= tol) mode = 2

	end subroutine rb_tr_rep


	function rb_element(kind1,kind2,h1,h2)

		! real matrix element between the real basis functions of kinds kind1 (row) and kind2 (column) from
		! h1 = H(a,b) and h2 = H(a,bbar) (h2 is used only when both are of a pair)

		implicit none

		integer,intent(in) :: kind1,kind2
		complex,intent(in) :: h1,h2
		real :: rb_element

		select case (10*kind1+kind2)
		case (0)
			rb_element = real(h1)
		case (1)
			rb_element = rb_s2*real(h1)
		case (2)
			rb_element = rb_s2*aimag(h1)
		case (10)
			rb_element = rb_s2*real(h1)
		case (11)
			rb_element = real(h1)+real(h2)
		case (12)
			rb_element = aimag(h1)-aimag(h2)
		case (20)
			rb_element = -rb_s2*aimag(h1)
		case (21)
			rb_element = -(aimag(h1)+aimag(h2))
		case (22)
			rb_element = real(h1)-real(h2)
		case default
			rb_element = 0.0
		end select

	end function rb_element


	subroutine rb_vector(h,hre,him)

		! h' = U^+ h of a vector of the complex basis (for the optical amplitudes sum_j conjg(A_j) h_j = sum_r a_r h'_r),
		! as its real and imaginary parts

		implicit none

		complex,dimension(rb_n),intent(in) :: h
		real,dimension(rb_n),intent(out) :: hre,him

		integer :: r
		complex :: c

		do r=1,rb_n
			select case (rb_kind(r))
			case (1)
				c = rb_r2*(h(rb_row(r))+real(rb_sgn(r))*h(rb_rowbar(r)))
			case (2)
				c = cmplx(0.0,rb_r2)*(h(rb_row(r))-real(rb_sgn(r))*h(rb_rowbar(r)))
			case default
				c = h(rb_row(r))
			end select
			hre(r) = real(c)
			him(r) = aimag(c)
		end do

	end subroutine rb_vector


	subroutine rb_expand(a,ac)

		! complex amplitudes A = U a on the states of the mesh, from the real eigenvector a

		implicit none

		real,dimension(rb_n),intent(in) :: a
		complex,dimension(rb_n),intent(out) :: ac

		integer :: r

		ac = cmplx(0.0,0.0)
		do r=1,rb_n
			select case (rb_kind(r))
			case (1)
				ac(rb_row(r)) = ac(rb_row(r))+rb_r2*a(r)
				ac(rb_rowbar(r)) = ac(rb_rowbar(r))+real(rb_sgn(r))*rb_r2*a(r)
			case (2)
				ac(rb_row(r)) = ac(rb_row(r))+cmplx(0.0,-rb_r2)*a(r)
				ac(rb_rowbar(r)) = ac(rb_rowbar(r))+cmplx(0.0,real(rb_sgn(r))*rb_r2)*a(r)
			case default
				ac(rb_row(r)) = ac(rb_row(r))+a(r)
			end select
		end do

	end subroutine rb_expand

end module bse_realbasis


subroutine rb_selfconj_vectors(dft,nbasis,nvec,rvec,hopmatrices,ihopmatrices,ffactor,ovp,kx,ky,kz,nb,first,w,vreal,info)

	! Real eigenvectors of the bands first, ..., first+nb-1 of a self-conjugate k point (k = -k + G): H(k) and S(k) are
	! real there for real H(R) and S(R), and the real generalized symmetric solver (SSYGV, or SSYEV for the orthogonal
	! DFT=W basis) gives real eigenvectors, also in a degenerate subspace; the complex solver would leave a free phase
	! (or a free unitary matrix) on them.  The eigenvalues w of these bands are returned (the shifts of eigsys are not
	! applied), info /= 0 when the solver fails.

	implicit none

	character(len=1),intent(in) :: dft
	integer,intent(in) :: nbasis,nvec,nb,first
	real,dimension(nvec,3),intent(in) :: rvec
	real,dimension(nvec,nbasis,nbasis),intent(in) :: hopmatrices,ihopmatrices,ovp
	integer,dimension(nvec),intent(in) :: ffactor
	real,intent(in) :: kx,ky,kz
	real,dimension(nb),intent(out) :: w
	complex,dimension(nbasis,nb),intent(out) :: vreal
	integer,intent(out) :: info

	complex,allocatable,dimension(:,:) :: htb,stb
	real,allocatable,dimension(:,:) :: hr,sr
	real,allocatable,dimension(:) :: wall,work
	integer :: lwork,i,j
	real :: wq(1)

	allocate(htb(nbasis,nbasis),stb(nbasis,nbasis),hr(nbasis,nbasis),sr(nbasis,nbasis),wall(nbasis))
	call hamiltonian(nbasis,nvec,rvec,hopmatrices,ihopmatrices,ffactor,kx,ky,kz,htb)
	hr = real(htb)
	sr = 0.0
	if (dft .eq. "S") then
		call overlap(nbasis,nvec,rvec,ovp,kx,ky,kz,stb)
		sr = real(stb)
	end if

	info = 0
	lwork = -1
	if (dft .eq. "S") then
		call SSYGV(1,'V','U',nbasis,hr,nbasis,sr,nbasis,wall,wq,lwork,info)
		lwork = max(1,int(wq(1)))
		allocate(work(lwork))
		call SSYGV(1,'V','U',nbasis,hr,nbasis,sr,nbasis,wall,work,lwork,info)
	else
		call SSYEV('V','U',nbasis,hr,nbasis,wall,wq,lwork,info)
		lwork = max(1,int(wq(1)))
		allocate(work(lwork))
		call SSYEV('V','U',nbasis,hr,nbasis,wall,work,lwork,info)
	end if
	if (info == 0) then
		do j=1,nb
			w(j) = wall(first+j-1)
			do i=1,nbasis
				vreal(i,j) = cmplx(hr(i,first+j-1),0.0)
			end do
		end do
	end if
	deallocate(htb,stb,hr,sr,wall,work)

end subroutine rb_selfconj_vectors


subroutine rb_kramers_gauge(dft,nbasis,nvec,rvec,ovp,kx,ky,kz,pmap,psgn,nb,w,vec,deve,devo)

	! Kramers gauge of the bands of a self-conjugate k point (k = -k + G) of a spin-orbit model (T^2 = -1): the bands of
	! the window come in doublets (2j-1, 2j), j = 1, ..., nb/2 (the window starts at an odd band index and has an even
	! number of bands: the caller checks it), equal in energy.  The first vector of each doublet is kept (any
	! eigenvector of the doublet is one), and the second is its time-reversed partner, vec(:,2j) = U conjg(vec(:,2j-1)),
	! so that T psi_{2j-1} = psi_{2j} and T psi_{2j} = -psi_{2j-1}.
	! deve: the largest energy difference inside a doublet (eV); devo: the largest deviation of the overlap matrix of the
	! window vectors, vec^+ S(k) vec, from the unit matrix (it fails for an accidental degeneracy of two doublets, where the
	! partner of a vector is not orthogonal to the other doublet).

	implicit none

	character(len=1),intent(in) :: dft
	integer,intent(in) :: nbasis,nvec,nb
	real,dimension(nvec,3),intent(in) :: rvec
	real,dimension(nvec,nbasis,nbasis),intent(in) :: ovp
	real,intent(in) :: kx,ky,kz
	integer,dimension(nbasis),intent(in) :: pmap
	real,dimension(nbasis),intent(in) :: psgn
	real,dimension(nb),intent(in) :: w
	complex,dimension(nbasis,nb),intent(inout) :: vec
	real,intent(out) :: deve,devo

	complex,allocatable,dimension(:,:) :: stb,sv,g
	integer :: i,j

	deve = 0.0
	do j=1,nb-1,2
		deve = max(deve,abs(w(j)-w(j+1)))
		do i=1,nbasis
			vec(i,j+1) = psgn(i)*conjg(vec(pmap(i),j))
		end do
	end do

	allocate(sv(nbasis,nb),g(nb,nb))
	if (dft .eq. "S") then
		allocate(stb(nbasis,nbasis))
		call overlap(nbasis,nvec,rvec,ovp,kx,ky,kz,stb)
		sv = matmul(stb,vec)
		deallocate(stb)
	else
		sv = vec
	end if
	g = matmul(conjg(transpose(vec)),sv)
	devo = 0.0
	do j=1,nb
		do i=1,nb
			if (i == j) then
				devo = max(devo,abs(g(i,j)-cmplx(1.0,0.0)))
			else
				devo = max(devo,abs(g(i,j)))
			end if
		end do
	end do
	deallocate(sv,g)

end subroutine rb_kramers_gauge
