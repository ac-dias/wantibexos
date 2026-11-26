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

	implicit none

	logical,save :: rb_active = .false.
	integer,save :: rb_nk = 0, rb_nt = 0, rb_n = 0, rb_nunit = 0
	integer,allocatable,save,dimension(:) :: rb_kbar
	integer,allocatable,save,dimension(:) :: rb_kind, rb_row, rb_rowbar
	! units of the real basis: one s,t couple or one single function; first function and length of each
	integer,allocatable,save,dimension(:) :: rb_ustart, rb_ulen
	real,parameter :: rb_r2 = 0.70710678118654752440
	real,parameter :: rb_s2 = 1.41421356237309504880

contains

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


	subroutine rb_setup(nk,nt,kbar)

		! tables of the real basis for nk k points of nt transitions each (row of transition it of k: (k-1) nt + it)

		implicit none

		integer,intent(in) :: nk,nt
		integer,dimension(nk),intent(in) :: kbar

		integer :: ik,it,r,u,a,ab

		call rb_free()
		rb_nk = nk
		rb_nt = nt
		rb_n = nk*nt
		allocate(rb_kbar(nk),rb_kind(rb_n),rb_row(rb_n),rb_rowbar(rb_n))
		rb_kbar = kbar

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
		do ik=1,nk
			if (kbar(ik) == ik) then
				do it=1,nt
					a = (ik-1)*nt+it
					r = r+1
					rb_kind(r) = 0
					rb_row(r) = a
					rb_rowbar(r) = a
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
		if (allocated(rb_ustart)) deallocate(rb_ustart)
		if (allocated(rb_ulen)) deallocate(rb_ulen)
		rb_active = .false.
		rb_n = 0
		rb_nunit = 0

	end subroutine rb_free


	subroutine rb_tr_gauge(nk,nb,norb,kbar,vector,eigt,nocpk,devmax,nmis)

		! time-reversal gauge: the band energies, the band index of the Fermi level and the eigenvectors of the k point
		! kbar > k are those of k, conjugated.  The vectors of the self-conjugate k points were made real by
		! rb_selfconj_vectors.  Done after the k points are shared, the same way on every rank.
		! devmax: the largest difference between the band energies of k and of kbar, as the eigensolver gave them (they
		! are equal for a time-reversal invariant Hamiltonian); nmis: the pairs with different Fermi level band index.

		implicit none

		integer,intent(in) :: nk,nb,norb
		integer,dimension(nk),intent(in) :: kbar
		complex,dimension(norb,nb,nk),intent(inout) :: vector
		real,dimension(nb,nk),intent(inout) :: eigt
		integer,dimension(nk),intent(inout) :: nocpk
		real,intent(out) :: devmax
		integer,intent(out) :: nmis

		integer :: ik,jk

		devmax = 0.0
		nmis = 0
		do ik=1,nk
			jk = kbar(ik)
			if (jk > ik) then
				devmax = max(devmax,maxval(abs(eigt(:,jk)-eigt(:,ik))))
				if (nocpk(jk) /= nocpk(ik)) nmis = nmis+1
				vector(:,:,jk) = conjg(vector(:,:,ik))
				eigt(:,jk) = eigt(:,ik)
				nocpk(jk) = nocpk(ik)
			end if
		end do

	end subroutine rb_tr_gauge


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
				c = rb_r2*(h(rb_row(r))+h(rb_rowbar(r)))
			case (2)
				c = cmplx(0.0,rb_r2)*(h(rb_row(r))-h(rb_rowbar(r)))
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
				ac(rb_rowbar(r)) = ac(rb_rowbar(r))+rb_r2*a(r)
			case (2)
				ac(rb_row(r)) = ac(rb_row(r))+cmplx(0.0,-rb_r2)*a(r)
				ac(rb_rowbar(r)) = ac(rb_rowbar(r))+cmplx(0.0,rb_r2)*a(r)
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
