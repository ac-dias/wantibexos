module bse_fastkernel

	! Direct-kernel vertices of the non-orthogonal DFT=S basis from per-state vectors.
	!
	! matrizelbse needs, for every kernel element, two vertices (electron: c1 k1 -> c2 k2, hole: v1 k1 -> v2 k2),
	! each the sandwich of sandwich_phase (general_subs.F90):
	!
	!   V = 1/2 sum_ij conjg(l_i) [ S1_ij ph_j + S2_ij ph_i ] r_j ,   ph_j = exp(i q.tau_j),  q = k1 - k2 - G ,
	!
	! S1 = S(k1), S2 = S(k2), l the state at k1, r the state at k2, G the image shift of bse_q_images. Done per element
	! it is a dense N x N double sum (N = number of basis functions, 3104 for the 4 x 4 supercell), repeated for every
	! (v, v') or (c, c') that shares the same (k1, k2). The phase factorizes, ph_j = e_j(k1) conjg(e_j(k2)) conjg(eG_j)
	! with e_j(k) = exp(i k.tau_j) and eG_j = exp(i G.tau_j), so with the per-state vectors
	!
	!   A_l = e(k1) * (S1^T conjg(l))    L_l = conjg(l) * e(k1)       (left state)
	!   R_r = conjg(e(k2)) * r           B_r = conjg(e(k2)) * (S2 r)  (right state)
	!
	! the vertex is   V = 1/2 sum_j conjg(eG_j) [ A_l,j R_r,j + L_l,j B_r,j ] ,
	! a product of small matrices (bands x N)(N x bands) for each (k1, k2, image). The table fk_vertex holds it for all
	! band pairs of the k pairs this rank needs, so a kernel element is two lookups. The overlap matrices are needed one
	! k at a time.
	! The arithmetic is the one of sandwich_phase in a different order (single precision rounding only).
	!
	! A rank evaluates the elements (i, j) of its own rows and columns of the Hamiltonian only, i = (k1, c1, v1),
	! j = (k2, c2, v2), so it needs the pairs (k1, k2) of its row k-points times its column k-points (and, when only
	! the upper triangle i <= j is built and the states are ordered by k, k1 <= k2). fk_setup stores those pairs
	! alone, with the images that bse_q_images returns for each (1 on most pairs, 2 or 3 on the zone boundary), so
	! the table of a run on P ranks is about 1/P of the table of all pairs, and it is freed by fk_free once the
	! Hamiltonian is built.
	!
	! BSE_QFOLD = M > 1 (an M x M supercell run): the images of a pair are those of bse_q_images_fold (the shortest image of each
	! of the M^2 classes q + g + Gamma_UC, up to 27 M^2 per pair, the numbering of matrizelbse); fk_setup takes M as mfold.
	! The table of that kernel is about M^2 times larger per pair than for M = 1, so the pairs a rank does not use matter most there.

	implicit none

	logical :: fk_active = .false.
	! fk_vertex(b1, b2, fk_slot(k1, k2, image)): band columns of vector(:, b, ik), image numbered as bse_q_images(k1, k2)
	complex,allocatable,dimension(:,:,:) :: fk_vertex
	! k point -> slot of the row (k1) and column (k2) k points of this rank, 0 when the rank does not need it
	integer,allocatable,dimension(:) :: fk_rmap,fk_cmap
	! fk_first(slot of k2, slot of k1): number of images stored before this pair, -1 for a pair that is not stored
	integer,allocatable,dimension(:,:) :: fk_first

	! BSE_QFOLD > 1: images, weights and potential W(q) of the last (k1, k2) pair a thread handled. They do not depend on the
	! band indices, so the transitions of one k pair (nc nv of them, consecutive in the loops of the solver) share them.
	! The key is (k1, k2, EDIEL, tolerance, potential); the rest of the arguments of W is constant during a run.
	logical,save :: fk_qc_valid = .false.
	integer,save :: fk_qc_nimg = 0
	real,save,dimension(3) :: fk_qc_k1 = 0.0, fk_qc_k2 = 0.0, fk_qc_ediel = 0.0
	real,save :: fk_qc_tolr = 0.0
	character(len=10),save :: fk_qc_coul = ''
	real,allocatable,save,dimension(:,:) :: fk_qc_shift
	real,allocatable,save,dimension(:) :: fk_qc_wimg, fk_qc_w
	!$omp threadprivate(fk_qc_valid,fk_qc_nimg,fk_qc_k1,fk_qc_k2,fk_qc_ediel,fk_qc_tolr,fk_qc_coul,fk_qc_shift,fk_qc_wimg,fk_qc_w)

contains

	subroutine fk_setup(n,nvec,nb,nk,rvec,ovp,kb,rlat,vector,tau,mfold,needr,needc,upper,nimgtot,t_states,t_table)

		use omp_lib

		implicit none

		integer,intent(in) :: n,nvec,nb,nk
		real,dimension(nvec,3),intent(in) :: rvec
		real,dimension(nvec,n,n),intent(in) :: ovp
		real,dimension(3,nk),intent(in) :: kb
		real,dimension(3,3),intent(in) :: rlat
		complex,dimension(n,nb,nk),intent(in) :: vector
		real,dimension(n,3),intent(in) :: tau
		! BSE_QFOLD: 1 for the images of bse_q_images, M > 1 for those of bse_q_images_fold
		integer,intent(in) :: mfold
		! k points of the rows (k1) and of the columns (k2) of the elements this rank evaluates
		logical,dimension(nk),intent(in) :: needr,needc
		! only the pairs k1 <= k2 are needed (upper triangle of a Hamiltonian ordered by k)
		logical,intent(in) :: upper
		integer,intent(out) :: nimgtot
		double precision,intent(out) :: t_states,t_table

		complex,allocatable,dimension(:,:,:) :: av,lv,rv,bv
		complex,allocatable,dimension(:,:) :: stb,tmp,tr,tb
		complex,allocatable,dimension(:) :: e1,g
		integer,allocatable,dimension(:) :: krow,kcol
		real,allocatable,dimension(:,:) :: gsh
		real,allocatable,dimension(:) :: wsh
		integer :: ik,ik1,ik2,s1,s2,n1,n2,img,nimg,j,b,nmax
		integer(kind=8) :: total
		double precision :: t0
		complex,parameter :: c0=(0.0,0.0),c1=(1.0,0.0),ch=(0.5,0.0)

		t0 = omp_get_wtime()

		call fk_free()

		! slots of the k points of this rank
		n1 = count(needr)
		n2 = count(needc)
		allocate(fk_rmap(nk),fk_cmap(nk),krow(max(1,n1)),kcol(max(1,n2)))
		fk_rmap = 0
		fk_cmap = 0
		s1 = 0
		s2 = 0
		do ik=1,nk
			if (needr(ik)) then
				s1 = s1 + 1
				fk_rmap(ik) = s1
				krow(s1) = ik
			end if
			if (needc(ik)) then
				s2 = s2 + 1
				fk_cmap(ik) = s2
				kcol(s2) = ik
			end if
		end do

		! images of the stored pairs: where the images of each pair start in fk_vertex
		nmax = 27*max(1,mfold)**2
		allocate(gsh(3,nmax),wsh(nmax))
		allocate(fk_first(max(1,n2),max(1,n1)))
		fk_first = -1
		total = 0
		do s1=1,n1
			do s2=1,n2
				ik1 = krow(s1)
				ik2 = kcol(s2)
				if (upper .and. ik1 > ik2) cycle
				call fk_images(kb(:,ik1),kb(:,ik2),rlat,mfold,nmax,nimg,gsh,wsh)
				fk_first(s2,s1) = int(total)
				total = total + nimg
				if (total >= huge(1)) stop 'fk_setup: more images than fit in a default integer; use more MPI ranks'
			end do
		end do
		nimgtot = int(total)
		deallocate(gsh,wsh)

		! per-state vectors of the k points of this rank, one k at a time (one S(k) in memory)
		allocate(av(n,nb,max(1,n1)),lv(n,nb,max(1,n1)),rv(n,nb,max(1,n2)),bv(n,nb,max(1,n2)))
		allocate(stb(n,n),tmp(n,nb),e1(n))
		do ik=1,nk
			if (.not. (needr(ik) .or. needc(ik))) cycle
			call overlap(n,nvec,rvec,ovp,kb(1,ik),kb(2,ik),kb(3,ik),stb)
			do j=1,n
				e1(j) = exp(cmplx(0.0,dot_product(kb(:,ik),tau(j,:))))
			end do
			if (needr(ik)) then
				s1 = fk_rmap(ik)
				! A: S^T conjg(l), then e(k)
				tmp = conjg(vector(:,:,ik))
				call cgemm('T','N',n,nb,n,c1,stb,n,tmp,n,c0,av(1,1,s1),n)
				do b=1,nb
					do j=1,n
						av(j,b,s1) = av(j,b,s1)*e1(j)
						lv(j,b,s1) = conjg(vector(j,b,ik))*e1(j)
					end do
				end do
			end if
			if (needc(ik)) then
				s2 = fk_cmap(ik)
				! B: S r, then conjg(e(k))
				call cgemm('N','N',n,nb,n,c1,stb,n,vector(1,1,ik),n,c0,bv(1,1,s2),n)
				do b=1,nb
					do j=1,n
						rv(j,b,s2) = conjg(e1(j))*vector(j,b,ik)
						bv(j,b,s2) = bv(j,b,s2)*conjg(e1(j))
					end do
				end do
			end if
		end do
		deallocate(stb,tmp,e1)
		t_states = omp_get_wtime() - t0

		! vertices of all band pairs for every stored (k1, k2, image)
		t0 = omp_get_wtime()
		allocate(fk_vertex(nb,nb,max(1,nimgtot)))
		fk_vertex = c0
		!$omp parallel default(shared) private(s1,s2,ik1,ik2,img,nimg,gsh,wsh,g,tr,tb,j,b)
		allocate(g(n),tr(n,nb),tb(n,nb),gsh(3,nmax),wsh(nmax))
		!$omp do collapse(2) schedule(dynamic)
		do s1=1,n1
			do s2=1,n2
				if (fk_first(s2,s1) >= 0) then
					ik1 = krow(s1)
					ik2 = kcol(s2)
					call fk_images(kb(:,ik1),kb(:,ik2),rlat,mfold,nmax,nimg,gsh,wsh)
					do img=1,nimg
						do j=1,n
							g(j) = exp(cmplx(0.0,-dot_product(gsh(:,img),tau(j,:))))
						end do
						do b=1,nb
							do j=1,n
								tr(j,b) = g(j)*rv(j,b,s2)
								tb(j,b) = g(j)*bv(j,b,s2)
							end do
						end do
						call cgemm('T','N',nb,nb,n,ch,av(1,1,s1),n,tr,n,c0,fk_vertex(1,1,fk_first(s2,s1)+img),nb)
						call cgemm('T','N',nb,nb,n,ch,lv(1,1,s1),n,tb,n,c1,fk_vertex(1,1,fk_first(s2,s1)+img),nb)
					end do
				end if
			end do
		end do
		!$omp end do
		deallocate(g,tr,tb,gsh,wsh)
		!$omp end parallel
		deallocate(av,lv,rv,bv,krow,kcol)
		t_table = omp_get_wtime() - t0

	end subroutine fk_setup

	! images of the pair (k1, k2): bse_q_images for mfold = 1, bse_q_images_fold for mfold > 1; shift and wimg hold nmax = 27 mfold^2
	subroutine fk_images(k1,k2,rlat,mfold,nmax,nimg,shift,wimg)

		implicit none

		real,dimension(3),intent(in) :: k1,k2
		real,dimension(3,3),intent(in) :: rlat
		integer,intent(in) :: mfold,nmax
		integer,intent(out) :: nimg
		real,dimension(3,nmax),intent(out) :: shift
		real,dimension(nmax),intent(out) :: wimg

		if (mfold > 1) then
			call bse_q_images_fold(k1,k2,rlat,mfold,nmax,nimg,shift,wimg)
		else
			! nmax = 27 here
			call bse_q_images(k1,k2,rlat,nimg,shift)
		end if

	end subroutine fk_images

	! third index of fk_vertex for the pair of k points (ik1 rows, ik2 columns) and the image img
	integer function fk_slot(ik1,ik2,img)

		implicit none

		integer,intent(in) :: ik1,ik2,img
		integer :: s1,s2

		s1 = fk_rmap(ik1)
		s2 = fk_cmap(ik2)
		if (s1 == 0 .or. s2 == 0) error stop 'fast kernel: k pair outside the vertex table of this MPI rank'
		if (fk_first(s2,s1) < 0) error stop 'fast kernel: k pair not stored in the vertex table'
		fk_slot = fk_first(s2,s1) + img

	end function fk_slot

	subroutine fk_free()

		implicit none

		if (allocated(fk_vertex)) deallocate(fk_vertex)
		if (allocated(fk_rmap)) deallocate(fk_rmap)
		if (allocated(fk_cmap)) deallocate(fk_cmap)
		if (allocated(fk_first)) deallocate(fk_first)
		fk_active = .false.

	end subroutine fk_free

end module bse_fastkernel
