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
	! band pairs of all k pairs, so a kernel element is two lookups. The overlap matrices are needed one k at a time.
	! The arithmetic is the one of sandwich_phase in a different order (single precision rounding only).

	implicit none

	logical :: fk_active = .false.
	! fk_vertex(b1, b2, image, ik2, ik1): band columns of vector(:, b, ik), image numbered as bse_q_images(k1, k2)
	complex,allocatable,dimension(:,:,:,:,:) :: fk_vertex

contains

	subroutine fk_setup(n,nvec,nb,nk,rvec,ovp,kb,rlat,vector,tau,maxnimg,t_states,t_table)

		use omp_lib

		implicit none

		integer,intent(in) :: n,nvec,nb,nk
		real,dimension(nvec,3),intent(in) :: rvec
		real,dimension(nvec,n,n),intent(in) :: ovp
		real,dimension(3,nk),intent(in) :: kb
		real,dimension(3,3),intent(in) :: rlat
		complex,dimension(n,nb,nk),intent(in) :: vector
		real,dimension(n,3),intent(in) :: tau
		integer,intent(out) :: maxnimg
		double precision,intent(out) :: t_states,t_table

		complex,allocatable,dimension(:,:,:) :: av,lv,rv,bv
		complex,allocatable,dimension(:,:) :: stb,tmp,tr,tb
		complex,allocatable,dimension(:) :: e1,g
		real,dimension(3,27) :: gsh
		integer :: ik,ik1,ik2,img,nimg,j,b
		double precision :: t0
		complex,parameter :: c0=(0.0,0.0),c1=(1.0,0.0),ch=(0.5,0.0)

		t0 = omp_get_wtime()

		! largest number of equally short images over all k pairs
		maxnimg = 1
		do ik1=1,nk
			do ik2=1,nk
				call bse_q_images(kb(:,ik1),kb(:,ik2),rlat,nimg,gsh)
				maxnimg = max(maxnimg,nimg)
			end do
		end do

		! per-state vectors, one k at a time (one S(k) in memory)
		allocate(av(n,nb,nk),lv(n,nb,nk),rv(n,nb,nk),bv(n,nb,nk))
		allocate(stb(n,n),tmp(n,nb),e1(n))
		do ik=1,nk
			call overlap(n,nvec,rvec,ovp,kb(1,ik),kb(2,ik),kb(3,ik),stb)
			do j=1,n
				e1(j) = exp(cmplx(0.0,dot_product(kb(:,ik),tau(j,:))))
			end do
			! A: S^T conjg(l), then e(k)
			tmp = conjg(vector(:,:,ik))
			call cgemm('T','N',n,nb,n,c1,stb,n,tmp,n,c0,av(1,1,ik),n)
			! B: S r, then conjg(e(k))
			call cgemm('N','N',n,nb,n,c1,stb,n,vector(1,1,ik),n,c0,bv(1,1,ik),n)
			do b=1,nb
				do j=1,n
					av(j,b,ik) = av(j,b,ik)*e1(j)
					lv(j,b,ik) = conjg(vector(j,b,ik))*e1(j)
					rv(j,b,ik) = conjg(e1(j))*vector(j,b,ik)
					bv(j,b,ik) = bv(j,b,ik)*conjg(e1(j))
				end do
			end do
		end do
		deallocate(stb,tmp,e1)
		t_states = omp_get_wtime() - t0

		! vertices of all band pairs for every (k1, k2, image)
		t0 = omp_get_wtime()
		if (allocated(fk_vertex)) deallocate(fk_vertex)
		allocate(fk_vertex(nb,nb,maxnimg,nk,nk))
		fk_vertex = c0
		!$omp parallel default(shared) private(ik1,ik2,img,nimg,gsh,g,tr,tb,j,b)
		allocate(g(n),tr(n,nb),tb(n,nb))
		!$omp do collapse(2) schedule(dynamic)
		do ik1=1,nk
			do ik2=1,nk
				call bse_q_images(kb(:,ik1),kb(:,ik2),rlat,nimg,gsh)
				do img=1,nimg
					do j=1,n
						g(j) = exp(cmplx(0.0,-dot_product(gsh(:,img),tau(j,:))))
					end do
					do b=1,nb
						do j=1,n
							tr(j,b) = g(j)*rv(j,b,ik2)
							tb(j,b) = g(j)*bv(j,b,ik2)
						end do
					end do
					call cgemm('T','N',nb,nb,n,ch,av(1,1,ik1),n,tr,n,c0,fk_vertex(1,1,img,ik2,ik1),nb)
					call cgemm('T','N',nb,nb,n,ch,lv(1,1,ik1),n,tb,n,c1,fk_vertex(1,1,img,ik2,ik1),nb)
				end do
			end do
		end do
		!$omp end do
		deallocate(g,tr,tb)
		!$omp end parallel
		deallocate(av,lv,rv,bv)
		t_table = omp_get_wtime() - t0

	end subroutine fk_setup

end module bse_fastkernel
