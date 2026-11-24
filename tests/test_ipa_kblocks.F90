! The independent-particle stage shared between MPI ranks (IPA stage of bsesolver, IPA= T of spoptics):
!   ipa_kblock           the block of consecutive k points of a rank: covers 1..nk once, sizes differ by at most one;
!   optdip_s_kall        the residual position term of all (c,v) pairs of a k point equals optdip_s of every pair;
!   ipa_share_kblocks_*  (MPI build) every rank fills its own block and receives the others, with blocks of 3 elements
!                        so that the transfers are split, and with fewer k points than ranks (empty blocks).
! make -j1 test-ipa-kblocks runs it on one rank; with an MPI build, mpiexec -n 3 build/tests/test_ipa_kblocks.x
program test_ipa_kblocks

#ifdef MPI
	use mpi
#endif
	implicit none

	integer :: failures,total_failures,rank,nodes,ierr

	rank = 0
	nodes = 1
	ierr = 0
#ifdef MPI
	call MPI_INIT(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,nodes,ierr)
#endif

	failures = 0
	call check_partition(failures)
	call check_dipole_kall(failures)
#ifdef MPI
	call check_sharing(rank,nodes,failures)
#endif

	total_failures = failures
#ifdef MPI
	call MPI_ALLREDUCE(failures,total_failures,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
	if (rank == 0) then
		if (total_failures == 0) then
			write(*,*) 'ipa_kblocks: PASS on ',nodes,' rank(s)'
		else
			write(*,*) 'ipa_kblocks: FAIL; mismatches: ',total_failures
		end if
	end if
#ifdef MPI
	call MPI_FINALIZE(ierr)
#endif
	if (total_failures /= 0) stop 1

contains

	subroutine check_partition(failures)

		integer,intent(inout) :: failures
		integer,parameter :: nks(7) = (/0,1,2,5,7,36,100/)
		integer :: nranks,ik,nk,r,k1,k2,expected_first,smallest,largest

		do ik=1,size(nks)
			nk = nks(ik)
			do nranks=1,9
				expected_first = 1
				smallest = nk+1
				largest = -1
				do r=0,nranks-1
					call ipa_kblock(nranks,r,nk,k1,k2)
					if (k1 /= expected_first) failures = failures+1
					if (k2 < k1-1) failures = failures+1
					smallest = min(smallest,k2-k1+1)
					largest = max(largest,k2-k1+1)
					expected_first = k2+1
				end do
				if (expected_first /= nk+1) failures = failures+1
				if (largest-smallest > 1) failures = failures+1
			end do
		end do

	end subroutine check_partition


	subroutine check_dipole_kall(failures)

		integer,intent(inout) :: failures
		integer,parameter :: n=7,nvec=5,nv=3,nc=2
		real,parameter :: scs=0.4,sme=0.05
		real :: rvec(nvec,3),dmat(nvec,n,n,3),evraw(nv),ecraw(nc)
		real :: kx,ky,kz,delta,scale,diff,u(2)
		complex :: vv(n,nv),vc(n,nc)
		complex :: hrx(nc,nv),hry(nc,nv),hrz(nc,nv),refx(nc,nv),refy(nc,nv),refz(nc,nv)
		integer :: c,v,a,b,seed_size
		integer,allocatable :: seed(:)

		call random_seed(size=seed_size)
		allocate(seed(seed_size))
		seed = 20261124
		call random_seed(put=seed)
		deallocate(seed)

		call random_number(rvec)
		call random_number(dmat)
		rvec = 4.0*(rvec-0.5)
		dmat = dmat-0.5
		do b=1,nv
			do a=1,n
				call random_number(u)
				vv(a,b) = cmplx(u(1)-0.5,u(2)-0.5)
			end do
		end do
		do b=1,nc
			do a=1,n
				call random_number(u)
				vc(a,b) = cmplx(u(1)-0.5,u(2)-0.5)
			end do
		end do
		! valence energies below the conduction ones (the gap is the DFT one, ecraw has the scissor in it)
		do v=1,nv
			evraw(v) = -3.0+0.5*real(v)
		end do
		do c=1,nc
			ecraw(c) = 1.0+scs+0.7*real(c)
		end do
		kx = 0.3
		ky = -0.7
		kz = 0.1

		! a nonzero vertex to add to, as after optsp_s_kall
		do v=1,nv
			do c=1,nc
				call random_number(u)
				hrx(c,v) = cmplx(u(1),u(2))
				hry(c,v) = cmplx(u(2),-u(1))
				hrz(c,v) = cmplx(-u(1),u(2))
			end do
		end do
		refx = hrx
		refy = hry
		refz = hrz
		do v=1,nv
			do c=1,nc
				delta = ecraw(c)-scs-evraw(v)
				call optdip_s(delta,sme,vv(:,v),vc(:,c),kx,ky,kz,n,nvec,rvec,dmat,refx(c,v),refy(c,v),refz(c,v))
			end do
		end do

		call optdip_s_kall(nv,nc,evraw,ecraw,scs,sme,vv,vc,kx,ky,kz,n,nvec,rvec,dmat,hrx,hry,hrz)

		scale = 1.0
		do v=1,nv
			do c=1,nc
				scale = max(scale,abs(refx(c,v)),abs(refy(c,v)),abs(refz(c,v)))
			end do
		end do
		do v=1,nv
			do c=1,nc
				diff = max(abs(hrx(c,v)-refx(c,v)),abs(hry(c,v)-refy(c,v)),abs(hrz(c,v)-refz(c,v)))
				if (diff > 2.0e-5*scale) then
					failures = failures+1
					write(*,*) 'optdip_s_kall differs from optdip_s at c,v = ',c,v,': ',diff,' of ',scale
				end if
			end do
		end do
		! the term itself is not negligible, so that the check means something
		if (abs(refx(1,1)-cmplx(0.0,0.0)) < 1.0e-3) failures = failures+1

	end subroutine check_dipole_kall


#ifdef MPI
	subroutine check_sharing(rank,nodes,failures)

		integer,intent(in) :: rank,nodes
		integer,intent(inout) :: failures
		integer,parameter :: nks(5) = (/1,2,3,5,11/)
		integer,parameter :: per_k=4,block_elements=3
		integer :: ik,nk,k,j,k1,k2,ierr
		complex,allocatable :: zc(:,:)
		real,allocatable :: zr(:,:)
		integer,allocatable :: zi(:)

		do ik=1,size(nks)
			nk = nks(ik)
			allocate(zc(per_k,nk),zr(per_k,nk),zi(nk))
			zc = cmplx(-999.0,-999.0)
			zr = -999.0
			zi = -1
			call ipa_kblock(nodes,rank,nk,k1,k2)
			do k=k1,k2
				do j=1,per_k
					zc(j,k) = cmplx(real(10*k+j),-real(k))
					zr(j,k) = real(100*k+j)+0.5
				end do
				zi(k) = 1000+k
			end do
			call ipa_share_kblocks_complex(zc,per_k,nk,MPI_COMM_WORLD,ierr,block_elements)
			if (ierr /= MPI_SUCCESS) failures = failures+1
			call ipa_share_kblocks_real(zr,per_k,nk,MPI_COMM_WORLD,ierr,block_elements)
			if (ierr /= MPI_SUCCESS) failures = failures+1
			call ipa_share_kblocks_integer(zi,nk,MPI_COMM_WORLD,ierr)
			if (ierr /= MPI_SUCCESS) failures = failures+1
			do k=1,nk
				do j=1,per_k
					if (abs(zc(j,k)-cmplx(real(10*k+j),-real(k))) > 1.0e-6) failures = failures+1
					if (abs(zr(j,k)-(real(100*k+j)+0.5)) > 1.0e-6) failures = failures+1
				end do
				if (zi(k) /= 1000+k) failures = failures+1
			end do
			deallocate(zc,zr,zi)
		end do

	end subroutine check_sharing
#endif

end program test_ipa_kblocks
