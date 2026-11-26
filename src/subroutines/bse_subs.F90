
function matrizelbse(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,rlat,est1,ec1,ev1,vbc1 &
	         ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,dft,nvec,rvec,sk,skp, &
	         use_center_phase,center_phase1,center_phase2, &
	         use_center_grad,w90dat,centers,ffdat) !funcao para calcular o elemento de matriz da matriz bse

	use bse_q_optics, only: rmn_data, rmn_bloch, center_phase_direct_vertices, center_phase_gradient_correction
	use bse_formfactor, only: ff_data, ff_direct_vertices, ff_exchange_element
	use bse_fastkernel, only: fk_active, fk_vertex, fk_slot, fk_qc_valid, fk_qc_nimg, fk_qc_k1, fk_qc_k2, fk_qc_ediel, &
	                          fk_qc_tolr, fk_qc_coul, fk_qc_shift, fk_qc_wimg, fk_qc_w
	use input_variables, only: bse_qfold

	implicit none

	character(len=10) :: coultype
	character(len=1) :: dft

	integer,dimension(3) :: ngrid

	integer :: w90basis,nvec
	real,dimension(w90basis,3) :: tau
	logical :: use_center_phase
	real :: a,vcell1
	real :: ez,w
	
	real,dimension(nvec,3) :: rvec
	!real,dimension(nvec,w90basis,w90basis) :: ovp
	
	complex,dimension(w90basis,w90basis) :: sk,skp
	complex,dimension(w90basis) :: center_phase1,center_phase2
	logical :: use_center_grad
	type(rmn_data),intent(in) :: w90dat
	type(ff_data),intent(in) :: ffdat

	integer, dimension(4) :: est1,est2
	real :: ec1,ec2,ev1,ev2
	real,dimension(3) :: kpt1,kpt2
	complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2

	complex, dimension(w90basis) :: vbc,vbv

	real,dimension(3,3) :: rlat

	real :: tolr
	integer :: ktol
	
	complex:: matrizelbse

	real,parameter:: pi=acos(-1.)

	real :: auxi

	real :: modk

	real,dimension(3) :: ediel

	real :: lc

	complex :: vc,vv,delta_k
	real,dimension(3,w90basis),intent(in) :: centers
	integer :: nimg,img,m,n,nmax
	real,dimension(3,27*max(1,bse_qfold)**2) :: gshift
	real,dimension(27*max(1,bse_qfold)**2) :: wimg
	logical :: qhit
	real,dimension(3) :: kpt2i
	complex,dimension(w90basis) :: phase2i
	complex :: melem
	real :: ang

	complex :: a_od_p(3,w90basis,w90basis),a_od_m(3,w90basis,w90basis)
	logical :: bloch_ok
	character(len=256) :: bloch_msg
	integer :: idiag

	real :: vcoul1

	real :: vcoul,v2dk,v3diel,v3davg,v2dt,v2dtavg,v0dt,v2dt2
	real :: v2dohono,v2drk,v1dt,v2d,v2diel,v2dlavg,v2davgtpv,v2dtpv
	real :: v1d,v1diel

	real :: r0


	! q = kpt1-kpt2 enters through V(q) and the Wannier-centre phases,
	! neither periodic in q: use the shortest images q+G, averaged over the
	! equally short ones on the zone boundary, so that the kernel keeps the
	! lattice symmetry and does not depend on how the k-grid cell is chosen.
	! BSE_QFOLD = M > 1 (an M x M supercell run): the images of the M^2 classes q+g+Gamma_UC, weight 1 each (bse_q_images_fold),
	! with W of every image kept for the k pair, which the nc nv transitions of the pair share
	nmax=27*max(1,bse_qfold)**2
	qhit=.false.
	if (bse_qfold .gt. 1) then
		if (.not. allocated(fk_qc_w)) allocate(fk_qc_shift(3,nmax),fk_qc_wimg(nmax),fk_qc_w(nmax))
		qhit=fk_qc_valid .and. all(fk_qc_k1 .eq. kpt1) .and. all(fk_qc_k2 .eq. kpt2) .and. all(fk_qc_ediel .eq. ediel) &
		     .and. fk_qc_tolr .eq. tolr .and. fk_qc_coul .eq. coultype
		if (qhit) then
			nimg=fk_qc_nimg
			gshift(:,1:nimg)=fk_qc_shift(:,1:nimg)
			wimg(1:nimg)=fk_qc_wimg(1:nimg)
		else
			call bse_q_images_fold(kpt1,kpt2,rlat,bse_qfold,nmax,nimg,gshift,wimg)
		end if
	else
		call bse_q_images(kpt1,kpt2,rlat,nimg,gshift)
	end if
	matrizelbse=cmplx(0.0,0.0)
	do img=1,nimg

	kpt2i=kpt2+gshift(:,img)

	if (qhit) then
		vcoul1=fk_qc_w(img)
	else
	select case (coultype)

	case("V2DK")

		vcoul1= v2dk(kpt1,kpt2i,ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vcoul1= vcoul(kpt1,kpt2i,rlat,ngrid,tolr)

	case("V3DL")

		vcoul1= v3diel(kpt1,kpt2i,ediel,rlat,ngrid,tolr)

	case("V3DAVG")

		vcoul1= v3davg(kpt1,kpt2i,1.0,rlat,ngrid,tolr)

	case("V3DLAVG")

		vcoul1= v3davg(kpt1,kpt2i,ediel(2),rlat,ngrid,tolr)
		
	case("V2D")

		vcoul1= v2d(kpt1,kpt2i,rlat,ngrid,tolr)

	case("V2DL")

		vcoul1= v2diel(kpt1,kpt2i,ediel,rlat,ngrid,tolr)		

	case("V2DLAVG")

		vcoul1= v2dlavg(kpt1,kpt2i,ediel,ngrid,rlat,tolr)

	case("V2DT")

		vcoul1= v2dt(kpt1,kpt2i,ngrid,rlat,tolr)

	case("V2DTAVG")

		vcoul1= v2dtavg(kpt1,kpt2i,ediel,ngrid,rlat,tolr)

	case("V2DAVGTPV")

		vcoul1= v2davgtpv(kpt1,kpt2i,ediel,ngrid,rlat,tolr)

	case("V2DTPV")

		vcoul1= v2dtpv(kpt1,kpt2i,ediel,ngrid,rlat,tolr)

	case("V2DT2")

		vcoul1= v2dt2(kpt1,kpt2i,ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vcoul1= v2dohono(kpt1,kpt2i,ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vcoul1= v2drk(kpt1,kpt2i,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		
	case("V1D")

		vcoul1= v1d(kpt1,kpt2i,ngrid,rlat,lc,tolr)

	case("V1DL")

		vcoul1= v1diel(kpt1,kpt2i,ngrid,rlat,lc,ediel,tolr)		
		
	case("V1DT")
	
		vcoul1= v1dt(kpt1,kpt2i,ngrid,rlat,tolr,lc)	
				

	case("V0DT")

		vcoul1= v0dt(kpt1,kpt2i,ngrid,rlat,tolr)

	case default

		write(*,*) "Wrong Coulomb Potential"
		STOP

	end select
	if (bse_qfold .gt. 1) fk_qc_w(img)=vcoul1
	end if

	



	! same transition: <c|c> = <v|v> = 1 at q = 0; with BSE_QFOLD > 1 the other classes (q = -g != 0) go through the
	! vertices below, which vanish for a crystal that has the unit-cell periodicity
	if (est1(1) .eq. est2(1) .and. (bse_qfold .le. 1 .or. all(gshift(:,img) .eq. 0.0))) then


		melem= (ec1-ev1) + vcoul1


	else

		if (ffdat%direct) then

		! BSE_FF: the pair densities of the orbitals of the basis,
		! <c1|exp(iq.r)|c2> <v2|exp(-iq.r)|v1> at this image q, for DFT=W
		! and DFT=S alike (F(R; 0) is the overlap)
		call ff_direct_vertices(ffdat,kpt1-kpt2i,kpt2,est2(4),rlat, &
			vbc1,vbc2,vbv1,vbv2,vc,vv)
		melem=vcoul1*vc*vv

		else

			
		select case (dft)
		
		case ("S")
		

		 !call overlap(w90basis,nvec,rvec,ovp,kpt1(1),kpt1(2),kpt1(3),sk)
		 !call overlap(w90basis,nvec,rvec,ovp,kpt2(1),kpt2(2),kpt2(3),skp)
		 
		 ! atomic gauge: the orbital centres enter as exp(i q.tau), q the
		 ! transfer of this image, as the Wannier centres do below; without
		 ! them C3 is broken and the E' exciton doublet of h-BN splits
		 if (fk_active) then
			! the same two sandwiches from the table of bse_fastkernel (band columns est(3) = c, est(2) = v)
			vc = fk_vertex(est1(3),est2(3),fk_slot(est1(4),est2(4),img))
			vv = fk_vertex(est1(2),est2(2),fk_slot(est1(4),est2(4),img))
		 else
		 call orbital_centres(w90basis,tau)
		 call sandwich_phase(w90basis,vbc1,sk,skp,vbc2,kpt1-kpt2i,tau,vc)
		 ! hole vertex <v k2|v k1> = conjg(<v k1|v k2>)
		 call sandwich_phase(w90basis,vbv1,sk,skp,vbv2,kpt1-kpt2i,tau,vv)
		 end if

		 melem=  vcoul1*vc*conjg(vv)		
		
		case default
	
		 if (use_center_phase) then

			! the image kpt2+G carries the phases exp(i(kpt2+G).t_m)
			do m=1,w90basis
				ang=dot_product(gshift(:,img),centers(:,m))
				phase2i(m)=center_phase2(m)*cmplx(cos(ang),sin(ang))
			end do
			call center_phase_direct_vertices(vbc1,vbc2,vbv1,vbv2, &
				center_phase1,phase2i,vc,vv)

			if (use_center_grad) then
				! first order of <m,0|exp(iq.r)|n,R> about the pair midpoint
				! (t_m+t_n+R)/2: exp(iq.(t_m+t_n)/2) iq.[A(kb)-diag(t)]_mn, with
				! A(k) = sum_R exp(ik.R) r(R) at the states' mean kb = (kpt1+kpt2i)/2
				! (the Bloch phases of R, not q), the diagonal less the centres, and
				! the phase conjugated for the hole vertex
				call rmn_bloch(w90dat,0.5*(kpt1+kpt2i),rlat,a_od_p,bloch_ok,bloch_msg)
				do idiag=1,w90basis
					a_od_p(:,idiag,idiag)=a_od_p(:,idiag,idiag)-centers(:,idiag)
				end do
				do n=1,w90basis
					do m=1,w90basis
						ang=0.5*dot_product(kpt1-kpt2i,centers(:,m)+centers(:,n))
						a_od_m(:,m,n)=a_od_p(:,m,n)*cmplx(cos(ang),-sin(ang))
						a_od_p(:,m,n)=a_od_p(:,m,n)*cmplx(cos(ang),sin(ang))
					end do
				end do
				call center_phase_gradient_correction(vbc1,vbc2,vbv1,vbv2, &
					a_od_p,a_od_m,w90basis,kpt1-kpt2i,vc,vv,delta_k)
				melem=vcoul1*(vc*vv+delta_k)
			else
				melem=vcoul1*vc*vv
			end if

		 else

			call vecconjg(vbc1,w90basis,vbc)

			call vecconjg(vbv2,w90basis,vbv)

			call prodintsq(vbc,vbc2,w90basis,vc)

			call prodintsq(vbv,vbv1,w90basis,vv)

			melem=vcoul1*vc*vv

		 end if
	         
		end select

		end if
	
		


	end if

	if (bse_qfold .gt. 1) then
		matrizelbse=matrizelbse+wimg(img)*melem
	else
		matrizelbse=matrizelbse+melem
	end if

	end do

	if (bse_qfold .gt. 1) then
		if (.not. qhit) then
			fk_qc_k1=kpt1
			fk_qc_k2=kpt2
			fk_qc_ediel=ediel
			fk_qc_tolr=tolr
			fk_qc_coul=coultype
			fk_qc_nimg=nimg
			fk_qc_shift(:,1:nimg)=gshift(:,1:nimg)
			fk_qc_wimg(1:nimg)=wimg(1:nimg)
			fk_qc_valid=.true.
		end if
	else
		matrizelbse=matrizelbse/real(nimg)
	end if
	! BSE_FF_EXCHANGE: the exchange term, independent of q
	if (ffdat%exchange) matrizelbse=matrizelbse+ff_exchange_element(ffdat,est1(1),est2(1))

end function matrizelbse




subroutine opticalactivity(dimse,excitonvec,hopt,activity,description)
! <0|r|S> = sum_j conjg(A_j^S) <c|r|v>_j, consistent with the <v2|v1> kernel

	implicit none


	integer :: dimse
	real,dimension(dimse) :: activity
	complex :: actaux
	real :: description(dimse)
	complex,dimension(dimse) :: hopt
	complex,dimension(dimse,dimse) :: excitonvec

	integer :: i,j
	

	activity=0.
	actaux=cmplx(0.,0.)

	do i=1,dimse

		do j=1,dimse

		
		actaux=actaux+(conjg(excitonvec(j,i))*hopt(j))

		end do

		activity(i)=actaux*conjg(actaux)

		if (activity(i) .gt. 0.1) then

			description(i)= 1.

		else if ( (activity(i) .lt. 0.1) .and. (activity(i) .gt. 1.0e-8)) then

			description(i)= 0.

		else

			description(i)= -1.


		end if
	
		actaux=cmplx(0.,0.)


	end do



end subroutine opticalactivity

subroutine excitonil(i,w90basis,nc,nv,ncount,stt,ndim,wf,lc,pinter,pintra)

	implicit none

	integer :: i !numero do estado excitonico

	integer :: j,k,nc,nv

	integer :: w90basis
	integer :: ncount,ndim
	integer,dimension(ncount*nc*nv,4) :: stt

	complex,dimension(ndim) :: wf

	real,dimension(ncount,2*w90basis) :: lc

	real :: pinter,pintra

	real :: inter,intra

	real :: p1c,p1v,p2c,p2v


	pinter = 0d0
	pintra = 0d0

	do j=1,ndim

		p1c = 1d0-lc(stt(j,4),stt(j,3))
		p1v = 1d0-lc(stt(j,4),stt(j,2))

		p2c = lc(stt(j,4),stt(j,3))
		p2v = lc(stt(j,4),stt(j,2))


		intra = (p1c*p1v)+(p2c*p2v)

		inter = (p1c*p2v)+(p2c*p1v)


		pinter = pinter + Real(conjg(wf(j))*wf(j)*inter)
 
		pintra = pintra + Real(conjg(wf(j))*wf(j)*intra)

	end do


end subroutine excitonil



subroutine dielbse(dimse,excitonvec,hopt1,hopt2,activity)
! <0|r|S> = sum_j conjg(A_j^S) <c|r|v>_j, consistent with the <v2|v1> kernel

	implicit none


	integer :: dimse
	real,dimension(dimse) :: activity, exciton
	real :: actaux
	complex,dimension(dimse) :: hopt1,hopt2
	complex,dimension(dimse,dimse) :: excitonvec

	integer :: i,j,k
	

	activity=0.0
	actaux=0.0

	do i=1,dimse



		do j=1,dimse
			do k=1,dimse


		actaux=actaux+(conjg(excitonvec(j,i))*hopt1(j)*excitonvec(k,i)*conjg(hopt2(k)))

		end do
			end do

		activity(i)=actaux

	
	
		actaux=0.0


	end do



end subroutine dielbse



subroutine dielbsep(nthread,dimse,excitonvec,hopt1,hopt2,activity)
! <0|r|S> = sum_j conjg(A_j^S) <c|r|v>_j, consistent with the <v2|v1> kernel

	use omp_lib
	implicit none


	integer :: dimse,nthread
	real,dimension(dimse) :: activity, exciton
	complex :: actaux,actaux2
	complex,dimension(dimse) :: hopt1,hopt2
	complex,dimension(dimse,dimse) :: excitonvec
	integer,dimension(dimse,4) :: stt

	integer :: i,j,k,no
	

	call OMP_SET_NUM_THREADS(nthread)

	activity=0.0
	!actaux=0.0


	! Each exciton state is independent.  This is one of the dominant
	! post-processing costs after the BSE diagonalization, so parallelize over
	! the output-state index rather than leaving the OpenMP directives commented.
	!$OMP PARALLEL DO DEFAULT(NONE) &
	!$OMP& SHARED(dimse, excitonvec, hopt1, hopt2, activity) &
	!$OMP& PRIVATE(i, j, actaux, actaux2) SCHEDULE(STATIC)
	do i=1,dimse


		actaux=0.0
		actaux2=0.0


		do j=1,dimse

			
			actaux=actaux+(conjg(excitonvec(j,i))*hopt1(j))
			actaux2=actaux2+(conjg(excitonvec(j,i))*hopt2(j))
			

			
		end do
			

		activity(i)=actaux*conjg(actaux2)

	
	
		


	end do
	!$OMP END PARALLEL DO



end subroutine dielbsep



subroutine dielbsev(nthread,dimse,excitonvec,hopt,activity)
! <0|r|S> = sum_j conjg(A_j^S) <c|r|v>_j, consistent with the <v2|v1> kernel

	use omp_lib
	implicit none


	integer :: dimse,nthread
	real,dimension(dimse) :: activity, exciton
	complex :: actaux
	complex,dimension(dimse) :: hopt
	complex,dimension(dimse,dimse) :: excitonvec

	integer :: i,j,k

	call OMP_SET_NUM_THREADS(nthread)	

	activity=0.0
	actaux=0.0

	!$OMP PARALLEL DO DEFAULT(NONE) &
	!$OMP& SHARED(dimse, excitonvec, hopt, activity) &
	!$OMP& PRIVATE(i, j, actaux) SCHEDULE(STATIC)
	do i=1,dimse

		actaux=0.0

		do j=1,dimse
			


		actaux=actaux+(conjg(excitonvec(j,i))*hopt(j))
	

		end do
		

		activity(i)=actaux*conjg(actaux)

	


	end do
	!$OMP END PARALLEL DO


end subroutine dielbsev



subroutine dielbsev_dist(dimse,excitonvec,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hopt,activity,MPIError)
! <0|r|S> = sum_j conjg(A_j^S) <c|r|v>_j, consistent with the <v2|v1> kernel

#ifdef MPI
	use mpi
#endif
	implicit none

	integer :: dimse,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,MPIError
	integer,parameter :: bse_mpi_block_elements=16777216
	integer :: li,lj,ig,jg
	integer, external :: indxl2g
	complex,dimension(lld,*) :: excitonvec
	complex,dimension(dimse) :: hopt
	real,dimension(dimse) :: activity
	complex,allocatable,dimension(:) :: amplitude

#ifdef MPI
	allocate(amplitude(dimse))
	amplitude = cmplx(0.0,0.0)

	do lj=1,locc
		jg = indxl2g(lj,nb,mycol,0,npcol)
		if (jg <= dimse) then
			do li=1,locr
				ig = indxl2g(li,mb,myrow,0,nprow)
				if (ig <= dimse) amplitude(jg) = amplitude(jg) + conjg(excitonvec(li,lj))*hopt(ig)
			end do
		end if
	end do

	call bse_mpi_allreduce_complex_sum_blocks(amplitude,int(dimse,kind=8),MPI_COMM_WORLD,MPIError,bse_mpi_block_elements)
	if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('distributed BSE optical amplitude',MPIError)
	activity = real(amplitude*conjg(amplitude))
	deallocate(amplitude)
#else
	activity = 0.0
	stop "Distributed BSE postprocessing requires MPI"
#endif

end subroutine dielbsev_dist


subroutine dielbsep_dist(dimse,excitonvec,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hopt1,hopt2,activity,MPIError)
! <0|r|S> = sum_j conjg(A_j^S) <c|r|v>_j, consistent with the <v2|v1> kernel

#ifdef MPI
	use mpi
#endif
	implicit none

	integer :: dimse,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,MPIError
	integer,parameter :: bse_mpi_block_elements=16777216
	integer :: li,lj,ig,jg
	integer, external :: indxl2g
	complex,dimension(lld,*) :: excitonvec
	complex,dimension(dimse) :: hopt1,hopt2
	real,dimension(dimse) :: activity
	complex,allocatable,dimension(:) :: amplitude1,amplitude2

#ifdef MPI
	allocate(amplitude1(dimse),amplitude2(dimse))
	amplitude1 = cmplx(0.0,0.0)
	amplitude2 = cmplx(0.0,0.0)

	do lj=1,locc
		jg = indxl2g(lj,nb,mycol,0,npcol)
		if (jg <= dimse) then
			do li=1,locr
				ig = indxl2g(li,mb,myrow,0,nprow)
				if (ig <= dimse) then
					amplitude1(jg) = amplitude1(jg) + conjg(excitonvec(li,lj))*hopt1(ig)
					amplitude2(jg) = amplitude2(jg) + conjg(excitonvec(li,lj))*hopt2(ig)
				end if
			end do
		end if
	end do

	call bse_mpi_allreduce_complex_sum_blocks(amplitude1,int(dimse,kind=8),MPI_COMM_WORLD,MPIError,bse_mpi_block_elements)
	if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('distributed BSE optical amplitude',MPIError)
	call bse_mpi_allreduce_complex_sum_blocks(amplitude2,int(dimse,kind=8),MPI_COMM_WORLD,MPIError,bse_mpi_block_elements)
	if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('distributed BSE optical amplitude',MPIError)
	activity = real(amplitude1*conjg(amplitude2))
	deallocate(amplitude1,amplitude2)
#else
	activity = 0.0
	stop "Distributed BSE postprocessing requires MPI"
#endif

end subroutine dielbsep_dist


subroutine bse_eigenvector_column_dist(dimse,excitonvec,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,icol,eigvec,MPIError)

#ifdef MPI
	use mpi
#endif
	implicit none

	integer :: dimse,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,icol,MPIError
	integer,parameter :: bse_mpi_block_elements=16777216
	integer :: li,lj,ig,jg
	integer, external :: indxl2g
	complex,dimension(lld,*) :: excitonvec
	complex,dimension(dimse) :: eigvec

#ifdef MPI
	eigvec = cmplx(0.0,0.0)

	do lj=1,locc
		jg = indxl2g(lj,nb,mycol,0,npcol)
		if (jg == icol) then
			do li=1,locr
				ig = indxl2g(li,mb,myrow,0,nprow)
				if (ig <= dimse) eigvec(ig) = excitonvec(li,lj)
			end do
		end if
	end do

	call bse_mpi_allreduce_complex_sum_blocks(eigvec,int(dimse,kind=8),MPI_COMM_WORLD,MPIError,bse_mpi_block_elements)
	if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('distributed BSE eigenvector',MPIError)
#else
	eigvec = cmplx(0.0,0.0)
	stop "Distributed BSE postprocessing requires MPI"
#endif

end subroutine bse_eigenvector_column_dist


! ---- BSE_REAL= T: the same post-processing for the real eigenvectors of the real basis (bse_realbasis) ----
! The optical amplitude is sum_r a_r h'_r, with h' = U^+ h (rb_vector) given by its real and imaginary parts.

subroutine dielbsev_real(nthread,dimse,excitonvec,hre,him,activity)

	use omp_lib
	implicit none

	integer :: dimse,nthread
	real,dimension(dimse) :: activity
	real,dimension(dimse) :: hre,him
	real,dimension(dimse,dimse) :: excitonvec

	integer :: i,j
	double precision :: ar,ai

	call OMP_SET_NUM_THREADS(nthread)

	activity=0.0

	!$OMP PARALLEL DO DEFAULT(NONE) &
	!$OMP& SHARED(dimse, excitonvec, hre, him, activity) &
	!$OMP& PRIVATE(i, j, ar, ai) SCHEDULE(STATIC)
	do i=1,dimse
		ar=0.0d0
		ai=0.0d0
		do j=1,dimse
			ar=ar+dble(excitonvec(j,i))*dble(hre(j))
			ai=ai+dble(excitonvec(j,i))*dble(him(j))
		end do
		activity(i)=real(ar*ar+ai*ai)
	end do
	!$OMP END PARALLEL DO

end subroutine dielbsev_real


subroutine dielbsep_real(nthread,dimse,excitonvec,h1re,h1im,h2re,h2im,activity)
! Re( amplitude1 * conjg(amplitude2) )

	use omp_lib
	implicit none

	integer :: dimse,nthread
	real,dimension(dimse) :: activity
	real,dimension(dimse) :: h1re,h1im,h2re,h2im
	real,dimension(dimse,dimse) :: excitonvec

	integer :: i,j
	double precision :: a1r,a1i,a2r,a2i

	call OMP_SET_NUM_THREADS(nthread)

	activity=0.0

	!$OMP PARALLEL DO DEFAULT(NONE) &
	!$OMP& SHARED(dimse, excitonvec, h1re, h1im, h2re, h2im, activity) &
	!$OMP& PRIVATE(i, j, a1r, a1i, a2r, a2i) SCHEDULE(STATIC)
	do i=1,dimse
		a1r=0.0d0
		a1i=0.0d0
		a2r=0.0d0
		a2i=0.0d0
		do j=1,dimse
			a1r=a1r+dble(excitonvec(j,i))*dble(h1re(j))
			a1i=a1i+dble(excitonvec(j,i))*dble(h1im(j))
			a2r=a2r+dble(excitonvec(j,i))*dble(h2re(j))
			a2i=a2i+dble(excitonvec(j,i))*dble(h2im(j))
		end do
		activity(i)=real(a1r*a2r+a1i*a2i)
	end do
	!$OMP END PARALLEL DO

end subroutine dielbsep_real


subroutine dielbsev_dist_real(dimse,excitonvec,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,hre,him,activity,MPIError)

#ifdef MPI
	use mpi
#endif
	implicit none

	integer :: dimse,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,MPIError
	integer :: li,lj,ig,jg
	integer, external :: indxl2g
	real,dimension(lld,*) :: excitonvec
	real,dimension(dimse) :: hre,him
	real,dimension(dimse) :: activity
	double precision,allocatable,dimension(:,:) :: amplitude

#ifdef MPI
	allocate(amplitude(dimse,2))
	amplitude = 0.0d0

	do lj=1,locc
		jg = indxl2g(lj,nb,mycol,0,npcol)
		if (jg <= dimse) then
			do li=1,locr
				ig = indxl2g(li,mb,myrow,0,nprow)
				if (ig <= dimse) then
					amplitude(jg,1) = amplitude(jg,1) + dble(excitonvec(li,lj))*dble(hre(ig))
					amplitude(jg,2) = amplitude(jg,2) + dble(excitonvec(li,lj))*dble(him(ig))
				end if
			end do
		end if
	end do

	call MPI_ALLREDUCE(MPI_IN_PLACE,amplitude,2*dimse,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIError)
	if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('distributed BSE optical amplitude',MPIError)
	activity = real(amplitude(:,1)*amplitude(:,1)+amplitude(:,2)*amplitude(:,2))
	deallocate(amplitude)
#else
	activity = 0.0
	stop "Distributed BSE postprocessing requires MPI"
#endif

end subroutine dielbsev_dist_real


subroutine dielbsep_dist_real(dimse,excitonvec,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,h1re,h1im,h2re,h2im,activity,MPIError)

#ifdef MPI
	use mpi
#endif
	implicit none

	integer :: dimse,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,MPIError
	integer :: li,lj,ig,jg
	integer, external :: indxl2g
	real,dimension(lld,*) :: excitonvec
	real,dimension(dimse) :: h1re,h1im,h2re,h2im
	real,dimension(dimse) :: activity
	double precision,allocatable,dimension(:,:) :: amplitude

#ifdef MPI
	allocate(amplitude(dimse,4))
	amplitude = 0.0d0

	do lj=1,locc
		jg = indxl2g(lj,nb,mycol,0,npcol)
		if (jg <= dimse) then
			do li=1,locr
				ig = indxl2g(li,mb,myrow,0,nprow)
				if (ig <= dimse) then
					amplitude(jg,1) = amplitude(jg,1) + dble(excitonvec(li,lj))*dble(h1re(ig))
					amplitude(jg,2) = amplitude(jg,2) + dble(excitonvec(li,lj))*dble(h1im(ig))
					amplitude(jg,3) = amplitude(jg,3) + dble(excitonvec(li,lj))*dble(h2re(ig))
					amplitude(jg,4) = amplitude(jg,4) + dble(excitonvec(li,lj))*dble(h2im(ig))
				end if
			end do
		end if
	end do

	call MPI_ALLREDUCE(MPI_IN_PLACE,amplitude,4*dimse,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIError)
	if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('distributed BSE optical amplitude',MPIError)
	activity = real(amplitude(:,1)*amplitude(:,3)+amplitude(:,2)*amplitude(:,4))
	deallocate(amplitude)
#else
	activity = 0.0
	stop "Distributed BSE postprocessing requires MPI"
#endif

end subroutine dielbsep_dist_real


subroutine bse_eigenvector_column_dist_real(dimse,excitonvec,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,icol,eigvec,MPIError)

#ifdef MPI
	use mpi
#endif
	implicit none

	integer :: dimse,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,icol,MPIError
	integer :: li,lj,ig,jg
	integer, external :: indxl2g
	real,dimension(lld,*) :: excitonvec
	real,dimension(dimse) :: eigvec

#ifdef MPI
	eigvec = 0.0

	do lj=1,locc
		jg = indxl2g(lj,nb,mycol,0,npcol)
		if (jg == icol) then
			do li=1,locr
				ig = indxl2g(li,mb,myrow,0,nprow)
				if (ig <= dimse) eigvec(ig) = excitonvec(li,lj)
			end do
		end if
	end do

	call MPI_ALLREDUCE(MPI_IN_PLACE,eigvec,dimse,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,MPIError)
	if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('distributed BSE eigenvector',MPIError)
#else
	eigvec = 0.0
	stop "Distributed BSE postprocessing requires MPI"
#endif

end subroutine bse_eigenvector_column_dist_real


subroutine bse_elpa_bounds_real(dimse,excitonvec,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,emin,emax,etrace,MPIError)
! bse_elpa_bounds for the real symmetric matrix of BSE_REAL= T

#ifdef MPI
	use mpi
#endif
	implicit none

	integer :: dimse,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,MPIError
	integer :: li,lj,ig,jg
	integer, external :: indxl2g
	real,dimension(lld,*) :: excitonvec
	real(kind=8) :: emin,emax,etrace
	real(kind=8),allocatable,dimension(:,:) :: acc

	emin = 0.0d0
	emax = 0.0d0
	etrace = 0.0d0
#ifdef MPI
	allocate(acc(dimse,2))
	acc = 0.0d0

	do lj=1,locc
		jg = indxl2g(lj,nb,mycol,0,npcol)
		do li=1,locr
			ig = indxl2g(li,mb,myrow,0,nprow)
			if (ig <= dimse .and. jg <= dimse) then
				if (ig == jg) then
					acc(ig,2) = real(excitonvec(li,lj),kind=8)
				else
					acc(ig,1) = acc(ig,1) + abs(real(excitonvec(li,lj),kind=8))
				end if
			end if
		end do
	end do

	call MPI_ALLREDUCE(MPI_IN_PLACE,acc,2*dimse,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIError)
	if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('BSE Hamiltonian bounds',MPIError)

	emin = minval(acc(:,2)-acc(:,1))
	emax = maxval(acc(:,2)+acc(:,1))
	etrace = sum(acc(:,2))
	deallocate(acc)
#else
	stop "Distributed BSE postprocessing requires MPI"
#endif

end subroutine bse_elpa_bounds_real


subroutine bse_elpa_bounds(dimse,excitonvec,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,emin,emax,etrace,MPIError)
! Gershgorin bounds [emin, emax] and trace of the distributed Hermitian BSE matrix (both triangles stored),
! taken before ELPA overwrites it: bse_elpa_verify uses them to catch a solver that returns wrong eigenvalues

#ifdef MPI
	use mpi
#endif
	implicit none

	integer :: dimse,lld,locr,locc,mb,nb,myrow,mycol,nprow,npcol,MPIError
	integer :: li,lj,ig,jg
	integer, external :: indxl2g
	complex,dimension(lld,*) :: excitonvec
	real(kind=8) :: emin,emax,etrace
	real(kind=8),allocatable,dimension(:,:) :: acc

	emin = 0.0d0
	emax = 0.0d0
	etrace = 0.0d0
#ifdef MPI
	! acc(:,1): sum of |H_ij|, j /= i, over the row i; acc(:,2): H_ii
	allocate(acc(dimse,2))
	acc = 0.0d0

	do lj=1,locc
		jg = indxl2g(lj,nb,mycol,0,npcol)
		do li=1,locr
			ig = indxl2g(li,mb,myrow,0,nprow)
			if (ig <= dimse .and. jg <= dimse) then
				if (ig == jg) then
					acc(ig,2) = real(excitonvec(li,lj),kind=8)
				else
					acc(ig,1) = acc(ig,1) + real(abs(excitonvec(li,lj)),kind=8)
				end if
			end if
		end do
	end do

	call MPI_ALLREDUCE(MPI_IN_PLACE,acc,2*dimse,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPIError)
	if (MPIError /= MPI_SUCCESS) call bse_mpi_collective_abort('BSE Hamiltonian bounds',MPIError)

	emin = minval(acc(:,2)-acc(:,1))
	emax = maxval(acc(:,2)+acc(:,1))
	etrace = sum(acc(:,2))
	deallocate(acc)
#else
	stop "Distributed BSE postprocessing requires MPI"
#endif

end subroutine bse_elpa_bounds


subroutine bse_elpa_verify(dimse,nstates,w,emin,emax,etrace,iam,ok)
! The ELPA result must lie inside the Gershgorin bounds of the Hamiltonian and, with every state computed,
! add up to its trace.  Only a gross failure is caught (a library whose 2-stage solver was miscompiled returns
! eigenvalues of the order of 1e4 times the spectrum, with a zero status), a small error is not.

	implicit none

	integer :: dimse,nstates,iam
	real,dimension(dimse) :: w
	real(kind=8) :: emin,emax,etrace
	logical :: ok
	real(kind=8) :: tol,esum,eabs

	ok = .true.
	tol = 1.0d-3*(emax-emin) + 1.0d-3
	if (.not. all(real(w(1:nstates),kind=8) >= emin-tol .and. real(w(1:nstates),kind=8) <= emax+tol)) then
		ok = .false.
		if (iam == 0) then
			write(*,*) 'ELPA returned eigenvalues outside the bounds of the BSE Hamiltonian (Gershgorin): [',emin,',',emax,'] eV'
			write(*,*) 'lowest and highest returned:',minval(w(1:nstates)),maxval(w(1:nstates))
		end if
	end if
	if (ok .and. nstates == dimse) then
		esum = sum(real(w(1:dimse),kind=8))
		eabs = sum(abs(real(w(1:dimse),kind=8)))
		if (.not. abs(esum-etrace) <= 1.0d-4*max(eabs,abs(etrace)) + 1.0d-3) then
			ok = .false.
			if (iam == 0) write(*,*) 'ELPA eigenvalues add up to',esum,'eV, the trace of the BSE Hamiltonian is',etrace,'eV'
		end if
	end if
	if (.not. ok .and. iam == 0) then
		write(*,*) 'The ELPA solver gave wrong eigenvalues.  The 2-stage solver (BSE_ALGO= elpa2) did so with an ELPA 2025.06'
		write(*,*) 'built by gfortran 15 at -O2 (its vectorizer): rebuild ELPA with -fno-tree-vectorize, or use BSE_ALGO= elpa1.'
	end if

end subroutine bse_elpa_verify


subroutine exclft(sysdim,ngrid,rlat,fosc,enexc,lft) !tempo de vida exciton

	implicit none
	integer,dimension(3) :: ngrid
	real :: lft,fosc,refr,enexc
	real,dimension(3,3) :: rlat
	
	real,parameter :: cfe= 7.297E-3 !constante estrutura fina
	real,parameter :: lspd= 299.792E16 !vel luz Ang/s^2
	real,parameter :: hbar= 6.582E-16 !const planck eV.s
	
	real,parameter :: cic= -(0.0904756)*10**(3)
	
	real,parameter:: pi=acos(-1.)	
	
	real :: aux1,aux2,aux3,vc
	
	character(len=5) :: sysdim
	

	
	select case (sysdim)
	
	 case("3D")
	 
	 	!for 3D systems- DOI: 10.1103/PhysRevLett.95.247402
		aux1 = (4.0*cfe)*(enexc*enexc*enexc)
		aux2 = fosc/dble(ngrid(1)*ngrid(2)*ngrid(3))
		aux3 = 3.0*(lspd*lspd)*(hbar*hbar*hbar)
		lft = (aux1*aux2)/aux3
		lft = (1.0/lft)
	 
	 case("2D")
	 
	 	 !for 2D systems- DOI:10.1021/nl503799t
	 	call vcell2D(rlat,vc)
		aux1 = dble(ngrid(1)*ngrid(2)*ngrid(3))*vc*hbar
		aux2 = (8.0*pi*cfe*enexc)*(fosc)	
		lft = aux1/aux2
	 
	 case("1D")
	 
	 	!for 1D systems- DOI: 10.1103/PhysRevLett.95.247402
		aux2 = (2.0*pi*cfe*enexc*enexc)*fosc
		aux1 = dble(ngrid(1)*ngrid(2)*ngrid(3))*rlat(1,1)*hbar*hbar*lspd
		lft = aux1/aux2
	 
	 case default

		write(*,*) "Wrong System dimension"
		STOP

	end select 
	

	

	
end subroutine exclft



subroutine excwf(outputfolder,ngkpt,kpt,nc,nv,nocp,stt,excenergy,excnum,ewf) !exciton wf direto

	implicit none
	
	integer :: i,erro,counter,j
	
	character(len=70) :: outputfolder    !pasta saida
	character(len=200) :: file1,file2
	CHARACTER(LEN=30) :: Format,Format2
	integer :: ngkpt,nc,nv,excnum
	integer,dimension(ngkpt) :: nocp
	integer,dimension(ngkpt*nc*nv,4) :: stt
	real :: excenergy
	real,dimension(ngkpt,3) :: kpt
	complex,dimension(ngkpt*nc*nv) :: ewf
	real,dimension(ngkpt) :: den
	
	
	WRITE (file1,"(a7,I0,a4)") 'exc_wf_',excnum,'.dat'
	Format = "(I6,3F15.4,2I6,2E15.4)"   ! the integer widths are required by gfortran
	
	WRITE (file2,"(a8,I0,a4)") 'exc_den_',excnum,'.dat'
	Format2 = "(3F15.4,1E15.4)"	
	
	OPEN(UNIT=700+excnum, FILE=trim(outputfolder)//trim(file1),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wf BSE output file"
    	
    	OPEN(UNIT=800+excnum, FILE=trim(outputfolder)//trim(file2),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening den BSE output file"
    	
    	write(700+excnum,*) "#"," ","excitonic state",excnum
    	write(700+excnum,*) "#"," ","exciton energy:",excenergy
    	write(700+excnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(700+excnum,*) "#"," ","Number of conduction states","",nc
    	write(700+excnum,*) "#"," ","Number of valence states","",nv
    	write(700+excnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(700+excnum,*) "#"," ","nocpk"," ","kx"," ","ky"," ","kz"," ","nc"," ","nv"," ","re_wf"," ","imag_wf"
    	

	do i=1,ngkpt*nc*nv
		
 	write(700+excnum,Format) nocp(stt(i,4)),kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),&
 			   	  nocp(stt(i,4))+stt(i,3)-nv,nocp(stt(i,4))-nv+stt(i,2),real(ewf(i)),aimag(ewf(i))
 
	end do
	
	write(800+excnum,*) "#"," ","excitonic state",excnum
    	write(800+excnum,*) "#"," ","exciton energy:",excenergy
    	write(800+excnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(800+excnum,*) "#"," ","Number of conduction states","",nc
    	write(800+excnum,*) "#"," ","Number of valence states","",nv
    	write(800+excnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(800+excnum,*) "#"," ","kx"," ","ky"," ","kz"," ","den"
    	
    	counter=0
    	
    	do i=1,ngkpt
    	
    		den(i) = 0.0
    	
    		do j=1,nc*nv
    		
    			counter=counter+1
    			
    			den(i) = den(i)+ewf(counter)*conjg(ewf(counter)) 
    		
    			
    		
    		end do
    		
    		 write(800+excnum,Format2) kpt(i,1),kpt(i,2),kpt(i,3),den(i)
    	
    	end do
    	
    
	close(700+excnum)
	close(800+excnum)

end subroutine excwf

subroutine excwfi(outputfolder,ngkpt,kpt,qpt,nc,nv,nocp,stt,excenergy,excnum,qptnum,ewf) !exciton wf indireto

	implicit none
	
	integer :: i,erro,counter,j
	
	character(len=70) :: outputfolder    !pasta saida
	character(len=200) :: file1,file2
	CHARACTER(LEN=30) :: Format,Format2
	integer :: ngkpt,nc,nv,excnum,qptnum
	integer,dimension(ngkpt) :: nocp
	integer,dimension(ngkpt*nc*nv,4) :: stt
	real :: excenergy
	real,dimension(ngkpt,3) :: kpt
	real,dimension(4) :: qpt
	complex,dimension(ngkpt*nc*nv) :: ewf
	real,dimension(ngkpt) :: den
	
	
	WRITE (file1,"(a7,I0,a1,I0,a4)") 'exc_wf_',excnum,"_",qptnum,'.dat'
	Format = "(I6,3F15.4,2I6,2E15.4)"   ! the integer widths are required by gfortran
	
	WRITE (file2,"(a8,I0,a1,I0,a4)") 'exc_den_',excnum,"_",qptnum,'.dat'
	Format2 = "(3F15.4,1E15.4)"
	
	OPEN(UNIT=700+excnum*qptnum, FILE=trim(outputfolder)//trim(file1),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening wf BSE output file"
    	
    	OPEN(UNIT=800+excnum*qptnum, FILE=trim(outputfolder)//trim(file2),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening den BSE output file"    	
    	
    	write(700+excnum*qptnum,*) "#"," ","excitonic momentum",real(qpt(2)),real(qpt(3)),real(qpt(4))
    	write(700+excnum*qptnum,*) "#"," ","excitonic state",excnum
    	write(700+excnum*qptnum,*) "#"," ","exciton energy:",excenergy
    	write(700+excnum*qptnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(700+excnum*qptnum,*) "#"," ","Number of conduction states","",nc
    	write(700+excnum*qptnum,*) "#"," ","Number of valence states","",nv
    	write(700+excnum*qptnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(700+excnum*qptnum,*) "#"," ","nocpk"," ","kx"," ","ky"," ","kz"," ","nc"," ","nv"," ","re_wf"," ","imag_wf"
    	

	do i=1,ngkpt*nc*nv
		
 	write(700+excnum*qptnum,Format) nocp(stt(i,4)),kpt(stt(i,4),1),kpt(stt(i,4),2),kpt(stt(i,4),3),&
 			   	  nocp(stt(i,4))+stt(i,3)-nv,nocp(stt(i,4))-nv+stt(i,2),real(ewf(i)),aimag(ewf(i))
 
	end do
	
	write(800+excnum*qptnum,*) "#"," ","excitonic momentum",real(qpt(2)),real(qpt(3)),real(qpt(4))
    	write(800+excnum*qptnum,*) "#"," ","excitonic state",excnum
    	write(800+excnum*qptnum,*) "#"," ","exciton energy:",excenergy
    	write(800+excnum*qptnum,*) "#"," ","number occupied states",nocp(1),"(not valid for metallic systems)"
    	write(800+excnum*qptnum,*) "#"," ","Number of conduction states","",nc
    	write(800+excnum*qptnum,*) "#"," ","Number of valence states","",nv
    	write(800+excnum*qptnum,*) "#"," ","Number of kpoints","",ngkpt    	
    	write(800+excnum*qptnum,*) "#"," ","kx"," ","ky"," ","kz"," ","den"
    	
    	counter=0
    	
    	do i=1,ngkpt
    	
    		den(i) = 0.0
    	
    		do j=1,nc*nv
    		
    			counter=counter+1
    			
    			den(i) = den(i)+ewf(counter)*conjg(ewf(counter)) 
    		
    			
    		
    		end do
    		
    		 write(800+excnum*qptnum,Format2) kpt(i,1),kpt(i,2),kpt(i,3),den(i)
    	
    	end do    	
    	
    
	close(700+excnum*qptnum)
	close(800+excnum*qptnum)

end subroutine excwfi
















subroutine bse_q_images(kpt1,kpt2,rlat,nimg,shift)
	! Periodic images of kpt2 closest to kpt1. A BSE kernel element
	! depends on q = kpt1-kpt2 through V(q) (and the Wannier-centre phases),
	! which is not periodic in q, so q must be the shortest vector of its
	! class q+G (likewise the exciton momentum Q in the exchange term). On the zone boundary two or three images are equally short
	! (K-type q on a hexagonal lattice, for instance): all are returned so
	! the caller can average over them. shift(:,i) is the reciprocal lattice
	! vector to add to kpt2; nimg is 1 away from the zone boundary.
	implicit none
	real,dimension(3),intent(in) :: kpt1,kpt2
	real,dimension(3,3),intent(in) :: rlat
	integer,intent(out) :: nimg
	real,dimension(3,27),intent(out) :: shift
	real,parameter :: pi=acos(-1.)
	real,parameter :: tol=1.0e-4
	real,dimension(3) :: blat1,blat2,blat3,q,g,qi
	real,dimension(3,27) :: gall
	real,dimension(27) :: qlen
	integer,dimension(3) :: n0
	integer :: i1,i2,i3,ic
	real :: qmin

	call recvec(rlat(1,:),rlat(2,:),rlat(3,:),blat1,blat2,blat3)
	q=kpt1-kpt2
	! crystal coordinates of q (a_i.b_j = 2 pi delta_ij), rounded
	n0(1)=nint(dot_product(q,rlat(1,:))/(2.0*pi))
	n0(2)=nint(dot_product(q,rlat(2,:))/(2.0*pi))
	n0(3)=nint(dot_product(q,rlat(3,:))/(2.0*pi))
	ic=0
	do i1=-1,1
		do i2=-1,1
			do i3=-1,1
				ic=ic+1
				g=real(n0(1)+i1)*blat1+real(n0(2)+i2)*blat2+real(n0(3)+i3)*blat3
				gall(:,ic)=g
				qi=q-g
				qlen(ic)=sqrt(dot_product(qi,qi))
			end do
		end do
	end do
	qmin=minval(qlen)
	nimg=0
	do ic=1,27
		if (qlen(ic) .le. qmin*(1.0+tol)+1.0e-6) then
			nimg=nimg+1
			shift(:,nimg)=gall(:,ic)
		end if
	end do

end subroutine bse_q_images


subroutine bse_q_images_fold(kpt1,kpt2,rlat,mfold,nmax,nimg,shift,wimg)
	! Images of the transfer q = kpt1-kpt2 for the direct kernel of an M x M supercell run (BSE_QFOLD = M = mfold). rlat is the
	! supercell lattice (a_1, a_2 in plane, a_3 not folded): its reciprocal lattice Gamma_SC has the spacing |b|/M, b the one of
	! the unit cell, and Gamma_UC = M Gamma_SC in plane. A pair of supercell k points stands for the M^2 transfers q+g+Gamma_UC
	! of the unit-cell kernel, g in Gamma_SC/Gamma_UC = {a B_1 + b B_2, a, b = 0..M-1}, B_i the reciprocal vectors of rlat: the
	! fold pair (g_c, g_c') of the unit-cell states is coupled by the class g = g_c - g_c' and by no other (translation symmetry),
	! with the transfer that is the SHORTEST over Gamma_UC of its class, as in the unit cell. So every class contributes its
	! own shortest image(s), weight 1 per class (equally short images of a class, within the 1e-4 of bse_q_images, share it:
	! wimg = 1/n_class), each with its own W(|q|) and vertices. bse_q_images alone keeps the globally shortest image, i.e. one
	! class, and drops the attraction at |q| > |b|/(2M).
	! shift(:,i) is added to kpt2 (q_i = q - shift), as in bse_q_images; shift and wimg hold nmax = 27 M^2 entries.
	! mfold = 1 is bse_q_images with wimg = 1/nimg.
	implicit none
	real,dimension(3),intent(in) :: kpt1,kpt2
	real,dimension(3,3),intent(in) :: rlat
	integer,intent(in) :: mfold,nmax
	integer,intent(out) :: nimg
	real,dimension(3,nmax),intent(out) :: shift
	real,dimension(nmax),intent(out) :: wimg
	real,parameter :: pi=acos(-1.)
	real,parameter :: tol=1.0e-4
	real,dimension(3) :: blat1,blat2,blat3,q,p,g,gc,qi
	real,dimension(3,27) :: gall
	real,dimension(27) :: qlen
	integer,dimension(3) :: n0
	integer :: i1,i2,i3,ic,ia,ib,ntie,imf
	real :: qmin

	if (mfold .le. 1) then
		call bse_q_images(kpt1,kpt2,rlat,nimg,shift)
		wimg(1:nimg)=1.0/real(nimg)
		return
	end if
	imf=mfold
	call recvec(rlat(1,:),rlat(2,:),rlat(3,:),blat1,blat2,blat3)
	q=kpt1-kpt2
	nimg=0
	do ia=0,mfold-1
		do ib=0,mfold-1
			gc=real(ia)*blat1+real(ib)*blat2
			p=q-gc
			! crystal coordinates of p in the lattice Gamma_UC (basis M B_1, M B_2, B_3), rounded
			n0(1)=nint(dot_product(p,rlat(1,:))/(2.0*pi*real(imf)))
			n0(2)=nint(dot_product(p,rlat(2,:))/(2.0*pi*real(imf)))
			n0(3)=nint(dot_product(p,rlat(3,:))/(2.0*pi))
			ic=0
			do i1=-1,1
				do i2=-1,1
					do i3=-1,1
						ic=ic+1
						g=gc+real(imf*(n0(1)+i1))*blat1+real(imf*(n0(2)+i2))*blat2+real(n0(3)+i3)*blat3
						gall(:,ic)=g
						qi=q-g
						qlen(ic)=sqrt(dot_product(qi,qi))
					end do
				end do
			end do
			qmin=minval(qlen)
			ntie=0
			do ic=1,27
				if (qlen(ic) .le. qmin*(1.0+tol)+1.0e-6) ntie=ntie+1
			end do
			if (nimg+ntie .gt. nmax) stop 'bse_q_images_fold: too many images (BSE_QFOLD)'
			do ic=1,27
				if (qlen(ic) .le. qmin*(1.0+tol)+1.0e-6) then
					nimg=nimg+1
					shift(:,nimg)=gall(:,ic)
					wimg(nimg)=1.0/real(ntie)
				end if
			end do
		end do
	end do

end subroutine bse_q_images_fold


subroutine bse_qfold_check(mfold,rlat,use_ff)
	! BSE_QFOLD = M > 1: the run is an M x M supercell whose in-plane lattice vectors are a_1, a_2 of rlat (a_3 along the
	! vacuum), and the pair densities come from the basis itself, not from a form-factor file, which holds the shortest
	! images of the run's lattice only.
	implicit none
	integer,intent(in) :: mfold
	real,dimension(3,3),intent(in) :: rlat
	logical,intent(in) :: use_ff
	real :: c13,c23,l1,l2,l3

	if (mfold .le. 1) return
	if (use_ff) stop 'BSE_QFOLD > 1 needs BSE_FF= F: the form-factor file (BSE_FF_FILE) holds only the shortest images of q'
	l1=sqrt(dot_product(rlat(1,:),rlat(1,:)))
	l2=sqrt(dot_product(rlat(2,:),rlat(2,:)))
	l3=sqrt(dot_product(rlat(3,:),rlat(3,:)))
	c13=abs(dot_product(rlat(1,:),rlat(3,:)))/(l1*l3)
	c23=abs(dot_product(rlat(2,:),rlat(3,:)))/(l2*l3)
	if (c13 .gt. 1.0e-3 .or. c23 .gt. 1.0e-3) &
		stop 'BSE_QFOLD > 1 folds the in-plane lattice vectors a_1, a_2: a_3 must be perpendicular to both'

end subroutine bse_qfold_check


subroutine orbital_centres(w90basis,tau)

	implicit none

	integer,intent(in) :: w90basis
	real,dimension(w90basis,3),intent(out) :: tau
	real,allocatable,dimension(:,:),save :: centres
	logical,save :: ready = .false.
	character(len=13),dimension(4),parameter :: files = (/ 'basis_set-NP ', &
		'basis_set-sp ','basis_set-nc ','basis_set-soc' /)
	character(len=16) :: species
	integer :: i,f,idx,erro
	real :: x,y,z

	! Orbital centres (Angstrom) of the DFT=S basis in tb-NP.dat order, for
	! the atomic-gauge vertices of matrizelbse and matrizelbsekq. Read once,
	! from the basis_set-* file utils/siesta2wtb.py writes with tb-*.dat,
	! in the run directory. Without one they are zero: the old
	! lattice-gauge vertices, with a warning.
	!$omp critical (orbital_centres_load)
	if (.not. ready) then
		allocate(centres(w90basis,3))
		centres = 0.0
		do f=1,size(files)
			open(unit=977,file=trim(files(f)),status='old',action='read',iostat=erro)
			if (erro /= 0) cycle
			read(977,*,iostat=erro)
			do i=1,w90basis
				if (erro == 0) read(977,*,iostat=erro) idx,species,x,y,z
				if (erro == 0) centres(i,:) = (/x,y,z/)
			end do
			close(977)
			if (erro == 0) then
				write(*,'(2A)') ' DFT=S BSE: orbital centres from ',trim(files(f))
				exit
			end if
			centres = 0.0
		end do
		if (erro /= 0) write(*,'(A)') ' WARNING: DFT=S BSE without a basis_set-* file:'// &
			' orbital centres set to zero, C3 symmetry is not preserved'
		ready = .true.
	end if
	!$omp end critical (orbital_centres_load)
	tau = centres

end subroutine orbital_centres


subroutine orbital_labels(norb,spec,lq,ok)

	! Species and shell l of each orbital of the DFT=S basis, from the same basis_set-* file as orbital_centres (columns
	! index species x y z [l m spin]); l = -1 when the file has no l column.  ok = .false. without the file.

	implicit none

	integer,intent(in) :: norb
	character(len=*),dimension(norb),intent(out) :: spec
	integer,dimension(norb),intent(out) :: lq
	logical,intent(out) :: ok
	character(len=13),dimension(4),parameter :: files = (/ 'basis_set-NP ', &
		'basis_set-sp ','basis_set-nc ','basis_set-soc' /)
	character(len=256) :: line
	character(len=16) :: species
	integer :: i,f,idx,erro,l,m,e2
	real :: x,y,z

	spec = ''
	lq = -1
	ok = .false.
	do f=1,size(files)
		open(unit=978,file=trim(files(f)),status='old',action='read',iostat=erro)
		if (erro /= 0) cycle
		read(978,'(A)',iostat=erro) line
		do i=1,norb
			if (erro == 0) read(978,'(A)',iostat=erro) line
			if (erro == 0) then
				read(line,*,iostat=e2) idx,species,x,y,z,l,m
				if (e2 /= 0) then
					l = -1
					read(line,*,iostat=erro) idx,species,x,y,z
				end if
			end if
			if (erro == 0) then
				spec(i) = species
				lq(i) = l
			end if
		end do
		close(978)
		if (erro == 0) then
			ok = .true.
			exit
		end if
		spec = ''
		lq = -1
	end do

end subroutine orbital_labels
