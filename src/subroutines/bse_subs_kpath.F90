function matrizelbsekq(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,q,rlat,est1,ec1,ev1,vbc1 &
                       ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,dft,nvec,rvec,sk,skp,skq,skpq, &
                       use_center_phase,centers) !funcao para calcular o elemento de matriz da matriz bse

	implicit none

	character(len=10) :: coultype
	character(len=1) :: dft

	integer,dimension(3) :: ngrid

	real,dimension(3,3) :: rlat
	real :: ez,w
	
	integer :: w90basis,nvec
	
	real,dimension(nvec,3) :: rvec
	!real,dimension(nvec,w90basis,w90basis) :: ovp
	
	complex,dimension(w90basis,w90basis) :: sk,skq,skp,skpq		


	real :: a,vcell1

	integer, dimension(4) :: est1,est2
	real :: ec1,ec2,ev1,ev2
	real,dimension(3) :: kpt1,kpt2
	real,dimension(3) :: kpt2f
	real,dimension(3,27) :: gshift
	integer :: nimg
	real,dimension(4) :: q
	real,dimension(3) :: vq,v0
	complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2

	complex, dimension(w90basis) :: vbc,vbv

	real :: tolr
	integer :: ktol
	
	complex:: matrizelbsekq

	real :: modk,modq

	real,dimension(3) :: ediel,ediel_bare

	real :: lc

	complex :: vc,vv
	complex :: vcv,vvc
	complex :: dirv,excv
	logical :: use_center_phase
	real,dimension(3,w90basis) :: centers
	real,dimension(3,27) :: gshiftq
	integer :: nimgq,img,m
	complex,dimension(w90basis) :: rphase
	real :: ang

	real :: vcoulk,vcoulq

	real :: v2dk,vcoul,v3diel,v3davg,v2dt,v2dtavg,v0dt,v2dt2
	real :: v2dohono,v2drk,v1dt,v2d,v2diel,v2davgtpv
	real :: v1d,v1diel

	real :: r0
	real,dimension(3) :: qimg
	real,dimension(w90basis,3) :: tau

	v0 = 0.0
	ediel_bare = 1.0

	vq = q(2:4)


	! Fold q = kpt1-kpt2 to its shortest image q+G: V(q) is not periodic
	! in q. Equally short boundary images give the same V; the centre-phase
	! vertices below are averaged over them.
	call bse_q_images(kpt1,kpt2,rlat,nimg,gshift)
	kpt2f=kpt2+gshift(:,1)
	! the exchange term needs Q itself as its shortest image too
	call bse_q_images(vq,v0,rlat,nimgq,gshiftq)
	vq=vq-gshiftq(:,1)
	modq=sqrt(dot_product(vq,vq))
	call modvec(kpt1,kpt2f,modk)

	select case (coultype)

	case("V2DK")

		vcoulk= v2dk(kpt1,kpt2f,ediel,rlat,ngrid,lc,tolr)
		vcoulq= v2dk(vq,v0,ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vcoulk= vcoul(kpt1,kpt2f,rlat,ngrid,tolr)
		vcoulq= vcoul(vq,v0,rlat,ngrid,tolr)

	case("V3DL")

		vcoulk= v3diel(kpt1,kpt2f,ediel,rlat,ngrid,tolr)
		vcoulq= v3diel(vq,v0,ediel,rlat,ngrid,tolr)

	case("V3DAVG")

		vcoulk= v3davg(kpt1,kpt2f,1.0,rlat,ngrid,tolr)
		vcoulq= v3davg(vq,v0,1.0,rlat,ngrid,tolr)

	case("V3DLAVG")

		vcoulk= v3davg(kpt1,kpt2f,ediel(2),rlat,ngrid,tolr)
		vcoulq= v3davg(vq,v0,ediel(2),rlat,ngrid,tolr)
		
	case("V2D")

		vcoulk= v2d(kpt1,kpt2f,rlat,ngrid,tolr)
		vcoulq= v2d(vq,v0,rlat,ngrid,tolr)		

	case("V2DL")

		vcoulk= v2diel(kpt1,kpt2f,ediel,rlat,ngrid,tolr)
		vcoulq= v2diel(vq,v0,ediel,rlat,ngrid,tolr)				

	case("V2DT")

		vcoulk= v2dt(kpt1,kpt2f,ngrid,rlat,tolr)
		vcoulq= v2dt(vq,v0,ngrid,rlat,tolr)

	case("V2DTAVG")

		vcoulk= v2dtavg(kpt1,kpt2f,ediel,ngrid,rlat,tolr)
		vcoulq= v2dtavg(vq,v0,ediel_bare,ngrid,rlat,tolr)

	case("V2DAVGTPV")

		vcoulk= v2davgtpv(kpt1,kpt2f,ediel,ngrid,rlat,tolr)
		vcoulq= v2davgtpv(vq,v0,ediel_bare,ngrid,rlat,tolr)

	case("V2DT2")

		vcoulk= v2dt2(kpt1,kpt2f,ngrid,rlat,lc,tolr)
		vcoulq= v2dt2(vq,v0,ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vcoulk= v2dohono(kpt1,kpt2f,ngrid,rlat,ediel,w,ez,tolr)
		vcoulq= v2dohono(vq,v0,ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vcoulk= v2drk(kpt1,kpt2f,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		vcoulq= v2drk(vq,v0,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		
	case("V1D")

		vcoulk= v1d(kpt1,kpt2f,ngrid,rlat,lc,tolr)
		vcoulq= v1d(vq,v0,ngrid,rlat,lc,tolr)

	case("V1DL")

		vcoulk= v1diel(kpt1,kpt2f,ngrid,rlat,lc,ediel,tolr)
		vcoulq= v1diel(vq,v0,ngrid,rlat,lc,ediel,tolr)		
		
	case("V1DT")
	
		vcoulk= v1dt(kpt1,kpt2f,ngrid,rlat,tolr,lc)	
		vcoulq= v1dt(vq,v0,ngrid,rlat,tolr,lc)		

	case("V0DT")

		vcoulk= v0dt(kpt1,kpt2f,ngrid,rlat,tolr)
		vcoulq= v0dt(vq,v0,ngrid,rlat,tolr)

	case default

		write(*,*) "Wrong Coulomb Potential"
		STOP

	end select

	




	! Direct vertices <c1|c2><v2|v1> (q = kpt1-kpt2) and exchange vertices
	! <c1|v1><v2|c2> (exciton momentum Q, none at Q=0). With the Wannier-centre
	! phases exp(+i q.t_m) and exp(+i Q.t_m) of <c1|c2> and <c1|v1> they
	! depend on the image q+G (Q+G), so, as in matrizelbse, they are averaged
	! over the equally short images.

	dirv=cmplx(0.0,0.0)
	excv=cmplx(0.0,0.0)

	select case (dft)

	case ("S")

		!call overlap(w90basis,nvec,rvec,ovp,kpt1(1),kpt1(2),kpt1(3),sk)
		!call overlap(w90basis,nvec,rvec,ovp,kpt2(1),kpt2(2),kpt2(3),skp)

		!call overlap(w90basis,nvec,rvec,ovp,kpt1(1)+q(2),kpt1(2)+q(3),kpt1(3)+q(4),skq)
		!call overlap(w90basis,nvec,rvec,ovp,kpt2(1)+q(2),kpt2(2)+q(3),kpt2(3)+q(4),skpq)

		! atomic gauge, as the centre-phase branch below: orbital centres
		! from orbital_centres, each vertex averaged over the equally short
		! images of its momentum transfer (q for the direct, Q for the
		! exchange term); the hole vertex is <v k2|v k1>
		call orbital_centres(w90basis,tau)
		if (est1(1) .ne. est2(1)) then
			do img=1,nimg
				qimg=kpt1-kpt2-gshift(:,img)
				call sandwich_phase(w90basis,vbc1,skq,skpq,vbc2,qimg,tau,vc)
				call sandwich_phase(w90basis,vbv1,sk,skp,vbv2,qimg,tau,vv)
				dirv=dirv+vc*conjg(vv)
			end do
			dirv=dirv/real(nimg)
		end if
		if (modq .ne. 0.) then
			do img=1,nimgq
				qimg=q(2:4)-gshiftq(:,img)
				call sandwich_phase(w90basis,vbc1,skq,sk,vbv1,qimg,tau,vcv)
				call sandwich_phase(w90basis,vbv2,skp,skpq,vbc2,-qimg,tau,vvc)
				excv=excv+vcv*vvc
			end do
			excv=excv/real(nimgq)
		end if

	case default

		if (use_center_phase) then

			if (est1(1) .ne. est2(1)) then
				do img=1,nimg
					do m=1,w90basis
						ang=dot_product(kpt1-kpt2-gshift(:,img),centers(:,m))
						rphase(m)=cmplx(cos(ang),sin(ang))
					end do
					vc=sum(conjg(vbc1)*vbc2*rphase)
					vv=sum(conjg(vbv2)*vbv1*conjg(rphase))
					dirv=dirv+vc*vv
				end do
				dirv=dirv/real(nimg)
			end if
			if (modq .ne. 0.) then
				do img=1,nimgq
					do m=1,w90basis
						ang=dot_product(q(2:4)-gshiftq(:,img),centers(:,m))
						rphase(m)=cmplx(cos(ang),sin(ang))
					end do
					vcv=sum(conjg(vbc1)*vbv1*rphase)
					vvc=sum(conjg(vbv2)*vbc2*conjg(rphase))
					excv=excv+vcv*vvc
				end do
				excv=excv/real(nimgq)
			end if

		else

			call vecconjg(vbc1,w90basis,vbc)

			call vecconjg(vbv2,w90basis,vbv)

			if (est1(1) .ne. est2(1)) then
				call prodintsq(vbc,vbc2,w90basis,vc)
				call prodintsq(vbv,vbv1,w90basis,vv)
				dirv=vc*vv
			end if
			if (modq .ne. 0.) then
				call prodintsq(vbc,vbv1,w90basis,vcv)
				call prodintsq(vbv,vbc2,w90basis,vvc)
				excv=vcv*vvc
			end if

		end if

	end select


	if (est1(1) .eq. est2(1)) then

		matrizelbsekq= (ec1-ev1) + vcoulk

	else

		matrizelbsekq= vcoulk*dirv

	end if

	if (modq .ne. 0.) matrizelbsekq= matrizelbsekq - vcoulq*excv



end function matrizelbsekq
