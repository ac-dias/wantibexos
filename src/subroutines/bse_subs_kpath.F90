function matrizelbsekq(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,q,rlat,est1,ec1,ev1,vbc1 &
                       ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,dft,nvec,rvec,sk,skp,skq,skpq) !funcao para calcular o elemento de matriz da matriz bse

	implicit none

	character(len=7) :: coultype
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
	real,dimension(4) :: q
	real,dimension(3) :: vq,v0
	complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2

	complex, dimension(w90basis) :: vbc,vbv,vbvkp

	real :: tolr
	integer :: ktol
	
	complex:: matrizelbsekq

	real :: modk,modq

	real,dimension(3) :: ediel,ediel_bare

	real :: lc

	complex :: vc,vv
	complex :: vcv,vvc

	real :: vcoulk,vcoulq

	real :: v2dk,vcoul,v3diel,v3davg,v2dt,v2dtavg,v0dt,v2dt2
	real :: v2dohono,v2drk,v1dt,v2d,v2diel
	real :: v1d,v1diel

	real :: r0

	v0 = 0.0
	ediel_bare = 1.0

	vq = q(2:4)


	call modvec(kpt1,kpt2,modk)

	call modvecq(q,modq)



	select case (coultype)

	case("V2DK")

		vcoulk= v2dk(kpt1,kpt2,ediel,rlat,ngrid,lc,tolr)
		vcoulq= v2dk(vq,v0,ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vcoulk= vcoul(kpt1,kpt2,rlat,ngrid,tolr)
		vcoulq= vcoul(vq,v0,rlat,ngrid,tolr)

	case("V3DL")

		vcoulk= v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)
		vcoulq= v3diel(vq,v0,ediel,rlat,ngrid,tolr)

	case("V3DAVG")

		vcoulk= v3davg(kpt1,kpt2,1.0,rlat,ngrid,tolr)
		vcoulq= v3davg(vq,v0,1.0,rlat,ngrid,tolr)

	case("V3DLAVG")

		vcoulk= v3davg(kpt1,kpt2,ediel(2),rlat,ngrid,tolr)
		vcoulq= v3davg(vq,v0,ediel(2),rlat,ngrid,tolr)
		
	case("V2D")

		vcoulk= v2d(kpt1,kpt2,rlat,ngrid,tolr)
		vcoulq= v2d(vq,v0,rlat,ngrid,tolr)		

	case("V2DL")

		vcoulk= v2diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)
		vcoulq= v2diel(vq,v0,ediel,rlat,ngrid,tolr)				

	case("V2DT")

		vcoulk= v2dt(kpt1,kpt2,ngrid,rlat,tolr)
		vcoulq= v2dt(vq,v0,ngrid,rlat,tolr)

	case("V2DTAVG")

		vcoulk= v2dtavg(kpt1,kpt2,ediel,ngrid,rlat,tolr)
		vcoulq= v2dtavg(vq,v0,ediel_bare,ngrid,rlat,tolr)

	case("V2DT2")

		vcoulk= v2dt2(kpt1,kpt2,ngrid,rlat,lc,tolr)
		vcoulq= v2dt2(vq,v0,ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vcoulk= v2dohono(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)
		vcoulq= v2dohono(vq,v0,ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vcoulk= v2drk(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		vcoulq= v2drk(vq,v0,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		
	case("V1D")

		vcoulk= v1d(kpt1,kpt2,ngrid,rlat,lc,tolr)
		vcoulq= v1d(vq,v0,ngrid,rlat,lc,tolr)

	case("V1DL")

		vcoulk= v1diel(kpt1,kpt2,ngrid,rlat,lc,ediel,tolr)
		vcoulq= v1diel(vq,v0,ngrid,rlat,lc,ediel,tolr)		
		
	case("V1DT")
	
		vcoulk= v1dt(kpt1,kpt2,ngrid,rlat,tolr,lc)	
		vcoulq= v1dt(vq,v0,ngrid,rlat,tolr,lc)		

	case("V0DT")

		vcoulk= v0dt(kpt1,kpt2,ngrid,rlat,tolr)
		vcoulq= v0dt(vq,v0,ngrid,rlat,tolr)

	case default

		write(*,*) "Wrong Coulomb Potential"
		STOP

	end select

	




if (modq .eq. 0.) then


	if (est1(1) .eq. est2(1)) then

		matrizelbsekq= (ec1-ev1) + vcoulk



	else
	
		select case (dft)
		
		case ("S")
		

		  !call overlap(w90basis,nvec,rvec,ovp,kpt1(1),kpt1(2),kpt1(3),sk)
		  !call overlap(w90basis,nvec,rvec,ovp,kpt2(1),kpt2(2),kpt2(3),skp)
		 
		  call sandwich_average(w90basis,vbc1,sk,skp,vbc2,vc)
		  call sandwich_average(w90basis,vbv1,sk,skp,vbv2,vv)
		 
		  matrizelbsekq=  vcoulk*vc*vv		
		
		case default	

	
		 call vecconjg(vbc1,w90basis,vbc)

		 call vecconjg(vbv1,w90basis,vbv)

		 call prodintsq(vbc,vbc2,w90basis,vc)

		 call prodintsq(vbv,vbv2,w90basis,vv)


		 matrizelbsekq= vcoulk*vc*vv

		end select

	end if
		


else 


	if (est1(1) .eq. est2(1)) then
	
		
		select case (dft)
		
		case ("S")


		  !call overlap(w90basis,nvec,rvec,ovp,kpt1(1),kpt1(2),kpt1(3),sk)
		  !call overlap(w90basis,nvec,rvec,ovp,kpt2(1),kpt2(2),kpt2(3),skp)
		  
		  !call overlap(w90basis,nvec,rvec,ovp,kpt1(1)+q(2),kpt1(2)+q(3),kpt1(3)+q(4),skq)
		  !call overlap(w90basis,nvec,rvec,ovp,kpt2(1)+q(2),kpt2(2)+q(3),kpt2(3)+q(4),skpq)
		  
		  call sandwich_average(w90basis,vbc1,skq,sk,vbv1,vcv)
		  call sandwich_average(w90basis,vbv2,skp,skpq,vbc2,vvc)

		 matrizelbsekq= (ec1-ev1) + vcoulk &
			     - vcoulq*vcv*vvc
		
		case default


		 call vecconjg(vbc1,w90basis,vbc)

		 call vecconjg(vbv1,w90basis,vbv)

		 call vecconjg(vbv2,w90basis,vbvkp)


		 call prodintsq(vbc,vbc2,w90basis,vc)

		 call prodintsq(vbv,vbv2,w90basis,vv)


		 call prodintsq(vbc,vbv1,w90basis,vcv)

		 call prodintsq(vbvkp,vbc2,w90basis,vvc)


		 matrizelbsekq= (ec1-ev1) + vcoulk &
			     - vcoulq*vcv*vvc

	
		end select


	else

		select case (dft)
		
		case ("S")
		
		  !call overlap(w90basis,nvec,rvec,ovp,kpt1(1),kpt1(2),kpt1(3),sk)
		  !call overlap(w90basis,nvec,rvec,ovp,kpt2(1),kpt2(2),kpt2(3),skp)
		  
		  !call overlap(w90basis,nvec,rvec,ovp,kpt1(1)+q(2),kpt1(2)+q(3),kpt1(3)+q(4),skq)
		  !call overlap(w90basis,nvec,rvec,ovp,kpt2(1)+q(2),kpt2(2)+q(3),kpt2(3)+q(4),skpq)
		  
		  call sandwich_average(w90basis,vbc1,skq,skpq,vbc2,vc)
		  call sandwich_average(w90basis,vbv1,sk,skp,vbv2,vv)
		  
		  call sandwich_average(w90basis,vbc1,skq,sk,vbv1,vcv)
		  call sandwich_average(w90basis,vbv2,skp,skpq,vbc2,vvc)

		 matrizelbsekq= vcoulk*vc*vv&
				- vcoulq*vcv*vvc		
		
		case default

		 call vecconjg(vbc1,w90basis,vbc)

		 call vecconjg(vbv1,w90basis,vbv)

		 call vecconjg(vbv2,w90basis,vbvkp)


		 call prodintsq(vbc,vbc2,w90basis,vc)

		 call prodintsq(vbv,vbv2,w90basis,vv)


		 call prodintsq(vbc,vbv1,w90basis,vcv)

		 call prodintsq(vbvkp,vbc2,w90basis,vvc)

		 matrizelbsekq= vcoulk*vc*vv&
				- vcoulq*vcv*vvc

		end select

	end if





end if





end function matrizelbsekq
