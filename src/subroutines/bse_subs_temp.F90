
function matrizelbsetemp(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,rlat,est1,ec1,ev1,vbc1 &
         ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,temp,dft,nvec,rvec,sk,skp) !funcao para calcular o elemento de matriz da matriz bse

	implicit none

	character(len=5) :: coultype
	character(len=1) :: dft

	integer,dimension(3) :: ngrid

	integer :: w90basis,nvec
	real :: a,vcell1
	real :: ez,w
	
	real,dimension(nvec,3) :: rvec
	!real,dimension(nvec,w90basis,w90basis) :: ovp
	
	complex,dimension(w90basis,w90basis) :: sk,skp

	integer, dimension(4) :: est1,est2
	real :: ec1,ec2,ev1,ev2
	real,dimension(3) :: kpt1,kpt2
	complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2

	complex, dimension(w90basis) :: vbc,vbv

	real,dimension(3,3) :: rlat

	real :: tolr
	integer :: ktol
	
	complex:: matrizelbsetemp

	real,parameter:: pi=acos(-1.)

	real :: auxi

	real :: modk

	real,dimension(3) :: ediel

	real :: lc

	complex :: vc,vv

	real :: vcoul1

	real :: vcoul,v2dk,v3diel,v2dt,v0dt,v2dt2
	real :: v2dohono,v2drk,v1dt,v2d,v2diel
	real :: v1d,v1diel,v3davg,v2dtavg
	
	real :: r0
	
	real :: temp, fermidisteh


	select case (coultype)

	case("V2DK")

		vcoul1= v2dk(kpt1,kpt2,ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vcoul1= vcoul(kpt1,kpt2,rlat,ngrid,tolr)

	case("V3DA")

		vcoul1= v3davg(kpt1,kpt2,ngrid,rlat,tolr)

	case("V3DL")

		vcoul1= v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)
		
	case("V2D")

		vcoul1= v2d(kpt1,kpt2,rlat,ngrid,tolr)

	case("V2DL")

		vcoul1= v2diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)		

	case("V2DT")

		vcoul1= v2dt(kpt1,kpt2,ngrid,rlat,tolr)

	case("V2DTA")

		vcoul1= v2dtavg(kpt1,kpt2,ngrid,rlat,tolr)

	case("V2DT2")

		vcoul1= v2dt2(kpt1,kpt2,ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vcoul1= v2dohono(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vcoul1= v2drk(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		
	case("V1D")

		vcoul1= v1d(kpt1,kpt2,ngrid,rlat,lc,tolr)

	case("V1DL")

		vcoul1= v1diel(kpt1,kpt2,ngrid,rlat,lc,ediel,tolr)		
		
	case("V1DT")
	
		vcoul1= v1dt(kpt1,kpt2,ngrid,rlat,tolr,lc)	
		
		
	case("V0DT")

		vcoul1= v0dt(kpt1,kpt2,ngrid,rlat,tolr)

	case default

		write(*,*) "Wrong Coulomb Potential"
		STOP

	end select

	



	if (est1(1) .eq. est2(1)) then


		matrizelbsetemp= (ec1-ev1) + vcoul1*fermidisteh(ec1,ev1,temp)


	else

	
		if (dft .eq. "S") then

		 !call overlap(w90basis,nvec,rvec,ovp,kpt1(1),kpt1(2),kpt1(3),sk)
		 !call overlap(w90basis,nvec,rvec,ovp,kpt2(1),kpt2(2),kpt2(3),skp)
		 
		 call sandwich(w90basis,vbc1,0.5*(sk+skp),vbc2,vc)
		 call sandwich(w90basis,vbv1,0.5*(sk+skp),vbv2,vv)		
		
		else
	
		 call vecconjg(vbc1,w90basis,vbc)

		 call vecconjg(vbv1,w90basis,vbv)

		 call prodintsq(vbc,vbc2,w90basis,vc)

		 call prodintsq(vbv,vbv2,w90basis,vv)
	
	
		end if
	
		matrizelbsetemp=  vcoul1*vc*vv*fermidisteh(ec1,ev1,temp)


	end if
		



end function matrizelbsetemp

function matrizelbsekqtemp(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,q,rlat,est1,ec1,ev1,vbc1 &
                       ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,temp,dft,nvec,rvec,sk,skp,skq,skpq) !funcao para calcular o elemento de matriz da matriz bse

	implicit none

	character(len=5) :: coultype
	character(len=1) :: dft

	integer,dimension(3) :: ngrid

	real,dimension(3,3) :: rlat
	real :: ez,w	

	integer :: w90basis,nvec
	real :: a,vcell1
	
	real,dimension(nvec,3) :: rvec
	!real,dimension(nvec,w90basis,w90basis) :: ovp
	
	complex,dimension(w90basis,w90basis) :: sk,skq,skp,skpq	

	integer, dimension(4) :: est1,est2
	real :: ec1,ec2,ev1,ev2
	real,dimension(3) :: kpt1,kpt2
	real,dimension(4) :: q
	real,dimension(3) :: vq,v0
	complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2

	complex, dimension(w90basis) :: vbc,vbv,vbvkp

	real :: tolr
	integer :: ktol
	
	complex:: matrizelbsekqtemp

	real :: modk,modq

	real,dimension(3) :: ediel

	real :: lc

	complex :: vc,vv
	complex :: vcv,vvc

	real :: vcoulk,vcoulq

	real :: v2dk,vcoul,v3diel,v2dt,v0dt,v2dt2
	real :: v2dohono,v2drk,v1dt,v2d,v2diel
	real :: v1d,v1diel,v3davg,v2dtavg

	real :: r0
	
	real :: temp, fermidisteh

	v0 = 0.0

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
		
	case("V3DA")

		vcoulk= v3davg(kpt1,kpt2,ngrid,rlat,tolr)
		vcoulq= v3davg(vq,v0,ngrid,rlat,tolr)		

	case("V3DL")

		vcoulk= v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)
		vcoulq= v3diel(vq,v0,ediel,rlat,ngrid,tolr)
		
	case("V2D")

		vcoulk= v2d(kpt1,kpt2,rlat,ngrid,tolr)
		vcoulq= v2d(vq,v0,rlat,ngrid,tolr)		

	case("V2DL")

		vcoulk= v2diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)
		vcoulq= v2diel(vq,v0,ediel,rlat,ngrid,tolr)			

	case("V2DT")

		vcoulk= v2dt(kpt1,kpt2,ngrid,rlat,tolr)
		vcoulq= v2dt(vq,v0,ngrid,rlat,tolr)
		
	case("V2DTA")

		vcoulk= v2dtavg(kpt1,kpt2,ngrid,rlat,tolr)
		vcoulq= v2dtavg(vq,v0,ngrid,rlat,tolr)		

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

		matrizelbsekqtemp= (ec1-ev1) + vcoulk*fermidisteh(ec1,ev1,temp)



	else

		select case (dft)

		case ("S")
		

		  !call overlap(w90basis,nvec,rvec,ovp,kpt1(1),kpt1(2),kpt1(3),sk)
		  !call overlap(w90basis,nvec,rvec,ovp,kpt2(1),kpt2(2),kpt2(3),skp)
		 
		  call sandwich(w90basis,vbc1,0.5*(sk+skp),vbc2,vc)
		  call sandwich(w90basis,vbv1,0.5*(sk+skp),vbv2,vv)
		 
		  matrizelbsekqtemp= vcoulk*vc*vv*fermidisteh(ec1,ev1,temp)		
		
		case default	

	
		 call vecconjg(vbc1,w90basis,vbc)

		 call vecconjg(vbv1,w90basis,vbv)

		 call prodintsq(vbc,vbc2,w90basis,vc)

		 call prodintsq(vbv,vbv2,w90basis,vv)


		 matrizelbsekqtemp= vcoulk*vc*vv*fermidisteh(ec1,ev1,temp)
		 
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
		  
		  call sandwich(w90basis,vbc1,0.5*(skq+sk),vbv1,vcv)		  
		  call sandwich(w90basis,vbv2,0.5*(skp+skpq),vbc2,vvc)
		  
		  matrizelbsekqtemp= (ec1-ev1) + (vcoulk - vcoulq*vcv*vvc)*fermidisteh(ec1,ev1,temp)
		  
		case default  		  

		 call vecconjg(vbc1,w90basis,vbc)

		 call vecconjg(vbv1,w90basis,vbv)

		 call vecconjg(vbv2,w90basis,vbvkp)


		 call prodintsq(vbc,vbc2,w90basis,vc)

		 call prodintsq(vbv,vbv2,w90basis,vv)


		 call prodintsq(vbc,vbv1,w90basis,vcv)

		 call prodintsq(vbvkp,vbc2,w90basis,vvc)


		 matrizelbsekqtemp= (ec1-ev1) + (vcoulk - vcoulq*vcv*vvc)*fermidisteh(ec1,ev1,temp)
		
		end select 
			     

	


	else

		select case (dft)
		
		case ("S")
		
		  !call overlap(w90basis,nvec,rvec,ovp,kpt1(1),kpt1(2),kpt1(3),sk)
		  !call overlap(w90basis,nvec,rvec,ovp,kpt2(1),kpt2(2),kpt2(3),skp)
		  
		  !call overlap(w90basis,nvec,rvec,ovp,kpt1(1)+q(2),kpt1(2)+q(3),kpt1(3)+q(4),skq)
		  !call overlap(w90basis,nvec,rvec,ovp,kpt2(1)+q(2),kpt2(2)+q(3),kpt2(3)+q(4),skpq)
		  
		  call sandwich(w90basis,vbc1,0.5*(skq+skpq),vbc2,vc)
		  call sandwich(w90basis,vbv1,0.5*(sk+skp),vbv2,vv)
		  
		  call sandwich(w90basis,vbc1,0.5*(skq+sk),vbv1,vcv)		  
		  call sandwich(w90basis,vbv2,0.5*(skp+skpq),vbc2,vvc)	
	

		   matrizelbsekqtemp= (vcoulk*vc*vv- vcoulq*vcv*vvc)*fermidisteh(ec1,ev1,temp)
	
		case default

		 call vecconjg(vbc1,w90basis,vbc)

		 call vecconjg(vbv1,w90basis,vbv)

		 call vecconjg(vbv2,w90basis,vbvkp)


		 call prodintsq(vbc,vbc2,w90basis,vc)

		 call prodintsq(vbv,vbv2,w90basis,vv)


		 call prodintsq(vbc,vbv1,w90basis,vcv)

		 call prodintsq(vbvkp,vbc2,w90basis,vvc)

		 matrizelbsekqtemp= (vcoulk*vc*vv- vcoulq*vcv*vvc)*fermidisteh(ec1,ev1,temp)
		
		end select		


	end if





end if





end function matrizelbsekqtemp


subroutine dielbseptemp(nthread,dimse,excitonvec,hopt1,hopt2,fdeh,activity)

	use omp_lib
	implicit none


	integer :: dimse,nthread
	real,dimension(dimse) :: activity, exciton
	complex :: actaux,actaux2
	complex,dimension(dimse) :: hopt1,hopt2
	real,dimension(dimse) :: fdeh
	complex,dimension(dimse,dimse) :: excitonvec


	integer :: i,j,k,no
	

	call OMP_SET_NUM_THREADS(nthread)

	activity=0.0
	!actaux=0.0


	! $OMP DO PRIVATE(actaux)
	!$OMP PARALLEL DO PRIVATE(actaux,actaux2)
	do i=1,dimse


		actaux=0.0
		actaux2=0.0


		do j=1,dimse

			
			actaux=actaux+(excitonvec(j,i)*hopt1(j))
			actaux2=actaux2+(excitonvec(j,i)*hopt2(j)*fdeh(j))
			

			
		end do
			

		activity(i)=actaux*conjg(actaux2)

	
	
		


	end do
	!$OMP END PARALLEL DO



end subroutine dielbseptemp






















