
function matrizelbsetemp(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,rlat,est1,ec1,ev1,vbc1 &
         ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,temp) !funcao para calcular o elemento de matriz da matriz bse

	implicit none

	character(len=5) :: coultype

	integer,dimension(3) :: ngrid

	integer :: w90basis
	double precision :: a,vcell1
	double precision :: ez,w

	integer, dimension(4) :: est1,est2
	double precision :: ec1,ec2,ev1,ev2
	double precision,dimension(3) :: kpt1,kpt2
	double complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2

	double complex, dimension(w90basis) :: vbc,vbv

	double precision,dimension(3,3) :: rlat

	double precision :: tolr
	integer :: ktol
	
	double complex:: matrizelbsetemp

	double precision,parameter:: pi=acos(-1.)

	double precision :: auxi

	double precision :: modk

	double precision,dimension(3) :: ediel

	double precision :: lc

	double complex :: vc,vv

	double precision :: vcoul1

	double precision :: vcoul,v2dk,v3diel,v2dt,v0dt,v2dt2
	double precision :: v2dohono,v2drk
	
	double precision :: r0
	
	double precision :: temp, fermidisteh


	select case (coultype)

	case("V2DK")

		vcoul1= v2dk(kpt1,kpt2,ediel,rlat,ngrid,lc,tolr)

	case("V3D")

		vcoul1= vcoul(kpt1,kpt2,rlat,ngrid,tolr)

	case("V3DL")

		vcoul1= v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)

	case("V2DT")

		vcoul1= v2dt(kpt1,kpt2,ngrid,rlat,tolr)

	case("V2DT2")

		vcoul1= v2dt2(kpt1,kpt2,ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vcoul1= v2dohono(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vcoul1= v2drk(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)

	case("V0DT")

		vcoul1= v0dt(kpt1,kpt2,ngrid,rlat,tolr)

	case default

		write(*,*) "Wrong Coulomb Potential"
		STOP

	end select

	



	if (est1(1) .eq. est2(1)) then


		matrizelbsetemp= (ec1-ev1) + vcoul1*fermidisteh(ec1,ev1,temp)


	else

	
		call vecconjg(vbc1,w90basis,vbc)

		call vecconjg(vbv1,w90basis,vbv)

		call prodintsq(vbc,vbc2,w90basis,vc)

		call prodintsq(vbv,vbv2,w90basis,vv)
	
		matrizelbsetemp=  vcoul1*vc*vv*fermidisteh(ec1,ev1,temp)


	end if
		



end function matrizelbsetemp

function matrizelbsekqtemp(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,q,rlat,est1,ec1,ev1,vbc1 &
                       ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,temp) !funcao para calcular o elemento de matriz da matriz bse

	implicit none

	character(len=5) :: coultype

	integer,dimension(3) :: ngrid

	double precision,dimension(3,3) :: rlat
	double precision :: ez,w	

	integer :: w90basis
	double precision :: a,vcell1

	integer, dimension(4) :: est1,est2
	double precision :: ec1,ec2,ev1,ev2
	double precision,dimension(3) :: kpt1,kpt2
	double precision,dimension(4) :: q
	double precision,dimension(3) :: vq,v0
	double complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2

	double complex, dimension(w90basis) :: vbc,vbv,vbvkp

	double precision :: tolr
	integer :: ktol
	
	double complex:: matrizelbsekqtemp

	double precision :: modk,modq

	double precision,dimension(3) :: ediel

	double precision :: lc

	double complex :: vc,vv
	double complex :: vcv,vvc

	double precision :: vcoulk,vcoulq

	double precision :: v2dk,vcoul,v3diel,v2dt,v0dt,v2dt2
	double precision :: v2dohono,v2drk

	double precision :: r0
	
	double precision :: temp, fermidisteh

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

	case("V3DL")

		vcoulk= v3diel(kpt1,kpt2,ediel,rlat,ngrid,tolr)
		vcoulq= v3diel(vq,v0,ediel,rlat,ngrid,tolr)

	case("V2DT")

		vcoulk= v2dt(kpt1,kpt2,ngrid,rlat,tolr)
		vcoulq= v2dt(vq,v0,ngrid,rlat,tolr)

	case("V2DT2")

		vcoulk= v2dt2(kpt1,kpt2,ngrid,rlat,lc,tolr)
		vcoulq= v2dt2(vq,v0,ngrid,rlat,lc,tolr)
		
	case("V2DOH")

		vcoulk= v2dohono(kpt1,kpt2,ngrid,rlat,ediel,w,ez,tolr)
		vcoulq= v2dohono(vq,v0,ngrid,rlat,ediel,w,ez,tolr)
		
	case("V2DRK")

		vcoulk= v2drk(kpt1,kpt2,ngrid,rlat,ediel,lc,ez,w,r0,tolr)
		vcoulq= v2drk(vq,v0,ngrid,rlat,ediel,lc,ez,w,r0,tolr)

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

	
		call vecconjg(vbc1,w90basis,vbc)

		call vecconjg(vbv1,w90basis,vbv)

		call prodintsq(vbc,vbc2,w90basis,vc)

		call prodintsq(vbv,vbv2,w90basis,vv)


		matrizelbsekqtemp= vcoulk*vc*vv*fermidisteh(ec1,ev1,temp)



	end if
		


else 


	if (est1(1) .eq. est2(1)) then


		call vecconjg(vbc1,w90basis,vbc)

		call vecconjg(vbv1,w90basis,vbv)

		call vecconjg(vbv2,w90basis,vbvkp)


		call prodintsq(vbc,vbc2,w90basis,vc)

		call prodintsq(vbv,vbv2,w90basis,vv)


		call prodintsq(vbc,vbv1,w90basis,vcv)

		call prodintsq(vbvkp,vbc2,w90basis,vvc)


		matrizelbsekqtemp= (ec1-ev1) + (vcoulk - vcoulq*vcv*vvc)*fermidisteh(ec1,ev1,temp)
			     

	


	else

		call vecconjg(vbc1,w90basis,vbc)

		call vecconjg(vbv1,w90basis,vbv)

		call vecconjg(vbv2,w90basis,vbvkp)


		call prodintsq(vbc,vbc2,w90basis,vc)

		call prodintsq(vbv,vbv2,w90basis,vv)


		call prodintsq(vbc,vbv1,w90basis,vcv)

		call prodintsq(vbvkp,vbc2,w90basis,vvc)

		matrizelbsekqtemp= (vcoulk*vc*vv- vcoulq*vcv*vvc)*fermidisteh(ec1,ev1,temp)
				


	end if





end if





end function matrizelbsekqtemp


subroutine dielbseptemp(nthread,dimse,excitonvec,hopt1,hopt2,fdeh,activity)

	use omp_lib
	implicit none


	integer :: dimse,nthread
	double precision,dimension(dimse) :: activity, exciton
	double complex :: actaux,actaux2
	double complex,dimension(dimse) :: hopt1,hopt2
	double precision,dimension(dimse) :: fdeh
	double complex,dimension(dimse,dimse) :: excitonvec


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






















