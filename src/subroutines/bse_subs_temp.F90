
function matrizelbsetemp(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,rlat,est1,ec1,ev1,vbc1 &
         ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,temp,dft,nvec,rvec,sk,skp, &
         use_center_phase,center_phase1,center_phase2, &
         use_center_grad,w90dat,centers) !funcao para calcular o elemento de matriz da matriz bse

	! The zero-temperature kernel matrizelbse (Coulomb potentials, shortest
	! image q, Wannier-centre phases), with the Coulomb term weighted by the
	! occupation differences fv-fc of the two transitions.

	use bse_q_optics, only: rmn_data

	implicit none

	character(len=10) :: coultype
	character(len=1) :: dft

	integer,dimension(3) :: ngrid

	integer :: w90basis,nvec
	real :: ez,w

	real,dimension(nvec,3) :: rvec

	complex,dimension(w90basis,w90basis) :: sk,skp

	integer, dimension(4) :: est1,est2
	real :: ec1,ec2,ev1,ev2
	real,dimension(3) :: kpt1,kpt2
	complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2

	logical :: use_center_phase,use_center_grad
	complex,dimension(w90basis) :: center_phase1,center_phase2
	type(rmn_data),intent(in) :: w90dat
	real,dimension(3,w90basis) :: centers

	real,dimension(3,3) :: rlat

	real :: tolr

	complex:: matrizelbsetemp
	complex:: matrizelbse,kernel

	real,dimension(3) :: ediel

	real :: lc

	real :: r0

	real :: temp, fermidisteh


	kernel= matrizelbse(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,rlat,est1,ec1,ev1,vbc1, &
	         vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,dft,nvec,rvec,sk,skp, &
	         use_center_phase,center_phase1,center_phase2, &
	         use_center_grad,w90dat,centers)

	if (est1(1) .eq. est2(1)) then

		matrizelbsetemp= (ec1-ev1) + (kernel-(ec1-ev1))*fermidisteh(ec1,ev1,temp)

	else

		matrizelbsetemp= kernel*sqrt(fermidisteh(ec1,ev1,temp)*fermidisteh(ec2,ev2,temp))

	end if

end function matrizelbsetemp

function matrizelbsekqtemp(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,q,rlat,est1,ec1,ev1,vbc1 &
                       ,vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,temp,dft,nvec,rvec,sk,skp,skq,skpq, &
                       use_center_phase,centers) !funcao para calcular o elemento de matriz da matriz bse

	! The zero-temperature kernel matrizelbsekq, with the Coulomb terms
	! weighted by the occupation differences fv-fc of the two transitions.

	implicit none

	character(len=10) :: coultype
	character(len=1) :: dft

	integer,dimension(3) :: ngrid

	real,dimension(3,3) :: rlat
	real :: ez,w

	integer :: w90basis,nvec

	real,dimension(nvec,3) :: rvec

	complex,dimension(w90basis,w90basis) :: sk,skq,skp,skpq

	integer, dimension(4) :: est1,est2
	real :: ec1,ec2,ev1,ev2
	real,dimension(3) :: kpt1,kpt2
	real,dimension(4) :: q
	complex, dimension(w90basis) :: vbc1,vbc2,vbv1,vbv2

	logical :: use_center_phase
	real,dimension(3,w90basis) :: centers

	real :: tolr

	complex:: matrizelbsekqtemp
	complex:: matrizelbsekq,kernel

	real,dimension(3) :: ediel

	real :: lc

	real :: r0

	real :: temp, fermidisteh


	kernel= matrizelbsekq(coultype,tolr,w90basis,ediel,lc,ez,w,r0,ngrid,q,rlat,est1,ec1,ev1,vbc1, &
	         vbv1,kpt1,est2,ec2,ev2,vbc2,vbv2,kpt2,dft,nvec,rvec,sk,skp,skq,skpq, &
	         use_center_phase,centers)

	if (est1(1) .eq. est2(1)) then

		matrizelbsekqtemp= (ec1-ev1) + (kernel-(ec1-ev1))*fermidisteh(ec1,ev1,temp)

	else

		matrizelbsekqtemp= kernel*sqrt(fermidisteh(ec1,ev1,temp)*fermidisteh(ec2,ev2,temp))

	end if

end function matrizelbsekqtemp


subroutine dielbseptemp(nthread,dimse,excitonvec,hopt1,hopt2,fdeh,activity)
! <0|r|S> = sum_j conjg(A_j^S) <c|r|v>_j, consistent with the <v2|v1> kernel

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
	!$OMP PARALLEL DO DEFAULT(NONE) &
	!$OMP& SHARED(dimse,excitonvec,hopt1,hopt2,fdeh,activity) &
	!$OMP& PRIVATE(i,j,actaux,actaux2) SCHEDULE(STATIC)
	do i=1,dimse


		actaux=0.0
		actaux2=0.0


		do j=1,dimse

			
			actaux=actaux+(conjg(excitonvec(j,i))*hopt1(j))
			actaux2=actaux2+(conjg(excitonvec(j,i))*hopt2(j)*fdeh(j))
			

			
		end do
			

		activity(i)=actaux*conjg(actaux2)

	
	
		


	end do
	!$OMP END PARALLEL DO



end subroutine dielbseptemp















