subroutine spincor(exc,w90basis,scor)

	implicit none

	integer :: i,j
	integer :: w90basis	

	double precision,parameter :: pi=acos(-1.)
	double complex,parameter :: imag = cmplx(0d0,1d0)

	double precision:: exc !0.0045
	double precision, dimension(3) :: mag !direção do vetor magnetização do substrato 
	double precision,parameter :: gama = pi/4. !angulo no plano xy
	double precision,parameter :: beta = pi/2. !angulo entre o eixo z e o plano xy

	double complex,dimension(w90basis,w90basis) :: scor
	double complex,dimension(w90basis/2,w90basis/2) :: scoruu,scorud,scordu,scordd

    	mag(1)= cos(gama)*sin(beta)
	mag(2)= sin(gama)*sin(beta)
	mag(3)= cos(beta)
	
	scor = 0.0
	scoruu = 0.0
	scorud = 0.0
	scordu = 0.0
	scordd = 0.0

	do i=1,w90basis/2

	scoruu(i,i) = exc*mag(3)
	scordd(i,i) = -exc*mag(3)

	scorud(i,i) = exc*(mag(1)-imag*mag(2))
	scordu(i,i) = exc*(mag(1)+imag*mag(2))

	end do

	do i=1,w90basis/2

		do j=1,w90basis/2

		scor(i,j) = scoruu(i,j)
		scor(i,j+(w90basis/2)) = scorud(i,j)
		scor(i+(w90basis/2),j) = scordu(i,j)
		scor((w90basis/2)+i,(w90basis/2)+j) = scordd(i,j)

		end do

	end do


end subroutine spincor

subroutine syscor(w90basis,hcor)

	implicit none

	integer :: i,j

	double precision,parameter :: pi=acos(-1.)

	double complex,parameter :: imag = cmplx(0d0,1d0)
	double precision,parameter :: ct1 = 0.5d0
	double precision,parameter :: ct2 = sqrt(3d0)/(2d0)

	integer :: w90basis

	double complex,dimension(w90basis,w90basis) :: hcor

	double complex,dimension(3,3) :: xexcuu,xexcdd,xexcud,xexcdu 
	double complex,dimension(5,5) :: mexcuu,mexcdd,mexcud,mexcdu

	double complex,dimension(w90basis/2,w90basis/2) :: hexcuu,hexcud,hexcdu,hexcdd

	double complex :: axy,al,aj,ak,am

	!termos da correção de desvio de singularidade magnética

	double precision,parameter:: exc2=  1.1999999999999999E-003  !0.0045
	double precision,parameter :: excl2= 1.1999999999999999E-003
	double precision, dimension(3) :: mag2 !direção do vetor magnetização do substrato 

	double precision,parameter :: gama2 = pi/4. !angulo no plano xy
	double precision,parameter :: beta2 = pi/2. !angulo entre o eixo z e o plano xy

	double complex :: axy2

    	mag2(1)= cos(gama2)*sin(beta2)
	mag2(2)= sin(gama2)*sin(beta2)
	mag2(3)= cos(beta2)

	al = cmplx(1./sqrt(2.),1./sqrt(2.))
	aj = cmplx(-sqrt(3.)/sqrt(2.),sqrt(3.)/sqrt(2.))
	ak = cmplx(-sqrt(2.),sqrt(2.))
	am = cmplx(0.0,sqrt(3.))

	axy2 = cmplx(mag2(1),-mag2(2))

	xexcuu = 0.0

	xexcuu(1,1) = exc2*mag2(3)
	xexcuu(1,2) = ct1*excl2*imag*mag2(2)
	xexcuu(1,3) = ct1*excl2*imag*mag2(1)
	xexcuu(2,1) = -ct1*excl2*imag*mag2(2)
	xexcuu(2,2) = exc2*mag2(3)
	xexcuu(2,3) = ct1*excl2*imag*mag2(3)
	xexcuu(3,1) = -ct1*excl2*imag*mag2(1)
	xexcuu(3,2) = -ct1*excl2*imag*mag2(3)
	xexcuu(3,3) = exc2*mag2(3)

	xexcud = 0.0

	xexcud(1,1) = exc2*axy2
	xexcud(2,2) = exc2*axy2
	xexcud(3,3) = exc2*axy2

	xexcdu = 0.0

	xexcdu(1,1) = exc2*conjg(axy2)
	xexcdu(2,2) = exc2*conjg(axy2)
	xexcdu(3,3) = exc2*conjg(axy2)

	xexcdd = 0.0

	xexcdd(1,1) = -exc2*mag2(3)
	xexcdd(1,2) = ct1*excl2*imag*mag2(2)
	xexcdd(1,3) = ct1*excl2*imag*mag2(1)
	xexcdd(2,1) = -ct1*excl2*imag*mag2(2)
	xexcdd(2,2) = -exc2*mag2(3)
	xexcdd(2,3) = ct1*excl2*imag*mag2(3)
	xexcdd(3,1) = -ct1*excl2*imag*mag2(1)
	xexcdd(3,2) = -ct1*excl2*imag*mag2(3)
	xexcdd(3,3) = -exc2*mag2(3)


	mexcuu = 0.0

	mexcuu(2-1,2-1) = exc2*mag2(3)
	mexcuu(2-1,3-1) = ct1*excl2*conjg(am)*mag2(2)
	mexcuu(2-1,4-1) = ct1*excl2*aj*mag2(1)
	mexcuu(3-1,2-1) = ct1*excl2*am*mag2(2)
	mexcuu(3-1,3-1) = exc2*mag2(3)
	mexcuu(3-1,4-1) = ct1*excl2*conjg(al)*mag2(3)
	mexcuu(3-1,5-1) = -imag*ct1*excl2*mag2(2)
	mexcuu(3-1,6-1) = ct1*excl2*conjg(al)*mag2(1)
	mexcuu(4-1,2-1) = ct1*excl2*conjg(aj)*mag2(1)
	mexcuu(4-1,3-1) = ct1*excl2*al*mag2(3)
	mexcuu(4-1,4-1) = exc2*mag2(3)
	mexcuu(4-1,5-1) = ct1*excl2*al*mag2(3)
	mexcuu(4-1,6-1) = -imag*ct1*excl2*mag2(2)
	mexcuu(5-1,3-1) = imag*ct1*excl2*mag2(2)
	mexcuu(5-1,4-1) = ct1*excl2*conjg(al)*mag2(1)
	mexcuu(5-1,5-1) = exc2*mag2(3)
	mexcuu(5-1,6-1) = ct1*excl2*ak*mag2(3)
	mexcuu(6-1,3-1) = ct1*excl2*al*mag2(1)
	mexcuu(6-1,4-1) = imag*ct1*excl2*mag2(2)
	mexcuu(6-1,5-1) = ct1*excl2*conjg(ak)*mag2(3)
	mexcuu(6-1,6-1) = exc2*mag2(3)


	mexcud = 0.0

	mexcud(2-1,2-1) = exc2*axy2
	mexcud(3-1,3-1) = exc2*axy2
	mexcud(4-1,4-1) = exc2*axy2
	mexcud(5-1,5-1) = exc2*axy2
	mexcud(6-1,6-1) = exc2*axy2

	mexcdu = 0.0

	mexcdu(2-1,2-1) = exc2*conjg(axy2)
	mexcdu(3-1,3-1) = exc2*conjg(axy2)
	mexcdu(4-1,4-1) = exc2*conjg(axy2)
	mexcdu(5-1,5-1) = exc2*conjg(axy2)
	mexcdu(6-1,6-1) = exc2*conjg(axy2)


	mexcdd = 0.0

	mexcdd(2-1,2-1) = -exc2*mag2(3)
	mexcdd(2-1,3-1) = ct1*excl2*conjg(am)*mag2(2)
	mexcdd(2-1,4-1) = ct1*excl2*aj*mag2(1)
	mexcdd(3-1,2-1) = ct1*excl2*am*mag2(2)
	mexcdd(3-1,3-1) = -exc2*mag2(3)
	mexcdd(3-1,4-1) = ct1*excl2*conjg(al)*mag2(3)
	mexcdd(3-1,5-1) = -imag*ct1*excl2*mag2(2)
	mexcdd(3-1,6-1) = ct1*excl2*conjg(al)*mag2(1)
	mexcdd(4-1,2-1) = ct1*excl2*conjg(aj)*mag2(1)
	mexcdd(4-1,3-1) = ct1*excl2*al*mag2(3)
	mexcdd(4-1,4-1) = -exc2*mag2(3)
	mexcdd(4-1,5-1) = ct1*excl2*al*mag2(3)
	mexcdd(4-1,6-1) = -imag*ct1*excl2*mag2(2)
	mexcdd(5-1,3-1) = imag*ct1*excl2*mag2(2)
	mexcdd(5-1,4-1) = ct1*excl2*conjg(al)*mag2(1)
	mexcdd(5-1,5-1) = -exc2*mag2(3)
	mexcdd(5-1,6-1) = ct1*excl2*ak*mag2(3)
	mexcdd(6-1,3-1) = ct1*excl2*al*mag2(1)
	mexcdd(6-1,4-1) = imag*ct1*excl2*mag2(2)
	mexcdd(6-1,5-1) = ct1*excl2*conjg(ak)*mag2(3)
	mexcdd(6-1,6-1) = -exc2*mag2(3)


	hexcuu = 0.0

	do i=1,5
	
		do j=1,5

			hexcuu(i,j) = mexcuu(i,j)	

		end do


	end do


	do i=1,3
	
		do j=1,3

			hexcuu(5+i,5+j) = xexcuu(i,j)

			hexcuu(8+i,8+j) = xexcuu(i,j)


		end do


	end do

	hexcud = 0.0

	do i=1,5
	
		do j=1,5

			hexcud(i,j) = mexcud(i,j)	

		end do


	end do


	do i=1,3
	
		do j=1,3

			hexcud(5+i,5+j) = xexcud(i,j)

			hexcud(8+i,8+j) = xexcud(i,j)
		

		end do


	end do

	hexcdu = 0.0

	do i=1,5
	
		do j=1,5

			hexcdu(i,j) = mexcdu(i,j)	

		end do


	end do


	do i=1,3
	
		do j=1,3

			hexcdu(5+i,5+j) = xexcdu(i,j)

			hexcdu(8+i,8+j) = xexcdu(i,j)


		end do


	end do


	hexcdd = 0.0

	do i=1,5
	
		do j=1,5

			hexcdd(i,j) = mexcdd(i,j)
	

		end do


	end do


	do i=1,3
	
		do j=1,3

			hexcdd(5+i,5+j) = xexcdd(i,j)

			hexcdd(8+i,8+j) = xexcdd(i,j)

		end do


	end do

	hcor = 0.0


		do i=1,w90basis/2

			do j=1,w90basis/2

			hcor(i,j) = hexcuu(i,j)
			hcor(i,j+(w90basis/2)) = hexcud(i,j)
			hcor(i+(w90basis/2),j) = hexcdu(i,j)
			hcor((w90basis/2)+i,(w90basis/2)+j) = hexcdd(i,j)

			end do

		end do




end subroutine syscor
