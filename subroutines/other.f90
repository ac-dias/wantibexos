subroutine lvl(w90basis,vec,lz)

	implicit none

	integer :: w90basis
	integer :: i,j
	double precision,parameter :: pi=acos(-1.)
	double complex,parameter :: imag = cmplx(0d0,1d0)

	double precision, dimension(3) :: mag !direção do vetor magnetização do substrato

	double complex,dimension(w90basis) :: vec
	double complex,dimension(w90basis) :: vcconj

	double complex,dimension(w90basis) :: zvvaux

	double complex,dimension(w90basis,w90basis) :: splz

	double precision,parameter :: ct1= 1.0
	double precision,parameter :: exc= 1.0

	double complex :: al,aj,ak,am !termos auxiliares

	double complex,dimension(5,5) :: dcor !correção orbitais d
	double complex,dimension(3,3) :: pcor !correção orbitais p
	double complex,dimension(1,1) :: scor !correção obitais s

	double complex :: lz


    	mag(1)= 0.0
	mag(2)= 0.0
	mag(3)= 1.0

	al = (1.0+imag)/sqrt(2.)
	aj = ((-1.0+imag)*(sqrt(3.)))/sqrt(2.)
	ak = (-1+imag)*sqrt(2.)
	am = imag*sqrt(3.)

	!escrevendo dcor

	dcor = 0.0
	
	dcor(1,2) = ct1*exc*dconjg(am)*mag(2)
	dcor(1,3) = ct1*exc*aj*mag(1)

	dcor(2,1) = ct1*exc*am*mag(2)
	dcor(2,3) = ct1*exc*dconjg(al)*mag(3)
	dcor(2,4) = ct1*exc*(-imag)*mag(2)
	dcor(2,5) = ct1*exc*dconjg(al)*mag(1)

	dcor(3,1) = ct1*exc*dconjg(aj)*mag(1)
	dcor(3,2) = ct1*exc*al*mag(3)
	dcor(3,4) = ct1*exc*al*mag(1)
	dcor(3,5) = ct1*exc*(-imag)*mag(2)

	dcor(4,2) = ct1*exc*(imag)*mag(2)
	dcor(4,3) = ct1*exc*dconjg(al)*mag(1)
	dcor(4,5) = ct1*exc*ak*mag(3)

	dcor(5,2) = ct1*exc*al*mag(1)
	dcor(5,3) = ct1*exc*(imag)*mag(2)
	dcor(5,4) = ct1*exc*dconjg(ak)*mag(3)

	!escrevendo pcor

	pcor = 0.0

	pcor(1,2) = ct1*exc*imag*mag(2)
	pcor(1,3) = ct1*exc*imag*mag(1)

	pcor(2,1) = ct1*exc*(-imag)*mag(2)
	pcor(2,3) = ct1*exc*imag*mag(3)

	pcor(3,1) = ct1*exc*(-imag)*mag(1)
	pcor(3,2) = ct1*exc*(-imag)*mag(2)

	!escrevendo scor

	scor = 0.0

	!escrevendo splz

	splz = 0.0

	do i=1,5

	 do j=1,5

	   splz(i,j) = dcor(i,j)
	   splz((w90basis/2)+i,(w90basis/2)+j) = dcor(i,j)

	 end do

	end do

	do i=1,3

	 do j=1,3

	   splz(5+i,5+j) = pcor(i,j)
	   splz(5+3+i,5+3+j) = pcor(i,j)

	   splz((w90basis/2)+5+i,(w90basis/2)+5+j) = pcor(i,j)
	   splz((w90basis/2)+5+3+i,(w90basis/2)+5+3+j) = pcor(i,j)

	 end do

	end do

	call vecconjg(vec,w90basis,vcconj)

	!calculando valor medio de lz

	call matvec(splz,vec,w90basis,zvvaux)

	call prodintsq(vcconj,zvvaux,w90basis,lz)




end subroutine lvl

subroutine exch(w90basis,exc,excl,mag,hexc)

	implicit none

	integer :: i,j

	double precision,parameter :: pi=acos(-1.)

	double complex,parameter :: imag = cmplx(0d0,1d0)
	double precision,parameter :: ct1 = 0.5d0
	double precision,parameter :: ct2 = sqrt(3d0)/(2d0)

	integer :: w90basis
	double precision :: exc,excl
	double precision,dimension(3) :: mag
	double complex,dimension(w90basis,w90basis) :: hexc

	double complex,dimension(3,3) :: xexcuu,xexcdd,xexcud,xexcdu 
	double complex,dimension(5,5) :: mexcuu,mexcdd,mexcud,mexcdu

	double complex,dimension(w90basis/2,w90basis/2) :: hexcuu,hexcud,hexcdu,hexcdd

	double complex :: axy,al,aj,ak,am


	axy = cmplx(mag(1),-mag(2))
	al = cmplx(1./sqrt(2.),1./sqrt(2.))
	aj = cmplx(-sqrt(3.)/sqrt(2.),sqrt(3.)/sqrt(2.))
	ak = cmplx(-sqrt(2.),sqrt(2.))
	am = cmplx(0.0,sqrt(3.))

	xexcuu = 0.0

	xexcuu(1,1) = exc*mag(3)
	xexcuu(1,2) = ct1*excl*imag*mag(2)
	xexcuu(1,3) = ct1*excl*imag*mag(1)
	xexcuu(2,1) = -ct1*excl*imag*mag(2)
	xexcuu(2,2) = exc*mag(3)
	xexcuu(2,3) = ct1*excl*imag*mag(3)
	xexcuu(3,1) = -ct1*excl*imag*mag(1)
	xexcuu(3,2) = -ct1*excl*imag*mag(3)
	xexcuu(3,3) = exc*mag(3)

	xexcud = 0.0

	xexcud(1,1) = exc*axy
	xexcud(2,2) = exc*axy
	xexcud(3,3) = exc*axy

	xexcdu = 0.0

	xexcdu(1,1) = exc*conjg(axy)
	xexcdu(2,2) = exc*conjg(axy)
	xexcdu(3,3) = exc*conjg(axy)

	xexcdd = 0.0

	xexcdd(1,1) = -exc*mag(3)
	xexcdd(1,2) = ct1*excl*imag*mag(2)
	xexcdd(1,3) = ct1*excl*imag*mag(1)
	xexcdd(2,1) = -ct1*excl*imag*mag(2)
	xexcdd(2,2) = -exc*mag(3)
	xexcdd(2,3) = ct1*excl*imag*mag(3)
	xexcdd(3,1) = -ct1*excl*imag*mag(1)
	xexcdd(3,2) = -ct1*excl*imag*mag(3)
	xexcdd(3,3) = -exc*mag(3)


	mexcuu = 0.0

	mexcuu(2-1,2-1) = exc*mag(3)
	mexcuu(2-1,3-1) = ct1*excl*conjg(am)*mag(2)
	mexcuu(2-1,4-1) = ct1*excl*aj*mag(1)
	mexcuu(3-1,2-1) = ct1*excl*am*mag(2)
	mexcuu(3-1,3-1) = exc*mag(3)
	mexcuu(3-1,4-1) = ct1*excl*conjg(al)*mag(3)
	mexcuu(3-1,5-1) = -imag*ct1*excl*mag(2)
	mexcuu(3-1,6-1) = ct1*excl*conjg(al)*mag(1)
	mexcuu(4-1,2-1) = ct1*excl*conjg(aj)*mag(1)
	mexcuu(4-1,3-1) = ct1*excl*al*mag(3)
	mexcuu(4-1,4-1) = exc*mag(3)
	mexcuu(4-1,5-1) = ct1*excl*al*mag(3)
	mexcuu(4-1,6-1) = -imag*ct1*excl*mag(2)
	mexcuu(5-1,3-1) = imag*ct1*excl*mag(2)
	mexcuu(5-1,4-1) = ct1*excl*conjg(al)*mag(1)
	mexcuu(5-1,5-1) = exc*mag(3)
	mexcuu(5-1,6-1) = ct1*excl*ak*mag(3)
	mexcuu(6-1,3-1) = ct1*excl*al*mag(1)
	mexcuu(6-1,4-1) = imag*ct1*excl*mag(2)
	mexcuu(6-1,5-1) = ct1*excl*conjg(ak)*mag(3)
	mexcuu(6-1,6-1) = exc*mag(3)


	mexcud = 0.0

	mexcud(2-1,2-1) = exc*axy
	mexcud(3-1,3-1) = exc*axy
	mexcud(4-1,4-1) = exc*axy
	mexcud(5-1,5-1) = exc*axy
	mexcud(6-1,6-1) = exc*axy

	mexcdu = 0.0

	mexcdu(2-1,2-1) = exc*conjg(axy)
	mexcdu(3-1,3-1) = exc*conjg(axy)
	mexcdu(4-1,4-1) = exc*conjg(axy)
	mexcdu(5-1,5-1) = exc*conjg(axy)
	mexcdu(6-1,6-1) = exc*conjg(axy)


	mexcdd = 0.0

	mexcdd(2-1,2-1) = -exc*mag(3)
	mexcdd(2-1,3-1) = ct1*excl*conjg(am)*mag(2)
	mexcdd(2-1,4-1) = ct1*excl*aj*mag(1)
	mexcdd(3-1,2-1) = ct1*excl*am*mag(2)
	mexcdd(3-1,3-1) = -exc*mag(3)
	mexcdd(3-1,4-1) = ct1*excl*conjg(al)*mag(3)
	mexcdd(3-1,5-1) = -imag*ct1*excl*mag(2)
	mexcdd(3-1,6-1) = ct1*excl*conjg(al)*mag(1)
	mexcdd(4-1,2-1) = ct1*excl*conjg(aj)*mag(1)
	mexcdd(4-1,3-1) = ct1*excl*al*mag(3)
	mexcdd(4-1,4-1) = -exc*mag(3)
	mexcdd(4-1,5-1) = ct1*excl*al*mag(3)
	mexcdd(4-1,6-1) = -imag*ct1*excl*mag(2)
	mexcdd(5-1,3-1) = imag*ct1*excl*mag(2)
	mexcdd(5-1,4-1) = ct1*excl*conjg(al)*mag(1)
	mexcdd(5-1,5-1) = -exc*mag(3)
	mexcdd(5-1,6-1) = ct1*excl*ak*mag(3)
	mexcdd(6-1,3-1) = ct1*excl*al*mag(1)
	mexcdd(6-1,4-1) = imag*ct1*excl*mag(2)
	mexcdd(6-1,5-1) = ct1*excl*conjg(ak)*mag(3)
	mexcdd(6-1,6-1) = -exc*mag(3)


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

	hexc = 0.0


		do i=1,w90basis/2

			do j=1,w90basis/2

			hexc(i,j) = hexcuu(i,j)
			hexc(i,j+(w90basis/2)) = hexcud(i,j)
			hexc(i+(w90basis/2),j) = hexcdu(i,j)
			hexc((w90basis/2)+i,(w90basis/2)+j) = hexcdd(i,j)

			end do

		end do




end subroutine exch

subroutine optdiel_new(rtype,vc,dimse,ngrid,elux,exciton,fosc,sme,rpart,ipart)


	implicit none

	double precision :: elux
	integer :: dimse
	double precision :: sme
	double precision :: rpart,ipart
	integer :: i,j

	integer :: rtype
	double precision :: delta

	integer,dimension(3) :: ngrid

	double precision :: vc

	double precision :: raux,iaux

	double precision,dimension(dimse) :: exciton
	double precision,dimension(dimse) :: fosc
	
	double precision :: caux,eaux

	double precision,parameter :: pi=acos(-1.)

	double precision,parameter :: gama= 180.90!e^2/e0

	double complex,parameter :: imag= cmplx(0.0,1.0)

	double precision :: baux,auxf,aux2,eaux2

	auxf = 1.0 !4.0*pi
	
	aux2 = sme

	baux = real(ngrid(1)*ngrid(2)*ngrid(3))*vc

	caux=gama/(baux*auxf)


	if ( rtype .eq. 1) then

		delta = 1.0
	else 

		delta = 0.0
	end if 


	rpart= 0.0
	ipart= 0.0


	do i=1,dimse

		eaux = ((elux-exciton(i))**2)+(sme**2)
		
		eaux2 = ((elux-exciton(i))**2)
		
		if (eaux2 .eq. 0.0) then
		
		raux = 0.0
				
		else

		raux = fosc(i)*(((exciton(i)-elux))/(eaux))*(sme)
		
		end if


		iaux = fosc(i)*((sme*aux2)/(eaux))		


	

		rpart=rpart+raux

		ipart=ipart+iaux


	end do

	!write(*,*) caux

	ipart= ipart*caux

	rpart= delta+(rpart*caux)


end subroutine optdiel_new

