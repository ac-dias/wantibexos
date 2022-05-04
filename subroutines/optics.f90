subroutine hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,hx,hy,hz) !light-matter interaction

	implicit none

	
	integer :: i,j,k

	integer :: w90basis,nvec

	double precision :: kx,ky,kz

	integer,dimension(w90basis) :: ffactor

	double precision,dimension(3,3) :: rlat

	integer,dimension(nvec,3) :: rvec

	double complex,dimension(w90basis,w90basis) :: htbx,htby,htbz

	double complex,dimension(w90basis,w90basis) :: hx,hy,hz

	double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	double precision,dimension(w90basis,w90basis) :: hopmatrices2,ihopmatrices2

	double complex,dimension(nvec) :: kvecs

	double complex,dimension(nvec) :: kvecsx,kvecsy,kvecsz

	double complex,dimension(w90basis,w90basis) :: basisc,ibasisc !matrizes de mudanca de base

	double complex,parameter :: imag=cmplx(0.0,1.0)



		do i=1,nvec

		kvecs(i) = cmplx(0.0,kx*(rvec(i,1)*rlat(1,1)+rvec(i,2)*rlat(2,1)+rvec(i,3)*rlat(3,1))+&
			   ky*(rvec(i,1)*rlat(1,2)+rvec(i,2)*rlat(2,2)+rvec(i,3)*rlat(3,2))+&
			   kz*(rvec(i,1)*rlat(1,3)+rvec(i,2)*rlat(2,3)+rvec(i,3)*rlat(3,3)))

		kvecsx(i) = cmplx(0.0,(rvec(i,1)*rlat(1,1)+rvec(i,2)*rlat(2,1)+rvec(i,3)*rlat(3,1)))

		kvecsy(i) = cmplx(0.0,(rvec(i,1)*rlat(1,2)+rvec(i,2)*rlat(2,2)+rvec(i,3)*rlat(3,2)))

		kvecsz(i) = cmplx(0.0,(rvec(i,1)*rlat(1,3)+rvec(i,2)*rlat(2,3)+rvec(i,3)*rlat(3,3)))


		end do


		htbx = 0.0
		htby = 0.0
		htbz = 0.0

		do i=1,nvec

			if ( (rvec(i,1) .eq. 0) .and. (rvec(i,2) .eq. 0) .and. (rvec(i,3) .eq. 0)) then

				hopmatrices2 = hopmatrices(i,:,:)
				ihopmatrices2 = ihopmatrices(i,:,:)

				do j=1,w90basis

					hopmatrices2(j,j) = 0.0
					ihopmatrices2(j,j) = 0.0
				end do

				htbx=htbx+kvecsx(i)*cdexp(kvecs(i))*(hopmatrices2+imag*ihopmatrices2)*(1.0/dble(ffactor(i)))

				htby=htby+kvecsy(i)*cdexp(kvecs(i))*(hopmatrices2+imag*ihopmatrices2)*(1.0/dble(ffactor(i)))

				htbz=htbz+kvecsz(i)*cdexp(kvecs(i))*(hopmatrices2+imag*ihopmatrices2)*(1.0/dble(ffactor(i)))
				

			else

				htbx=htbx+kvecsx(i)*cdexp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/dble(ffactor(i)))

				htby=htby+kvecsy(i)*cdexp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/dble(ffactor(i)))

				htbz=htbz+kvecsz(i)*cdexp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/dble(ffactor(i)))


			end if

		end do


		hx=htbx
		hy=htby
		hz=htbz

	

end subroutine hlm




subroutine optsp(ev,vv,ec,vc,kx,ky,kz,ffactor,sme,&
		  w90basis,nvec,rlat,rvec,hopmatrices,ihopmatrices,&
		    hrx,hry,hrz)

	implicit none

	integer :: i,j,k

	integer :: w90basis,nvec

	double precision :: kx,ky,kz

	double precision :: ev,ec,sme,factor

	integer,dimension(w90basis) :: ffactor

	double complex :: esp

	double precision,dimension(3,3) :: rlat

	integer,dimension(nvec,3) :: rvec

	double complex,dimension(w90basis,w90basis) :: hx,hy,hz

	double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	double complex,dimension(w90basis) :: vc,vv

	double complex :: hrx,hry,hrz

	double precision :: exx,exy,exz,eyy,eyz,ezz



	call hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,hx,hy,hz)

	esp = ec - ev + cmplx(0.0,sme)

	
	call sandwich(w90basis,vc,hx,vv,hrx)

	call sandwich(w90basis,vc,hy,vv,hry)

	call sandwich(w90basis,vc,hz,vv,hrz)

	factor = 1./esp

	hrx = factor*hrx
	hry = factor*hry
	hrz = factor*hrz



end subroutine optsp



subroutine optdiel(rtype,vc,dimse,ngrid,elux,exciton,fosc,sme,rpart,ipart)


	implicit none

	double precision :: elux
	double precision :: gaussian	
	integer :: dimse
	double precision :: sme,smeavg,smef
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

	double precision :: baux,auxf

	auxf = 1.0 !4.0*pi

	baux = real(ngrid(1)*ngrid(2)*ngrid(3))*vc

	caux=gama/(baux*auxf)


	if ( rtype .eq. 1) then

		delta = 1.0
	else 

		delta = 0.0
	end if 


	rpart= 0.0
	ipart= 0.0

	call avgsme(dimse,exciton,smeavg)

	smef = smeavg+sme
	
	!write(*,*) sme,smeavg,smef



	do i=1,dimse

		eaux = ((elux-exciton(i))**2)+(smef**2)


		raux = fosc(i)*((exciton(i)-elux)/(eaux))
								
		iaux = fosc(i)*((smef)/(eaux))
				
		rpart=rpart+raux

		ipart=ipart+iaux


	end do

	!write(*,*) caux

	ipart= ipart*caux

	rpart= delta+(rpart*caux)


end subroutine optdiel

subroutine optspbz(vv,vc,kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,ihopmatrices,&
		    hxsp,hysp,hzsp)

	implicit none

	integer :: i,j,k

	integer :: w90basis,nvec

	double precision :: kx,ky,kz

	double precision,dimension(3,3) :: rlat

	integer,dimension(nvec,3) :: rvec

	integer,dimension(w90basis) :: ffactor

	double complex,dimension(w90basis,w90basis) :: htbx,htby,htbz

	double complex,dimension(w90basis,w90basis) :: hx,hy,hz

	double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	double complex,dimension(w90basis) :: vc,vv

	double complex,dimension(w90basis) :: vcconj

	double complex,dimension(w90basis) :: vvauxfx,vvauxfy,vvauxfz

	double complex :: hxsp,hysp,hzsp



	call hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,hx,hy,hz)


	call vecconjg(vc,w90basis,vcconj)

	!parte total

	call matvec(hx,vv,w90basis,vvauxfx)
	call prodintsq(vcconj,vvauxfx,w90basis,hxsp)

	call matvec(hy,vv,w90basis,vvauxfy)
	call prodintsq(vcconj,vvauxfy,w90basis,hysp)

	call matvec(hz,vv,w90basis,vvauxfz)
	call prodintsq(vcconj,vvauxfz,w90basis,hzsp)

	!termino parte total



end subroutine optspbz

subroutine avgsme(ndim,energyvec,smeavg)

	implicit none
	
	integer :: ndim,i
	double precision :: smeavg,factor
	double precision,dimension(ndim) :: energyvec	
	
	factor = 10.0
	
	smeavg = 0.0
	
	do i=2,ndim
	
		smeavg = smeavg+ energyvec(i)-energyvec(i-1)
	
	end do
	
		smeavg = factor*(smeavg/dble(ndim))
		
		!write(*,*) smeavg

end subroutine avgsme


subroutine imagrenorm(ndim,nelux,exciton,dielf)

	implicit none
	
	integer :: ndim,nelux,i
	double precision,dimension(ndim) :: exciton
	double precision,dimension(nelux,3) :: dielf
	
	double precision,parameter :: factor=0.15
	integer :: ibegin
	double precision :: aux,aux2,aux3


	aux = exciton(1)-factor

	if (aux .le. 0.000) go to 131
	
	ibegin = 1
		
	do i=1,nelux
	
	if (dielf(i,1) .le. aux ) then
	
	 ibegin = i
	
	else
	
	go to 130
	
	end if
	
	end do	
	
130 	continue	
		
	aux3 = dielf(ibegin+1,3)
		
	do i=1,ibegin
	
		dielf(i,3) = 0.000
	
	end do

	!write(*,*) ibegin,dielf(ibegin+1,3)	
	
	do i=ibegin+1,nelux
	
			aux2 = dielf(i,3) - aux3
			dielf(i,3) = aux2
			
			if (dielf(i,3) .lt. 0d0 ) then
				dielf(i,3) = 0d0
			else 
				continue
			end if
	
	end do
	

	!write(*,*) ibegin,dielf(ibegin+1,3)

131	continue	

end subroutine imagrenorm

