subroutine hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,hx,hy,hz) !light-matter interaction

	implicit none

	
	integer :: i,j,k

	integer :: w90basis,nvec

	real :: kx,ky,kz

	integer,dimension(nvec) :: ffactor

	real,dimension(3,3) :: rlat

	real,dimension(nvec,3) :: rvec

	complex,dimension(w90basis,w90basis) :: htbx,htby,htbz

	complex,dimension(w90basis,w90basis) :: hx,hy,hz

	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	real,dimension(w90basis,w90basis) :: hopmatrices2,ihopmatrices2

	complex,dimension(nvec) :: kvecs

	complex,dimension(nvec) :: kvecsx,kvecsy,kvecsz

	complex,dimension(w90basis,w90basis) :: basisc,ibasisc !matrizes de mudanca de base

	complex,parameter :: imag=cmplx(0.0,1.0)



		do i=1,nvec

		kvecs(i) = cmplx(0.0,kx*rvec(i,1))+cmplx(0.0,ky*rvec(i,2))+cmplx(0.0,kz*rvec(i,3))

		kvecsx(i) = cmplx(0.0,rvec(i,1))

		kvecsy(i) =  cmplx(0.0,rvec(i,2))

		kvecsz(i) =  cmplx(0.0,rvec(i,3))


		end do


		htbx = 0.0
		htby = 0.0
		htbz = 0.0

		do i=1,nvec

			!if ( (rvec(i,1) .eq. 0) .and. (rvec(i,2) .eq. 0) .and. (rvec(i,3) .eq. 0)) then

			!	hopmatrices2 = hopmatrices(i,:,:)
			!	ihopmatrices2 = ihopmatrices(i,:,:)

			!	do j=1,w90basis

			!		hopmatrices2(j,j) = 0.0
			!		ihopmatrices2(j,j) = 0.0
			!	end do

			!	htbx=htbx+kvecsx(i)*exp(kvecs(i))*(hopmatrices2+imag*ihopmatrices2)*(1.0/real(ffactor(i)))

			!	htby=htby+kvecsy(i)*exp(kvecs(i))*(hopmatrices2+imag*ihopmatrices2)*(1.0/real(ffactor(i)))

			!	htbz=htbz+kvecsz(i)*exp(kvecs(i))*(hopmatrices2+imag*ihopmatrices2)*(1.0/real(ffactor(i)))
				

			!else

				htbx=htbx+kvecsx(i)*exp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/real(ffactor(i)))

				htby=htby+kvecsy(i)*exp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/real(ffactor(i)))

				htbz=htbz+kvecsz(i)*exp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/real(ffactor(i)))


			!end if

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

	real :: kx,ky,kz

	real :: ev,ec,sme,factor

	integer,dimension(nvec) :: ffactor

	complex :: esp

	real,dimension(3,3) :: rlat

	real,dimension(nvec,3) :: rvec

	complex,dimension(w90basis,w90basis) :: hx,hy,hz

	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	complex,dimension(w90basis) :: vc,vv

	complex :: hrx,hry,hrz

	real :: exx,exy,exz,eyy,eyz,ezz



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


subroutine hlm_s(kx,ky,kz,w90basis,nvec,rvec,hopmatrices,ihopmatrices,ovp,tau,&
		 hx,hy,hz,sx,sy,sz)

	implicit none

	integer :: i,a,b
	integer :: w90basis,nvec
	real :: kx,ky,kz
	real,dimension(nvec,3) :: rvec
	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices,ovp
	real,dimension(w90basis,3) :: tau
	complex,dimension(w90basis,w90basis) :: hx,hy,hz,sx,sy,sz
	complex :: ph,hab,dx,dy,dz

	! dH/dk and dS/dk of the non-orthogonal DFT=S basis in the atomic
	! gauge: element ab of the R term carries i(R+tau_b-tau_a), the vector
	! from orbital a to orbital b in cell R, where hlm has iR. Between the
	! lattice-gauge eigenvectors this is the derivative for the
	! coefficients c_a(k) exp(-i k.tau_a) of sandwich_phase.
	hx = 0.0
	hy = 0.0
	hz = 0.0
	sx = 0.0
	sy = 0.0
	sz = 0.0
	do i=1,nvec
		ph = exp(cmplx(0.0,kx*rvec(i,1)+ky*rvec(i,2)+kz*rvec(i,3)))
		do b=1,w90basis
			do a=1,w90basis
				dx = cmplx(0.0,rvec(i,1)+tau(b,1)-tau(a,1))*ph
				dy = cmplx(0.0,rvec(i,2)+tau(b,2)-tau(a,2))*ph
				dz = cmplx(0.0,rvec(i,3)+tau(b,3)-tau(a,3))*ph
				hab = cmplx(hopmatrices(i,a,b),ihopmatrices(i,a,b))
				hx(a,b) = hx(a,b)+dx*hab
				hy(a,b) = hy(a,b)+dy*hab
				hz(a,b) = hz(a,b)+dz*hab
				sx(a,b) = sx(a,b)+dx*ovp(i,a,b)
				sy(a,b) = sy(a,b)+dy*ovp(i,a,b)
				sz(a,b) = sz(a,b)+dz*ovp(i,a,b)
			end do
		end do
	end do

end subroutine hlm_s


subroutine optsp_s(ev,vv,ec,vc,kx,ky,kz,sme,w90basis,nvec,rvec,&
		   hopmatrices,ihopmatrices,ovp,tau,hrx,hry,hrz)

	implicit none

	integer :: w90basis,nvec
	real :: kx,ky,kz
	real :: ev,ec,sme
	real,dimension(nvec,3) :: rvec
	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices,ovp
	real,dimension(w90basis,3) :: tau
	complex,dimension(w90basis) :: vc,vv
	complex :: hrx,hry,hrz

	! optsp for the non-orthogonal DFT=S basis:
	!   hr = <c| dH/dk - (ec+ev)/2 dS/dk |v> / (ec - ev + i sme)
	! with both derivatives in the atomic gauge (hlm_s). For the
	! S-orthonormal eigenvectors of CHEGV, i<c|dH - (ec+ev)/2 dS|v> =
	! (ev-ec) A_cv, A_cv = i<c|S dv> + (i/2)<c|dS|v> the (Hermitian)
	! interband connection, up to the dipoles <a|r-(r_a+r_b)/2|b> between
	! basis orbitals, which optdip_s adds when BSE_CENTER_FILE gives the
	! position matrix of the basis. optsp's lattice-gauge dH/dk alone
	! breaks C3 (eps_xx /= eps_yy and eps_xy /= 0 for h-BN) and lacks dS/dk.
	call optspbz_s(ev,vv,ec,vc,kx,ky,kz,w90basis,nvec,rvec,hopmatrices,&
		       ihopmatrices,ovp,tau,hrx,hry,hrz)
	hrx = hrx/cmplx(ec-ev,sme)
	hry = hry/cmplx(ec-ev,sme)
	hrz = hrz/cmplx(ec-ev,sme)

end subroutine optsp_s


subroutine optdip_s(vv,vc,kx,ky,kz,w90basis,nvec,rvec,dmat,hrx,hry,hrz)

	implicit none

	integer :: i
	integer :: w90basis,nvec
	real :: kx,ky,kz
	real,dimension(nvec,3) :: rvec
	real,dimension(nvec,w90basis,w90basis,3) :: dmat
	complex,dimension(w90basis) :: vc,vv
	complex :: hrx,hry,hrz
	complex :: ph,drx,dry,drz
	complex,dimension(w90basis,w90basis) :: dx,dy,dz

	! Adds to the DFT=S optical vertex hr of optsp_s the dipoles between basis
	! orbitals, i<c|D(k)|v>, D(k) = sum_R exp(ik.R) D(R) with D(R) from
	! rmn_dfts_dipoles (between the lattice-gauge eigenvectors the orbital-
	! centre phases of the atomic gauge cancel). hr is then i r_cv, r_cv the
	! interband dipole of the basis without approximation.
	dx = 0.0
	dy = 0.0
	dz = 0.0
	do i=1,nvec
		ph = exp(cmplx(0.0,kx*rvec(i,1)+ky*rvec(i,2)+kz*rvec(i,3)))
		dx = dx+ph*dmat(i,:,:,1)
		dy = dy+ph*dmat(i,:,:,2)
		dz = dz+ph*dmat(i,:,:,3)
	end do
	call sandwich(w90basis,vc,dx,vv,drx)
	call sandwich(w90basis,vc,dy,vv,dry)
	call sandwich(w90basis,vc,dz,vv,drz)
	hrx = hrx+cmplx(0.0,1.0)*drx
	hry = hry+cmplx(0.0,1.0)*dry
	hrz = hrz+cmplx(0.0,1.0)*drz

end subroutine optdip_s



subroutine optdiel(rtype,vc,dimse,ngrid,elux,exciton,fosc,sme,rpart,ipart)


	implicit none

	real :: elux
	real :: gaussian	
	integer :: dimse
	real :: sme,smeavg,smef
	real :: rpart,ipart
	integer :: i,j

	integer :: rtype
	real :: delta

	integer,dimension(3) :: ngrid

	real :: vc

	real :: raux,iaux

	real,dimension(dimse) :: exciton
	real,dimension(dimse) :: fosc
	
	real :: caux,eaux

	real,parameter :: pi=acos(-1.)

	real,parameter :: gama= 180.90!e^2/e0

	complex,parameter :: imag= cmplx(0.0,1.0)

	real :: baux,auxf

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

subroutine jdoscalc(vc,dimbse,ngrid,elux,exciton,sme,jdos)

	implicit none
	
	real :: vc
	integer :: dimbse
	integer,dimension(3) :: ngrid
	real,dimension(dimbse) :: exciton
	real :: elux,sme,jdos,smeavg,smef
	
	integer :: i
	real :: gaussian2
	real :: baux
	
	call avgsme(dimbse,exciton,smeavg)
	
	smef = smeavg+sme
	
	baux = real(ngrid(1)*ngrid(2)*ngrid(3))
	
	jdos = 0.0
	
	do i=1,dimbse
	
	jdos = jdos + (gaussian2(elux-exciton(i),smef)/baux)
	
	end do
	

end subroutine jdoscalc

function gaussian2(deltaen,sme)

	implicit none
	
	real :: gaussian2
	real :: deltaen,sme
	real,parameter :: pi=acos(-1.)
	
	real :: norm,deno,aux1	
	
	norm = 1.0D0/(sme*sqrt(2.0D0*pi))
	deno = 2.0D0*sme*sme
	aux1 = -1.0D0*((deltaen)*(deltaen))/deno
	
	gaussian2 = norm*exp(aux1)	

end function gaussian2

subroutine optspbz(vv,vc,kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,ihopmatrices,&
		    hxsp,hysp,hzsp)

	implicit none

	integer :: i,j,k

	integer :: w90basis,nvec

	real :: kx,ky,kz

	real,dimension(3,3) :: rlat

	real,dimension(nvec,3) :: rvec

	integer,dimension(nvec) :: ffactor

	complex,dimension(w90basis,w90basis) :: htbx,htby,htbz

	complex,dimension(w90basis,w90basis) :: hx,hy,hz

	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	complex,dimension(w90basis) :: vc,vv

	complex,dimension(w90basis) :: vcconj

	complex,dimension(w90basis) :: vvauxfx,vvauxfy,vvauxfz

	complex :: hxsp,hysp,hzsp



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


subroutine optspbz_s(ev,vv,ec,vc,kx,ky,kz,w90basis,nvec,rvec,hopmatrices,&
		     ihopmatrices,ovp,tau,hxsp,hysp,hzsp)

	implicit none

	integer :: w90basis,nvec
	real :: kx,ky,kz
	real :: ev,ec,ebar
	real,dimension(nvec,3) :: rvec
	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices,ovp
	real,dimension(w90basis,3) :: tau
	complex,dimension(w90basis) :: vc,vv
	complex :: hxsp,hysp,hzsp,sxsp,sysp,szsp
	complex,dimension(w90basis,w90basis) :: hx,hy,hz,sx,sy,sz

	! optspbz for the non-orthogonal DFT=S basis (see optsp_s):
	! <c| dH/dk - (ec+ev)/2 dS/dk |v>, atomic gauge
	call hlm_s(kx,ky,kz,w90basis,nvec,rvec,hopmatrices,ihopmatrices,ovp,tau,&
		   hx,hy,hz,sx,sy,sz)
	ebar = 0.5*(ec+ev)
	call sandwich(w90basis,vc,hx,vv,hxsp)
	call sandwich(w90basis,vc,hy,vv,hysp)
	call sandwich(w90basis,vc,hz,vv,hzsp)
	call sandwich(w90basis,vc,sx,vv,sxsp)
	call sandwich(w90basis,vc,sy,vv,sysp)
	call sandwich(w90basis,vc,sz,vv,szsp)
	hxsp = hxsp-ebar*sxsp
	hysp = hysp-ebar*sysp
	hzsp = hzsp-ebar*szsp

end subroutine optspbz_s

subroutine avgsme(ndim,energyvec,smeavg)

	implicit none
	
	integer :: ndim,i
	real :: smeavg,factor
	real,dimension(ndim) :: energyvec	
	
	factor = 10.0
	
	smeavg = 0.0
	
	do i=2,ndim
	
		smeavg = smeavg+ energyvec(i)-energyvec(i-1)
	
	end do
	
		smeavg = factor*(smeavg/real(ndim))
		
		!write(*,*) smeavg

end subroutine avgsme


subroutine imagrenorm(ndim,nelux,exciton,dielf)

	implicit none
	
	integer :: ndim,nelux,i
	real,dimension(ndim) :: exciton
	real,dimension(nelux,3) :: dielf
	
	real,parameter :: factor=0.15
	integer :: ibegin
	real :: aux,aux2,aux3


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

