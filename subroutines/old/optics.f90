subroutine optrpa(ev,vv,ec,vc,kx,ky,kz,ffactor,sme,&
		  w90basis,nvec,rlat,rvec,hopmatrices,ihopmatrices,&
		    exx,exy,exz,eyy,eyz,ezz)

	implicit none

	integer :: i,j,k

	integer :: w90basis,nvec

	integer,dimension(w90basis) :: ffactor

	double precision :: kx,ky,kz

	double precision :: ev,ec,erpa,sme,factor

	double precision,dimension(3,3) :: rlat

	integer,dimension(nvec,3) :: rvec

	double complex,dimension(w90basis,w90basis) :: hx,hy,hz

	double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	double complex,dimension(w90basis) :: vc,vv

	double complex :: hrx,hry,hrz

	double precision :: exx,exy,exz,eyy,eyz,ezz



	call hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,hx,hy,hz)

	erpa = ec - ev

	
	call sandwich(w90basis,vc,hx,vv,hrx)

	call sandwich(w90basis,vc,hy,vv,hry)

	call sandwich(w90basis,vc,hz,vv,hrz)

	factor = 1./((erpa*erpa)+(sme*sme))

	exx = factor*(hrx*conjg(hrx))
	exy = factor*(hrx*conjg(hry))
	exz = factor*(hrx*conjg(hrz))

	eyy = factor*(hry*conjg(hry))
	eyz = factor*(hry*conjg(hrz))

	ezz = factor*(hrz*conjg(hrz))


end subroutine optrpa


subroutine optsp(ec,ev,sme,vv,vc,kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,ihopmatrices,&
		    hxsp,hysp,hzsp)

	implicit none

	integer :: i,j,k

	integer :: w90basis,nvec

	double precision :: kx,ky,kz
	
	double precision :: ec,ev,sme

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
	
	double complex :: factor,erpa
	




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

	erpa = ec - ev + cmplx(0.0,sme)

	factor = 1.0/(erpa)
	
	hxsp = factor*hxsp
	
	hysp = factor*hysp
	
	hzsp = factor*hzsp

end subroutine optsp


subroutine hlm(kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
		    ihopmatrices,hx,hy,hz)

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

