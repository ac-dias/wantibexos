
subroutine eigsys(nthread,scs,exc,nocp,ffactor,kx,ky,kz,w90basis,nvec,rlat,rvec,hopmatrices,&
		   ihopmatrices,efermi,energias,autovetores)

	use omp_lib
        !use f95_lapack

	implicit none

	integer :: i,j,k,nthread

	integer :: w90basis,nvec

	double precision :: kx,ky,kz

	!double precision :: exc,exc2

	double precision :: scs,exc

	integer :: nocp

	double precision,dimension(3,3) :: rlat

	integer,dimension(nvec,3) :: rvec

	integer,dimension(w90basis) :: ffactor

	!double precision,dimension(3) :: mag

	double complex,dimension(w90basis,w90basis) :: htb,fermimt,scor,hexc

	double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	double complex,dimension(nvec) :: kvecs

	double precision :: efermi

	double precision,dimension(w90basis) :: energias

	double complex,dimension(w90basis,w90basis) :: autovetores

	double complex,parameter :: imag=cmplx(0.0,1.0)

	!definicoes diagonalizacao
	INTEGER   ::       ifail
	double precision,parameter :: ABSTOL=1.0e-6
        INTEGER          LWMAX
	INTEGER ::         LIWORK, LRWORK
        INTEGER  ::        INFO, LWORK
        DOUBLE PRECISION RWORK( 1 + 5*w90basis + 2*w90basis**2 )
        double complex, dimension (2*w90basis+w90basis**2) :: WORK
	INTEGER,dimension(3 + 5*w90basis) :: IWORK

      	LWORK = 2*w90basis+w90basis**2
      	LIWORK = 3 + 5*w90basis
      	LRWORK = 1 + 5*w90basis + 2*w90basis**2
	LWMAX= 2*w90basis-1

	call OMP_SET_NUM_THREADS(nthread)


		do i=1,nvec

		kvecs(i) = cmplx(0.0,kx*(rvec(i,1)*rlat(1,1)+rvec(i,2)*rlat(2,1)+rvec(i,3)*rlat(3,1))+&
			   ky*(rvec(i,1)*rlat(1,2)+rvec(i,2)*rlat(2,2)+rvec(i,3)*rlat(3,2))+&
			   kz*(rvec(i,1)*rlat(1,3)+rvec(i,2)*rlat(2,3)+rvec(i,3)*rlat(3,3)))

		end do


		htb = 0.0

		do i=1,nvec

			htb=htb+cdexp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/dble(ffactor(i)))

		end do

		call fermilvl(w90basis,efermi,fermimt)

		!call syscor(w90basis,scor)

		!call exch(w90basis,exc,exc2,mag,hexc)

		!call spincor(exc,w90basis,scor)

		htb = htb+fermimt!+scor+hexc

       		!call LA_HEEVR( htb, energias, JOBZ='V', UPLO='U', ABSTOL=ABSTOL, INFO=INFO)

      		!WRITE(*,*)'ZHEEV Example Program Results'
      		LWORK = -1
      		CALL ZHEEV( 'Vectors', 'Lower', w90basis, htb, w90basis, energias, WORK, LWORK, RWORK, INFO )
      		LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      		CALL ZHEEV( 'Vectors', 'Lower', w90basis, htb, w90basis, energias, WORK, LWORK, RWORK,INFO )
		!CALL ZHEEVD('V','U',w90basis,htb,w90basis,energias,WORK, LWORK,&
		!		RWORK,LRWORK,IWORK,INFO)
      		IF( INFO.GT. 0 ) THEN
         	WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         	STOP
      		END IF

		do i=1,w90basis

			do j=1,w90basis

				autovetores(i,j) = htb(j,i)
			end do

		end do

		do i=1,w90basis

			if (energias(i) .gt. 0.0 ) then

				nocp = i-1
				EXIT

				

			else

			continue

			end if

		end do

		do i=nocp+1,w90basis

			energias(i)=energias(i)+scs

		end do

	

end subroutine eigsys

subroutine eigsys2(nthread,scs,exc,nocp,ffactor,kx,ky,kz,w90basis,nvec,rlat,rvec,hopmatrices,&
		   ihopmatrices,efermi,energias,autovetores)

	use omp_lib
        !use f95_lapack

	implicit none

	integer :: i,j,k,nthread

	integer :: w90basis,nvec

	double precision :: kx,ky,kz

	!double precision :: exc,exc2

	double precision :: scs,exc

	integer :: nocp

	double precision,dimension(3,3) :: rlat

	integer,dimension(nvec,3) :: rvec

	integer,dimension(w90basis) :: ffactor

	!double precision,dimension(3) :: mag

	double complex,dimension(w90basis,w90basis) :: htb,fermimt,scor,hexc

	double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices

	double complex,dimension(nvec) :: kvecs

	double precision :: efermi

	double precision,dimension(w90basis) :: energias

	double complex,dimension(w90basis,w90basis) :: autovetores

	double complex,parameter :: imag=cmplx(0.0,1.0)

	!definicoes diagonalizacao
	INTEGER   ::       ifail
	double precision,parameter :: ABSTOL=1.0e-6
        INTEGER          LWMAX
	INTEGER ::         LIWORK, LRWORK
        INTEGER  ::        INFO, LWORK
        DOUBLE PRECISION RWORK( 1 + 5*w90basis + 2*w90basis**2 )
        double complex, dimension (2*w90basis+w90basis**2) :: WORK
	INTEGER,dimension(3 + 5*w90basis) :: IWORK

      	LWORK = 2*w90basis+w90basis**2
      	LIWORK = 3 + 5*w90basis
      	LRWORK = 1 + 5*w90basis + 2*w90basis**2
	LWMAX= 2*w90basis-1

	!call OMP_SET_NUM_THREADS(nthread)


		do i=1,nvec

		kvecs(i) = cmplx(0.0,kx*(rvec(i,1)*rlat(1,1)+rvec(i,2)*rlat(2,1)+rvec(i,3)*rlat(3,1))+&
			   ky*(rvec(i,1)*rlat(1,2)+rvec(i,2)*rlat(2,2)+rvec(i,3)*rlat(3,2))+&
			   kz*(rvec(i,1)*rlat(1,3)+rvec(i,2)*rlat(2,3)+rvec(i,3)*rlat(3,3)))

		end do


		htb = 0.0

		do i=1,nvec

			htb=htb+cdexp(kvecs(i))*(hopmatrices(i,:,:)+imag*ihopmatrices(i,:,:))*(1.0/dble(ffactor(i)))

		end do

		call fermilvl(w90basis,efermi,fermimt)

		!call syscor(w90basis,scor)

		!call exch(w90basis,exc,exc2,mag,hexc)

		!call spincor(exc,w90basis,scor)

		htb = htb+fermimt!+scor+hexc

       		!call LA_HEEVR( htb, energias, JOBZ='V', UPLO='U', ABSTOL=ABSTOL, INFO=INFO)

      		!WRITE(*,*)'ZHEEV Example Program Results'
      		LWORK = -1
      		CALL ZHEEV( 'Vectors', 'Lower', w90basis, htb, w90basis, energias, WORK, LWORK, RWORK, INFO )
      		LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      		CALL ZHEEV( 'Vectors', 'Lower', w90basis, htb, w90basis, energias, WORK, LWORK, RWORK,INFO )
		!CALL ZHEEVD('V','U',w90basis,htb,w90basis,energias,WORK, LWORK,&
		!		RWORK,LRWORK,IWORK,INFO)
      		IF( INFO.GT. 0 ) THEN
         	WRITE(*,*)'The algorithm failed to compute eigenvalues.', INFO
         	!STOP
      		END IF

		do i=1,w90basis

			do j=1,w90basis

				autovetores(i,j) = htb(j,i)
			end do

		end do

		do i=1,w90basis

			if (energias(i) .gt. 0.0 ) then

				nocp = i-1
				EXIT

				

			else

			continue

			end if

		end do

		do i=nocp+1,w90basis

			energias(i)=energias(i)+scs

		end do

	

end subroutine eigsys2






subroutine spinvl(w90basis,vec,dft,systype,spz)


	implicit none

	integer :: w90basis,w2

	double complex,dimension(w90basis) :: vec
	double complex,dimension(w90basis) :: vcconj

	double complex,dimension(w90basis) :: zvvaux

	double complex,dimension(w90basis,w90basis) :: spmz

	double complex :: spz

	integer :: i,j
	
	character(len=1) :: dft
	character(len=4) :: systype	

	spmz=0.0

	w2=w90basis/2
	
	if (systype .eq. "NP") then 
	
	 do i=1,w90basis


		spmz(i,i) = cmplx(1.,0.)


	 end do
	
	else if (systype .eq. "SP") then
	
	 do i=1,w2


		spmz(i,i) = cmplx(1.,0.)
		spmz(w2+i,w2+i) = cmplx(-1.,0.)


	 end do
	
	
	else if ((systype .eq. "SOC") .and. (dft .eq. "V") ) then
	
	 do i=1,w2


		spmz(i,i) = cmplx(1.,0.)
		spmz(w2+i,w2+i) = cmplx(-1.,0.)


	 end do
	
	else
	
	 do i=1,w2


		spmz((2*i)-1,(2*i)-1) = cmplx(1.,0.)
		spmz(2*i,2*i) = cmplx(-1.,0.)


	 end do
	
	
	end if
	
	
	


	call vecconjg(vec,w90basis,vcconj)

	!calculando valor medio de sz

	call matvec(spmz,vec,w90basis,zvvaux)

	call prodintsq(vcconj,zvvaux,w90basis,spz)




end subroutine spinvl


subroutine fermilvl(w90basis,efermi,fermimt)

	implicit none

	integer :: w90basis,i
	double precision :: efermi
	double complex,dimension(w90basis,w90basis) :: fermimt

	fermimt = cmplx(0.0,0.0)

	do i=1,w90basis

		fermimt(i,i) = cmplx(-efermi,0.0)

	end do


end subroutine fermilvl



