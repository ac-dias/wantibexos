
subroutine eigsys(nthread,dft,systype,scs,exc,nocp,ffactor,kx,ky,kz,w90basis,nvec,rlat,rvec,hopmatrices,&
		   ihopmatrices,ovp,efermi,energias,autovetores,nocpf,fermishift,mag)

	use omp_lib
        !use f95_lapack

	implicit none

	integer :: i,j,k,nthread

	integer :: w90basis,nvec

	double precision :: kx,ky,kz

	double precision :: scs,exc
	
	double precision,dimension(3) :: mag

	integer :: nocp

	double precision,dimension(3,3) :: rlat

	double precision,dimension(nvec,3) :: rvec

	integer,dimension(nvec) :: ffactor


	double complex,dimension(w90basis,w90basis) :: htb,fermimt,ovptb
	double complex,dimension(w90basis,w90basis) :: hexc
	

	double precision,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices
	double precision,dimension(nvec,w90basis,w90basis) :: ovp


	double precision :: efermi,fermishift

	double precision,dimension(w90basis) :: energias

	double complex,dimension(w90basis,w90basis) :: autovetores

	double complex,parameter :: imag=cmplx(0.0,1.0)
	
	integer :: nocpf
	
	character(len=1) :: dft
	character(len=4) :: systype

	!definicoes diagonalizacao
	INTEGER   ::       ifail
	double precision,parameter :: ABSTOL=1.0e-6
	INTEGER ::         LIWORK, LRWORK
        INTEGER  ::        INFO, LWORK,LWMAX
        DOUBLE PRECISION,dimension(1 + 5*w90basis + 2*w90basis**2) :: RWORK
        DOUBLE COMPLEX,dimension(2*w90basis + w90basis**2) :: WORK
	INTEGER,dimension(3 + 5*w90basis) :: IWORK

	LWMAX = (2*w90basis-1)**2

	!call OMP_SET_NUM_THREADS(nthread)
	!call MKL_SET_NUM_THREADS(1)

		! $omp single
		select case (dft)
		
		case ("S")
		
		 
		 call hamiltonian(w90basis,nvec,rvec,hopmatrices,ihopmatrices,ffactor,kx,ky,kz,htb)
		 call overlap(w90basis,nvec,rvec,ovp,kx,ky,kz,ovptb)
		
		case default
		
		 call hamiltonian(w90basis,nvec,rvec,hopmatrices,ihopmatrices,ffactor,kx,ky,kz,htb)
		 

		
		end select
		
		

		call fermilvl(w90basis,efermi+fermishift,fermimt)

		
		
		if (abs(exc) .gt. 0.0 ) then
		
		!write(*,*) "before exc"
		
		 call spin_exchange(w90basis,dft,systype,exc,mag,ovptb,hexc)
		 
		 !write(*,*) "after exc" 
		 
		
		else
		
		 hexc = 0.0
		
		end if
		
	       

		htb = htb+hexc+fermimt

		
  
  		select case (dft)
  		
  		case ("S") 
  		
  		 !LWORK = -1
      		 !CALL ZHEGV( 1,'Vectors', 'Lower', w90basis, htb, w90basis,ovptb,w90basis, energias, WORK, LWORK, RWORK, INFO )
      		 !LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      		 !CALL ZHEGV( 1,'Vectors', 'Lower', w90basis, htb, w90basis,ovptb,w90basis, energias, WORK, LWORK, RWORK, INFO )

      		 !IF( INFO.GT. 0 ) THEN
         	 ! WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         	 ! STOP
      		 !END IF
  		
  		 LWORK = -1
      		 LIWORK = -1
      		 LRWORK = -1
  		 CALL ZHEGVD(1,'V','Lower',w90basis,htb,w90basis,ovptb,w90basis,energias,WORK, LWORK,&
				RWORK,LRWORK,IWORK,LIWORK,INFO)
		 LWORK = MIN( 2*w90basis + w90basis**2, INT( WORK( 1 ) ) )
      		 LRWORK = MIN( 1 + 5*w90basis + 2*w90basis**2, INT( RWORK( 1 ) ) )
      		 LIWORK = MIN( 3 + 5*w90basis, IWORK( 1 ) )
		 CALL ZHEGVD(1,'V','Lower',w90basis,htb,w90basis,ovptb,w90basis,energias,WORK, LWORK,&
				RWORK,LRWORK,IWORK,LIWORK,INFO) 
				
		 IF( INFO.GT. 0 ) THEN
         	  WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         	  STOP
      		 END IF		
  		
  		case default  
  		
  		  		
      		
      		 LWORK = -1
      		 LIWORK = -1
      		 LRWORK = -1
		 CALL ZHEEVD('V','Lower',w90basis,htb,w90basis,energias,WORK, LWORK,&
				RWORK,LRWORK,IWORK,LIWORK,INFO)
		 LWORK = MIN( 2*w90basis + w90basis**2, INT( WORK( 1 ) ) )
      		 LRWORK = MIN( 1 + 5*w90basis + 2*w90basis**2, INT( RWORK( 1 ) ) )
      		 LIWORK = MIN( 3 + 5*w90basis, IWORK( 1 ) )
		 CALL ZHEEVD('V','Lower',w90basis,htb,w90basis,energias,WORK, LWORK,&
				RWORK,LRWORK,IWORK,LIWORK,INFO)      				
      		 IF( INFO.GT. 0 ) THEN
         	  WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         	  STOP
      		 END IF
      		
      		 
      		
      		end select

		do i=1,w90basis

			do j=1,w90basis

				autovetores(i,j) = htb(j,i)
			end do

		end do
		
		
		if (nocpf .gt. 0) then
		
			nocp = nocpf
		
		else

		do i=1,w90basis

			if (energias(i) .gt. 0.0 ) then

				nocp = i-1
				EXIT

				

			else

			continue

			end if

		end do
		
		end if

		do i=nocp+1,w90basis

			energias(i)=energias(i)+scs

		end do
		! $omp end single
	

end subroutine eigsys










