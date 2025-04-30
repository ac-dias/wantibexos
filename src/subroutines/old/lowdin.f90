subroutine eigsysl(nthread,dft,systype,scs,exc,nocp,ffactor,kx,ky,kz,w90basis,nvec,rlat,rvec,hopmatrices,&
		   ihopmatrices,ovp,efermi,energias,autovetores,nocpf,fermishift,mag,lowdin,power,ovptb,smhalf)

	use omp_lib
        !use f95_lapack

	implicit none

	integer :: i,j,k,nthread

	integer :: w90basis,nvec

	real :: kx,ky,kz

	real :: scs,exc
	
	real,dimension(3) :: mag

	integer :: nocp

	real,dimension(3,3) :: rlat

	real,dimension(nvec,3) :: rvec

	integer,dimension(nvec) :: ffactor


	complex,dimension(w90basis,w90basis) :: htb,fermimt,ovptb
	complex,dimension(w90basis,w90basis) :: hexc
	
	complex,dimension(w90basis,w90basis) :: htbl,smhalf
	

	real,dimension(nvec,w90basis,w90basis) :: hopmatrices,ihopmatrices
	real,dimension(nvec,w90basis,w90basis) :: ovp


	real :: efermi,fermishift

	real,dimension(w90basis) :: energias

	complex,dimension(w90basis,w90basis) :: autovetores

	complex,parameter :: imag=cmplx(0.0,1.0)
	
	integer :: nocpf
	
	character(len=1) :: dft
	character(len=4) :: systype
	
	integer :: power
	logical :: lowdin

	!definicoes diagonalizacao
	INTEGER   ::       ifail
	real,parameter :: ABSTOL=1.0e-6
	INTEGER ::         LIWORK, LRWORK
        INTEGER  ::        INFO, LWORK,LWMAX
        real,dimension(1 + 5*w90basis + 2*w90basis**2) :: RWORK
        complex,dimension(2*w90basis + w90basis**2) :: WORK
	INTEGER,dimension(3 + 5*w90basis) :: IWORK

	LWMAX = (2*w90basis-1)**2

	!call OMP_SET_NUM_THREADS(nthread)
	!call MKL_SET_NUM_THREADS(1)

		! $omp single
		select case (dft)
		
		case ("S")
		
		 
		 call hamiltonian(w90basis,nvec,rvec,hopmatrices,ihopmatrices,ffactor,kx,ky,kz,htb)
		 call overlap(w90basis,nvec,rvec,ovp,kx,ky,kz,ovptb)
		 
		 write(*,*) "Calculou H e S"
		 
		 	if (lowdin) then
		 	
		 		call shalf(w90basis,power,ovptb,smhalf)
		 		
		 		!smhalf = ovptb
		 		
		 		write(*,*) "Calculou S^-1/2"
		 		
		 		!call lowdinhamiltonian(w90basis,htb,smhalf,htbl)
		 		
		 		htbl = matmul(smhalf,matmul(htb,smhalf))
		 		
		 		write(*,*) "Calculou HL"
		 	
		 	else
		 	
		 	 continue
		 	
		 	end if
		
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
  		
  		
  		 if (lowdin) then 
  		 
  		 LWORK = -1
      		 LIWORK = -1
      		 LRWORK = -1
		 CALL CHEEVD('V','U',w90basis,htbl,w90basis,energias,WORK, LWORK,&
				RWORK,LRWORK,IWORK,LIWORK,INFO)
		 LWORK = MIN( 2*w90basis + w90basis**2, INT( WORK( 1 ) ) )
      		 LRWORK = MIN( 1 + 5*w90basis + 2*w90basis**2, INT( RWORK( 1 ) ) )
      		 LIWORK = MIN( 3 + 5*w90basis, IWORK( 1 ) )
		 CALL CHEEVD('V','U',w90basis,htbl,w90basis,energias,WORK, LWORK,&
				RWORK,LRWORK,IWORK,LIWORK,INFO)      				
      		 IF( INFO.GT. 0 ) THEN
         	  WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         	  STOP
      		 END IF
  		 
  		 write(*,*) 'Diagonalizou o HL'
  		 
  		 else
  		
  		 !LWORK = -1
      		 !CALL CHEGV( 1,'Vectors', 'U', w90basis, htb, w90basis,ovptb,w90basis, energias, WORK, LWORK, RWORK, INFO )
      		 !LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      		 !CALL CHEGV( 1,'Vectors', 'U', w90basis, htb, w90basis,ovptb,w90basis, energias, WORK, LWORK, RWORK, INFO )

      		 !IF( INFO.GT. 0 ) THEN
         	 ! WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         	 ! STOP
      		 !END IF
  		
  		 LWORK = -1
      		 LIWORK = -1
      		 LRWORK = -1
  		 CALL CHEGVD(1,'V','U',w90basis,htb,w90basis,ovptb,w90basis,energias,WORK, LWORK,&
		   		RWORK,LRWORK,IWORK,LIWORK,INFO)
		 LWORK = MIN( 2*w90basis + w90basis**2, INT( WORK( 1 ) ) )
      		 LRWORK = MIN( 1 + 5*w90basis + 2*w90basis**2, INT( RWORK( 1 ) ) )
      		 LIWORK = MIN( 3 + 5*w90basis, IWORK( 1 ) )
		 CALL CHEGVD(1,'V','U',w90basis,htb,w90basis,ovptb,w90basis,energias,WORK, LWORK,&
		 		RWORK,LRWORK,IWORK,LIWORK,INFO) 
				
		 IF( INFO.GT. 0 ) THEN
         	  WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         	  STOP
      		 END IF
      		 
      		 end if		
  		
  		case default  
  		
  		  		
      		
      		 LWORK = -1
      		 LIWORK = -1
      		 LRWORK = -1
		 CALL CHEEVD('V','U',w90basis,htb,w90basis,energias,WORK, LWORK,&
				RWORK,LRWORK,IWORK,LIWORK,INFO)
		 LWORK = MIN( 2*w90basis + w90basis**2, INT( WORK( 1 ) ) )
      		 LRWORK = MIN( 1 + 5*w90basis + 2*w90basis**2, INT( RWORK( 1 ) ) )
      		 LIWORK = MIN( 3 + 5*w90basis, IWORK( 1 ) )
		 CALL CHEEVD('V','U',w90basis,htb,w90basis,energias,WORK, LWORK,&
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
		
		!do i = 1,w90basis
		
		 !call normalizewf(w90basis,autovetores(i,:))
		
		!end do
		
		
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
	

end subroutine eigsysl



subroutine basischange(ndim,smh,wf,wfno)

	implicit none
	
	integer :: ndim
	complex,dimension(ndim,ndim) :: smh
	complex,dimension(ndim) :: wf,wfno
	
	call matvec(smh,wf,ndim,wfno)


end subroutine basischange

subroutine lowdinhamiltonian(ndim,htb,smh,htbl)

	implicit none
	integer :: ndim,i,j,k
	complex,dimension(ndim,ndim) :: htb,smh,htblaux,htbl
	
	
	!multiplicando htb por smh
	do i=1,ndim
	 do j=1,ndim
	  
	  htblaux(i,j) = cmplx(0.0,0.0)
	   
	   do k=1,ndim
	   
	    htblaux(i,j) = htblaux(i,j) + htb(i,k)*smh(k,j)
	    
	   end do
	 
	 end do
	end do
	
	!multiplicando smh por htblaux
	do i=1,ndim
	 do j=1,ndim
	  
	  htbl(i,j) = cmplx(0.0,0.0)
	   
	   do k=1,ndim
	   
	    htbl(i,j) = htbl(i,j) + smh(i,k)*htblaux(k,j)
	    
	   end do
	 
	 end do
	end do	

end subroutine lowdinhamiltonian


subroutine shalf(ndim,power,ovp,smh) !calcula a matriz de overlap elevada a -1/2

	implicit none

	integer :: i,j,k,ndim,power
	complex,dimension(ndim,ndim) :: ovp, xmatrix, smh, idtm
	complex,dimension(power,ndim,ndim) :: xpmatrix
	real,dimension(power) :: sumindex
	
	call idtmatrix(ndim,idtm)
	
	xmatrix = ovp - idtm
	
	call binomialterms(-0.5,power,sumindex)
	
	call matrixpower(ndim,power,xmatrix,xpmatrix)
	
	smh = cmplx(0.0,0.0)
	
	do i=1,power
	
		do j=1,power
		 do k=1,power
	
		smh(j,k) = smh(j,k) + sumindex(i)*xpmatrix(i,j,k)
		 
	         end do
		end do
	
	end do
	
		smh = smh + idtm

end subroutine shalf





subroutine binomialterms(k,nterms,sumindex) !calcula os n primeiros termos (nterms) da série binomial (1+x)^{k}

	implicit none
	
	integer :: i
	real :: k
	integer :: nterms
	real,dimension(nterms) :: sumindex
	
	!sumindex(0) = 1.0
	sumindex(1) = k
	
	do i=2,nterms
	
		sumindex(i) = (sumindex(i-1)*(k-real(i)+1.0))/real(i)
	
	end do
	

end subroutine binomialterms


subroutine matrixpower(ndim,power,matrix,matrixp) !calcula uma matriz elevada da potencia 0 até a potencia "power", guarda todas essas matrizes

	implicit none
	
	integer :: ndim,power,i
	complex,dimension(ndim,ndim) :: matrix, idtm
	complex,dimension(power,ndim,ndim) :: matrixp
	
	!matrixp(0,:,:)= cmplx(0.0,0.0)
	
	!do i=1,ndim
	
	!	matrixp(0,i,i) = cmplx(1.0,0.0)
	
	!end do
	
	
	matrixp(1,:,:) = matrix
	
	do i=2,power
	
		matrixp(i,:,:)= matmul(matrix,matrixp(i-1,:,:))
	
	end do	

	


end subroutine matrixpower





subroutine idtmatrix(n,idtm) !construir matrix identidade

	implicit none
	
	integer :: i,n
	complex,dimension(n,n) :: idtm
	
	idtm= cmplx(0.0,0.0)
	
	do i=1,n
	
		idtm(i,i) = cmplx(1.0,0.0)
	
	end do

end subroutine
