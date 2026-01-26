subroutine gw_pi0_calc(nthreads,outputfolder,ngrid,smegw,params,&
                       exc,mshift,nocpf,fermishift,dft,mag,&
                       nomega,omegamax,sysdim)

	use omp_lib
	use hamiltonian_input_variables

	implicit none

	!subroutine input variables
	real,parameter:: pi=acos(-1.)
	integer :: i,j,k,l,h,m,n,erro
	complex,parameter :: imag=cmplx(0.0,1.0)
	
	integer :: nthreads
	character(len=70) :: outputfolder    !pasta saida
	integer,dimension(3) :: ngrid
	real :: smegw
	character(len=70) :: params   !parametros TB
	real :: exc
	real,dimension(3) :: mshift	
	integer :: nocpf
	real :: fermishift
	character(len=1) :: dft	
	real,dimension(3) :: mag
	integer :: nomega
	real :: omegamax
	
	character(len=5) :: sysdim	
	
	!variaveis relacionadas a marcacao do tempo
	real:: t0,tf
	integer,dimension(8) :: values,values2
	
	!internal variables
	real,allocatable,dimension(:) :: omega
		
	integer :: ngqpt
	real,allocatable,dimension(:,:) :: qpt
	
	real,allocatable,dimension(:,:,:) :: qptpq
	
	real,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores		
	
	
	real,allocatable,dimension(:,:) :: qeigv
	complex,allocatable,dimension(:,:,:) :: qvector
	integer,allocatable,dimension(:) :: nocpq
	complex,allocatable,dimension(:,:,:) :: qovp	
	
	
	integer :: nocpaux	
	
	real,allocatable,dimension(:,:,:) :: qpqeigv
	complex,allocatable,dimension(:,:,:,:) :: qpqvector
	complex,allocatable,dimension(:,:,:,:) :: qpqovp	
	
	real,allocatable,dimension(:,:,:,:) :: mmnqpq	
	
	complex,allocatable,dimension(:,:) :: pi0	
	
	real :: gs,wk,waux
	
	real :: aux1,enocp
	complex :: aux2
	

	
	
	!OUTOUT files	
	OPEN(UNIT=300, FILE=trim(outputfolder)//"log_gw_pi0_mesh.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_gw_pi0_mesh output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"gw_pi0_mesh.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_pi0_mesh output file"
	OPEN(UNIT=302, FILE=trim(outputfolder)//"gw_pi0_gamma.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_pi0_gamma output file"    	
    	
	call cpu_time(t0)
	call date_and_time(VALUES=values)
	call OMP_SET_NUM_THREADS(nthreads)    		
	
	
	!Hamiltonian read
	
	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	allocate(ovp(nvec,w90basis,w90basis))

	call hamiltonian_input_read(200,params)
	

	end select	
	
	write(300,*)
	write(300,*)
	write(300,*)
	write(300,*) 'threads:', nthreads
	write(300,*)
	write(300,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
	write(300,*)
	write(300,"(A13,3F15.4)") 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	write(300,*)
	write(300,*) 'NOMEGA:',nomega
	write(300,*) 'OMEGAMAX:',omegamax
	write(300,*) 'SIGMA_GW:',smegw		
	
	
	write(300,*)
	write(300,*) 'begin','  ','month',values(2),'day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(300,*)	
	call flush(300)	
	
	allocate(omega(nomega))	
	
	!define omega values
	
	do i=1,nomega
	
		omega(i) = omegamax*(real(i-1.0)/real(nomega-1.0))
		!write(*,*) omega(i) !for debug only	
	end do
	
	
	!calculate kmesh (q)
	
	ngqpt = ngrid(1)*ngrid(2)*ngrid(3)
	
	allocate(qpt(ngqpt,3))
	
	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),qpt)	
	
	
	!calculate eigenvalues and eigenvectors in kmesh (q)
	
	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	allocate(qeigv(ngqpt,w90basis),qvector(ngqpt,w90basis,w90basis))
	allocate(nocpq(ngqpt))
	
	if (dft .eq. "S") then
		 allocate(qovp(ngqpt,w90basis,w90basis))
	else
	 continue
	end if		
	
	!$omp parallel do default(shared) private(i,j,l,h,eaux,vaux)
	do i=1,ngqpt

#ifdef MKL
		call MKL_SET_NUM_THREADS(1)
#endif		


#ifdef AOCL		
		call bli_thread_set_num_threads(1)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(1)
		
#endif
	
		call eigsys(nthreads,dft,systype,scs,exc,nocpq(i),ffactor,qpt(i,1),qpt(i,2),qpt(i,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,ovp,efermi,eaux,vaux,nocpf,fermishift,mag)
#ifdef MKL
		call MKL_SET_NUM_THREADS(nthreads)
#endif		

#ifdef AOCL		
		call bli_thread_set_num_threads(nthreads)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(nthreads)
		
#endif	

			
		do j=1,w90basis
				qeigv(i,j)= eaux(j)
				!ebands(j,i) = qeigv(i,j)

		end do	
		
		do l=1,w90basis


			do h=1,w90basis

				qvector(i,l,h)=vaux(l,h)


			end do
			

		end do		
			

		if (dft .eq. "S") then
		 call overlap(w90basis,nvec,rvec,ovp,qpt(i,1),qpt(i,2),qpt(i,3),qovp(j,:,:))
		else
		 continue
		end if
	
	end do
	!$omp end parallel do	
	
	
	write(300,*)
	write(300,*)"Eigenvalues and Eigenvectors in q-mesh calculated"	
	write(300,*)	
	
	
	
	allocate(qptpq(ngqpt,ngqpt,3))
	
	do i=1,ngqpt
	 do j=1,ngqpt
	 
	 	qptpq(i,j,:) = qpt(i,:) + qpt(j,:)
	 
	 end do
	end do			
	
	
	!calculate eigenvalues and eigenvectors in kmesh (q+q')
	
	allocate(qpqeigv(ngqpt,ngqpt,w90basis),qpqvector(ngqpt,ngqpt,w90basis,w90basis))
	
	if (dft .eq. "S") then
		 allocate(qpqovp(ngqpt,ngqpt,w90basis,w90basis))
	else
	 continue
	end if	
	
	!$omp parallel do default(shared) private(i,j,k,l,h,eaux,vaux)
	do i=1,ngqpt
	 do j=1,ngqpt		

#ifdef MKL
		call MKL_SET_NUM_THREADS(1)
#endif		


#ifdef AOCL		
		call bli_thread_set_num_threads(1)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(1)
		
#endif

		call eigsys(nthreads,dft,systype,scs,exc,nocpaux,ffactor,qptpq(i,j,1),qptpq(i,j,2),qptpq(i,j,3),w90basis,nvec,&
			    rlat,rvec,hopmatrices,&
		             ihopmatrices,ovp,efermi,eaux,vaux,nocpf,fermishift,mag)

#ifdef MKL
		call MKL_SET_NUM_THREADS(nthreads)
#endif		

#ifdef AOCL		
		call bli_thread_set_num_threads(nthreads)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(nthreads)
		
#endif	

		if (dft .eq. "S") then
		 call overlap(w90basis,nvec,rvec,ovp,qptpq(i,j,1),qptpq(i,j,2),qptpq(i,j,3),qpqovp(i,j,:,:))
		else
		 continue
		end if
		
		
		do k=1,w90basis
				qpqeigv(i,j,k)= eaux(k)


		end do	
		
		do l=1,w90basis


			do h=1,w90basis

				qpqvector(i,j,l,h)=vaux(l,h)


			end do
			

		end do			

	 end do
	end do
	!$omp end parallel do
		
	write(300,*)
	write(300,*)"Eigenvalues and Eigenvectors in q+q'-mesh calculated"	
	write(300,*)	
	call flush(300)	
	
	!calculate Mmn(q,q+q')
	
	allocate(mmnqpq(w90basis,w90basis,ngqpt,ngqpt))
	
	!$omp parallel do collapse(4) default(shared) private(i,j,k,l)
	do i=1,ngqpt
	 do j=1,w90basis
	  do k=1,ngqpt
	   do l=1,w90basis
	
		call ovpsqr(dft,w90basis,qvector(i,j,:),qovp(i,:,:),qpqvector(i,k,l,:),qpqovp(i,k,:,:),mmnqpq(j,l,i,k))
	
	   end do
	  end do
	 end do  
	end do	
	!$omp end parallel do	
	
	
	deallocate(qpqvector,qvector)
		
	if (dft .eq. "S") then
		 deallocate(qpqovp,qovp)
	else
	 continue
	end if	
	
	write(300,*)
	write(300,*)"Mmn(q,q+q') calculated"	
	write(300,*)	
	call flush(300)	
	
	deallocate(eaux,vaux)
	
	!calculate Pi0(q,w)	
	
	
	allocate(pi0(ngqpt,nomega))	
	
	!$omp parallel do default(shared) private(i)
	do i=1,ngqpt
	 do j=1,nomega
	
		call polarization(fermishift,sysdim,ngrid,ngqpt,rlat,omega(j),w90basis,systype,&
			                  smegw,qeigv,qpqeigv(i,:,:),mmnqpq(:,:,i,:),pi0(i,j))
	
	
	 	!pol(i,:) = 0.0
	 end do	
	end do
	!$omp end parallel do
	
!	select case (systype)
	
!	case("NP")
	
!		gs = 2.0
	
!	case default
	
!		gs = 1.0	
!	end select
	
!	waux = wk(sysdim,ngrid,rlat)	
	
!	pi0 = 0.0
	
	! $omp parallel do default(shared) private(h)
!	do h=1,ngqpt
!	 do m=1,nomega
	 	
!	  do i=1,ngqpt
!	    do j=1,w90basis
!	      do k=1,w90basis
	      
	      
	    
!	    	aux1 = enocp(fermishift,qeigv(i,j)) - enocp(fermishift,qpqeigv(h,i,k))
	    	
	    	
	    	
!	    	aux2 = omega(m) - (qpqeigv(h,i,k)-qeigv(i,j)) + smegw*cmplx(0.0,1.0)
	    	
!	    	pi0(h,m) = pi0(h,m) + ((aux1/aux2)*mmnqpq(j,k,h,i))
	    
!	      end do
!	    end do
!	  end do		
	
	
!		pi0(h,m) = pi0(h,m)*gs*waux
	
!	 end do
!	end do
	! $omp end parallel do
	
	write(301,*) ngrid(1),ngrid(2),ngrid(3)
	write(301,*) 
	
	do i=1,ngqpt
	
		write(301,*) i
		write(301,*)		
	
	 do j=1,nomega
	
		write(301,*) omega(j),real(pi0(i,j)),aimag(pi0(i,j))
	
	 end do
	 
	 	write(301,*)
	end do
	

	do i=1,nomega
	
		write(302,*) omega(i),real(pi0(1,i)),aimag(pi0(1,i))
	
	end do
	
	deallocate(qpqeigv,qeigv)
	deallocate(pi0)
	deallocate(mmnqpq)
	
	write(300,*)
	write(300,*)"Pi0(q,w)  calculated"	
	write(300,*)	
	call flush(300)	
	
		
	call cpu_time(tf)
	call date_and_time(VALUES=values2)

	write(300,*)
	write(300,*) 'end','   ','month',values2(2),'day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(300,*)
	call flush(300)
		
	
	close(300)
	close(301)
	close(302)		
	
	

end subroutine gw_pi0_calc
