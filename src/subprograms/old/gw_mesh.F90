subroutine gwmeshcalc(nthreads,outputfolder,calcparams,ngrid,smegw,ktolgw,params,edielgw,&
                       exc,mshift,coultypegw,ezgw,wgw,r0gw,lcgw,rk,meshtype,nocpf,fermishift,dft,mag,&
                       nomega,omegamax,sysdim,selfxonly)

	use omp_lib
	use hamiltonian_input_variables

	implicit none
	
	!subroutine input variables
	real,parameter:: pi=acos(-1.)
	integer :: i,j,k,l,h,m,n,erro
	complex,parameter :: imag=cmplx(0.0,1.0)
	
	integer :: nthreads
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: calcparams
	integer,dimension(3) :: ngrid
	real :: smegw
	real :: ktolgw	
	character(len=70) :: params   !parametros TB
	real,dimension(3) :: edielgw
	real :: exc
	real,dimension(3) :: mshift	
	character(len=5) :: coultypegw
	real :: ezgw,wgw,r0gw,lcgw
	real :: rk
	character(len=70) :: meshtype
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
	real,allocatable,dimension(:,:,:) :: qptmq,qptpq
	
	real,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores	
		
	real,allocatable,dimension(:,:) :: qeigv,ebands
	complex,allocatable,dimension(:,:,:) :: qvector
	integer,allocatable,dimension(:) :: nocpq
	complex,allocatable,dimension(:,:,:) :: qovp
	
	integer :: nkc,nkv,nkgap
	real :: cbm,vbm,gap
	
	integer :: nocpaux
	real,allocatable,dimension(:,:,:) :: qpqeigv,qmqeigv
	complex,allocatable,dimension(:,:,:,:) :: qpqvector,qmqvector
	complex,allocatable,dimension(:,:,:,:) :: qpqovp,qmqovp
	
	real,allocatable,dimension(:,:,:,:) :: mmnqq,mmnqpq,mmnqmq
	
	complex,allocatable,dimension(:,:) :: pol,w0,selfc
	
	real,allocatable,dimension(:,:) :: vqa	
	
	real,allocatable,dimension(:,:) :: selfx,znk,eigcor
	
	complex,dimension(nomega) :: polaux
	
	logical :: selfxonly			
	
	!OUTOUT files	
	OPEN(UNIT=300, FILE=trim(outputfolder)//"log_gw_mesh.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_gw_mesh output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"gw_qp_energy_cor_mesh.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_qp_energy_cor_mesh output file"
	OPEN(UNIT=302, FILE=trim(outputfolder)//"gw_pi0_gamma.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_qp_gamma_pol output file"    	
	OPEN(UNIT=303, FILE=trim(outputfolder)//"gw_pi0_mesh.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_pi0_mesh output file"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"gw_w0_mesh.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_w0_mesh output file"
    	
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
	write(300,"(A14,1E15.4)") 'ktol-coulomb:', ktolgw
	write(300,*)
	write(300,"(A13,3F15.4)") 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	write(300,*)
	write(300,*) 'Coulomb Potential:',coultypegw
	write(300,*)
	write(300,*) 'NOMEGA:',nomega
	write(300,*) 'OMEGAMAX:',omegamax
	write(300,*) 'SIGMA_GW:',smegw
	write(300,*) 'SELF_X_ONLY:',selfxonly			
	
	
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
	allocate(ebands(w90basis,ngqpt))
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
				ebands(j,i) = qeigv(i,j)

		end do	
		
		do l=1,w90basis


			do h=1,w90basis

				qvector(i,l,h)=vaux(l,h)


			end do
			

		end do		
			

		if (dft .eq. "S") then
		 call overlap(w90basis,nvec,rvec,ovp,qpt(i,1),qpt(i,2),qpt(i,3),qovp(i,:,:))
		else
		 continue
		end if
	
	end do
	!$omp end parallel do	
	
	
	write(300,*)
	write(300,*)"Eigenvalues and Eigenvectors in q-mesh calculated"	
	write(300,*)
	
	call gapfinder(w90basis,ngqpt,nocpq,ebands,nkc,nkv,nkgap,gap,cbm,vbm)
	
	write(300,*) "TB fundamental band gap (eV):",cbm-vbm
	write(300,*) "TB direct band gap (eV):",gap
	write(300,*)
	write(300,*)
	write(300,*)"TB kpoint cbm"
	write(300,"(3E18.8)")qpt(nkc,1),qpt(nkc,2),qpt(nkc,3)
	write(300,*)	
	write(300,*)"TB kpoint vbm"
	write(300,"(3E18.8)")qpt(nkv,1),qpt(nkv,2),qpt(nkv,3)
	write(300,*)		 
	write(300,*)"TB kpoint direct band gap"
	write(300,"(3E18.8)")qpt(nkgap,1),qpt(nkgap,2),qpt(nkgap,3)
	call flush(300)
	
	!calculate Mmn(q,q')
	
	allocate(mmnqq(w90basis,w90basis,ngqpt,ngqpt))
	
	!$omp parallel do collapse(4) default(shared) private(i,j,k,l)
	do i=1,ngqpt
	 do j=1,w90basis
	  do k=1,ngqpt
	   do l=1,w90basis
	
		call ovpsqr(dft,w90basis,qvector(i,j,:),qovp(i,:,:),qvector(k,l,:),qovp(k,:,:),mmnqq(j,l,i,k))
	
	   end do
	  end do
	 end do  
	end do	
	!$omp end parallel do
		
	write(300,*)
	write(300,*)"Mmn(q,q') calculated"	
	write(300,*)	
	call flush(300)
	
	!calculate self-energy exchange part
	
	allocate(selfx(w90basis,ngqpt))
	
	!$omp parallel do default(shared) private(i,j)
	do i=1,ngqpt
	 do j=1,w90basis
	
	  call self_en_x(w90basis,sysdim,rlat,j,qpt(i,:),ngrid,ngqpt,qpt,nocpq,coultypegw,edielgw,&
	               lcgw,ezgw,wgw,r0gw,ktolgw,mmnqq(:,:,i,:),selfx(j,i))
	
	 end do
	end do
	!$omp end parallel do
	
	deallocate(mmnqq)
	
	write(300,*)
	write(300,*)"self_x(n,q) calculated"	
	write(300,*)	
	call flush(300)		
	
	if (selfxonly) then
	
	 allocate(selfc(w90basis,ngqpt),znk(w90basis,ngqpt))
	 
	 selfc = 0.0
	 
	 znk = 1.0
	 
	 allocate(pol(ngqpt,nomega))
	 
	 allocate(w0(ngqpt,nomega),vqa(ngqpt,nomega))
	 
	 go to 991
	
	else 
	
		continue
	
	end if
	
		
	!calculate kmesh (q+q')
	
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
	
	deallocate(qpqvector)
		
	if (dft .eq. "S") then
		 deallocate(qpqovp)
	else
	 continue
	end if	
	
	write(300,*)
	write(300,*)"Mmn(q,q+q') calculated"	
	write(300,*)	
	call flush(300)
	

	
	!calculate kmesh (q-q')
	
	allocate(qptmq(ngqpt,ngqpt,3))
	
	do i=1,ngqpt
	 do j=1,ngqpt
	 
	 	qptmq(i,j,:) = qpt(i,:) - qpt(j,:)
	 
	 end do
	end do		
	
	!calculate eigenvalues and eigenvectors in kmesh (q-q')
	
	allocate(qmqeigv(ngqpt,ngqpt,w90basis),qmqvector(ngqpt,ngqpt,w90basis,w90basis))
	
	if (dft .eq. "S") then
		 allocate(qmqovp(ngqpt,ngqpt,w90basis,w90basis))
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

		call eigsys(nthreads,dft,systype,scs,exc,nocpaux,ffactor,qptmq(i,j,1),qptmq(i,j,2),qptmq(i,j,3),w90basis,nvec,&
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
		 call overlap(w90basis,nvec,rvec,ovp,qptmq(i,j,1),qptmq(i,j,2),qptmq(i,j,3),qmqovp(i,j,:,:))
		else
		 continue
		end if
		
		
		do k=1,w90basis
				qmqeigv(i,j,k)= eaux(k)


		end do	
		
		do l=1,w90basis


			do h=1,w90basis

				qmqvector(i,j,l,h)=vaux(l,h)


			end do
			

		end do			

	 end do
	end do
	!$omp end parallel do
		
	write(300,*)
	write(300,*)"Eigenvalues and Eigenvectors in q-q'-mesh calculated"	
	write(300,*)
	call flush(300)
	
	deallocate(eaux,vaux)			

	!calculate Mmn(q,q-q')	
	
	allocate(mmnqmq(w90basis,w90basis,ngqpt,ngqpt))
	
	!$omp parallel do collapse(4) default(shared) private(i,j,k,l)
	do i=1,ngqpt
	 do j=1,w90basis
	  do k=1,ngqpt
	   do l=1,w90basis
	
		call ovpsqr(dft,w90basis,qvector(i,j,:),qovp(i,:,:),qmqvector(i,k,l,:),qmqovp(i,k,:,:),mmnqmq(j,l,i,k))
	
	   end do
	  end do
	 end do  
	end do	
	!$omp end parallel do
	
	deallocate(qmqvector)
		
	if (dft .eq. "S") then
		 deallocate(qpqovp)
	else
	 continue
	end if	
	
	write(300,*)
	write(300,*)"Mmn(q,q-q') calculated"	
	write(300,*)		
	call flush(300)
	
	!calculate Pi0(q,w) and W0(q,w)
	
	allocate(pol(ngqpt,nomega))
	
	
	
	!$omp parallel do default(shared) private(i)
	do i=1,ngqpt
	 do j=1,nomega
	
		call polarization(fermishift,sysdim,ngrid,ngqpt,rlat,omega(j),w90basis,systype,&
		                   smegw,qeigv,qpqeigv(i,:,:),mmnqpq(:,:,i,:),pol(i,j))
	
	
	 	!pol(i,:) = 0.0
	 end do	
	end do
	!$omp end parallel do
	
	!pol = 0.0
	
	do i=1,nomega
	
		write(302,*) omega(i),real(pol(1,i)),aimag(pol(1,i))
	
	end do
	
	deallocate(qpqeigv)
	
	write(303,*) ngrid(1),ngrid(2),ngrid(3)
	write(303,*) 
	
	do i=1,ngqpt
	
		write(303,*) i
		write(303,*)		
	
	 do j=1,nomega
	
		write(303,*) omega(j),real(pol(i,j)),aimag(pol(i,j))
	
	 end do
	 
	 	write(303,*)
	end do	
	
	write(300,*)
	write(300,*)"Pi0(q,w)  calculated"	
	write(300,*)	
	call flush(300)
	
	allocate(w0(ngqpt,nomega),vqa(ngqpt,nomega))	
	
	!$omp parallel do default(shared) private(i)
	do i=1,ngqpt
	
		call w0coul(qpt(i,:),nomega,rlat,ngrid,ngqpt,coultypegw,edielgw,lcgw,ezgw,wgw,r0gw,ktolgw,pol(i,:),vqa(i,:),w0(i,:))
	

	
	end do
	!$omp end parallel do
	
	write(304,*) ngrid(1),ngrid(2),ngrid(3)
	write(304,*) 
	
	do i=1,ngqpt
	
		write(304,*) i
		write(304,*)		
	
	 do j=1,nomega
	
		write(304,*) omega(j),vqa(i,j),real(w0(i,j)),aimag(w0(i,j))
	
	 end do
	 
	 	write(304,*)
	end do		
	
	write(300,*)
	write(300,*)"W0(q,w)  calculated"	
	write(300,*)	
	call flush(300)
	
	deallocate(mmnqpq)
	
	deallocate(qptpq)	
	

	
	
	!calculate self-energy correlation part
	
	allocate(selfc(w90basis,ngqpt),znk(w90basis,ngqpt))
	
	!$omp parallel do default(shared) private(i,j)
	do i=1,ngqpt
	 do j=1,w90basis
	


	  call self_c_znk(w90basis,sysdim,rlat,ngrid,ngqpt,fermishift,j,&
	                  qeigv(i,j),qmqeigv(i,:,:),qpt(i,:),omega,nomega,w0,mmnqmq(:,:,i,:),smegw,selfc(j,i),znk(j,i))
	
	 end do
	end do	
	!$omp end parallel do
	
	deallocate(mmnqmq,qmqeigv)
	
	write(300,*)
	write(300,*)"self_c(n,q) calculated"	
	write(300,*)	
	call flush(300)
	
	!eigenvalues correction

991     continue
	
	allocate(eigcor(w90basis,ngqpt))
	
	do i=1,ngqpt
	 do j=1,w90basis
	
		eigcor(j,i) = znk(j,i)*(selfx(j,i)+real(selfc(j,i)))
	
	 end do
	end do
	
	write(301,*) w90basis
	write(301,*) ngrid(1),ngrid(2),ngrid(3)
	write(301,*)
	
	do i=1,ngqpt
	
		write(301,*) '#qpt n:',i,qpt(i,1),qpt(i,2),qpt(i,3)
		write(301,*)
		
	 do j=1,w90basis
	 
	 	write(301,*) eigcor(j,i),selfx(j,i),real(selfc(j,i)),aimag(selfc(j,i)),znk(j,i)
	 
	 end do
	
		write(301,*)
	
	end do
	
	deallocate(znk,selfx,selfc)
	
	
	do i=1,ngqpt
	 do j=1,w90basis
	
		ebands(j,i) = ebands(j,i)+eigcor(j,i)
	
	 end do
	end do	
	
	!estimating GW gap
	

	write(300,*)
	write(300,*)"QP Eigenvalues and Eigenvectors in q-mesh calculated"	
	write(300,*)
	
	call gapfinder(w90basis,ngqpt,nocpq,ebands,nkc,nkv,nkgap,gap,cbm,vbm)
	
	write(300,*) "G0W0 QP fundamental band gap (eV):",cbm-vbm
	write(300,*) "G0W0 QP direct band gap (eV):",gap
	write(300,*)
	write(300,*)
	write(300,*)"G0W0 QP kpoint cbm"
	write(300,"(3E18.8)")qpt(nkc,1),qpt(nkc,2),qpt(nkc,3)
	write(300,*)	
	write(300,*)"G0W0 QP kpoint vbm"
	write(300,"(3E18.8)")qpt(nkv,1),qpt(nkv,2),qpt(nkv,3)
	write(300,*)		 
	write(300,*)"G0W0 QP kpoint direct band gap"
	write(300,"(3E18.8)")qpt(nkgap,1),qpt(nkgap,2),qpt(nkgap,3)	
	call flush(300)
	
	deallocate(omega)
	deallocate(ebands)
	deallocate(eigcor)
	deallocate(nocpq)
	deallocate(pol,w0,vqa)

	call cpu_time(tf)
	call date_and_time(VALUES=values2)

	write(300,*)
	write(300,*) 'end','   ','month',values2(2),'day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(300,*)
	call flush(300)
		
	
	close(300)
	close(301)
	close(302)
	close(303)
	close(304)		
		
end subroutine gwmeshcalc
