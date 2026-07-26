subroutine sgw_mesh(nthreads,outputfolder,ngrid,ifactor,smegw,ktolgw,params,edielgw,&
                       exc,mshift,coultypegw,ezgw,wgw,r0gw,lcgw,rk,meshtype,nocpf,fermishift,dft,mag,&
                       nomega,omegamax,sysdim,selfxonly,rnmgw)
                       
	use omp_lib
	use hamiltonian_input_variables

	implicit none                       
	
	!subroutine input variables
	real,parameter:: pi=acos(-1.)
	integer :: i,j,k,l,h,m,n,erro,erro1
	complex,parameter :: imag=cmplx(0.0,1.0)
	
	integer :: nthreads
	character(len=70) :: outputfolder    !pasta saida
	integer,dimension(3) :: ngrid,ngridint
	real :: ifactor
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
	character(len=70) :: kpaths    !kpath	input file
	integer :: nomega
	real :: omegamax
	
	character(len=5) :: sysdim
		
	
	!variaveis relacionadas a marcacao do tempo
	real:: t0,tf
	integer,dimension(8) :: values,values2
	
	!internal variables	
	
	real,allocatable,dimension(:) :: omega	
	
	integer :: ngqpt,nkpt	
	
	real,allocatable,dimension(:,:) :: qpt,kpt
	real,allocatable,dimension(:,:,:) :: kptmq,qptpq	
	
	real,allocatable,dimension(:) :: eaux !variavel auxiliar para energia
	complex,allocatable,dimension(:,:) :: vaux !variavel auxiliar para os autovetores	
	
	real,allocatable,dimension(:,:) :: keigv,qeigv,ebands
	complex,allocatable,dimension(:,:,:) :: kvector,qvector
	integer,allocatable,dimension(:) :: nocpq,nocpk
	complex,allocatable,dimension(:,:,:) :: qovp,kovp
	
	integer :: nkc,nkv,nkgap
	real :: cbm,vbm,gap		
	
	integer :: nocpaux	
		
	real,allocatable,dimension(:,:,:) :: qpqeigv,kmqeigv
	complex,allocatable,dimension(:,:,:,:) :: qpqvector,kmqvector
	complex,allocatable,dimension(:,:,:,:) :: qpqovp,kmqovp	
	
	real,allocatable,dimension(:,:,:,:) :: mmnqq,mmnkq	
	
	complex,allocatable,dimension(:,:) :: pol,w0,selfc
	
	real,allocatable,dimension(:,:) :: vqa	
	
	real,allocatable,dimension(:,:) :: selfx,znk,eigcor,eigcoravg	
	
	real,allocatable,dimension(:) :: selfxavg,znkavg
	
	complex,allocatable,dimension(:) :: selfcavg	
	
	logical :: selfxonly	
	
	real,allocatable,dimension(:,:) :: encor,encoravg
	
	real :: aux,aux1,aux2
	
	real :: rnmgw
	
	
	!OUTPUT files
	OPEN(UNIT=300, FILE=trim(outputfolder)//"log_gw_mesh.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening log_gw_mesh output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//"gw_qp_energy_cor_mesh.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_qp_energy_cor_mesh output file"
	OPEN(UNIT=302, FILE=trim(outputfolder)//"gw_pi0_gamma.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_qp_gamma_pol output file" 		
	OPEN(UNIT=303, FILE=trim(outputfolder)//"gw_pi0_w0_mesh.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_pi0_w0_mesh output file"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"gw_qp_energy_cor_avg.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_qp_energy_cor_avg output file"
	OPEN(UNIT=305, FILE=trim(outputfolder)//"gw_qp_energy_cor_mesh_renorm.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_qp_energy_cor_mesh_renorm output file"    	

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
	!write(300,*) 'integration grid:',ngridint(1),ngridint(2),ngridint(3)		
	write(300,*)
	!write(300,"(A14,1E15.4)") 'ktol-coulomb:', ktolgw
	write(300,*)
	write(300,"(A13,3F15.4)") 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	write(300,*)
	write(300,*) 'Coulomb Potential:',coultypegw
	write(300,*)
	write(300,*) 'NOMEGA:',nomega
	write(300,*) 'OMEGAMAX:',omegamax
	write(300,*) 'SIGMA_GW:',smegw
	write(300,*) 'SELF_X_ONLY:',selfxonly
	write(300,*) 'KTOL_GW:', ktolgw
	write(300,*) 'RNM_GW:', rnmgw	
	
	write(300,*)
	write(300,*) 'begin','  ','month',values(2),'day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(300,*)	
	call flush(300)	

	!calculate kmesh (k)
	
	nkpt = ngrid(1)*ngrid(2)*ngrid(3)
	
	allocate(kpt(nkpt,3))

	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpt)

	write(300,*)	"k-mesh points calculated"
	call flush(300)


	
	!define omega values
	
	allocate(omega(nomega))	
	
	
	do i=1,nomega
	
		omega(i) = omegamax*(real(i-1.0)/real(nomega-1.0))
		!write(*,*) omega(i) !for debug only	
	end do
	
	write(300,*)	"omega frequencies calculated"
	call flush(300)	
	
	allocate(eaux(w90basis),vaux(w90basis,w90basis))
	
	!calculate eigenvalues and eigenvectors in kpath (k)	
	
	allocate(ebands(w90basis,nkpt))
	allocate(keigv(nkpt,w90basis),kvector(nkpt,w90basis,w90basis))
	allocate(nocpk(nkpt))
	
	if (dft .eq. "S") then
		 allocate(kovp(nkpt,w90basis,w90basis))
	else
	 continue
	end if	
	
	!$omp parallel do default(shared) private(i,j,l,h,eaux,vaux)
	do i=1,nkpt

#ifdef MKL
		call MKL_SET_NUM_THREADS(1)
#endif		


#ifdef AOCL		
		call bli_thread_set_num_threads(1)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(1)
		
#endif
	
		call eigsys(nthreads,dft,systype,scs,exc,nocpk(i),ffactor,kpt(i,1),kpt(i,2),kpt(i,3),w90basis,nvec,&
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
				keigv(i,j)= eaux(j)
				ebands(j,i) = keigv(i,j)

		end do	
		
		do l=1,w90basis


			do h=1,w90basis

				kvector(i,l,h)=vaux(l,h)


			end do
			

		end do		
			

		if (dft .eq. "S") then
		 call overlap(w90basis,nvec,rvec,ovp,kpt(i,1),kpt(i,2),kpt(i,3),kovp(i,:,:))
		else
		 continue
		end if
	
	end do
	!$omp end parallel do	


	
	write(300,*)
	write(300,*)"Eigenvalues and Eigenvectors in k-mesh calculated"	
	write(300,*)
	
	call gapfinder(w90basis,nkpt,nocpk,ebands,nkc,nkv,nkgap,gap,cbm,vbm)
	
	write(300,*) "TB fundamental band gap (eV):",cbm-vbm
	write(300,*) "TB direct band gap (eV):",gap
	write(300,*)
	write(300,*)
	write(300,*)"TB kpoint cbm"
	write(300,"(3E18.8)")kpt(nkc,1),kpt(nkc,2),kpt(nkc,3)
	write(300,*)	
	write(300,*)"TB kpoint vbm"
	write(300,"(3E18.8)")kpt(nkv,1),kpt(nkv,2),kpt(nkv,3)
	write(300,*)		 
	write(300,*)"TB kpoint direct band gap"
	write(300,"(3E18.8)")kpt(nkgap,1),kpt(nkgap,2),kpt(nkgap,3)
	call flush(300)
	

										

	!calculate Mmn(k,q)
	
	allocate(mmnkq(w90basis,w90basis,nkpt,nkpt))
	
	!$omp parallel do collapse(4) default(shared) private(i,j,k,l)
	do i=1,nkpt
	 do j=1,w90basis
	  do k=1,nkpt
	   do l=1,w90basis
	
		call ovpsqr(dft,w90basis,kvector(i,j,:),kovp(i,:,:),kvector(k,l,:),kovp(k,:,:),mmnkq(j,l,i,k))
	
	   end do
	  end do
	 end do  
	end do	
	!$omp end parallel do
		
	!if (dft .eq. "S") then
	!	 deallocate(kovp)
	!else
	! continue
	!end if	
	
	!deallocate(kvector)	
	
	write(300,*)
	write(300,*)"Mmn(k,q) calculated"	
	write(300,*)	
	call flush(300)
	
	!calculate self-energy exchange part
	
	allocate(selfx(w90basis,nkpt))
	
	!$omp parallel do default(shared) private(i,j)
	do i=1,nkpt
	 do j=1,w90basis
	
	  call self_en_x(w90basis,sysdim,rlat,j,kpt(i,:),ngrid,nkpt,kpt,nocpk,coultypegw,edielgw,&
	               lcgw,ezgw,wgw,r0gw,ktolgw,mmnkq(:,:,i,:),selfx(j,i))
	
	 end do
	end do
	!$omp end parallel do
	

	write(300,*)
	write(300,*)"self_x(n,q) calculated"	
	write(300,*)	
	call flush(300)
	
if (selfxonly) then

	!if (dft .eq. "S") then
	!	 deallocate(qovp)
	!else
	! continue
	!end if	
	
	!deallocate(qvector)
	
	 allocate(selfc(w90basis,nkpt),znk(w90basis,nkpt))
	 
	 selfc = 0.0
	 
	 znk = 1.0
	 
	 allocate(pol(nkpt,nomega))
	 
	 allocate(w0(nkpt,nomega),vqa(nkpt,nomega))	


else

	!calculate Mmn(q,q')
	
	allocate(mmnqq(w90basis,w90basis,nkpt,nkpt))
	
	!$omp parallel do collapse(4) default(shared) private(i,j,k,l)
	do i=1,nkpt
	 do j=1,w90basis
	  do k=1,nkpt
	   do l=1,w90basis
	
		call ovpsqr(dft,w90basis,kvector(i,j,:),kovp(i,:,:),kvector(k,l,:),kovp(k,:,:),mmnqq(j,l,i,k))
	
	   end do
	  end do
	 end do  
	end do	
	!$omp end parallel do
		

	
	write(300,*)
	write(300,*)"Mmn(q,q') calculated"	
	write(300,*)	
	call flush(300)
	
	if (dft .eq. "S") then
		 deallocate(kovp)
	else
	 continue
	end if	
	
	deallocate(kvector)
	
	!calculate kmesh (q+q')
	
	allocate(qptpq(nkpt,nkpt,3))
	
	do i=1,nkpt
	 do j=1,nkpt
	 
	 	qptpq(i,j,:) = kpt(i,:) + kpt(j,:)
	 
	 end do
	end do							

	!calculate eigenvalues and eigenvectors in kmesh (q+q')
	
	allocate(qpqeigv(nkpt,nkpt,w90basis))
	
	!$omp parallel do default(shared) private(i,j,k,l,h,eaux,vaux)
	do i=1,nkpt
	 do j=1,nkpt		

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

		!if (dft .eq. "S") then
		 !call overlap(w90basis,nvec,rvec,ovp,qptpq(i,j,1),qptpq(i,j,2),qptpq(i,j,3),qpqovp(i,j,:,:))
		!else
		 !continue
		!end if
		
		
		do k=1,w90basis
				qpqeigv(i,j,k)= eaux(k)


		end do	
		
		!do l=1,w90basis


		!	do h=1,w90basis

		!		qpqvector(i,j,l,h)=vaux(l,h)


		!	end do
			

		!end do			

	 end do
	end do
	!$omp end parallel do	

	write(300,*)
	write(300,*)"Eigenvalues and Eigenvectors in q+q'-mesh calculated"	
	write(300,*)	
	call flush(300)
	
	!calculate kmesh (k-q)
	
	allocate(kptmq(nkpt,nkpt,3))
	
	do i=1,nkpt
	 do j=1,nkpt
	 
	 	kptmq(i,j,:) = kpt(i,:) - kpt(j,:)
	 
	 end do
	end do				

	!calculate eigenvalues and eigenvectors in kmesh (k-q)
	
	allocate(kmqeigv(nkpt,nkpt,w90basis))	
	

	
	!$omp parallel do default(shared) private(i,j,k,l,h,eaux,vaux)
	do i=1,nkpt
	 do j=1,nkpt		

#ifdef MKL
		call MKL_SET_NUM_THREADS(1)
#endif		


#ifdef AOCL		
		call bli_thread_set_num_threads(1)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(1)
		
#endif

		call eigsys(nthreads,dft,systype,scs,exc,nocpaux,ffactor,kptmq(i,j,1),kptmq(i,j,2),kptmq(i,j,3),w90basis,nvec,&
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

		!if (dft .eq. "S") then
		 !call overlap(w90basis,nvec,rvec,ovp,qptmq(i,j,1),qptmq(i,j,2),qptmq(i,j,3),qmqovp(i,j,:,:))
		!else
		 !continue
		!end if
		
		
		do k=1,w90basis
				kmqeigv(i,j,k)= eaux(k)


		end do	
		
		!do l=1,w90basis


		!	do h=1,w90basis

		!		qmqvector(i,j,l,h)=vaux(l,h)


		!	end do
			

		!end do			

	 end do
	end do
	!$omp end parallel do
		
	write(300,*)
	write(300,*)"Eigenvalues and Eigenvectors in k-q -mesh calculated"	
	write(300,*)
	call flush(300)
	
	!calculate Pi0(q,w) and W0(q,w)
	
	allocate(pol(nkpt,nomega))	
	
	!$omp parallel do default(shared) private(i)
	do i=1,nkpt
	 do j=1,nomega
	
		call polarization(fermishift,sysdim,ngrid,nkpt,rlat,omega(j),w90basis,systype,&
		                   smegw,keigv,qpqeigv(i,:,:),mmnqq(:,:,i,:),pol(i,j))
	
	
	 	!pol(i,:) = 0.0
	 end do	
	end do
	!$omp end parallel do	
	
	deallocate(qpqeigv)	
	
	write(300,*)
	write(300,*)"Pi0(q,w)  calculated"	
	write(300,*)	
	call flush(300)
	
	allocate(w0(nkpt,nomega),vqa(nkpt,nomega))	
	
	!$omp parallel do default(shared) private(i)
	do i=1,nkpt
	
		call w0coul(kpt(i,:),nomega,rlat,ngrid,nkpt,coultypegw,edielgw,lcgw,ezgw,wgw,r0gw,ktolgw,pol(i,:),vqa(i,:),w0(i,:))
	

	
	end do
	!$omp end parallel do	
	
	write(300,*)
	write(300,*)"W0(q,w)  calculated"	
	write(300,*)	
	call flush(300)
	
	do i=1,nomega
	
		write(302,*) omega(i),real(pol(1,i)),aimag(pol(1,i))
	
	end do	



	write(303,*) ngrid(1),ngrid(2),ngrid(3)
	write(303,*) 
	
	do i=1,nkpt
	
		write(303,*) '#kpt n:',i,kpt(i,1),kpt(i,2),kpt(i,3)
		write(303,*) '#Omega Re(Pi0(q,omega)) Im(Pi0(q,omega)) V(q) Re(W0(q,omega)) Im(W0(q,omega))'
		write(303,*)		
	
	 do j=1,nomega
	
		write(303,"(6E18.8)") omega(j),real(pol(i,j)),aimag(pol(i,j)),vqa(i,j),real(w0(i,j)),aimag(w0(i,j))
	
	 end do
	 
	 	write(303,*)
	end do		
	
	deallocate(qptpq,kptmq)
	
	!calculate self-energy correlation part
	
	allocate(selfc(w90basis,nkpt),znk(w90basis,nkpt))
	
	!$omp parallel do default(shared) private(i,j)
	do i=1,nkpt
	 do j=1,w90basis
	


	  call self_c_znk(w90basis,sysdim,rlat,ngrid,nkpt,fermishift,j,&
	                  keigv(i,j),kmqeigv(i,:,:),kpt(i,:),omega,nomega,w0,mmnkq(:,:,i,:),smegw,selfc(j,i),znk(j,i))
	
	 end do
	end do	
	!$omp end parallel do
	
	write(300,*)
	write(300,*)"self_c(n,q) calculated"	
	write(300,*)	
	call flush(300)			

end if	

	

	deallocate(ovp)
	
	deallocate(mmnkq)

	deallocate(pol,w0,vqa)	
	
	!correct eigenvalues

	allocate(eigcor(w90basis,nkpt),eigcoravg(w90basis,nkpt))
	
	allocate(encor(w90basis,nkpt),encoravg(w90basis,nkpt))
	
	do i=1,nkpt
	 do j=1,w90basis
	
		eigcor(j,i) = znk(j,i)*(selfx(j,i)+real(selfc(j,i)))
	
	 end do
	end do	
	
	write(301,*) w90basis
	write(301,*) nkpt	
	write(301,*) ngrid(1),ngrid(2),ngrid(3)
	write(301,*)
	
	do i=1,nkpt
	
		write(301,*) '#kpt n:',i,kpt(i,1),kpt(i,2),kpt(i,3)
		write(301,*) '#EnTB G0W0cor selfx Re(selfc) Im(selfc) Znk'
		write(301,*)
		
	 do j=1,w90basis
	 
	 	write(301,'(6E18.8)') ebands(j,i),eigcor(j,i),selfx(j,i),real(selfc(j,i)),aimag(selfc(j,i)),znk(j,i)
	 
	 end do
	
		write(301,*)
	
	end do
	
	allocate(selfcavg(w90basis),selfxavg(w90basis),znkavg(w90basis))
	
	call selfavg(nkpt,w90basis,selfx,selfc,znk,selfxavg,selfcavg,znkavg)
	
	do i=1,nkpt
	 do j=1,w90basis
	
		eigcoravg(j,i) = znkavg(j)*(selfxavg(j)+real(selfcavg(j)))
	
	 end do
	end do	
	
	write(304,*) '#G0W0cor-avg selfx-avg Re(selfc-avg) Im(selfc-avg) Znk-avg'
	
	do i=1,w90basis
	
	  write(304,*) eigcoravg(i,1),selfxavg(i),real(selfcavg(i)),aimag(selfcavg(i)),znkavg(i)
	
	end do				
	

	
	
	write(305,*) w90basis
	write(305,*) nkpt	
	write(305,*) ngrid(1),ngrid(2),ngrid(3)
	write(305,*)
	
	do i=1,nkpt
	
		write(305,*) '#kpt n:',i,kpt(i,1),kpt(i,2),kpt(i,3)
		write(305,*) '#EnTB G0W0cor selfx Re(selfc) Im(selfc) Znk'
		write(305,*)
		
	 do j=1,w90basis
	 
	 	if (abs(eigcor(j,i)-eigcoravg(j,i)) .gt. rnmgw ) then
	 	
	 		eigcor(j,i)= eigcoravg(j,i)
	 		selfx(j,i) = selfxavg(j)
	 		selfc(j,i) = selfcavg(j)
	 		znk(j,i) = znkavg(j)
	 	
	 	end if
	 
	 	write(305,'(6E18.8)') ebands(j,i),eigcor(j,i),selfx(j,i),real(selfc(j,i)),aimag(selfc(j,i)),znk(j,i)
	 
	 end do
	
		write(305,*)
	
	end do	
	
	do i=1,w90basis
	
		
	 do j=1,nkpt
	 
	 	encor(i,j) = ebands(i,j)+eigcor(i,j)
	 	encoravg(i,j) = ebands(i,j)+eigcoravg(i,j)	 
	 
	 end do	
	
	end do		

	!estimating GW gap
	
	write(300,*)
	write(300,*)"QP Eigenvalues and Eigenvectors in k-mesh calculated"	
	write(300,*)
	
	call gapfinder(w90basis,nkpt,nocpk,encor,nkc,nkv,nkgap,gap,cbm,vbm)
	
	write(300,*) "G0W0 QP fundamental band gap (eV):",cbm-vbm
	write(300,*) "G0W0 QP direct band gap (eV):",gap
	write(300,*)
	write(300,*)
	write(300,*)"G0W0 QP kpoint cbm"
	write(300,"(3E18.8)")kpt(nkc,1),kpt(nkc,2),kpt(nkc,3)
	write(300,*)	
	write(300,*)"G0W0 QP kpoint vbm"
	write(300,"(3E18.8)")kpt(nkv,1),kpt(nkv,2),kpt(nkv,3)
	write(300,*)		 
	write(300,*)"G0W0 QP kpoint direct band gap"
	write(300,"(3E18.8)")kpt(nkgap,1),kpt(nkgap,2),kpt(nkgap,3)	
	call flush(300)	
							
	deallocate(encor)
	
	write(300,*)
	write(300,*)"QP Eigenvalues and Eigenvectors avg in k-mesh calculated"	
	write(300,*)
	
	call gapfinder(w90basis,nkpt,nocpk,encoravg,nkc,nkv,nkgap,gap,cbm,vbm)
	
	write(300,*) "G0W0 QP avg fundamental band gap (eV):",cbm-vbm
	write(300,*) "G0W0 QP avg direct band gap (eV):",gap
	write(300,*)
	write(300,*)
	write(300,*)"G0W0 QP avg kpoint cbm"
	write(300,"(3E18.8)")kpt(nkc,1),kpt(nkc,2),kpt(nkc,3)
	write(300,*)	
	write(300,*)"G0W0 QP avg kpoint vbm"
	write(300,"(3E18.8)")kpt(nkv,1),kpt(nkv,2),kpt(nkv,3)
	write(300,*)		 
	write(300,*)"G0W0 QP avg kpoint direct band gap"
	write(300,"(3E18.8)")kpt(nkgap,1),kpt(nkgap,2),kpt(nkgap,3)	
	call flush(300)
	
	deallocate(selfcavg,selfxavg,znkavg,encoravg,ebands)
	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)
	deallocate(kpt)
	
	call cpu_time(tf)
	call date_and_time(VALUES=values2)

	write(300,*)
	write(300,*) 'end','   ','month',values2(2),'day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(300,*)
	call flush(300)		
	
	close(200)
	close(300)
	close(301)
	close(302)
	close(303)
	close(304)
	close(305)
	
end subroutine sgw_mesh 	
	
