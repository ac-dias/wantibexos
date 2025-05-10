subroutine boltztransport(nthreads,outputfolder,ngrid,nsteps,smeboltz,params,exc,mag,mshift,dft,nocpf,fermishift, &
                           mu0,muf,nmu,btemp,klat,elft,hlft)

	use omp_lib
	use hamiltonian_input_variables
	use constantsboltz
	
	implicit none
	
	integer :: i,j,k1,erro
	
	integer :: nthreads
	integer,dimension(3) :: ngrid
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: params   !parametros TB	
	character(len=1) :: dft
	real,dimension(3) :: mag
	real,dimension(3) :: mshift	
	real :: smeboltz,exc
	integer :: nocpf
	real :: fermishift	
	
	integer :: nsteps,nmu
	real :: mu0,muf,klat
	real,dimension(3) :: elft,hlft
	real :: btemp
	real :: e0,ef
	
	integer :: ngkpt
	
	real,allocatable,dimension(:) :: eaux
	complex,allocatable,dimension(:,:) :: vaux 
	
	real,allocatable,dimension(:,:) :: autovalores
	complex,allocatable,dimension(:,:,:) :: autovetores
	
	real,allocatable,dimension(:,:) :: kpts
	integer,allocatable,dimension(:) :: nocpj
	
	real,allocatable,dimension(:) :: enint
	real,allocatable,dimension(:) :: mu
	real,allocatable,dimension(:,:,:)  :: elvel
	real,allocatable,dimension(:,:) :: tdfres
	real,allocatable,dimension(:,:) :: dfermi
	
	real,allocatable,dimension(:,:) :: kij,seij,sij,ztij,pfij,sigsij
	real :: auxkij,auxsij,auxseij,auxpfij
	real :: dfermidist
	
	real:: t0,tf
	integer,dimension(8) :: values,values2
	
	character(len=200) :: file1,file2,file3,file4,file5,file6
	
	!real,parameter :: a2m = 1.0E-10 !convert angstrom to meter
	!real,parameter :: echarge = 1.602176620898E-19 !electron charge in coulomb

	
	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=400, FILE=trim(outputfolder)//"log_boltzmann_transport.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening boltzmann transport log file"
    	
    	write (file1,"(a17,F7.3,a6)") "elec_condutivity_",btemp,"_K.dat"
	write (file2,"(a8,F7.3,a6)") "seebeck_",btemp,"_K.dat"
	write (file3,"(a30,F7.3,a6)") "electron_thermal_conductivity_",btemp,"_K.dat"
	write (file4,"(a13,F7.3,a6)") "power_factor_",btemp,"_K.dat"
	write (file5,"(a3,F7.3,a6)") "ZT_",btemp,"_K.dat"
	write (file6,"(a4,F7.3,a6)") "TDF_",btemp,"_K.dat"    	
	
	OPEN(UNIT=300, FILE=trim(outputfolder)//trim(file1),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening electrical conductivity output file"
	OPEN(UNIT=301, FILE=trim(outputfolder)//trim(file2),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening seebeck coeficient output file"    	
	OPEN(UNIT=302, FILE=trim(outputfolder)//trim(file3),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening electron thermal conductivity output file"     	
	OPEN(UNIT=303, FILE=trim(outputfolder)//trim(file4),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening power factor output file"  
	OPEN(UNIT=304, FILE=trim(outputfolder)//trim(file5),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening ZT - figure of merit output file" 
	OPEN(UNIT=305, FILE=trim(outputfolder)//trim(file6),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening TDF - transport distribution function output file"  
	!OPEN(UNIT=306, FILE=trim(outputfolder)//"DFERMI.dat",STATUS='unknown', IOSTAT=erro)
    	!if (erro/=0) stop "Error opening derivative fermi distribution output file"     	  	   	  	
    	
	call cpu_time(t0)
	call date_and_time(VALUES=values)
	    	
 	call OMP_SET_NUM_THREADS(nthreads)
 	
	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	allocate(ovp(nvec,w90basis,w90basis))

	call hamiltonian_input_read(200,params)

	end select 	   	
    	
   	ngkpt = ngrid(1)*ngrid(2)*ngrid(3) 	
   	allocate(kpts(ngkpt,3))
   	
   	call monhkhorst_pack(ngrid(1),ngrid(2),ngrid(3),mshift,rlat(1,:),rlat(2,:),rlat(3,:),kpts)
   	
	write(400,*)
	write(400,*)
	write(400,*)
	write(400,*) 'threads:', nthreads
	write(400,*)
	write(400,*) 'grid:',ngrid(1),ngrid(2),ngrid(3)
	write(400,*)
	write(400,"(A13,3F15.4)") 'kmesh shift:',mshift(1),mshift(2),mshift(3)
	write(400,*)
	write(400,*) 'temperature:',btemp
	write(400,*) 'integration steps:',nsteps
	write(400,*) 'smearing:',smeboltz
	write(400,*) 		
	write(400,*) 'electron lifetime - x:',elft(1)
	write(400,*) 'electron lifetime - y:',elft(2)
	write(400,*) 'electron lifetime - z:',elft(3)		
	write(400,*) 'hole lifetime - x:',hlft(1)
	write(400,*) 'hole lifetime - y:',hlft(2)
	write(400,*) 'hole lifetime - z:',hlft(3)		
	write(400,*) 
	write(400,*) 'lattice thermal conductivity:',klat					
	write(400,*)
	write(400,*) 'mu steps:',nmu
	write(400,*) 'mu initial:', mu0,' ','mu final:',muf
	write(400,*)				
	write(400,*) 'begin','  ','day',values(3),'',values(5),'hours',values(6),'min',values(7),'seg'
	write(400,*) 

	call flush(400)	   	
   	
   	
   	allocate(nocpj(ngkpt))
   	
   	allocate(autovalores(ngkpt,w90basis))
   	allocate(autovetores(ngkpt,w90basis,w90basis))
   	
   	allocate(vaux(w90basis,w90basis),eaux(w90basis))
   	
   	
   	e0= 0.0
   	ef= 0.0
   	
   	!$omp parallel do default(shared) private(j,i,k1,eaux,vaux)
   	do j=1,ngkpt
   	
#ifdef MKL
		call MKL_SET_NUM_THREADS(1)
#endif		


#ifdef AOCL		
		call bli_thread_set_num_threads(1)
#endif

#ifdef OPENBLAS
		call OPENBLAS_SET_NUM_THREADS(1)
		
#endif		
   		
		call eigsys(nthreads,dft,systype,scs,exc,nocpj(j),ffactor,kpts(j,1),kpts(j,2),kpts(j,3),w90basis,nvec,rlat,rvec,hopmatrices,&
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
		
		
		
		do i=1,w90basis
		
			autovalores(j,i) = eaux(i)
		
		end do
		
		do i=1,w90basis
		 do k1=1,w90basis
		
			autovetores(j,i,k1) = vaux(i,k1)
		
		 end do
		end do
		
		if (eaux(1) .lt. e0 ) then
		
			e0 = eaux(1)
			
		end if
		
		
		if (eaux(w90basis) .gt. ef ) then
		
			ef = eaux(w90basis)
			
		end if				
   	
   	end do
   	!$omp end parallel do
   	
	write(400,*) 'eigenvalues and eigenvectors calculated'
	call flush(400)   	
   	
   	deallocate(eaux,vaux)
   	
   	allocate(elvel(w90basis,ngkpt,3))
   	
   	!$omp parallel do default(shared) private(i,j)
   	do i=1,ngkpt
   	 do j=1,w90basis
   	 
   	 	call bndvel(kpts(i,1),kpts(i,2),kpts(i,3),ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
	                     ihopmatrices,autovetores(i,j,:),elvel(j,i,1),elvel(j,i,2),elvel(j,i,3))
   	 
   	 	!write(401,*) elvel(j,i,1),elvel(j,i,2),elvel(j,i,3)
   	 
   	 end do
   	end do
	!$omp end parallel do
	
	write(400,*) 'band velocities calculated'
	call flush(400)	
   	
   	allocate(enint(nsteps))
   	allocate(tdfres(nsteps,3))
   	
   	!write(*,*) e0,ef
   	
   	!$omp parallel do default(shared) private(i)
   	do i=1,nsteps
   	
   	 enint(i) = e0 + (ef-e0)*((real(i)-1.0)/(real(nsteps)-1.0))   
   	
   	call tdf(rlat,w90basis,ngkpt,autovalores,nocpj(i),smeboltz,enint(i),elft,hlft,elvel(:,:,1),&
   	          elvel(:,:,2),elvel(:,:,3),tdfres(i,:))
   	 
   	
   	end do
   	!$omp end parallel do
   	
   	write(305,*) "#","  ", "energy","  ","xx","  ","yy","  ","zz"
   	
   	do i=1,nsteps
   	   write(305,"(1F18.6,3E18.8)") enint(i),tdfres(i,1),tdfres(i,2),tdfres(i,3)
   	end do
   	
	write(400,*) 'transport distribution function calculated'
	call flush(400)   	
   	
   	deallocate(nocpj)
   	deallocate(kpts)
   	deallocate(autovetores,autovalores)
   	deallocate(elvel)
   	deallocate(ovp)
   	
   	allocate(mu(nmu))
   	allocate(sij(nmu,3))
   	allocate(seij(nmu,3))
   	allocate(kij(nmu,3))
   	allocate(sigsij(nmu,3))
   	
   	!allocate(dfermi(6,nsteps))
   	
   	!write(306,*) "#","  ", "mu","  ","energy","  ","dfermi"
   	
   	!do i=1,nmu
   	!        mu(i) = mu0 + (muf-mu0)*((real(i)-1.0)/(real(nmu)-1.0))
   	! do j=1,nsteps
   	 
   	! 	dfermi(i,j) = dfermidist(mu(i),btemp,enint(j))
   	 	
   	! 	write(306,"(2F18.6,1E18.8)") mu(i),enint(j),dfermi(i,j)
   	 	
   	! end do
   	 
   	!end do
   	
   	!write(400,*) 'Fermi Distribution derivative calculated'
	!call flush(400) 
  	 
        !write(*,*) mu0,muf  	   	

   	!$omp parallel do default(shared) private(i)   	
   	do i=1,nmu
   	
   		mu(i) = mu0 + (muf-mu0)*((real(i)-1.0)/(real(nmu)-1.0))
   	
   		call eleccond(nsteps,enint,tdfres,mu(i),btemp,sij(i,:))
   		!call seebeck(nsteps,enint,tdfres,mu(i),btemp,sij(i,:),seij(i,:))
		call sigmas(nsteps,enint,tdfres,mu(i),btemp,sigsij(i,:))   		
   		call electermcond(nsteps,enint,tdfres,mu(i),btemp,kij(i,:))
   		
   		if ( mu(i) .gt. 0.0) then
   		
   			sij(i,1) = sij(i,1)*elft(1)
   			sij(i,2) = sij(i,2)*elft(2)
   			sij(i,3) = sij(i,3)*elft(3)   
   			
   			kij(i,1) = kij(i,1)*elft(1)
   			kij(i,2) = kij(i,2)*elft(2)
   			kij(i,3) = kij(i,3)*elft(3)
   			
   			seij(i,1) = (sigsij(i,1)*elft(1))/sij(i,1)
   			seij(i,2) = (sigsij(i,2)*elft(2))/sij(i,2)
   			seij(i,3) = (sigsij(i,3)*elft(3))/sij(i,3)   			   						   						   			
   		
   		else
   		
   			sij(i,1) = sij(i,1)*hlft(1)
   			sij(i,2) = sij(i,2)*hlft(2)
   			sij(i,3) = sij(i,3)*hlft(3)
   			
   			kij(i,1) = kij(i,1)*hlft(1)
   			kij(i,2) = kij(i,2)*hlft(2)
   			kij(i,3) = kij(i,3)*hlft(3) 
   			
   			seij(i,1) = (sigsij(i,1)*hlft(1))/sij(i,1)
   			seij(i,2) = (sigsij(i,2)*hlft(2))/sij(i,2)
   			seij(i,3) = (sigsij(i,3)*hlft(3))/sij(i,3)    			   		
   		
   		end if
   		
   		!write(403,*) mu(i)
 
   	end do
   	!$omp end parallel do   	
   	
   	allocate(pfij(nmu,3))
   	allocate(ztij(nmu,3)) 
   	
   	deallocate(enint,tdfres)
   	
   	auxkij = echarge/a2m
   	auxsij = 1.0/(echarge*a2m)
   	auxseij = echarge
   	auxpfij = echarge/a2m
   	
   	kij  = kij*auxkij
   	sij  = sij*auxsij
   	seij = seij*auxseij
   	!pfij = pfij*auxpfij
   	
   	!$omp parallel do default(shared) private(i) 
 	do i=1,nmu
 	 do j=1,3
 	 
 	 	pfij(i,j) = sij(i,j)*((seij(i,j))**2)
 	 	
 	 	ztij(i,j) = (pfij(i,j)*btemp)/(kij(i,j)+klat)
 	 
 	 end do
 	end do  	
   	!$omp end parallel do
   	
   	


	write(300,*) "#","  ", "mu","  ","xx","  ","yy","  ","zz"
	write(301,*) "#","  ", "mu","  ","xx","  ","yy","  ","zz"
	write(302,*) "#","  ", "mu","  ","xx","  ","yy","  ","zz"
	write(303,*) "#","  ", "mu","  ","xx","  ","yy","  ","zz"
	write(304,*) "#","  ", "mu","  ","xx","  ","yy","  ","zz"				

   	!"(7F15.8)"
   	do i=1,nmu
   	
   		write(300,"(1F18.6,3E18.8)") mu(i),sij(i,1),sij(i,2),sij(i,3)
     		write(301,"(1F18.6,3E18.8)") mu(i),seij(i,1),seij(i,2),seij(i,3)
     		write(302,"(1F18.6,3E18.8)") mu(i),kij(i,1),kij(i,2),kij(i,3)		
   		write(303,"(1F18.6,3E18.8)") mu(i),pfij(i,1),pfij(i,2),pfij(i,3)
   		write(304,"(1F18.6,3E18.8)") mu(i),ztij(i,1),ztij(i,2),ztij(i,3)
   		  		   		   	
   	end do
   	
	write(400,*) 'Termoelectric properties finished'
	call flush(400)   
	
	call cpu_time(tf)
	call date_and_time(VALUES=values2)

	write(400,*)
	write(400,*) 'end','   ','day',values2(3),'',values2(5),'hours',values2(6),'min',values2(7),'seg'
	write(400,*)	
	

   	deallocate(mu)
   	deallocate(sij)
   	deallocate(seij)
   	deallocate(kij)		
   	deallocate(pfij)
   	deallocate(ztij)
   	deallocate(sigsij)  
   	!deallocate(dfermi) 	
   	
   	close(400)
   	close(300)
   	close(301)
   	close(302)
   	close(303)
   	close(304)
   	close(305)
   	!close(306)   	
   	
   	
    		

end subroutine boltztransport
