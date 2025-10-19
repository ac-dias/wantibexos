subroutine efmass_num(nthreads,dft,outputfolder,params,emfile,h,nocpf,fermishift,exc,mag,sysdim)


	use omp_lib
	use hamiltonian_input_variables
	implicit none
	
	integer :: i,j,k,erro,j2,n
	integer :: nthreads
	character(len=70) :: params   !parametros TB
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: emfile   !parametros TB
	character(len=5) :: sysdim   !parametros TB
	real,parameter :: pi=acos(-1.)
	real,allocatable,dimension(:) :: eigv,eigvm
	complex,allocatable,dimension(:,:) :: autovetores
	real :: exc,sme
	real,dimension(3) :: mag

	
	real,allocatable,dimension(:,:) :: kpt
	integer,allocatable,dimension(:) :: nbnd
	real :: flag
	real,dimension(3,3) :: emtensor
	character(len=1) :: coordtype,dft
	real :: kx,ky,kz,h
	
	integer :: nocpf
	real :: fermishift
	
	character(len=200) outputfile
	
	write (outputfile,"(a10,F6.4,a7)") "em_tensor_",h,"_dK.dat"
	
	OPEN(UNIT=202, FILE= emfile,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening effective mass input file"
    	
	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//trim(outputfile),STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening effective mass tensor output file"    	
	
	call OMP_SET_NUM_THREADS(nthreads)
	
	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	allocate(ovp(nvec,w90basis,w90basis))

	call hamiltonian_input_read(200,params)

	end select	
	
		!reading efective mass input file
	read(202,*) n
	read(202,*) coordtype
	
	allocate(kpt(n,3))
	allocate(nbnd(n))
	allocate(eigvm(w90basis))
	allocate(autovetores(w90basis,w90basis),eigv(w90basis))
	
	do i=1,n
	
		read(202,*) kpt(i,1),kpt(i,2),kpt(i,3),nbnd(i),flag
	
	end do
	
	if (coordtype .eq. "D") then
	
		do i=1,n
		

		    
			call conv_kpt_d2c(rlat(1,:),rlat(2,:),rlat(3,:),kpt(i,:))
		
		end do
	
	else
		continue
	end if
	
	do i=1,n
	
		kx= kpt(i,1)
		ky= kpt(i,2)
		kz= kpt(i,3)
		
		
		call num_effective_mass_calculator(dft,systype,nthreads,scs,exc,kx,ky,kz,ffactor,w90basis,nvec,rlat,rvec,hopmatrices,&
                                     ihopmatrices,ovp,efermi,eigvm,nbnd(i),h,emtensor,nocpf,fermishift,mag)
		
	        write(300,*) "##################################################################" 
                write(300,*) "kx,ky,kz,band number,energy" 
                write(300,"(3F15.4,I0,1F15.4)") kx,ky,kz,nbnd(i),eigvm(nbnd(i))
                write(300,*)               
                write(300,*) "Effective Mass Tensor "
                if (trim(sysdim) .eq. "2D") then
                  write(300,"(3F15.4)") real(emtensor(1,1)), real(emtensor(1,2))
                  write(300,"(3F15.4)") real(emtensor(2,1)), real(emtensor(2,2))               
                else if (trim(sysdim) .eq. "1D") then
                  write(300,"(3F15.4)") real(emtensor(1,1))            
                else
                  write(300,"(3F15.4)") real(emtensor(1,1)), real(emtensor(1,2)),real(emtensor(1,3))
                  write(300,"(3F15.4)") real(emtensor(2,1)), real(emtensor(2,2)),real(emtensor(2,3))               
                  write(300,"(3F15.4)") real(emtensor(3,1)), real(emtensor(3,2)),real(emtensor(3,3))
                end if
                !write(300,*)               
                !write(300,*) "Effective Mass Tensor - imaginary part"
                !write(300,"(3F15.4)") aimag(emtensor(1,1)), aimag(emtensor(1,2)),aimag(emtensor(1,3))
                !write(300,"(3F15.4)") aimag(emtensor(2,1)), aimag(emtensor(2,2)),aimag(emtensor(2,3))               
                !write(300,"(3F15.4)") aimag(emtensor(3,1)), aimag(emtensor(3,2)),aimag(emtensor(3,3))                                
                write(300,*) "##################################################################"	
		

	end do		

	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)
	deallocate(eigv,eigvm)
	deallocate(autovetores)
	deallocate(kpt)
	deallocate(nbnd)
	deallocate(ovp)		
		
	close(200)
	close(202)
	close(300)		
		
end subroutine efmass_num		
		
		
subroutine efmass(nthreads,dft,outputfolder,params,emfile,nocpf,fermishift,exc,mag,sysdim)

	use omp_lib
	use hamiltonian_input_variables
	implicit none

	integer :: i,j,k,erro,j2,n
	integer :: nthreads
	character(len=70) :: params   !parametros TB
	character(len=70) :: outputfolder    !pasta saida
	character(len=70) :: emfile   !parametros TB
	real,parameter :: pi=acos(-1.)
	real,allocatable,dimension(:) :: eigv
	complex,allocatable,dimension(:,:) :: autovetores
	real :: exc,sme
	
	integer :: nocpf
	real :: fermishift

	
	real,allocatable,dimension(:,:) :: kpt
	integer,allocatable,dimension(:) :: nbnd
	real :: flag
	complex,dimension(3,3) :: emtensor
	character(len=1) :: coordtype,dft
	real :: kx,ky,kz
	real,dimension(3) :: mag
	character(len=5) :: sysdim
	
	OPEN(UNIT=202, FILE= emfile,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening effective mass input file"
    	
	!OUTPUT : criando arquivos de saida
	OPEN(UNIT=300, FILE=trim(outputfolder)//"em_tensor.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening effective mass tensor output file"    	
	
	call OMP_SET_NUM_THREADS(nthreads)
	
	select case (dft)
	
	case ("S")
		
	call hamiltonian_nort_input_read(200,params)
	
	
	case default
	
	allocate(ovp(nvec,w90basis,w90basis))

	call hamiltonian_input_read(200,params)

	end select	
	
	allocate(autovetores(w90basis,w90basis),eigv(w90basis))
	
	!reading efective mass input file
	read(202,*) n
	read(202,*) coordtype
	
	allocate(kpt(n,3))
	allocate(nbnd(n))
	
	do i=1,n
	
		read(202,*) kpt(i,1),kpt(i,2),kpt(i,3),nbnd(i),flag
	
	end do
	
	if (coordtype .eq. "D") then
	
		do i=1,n
		

		    
			call conv_kpt_d2c(rlat(1,:),rlat(2,:),rlat(3,:),kpt(i,:))
		
		end do
	
	else
		continue
	end if
	
	write(*,*) sysdim
	write(*,*) coordtype
	
	do i=1,n
	
		kx= kpt(i,1)
		ky= kpt(i,2)
		kz= kpt(i,3)
		
		call eigsys(nthreads,dft,systype,scs,exc,nocp,ffactor,kx,ky,kz,w90basis,nvec,rlat,rvec,hopmatrices,&
		            ihopmatrices,ovp,efermi,eigv,autovetores,nocpf,fermishift,mag)
		
		call effective_mass_calculator(dft,systype,kx,ky,kz,ffactor,w90basis,autovetores(nbnd(i),:),nvec,&
                                               rlat,rvec,hopmatrices,ihopmatrices,ovp,exc,mag,emtensor)
                                               
                write(300,*) "##################################################################" 
                write(300,*) "kx,ky,kz,band number,energy" 
                write(300,"(3F15.4,I8,1F15.4)") kx,ky,kz,nbnd(i),eigv(nbnd(i))
                write(300,*)               
                write(300,*) "Effective Mass Tensor "
                if (trim(sysdim) .eq. "2D") then
                  write(300,"(3F15.4)") real(emtensor(1,1)), real(emtensor(1,2))
                  write(300,"(3F15.4)") real(emtensor(2,1)), real(emtensor(2,2))               
                else if (trim(sysdim) .eq. "1D") then
                  write(300,"(3F15.4)") real(emtensor(1,1))            
                else
                  write(300,"(3F15.4)") real(emtensor(1,1)), real(emtensor(1,2)),real(emtensor(1,3))
                  write(300,"(3F15.4)") real(emtensor(2,1)), real(emtensor(2,2)),real(emtensor(2,3))               
                  write(300,"(3F15.4)") real(emtensor(3,1)), real(emtensor(3,2)),real(emtensor(3,3))
                end if
                !write(300,*)               
                !write(300,*) "Effective Mass Tensor - imaginary part"
                !write(300,"(3F15.4)") aimag(emtensor(1,1)), aimag(emtensor(1,2)),aimag(emtensor(1,3))
                !write(300,"(3F15.4)") aimag(emtensor(2,1)), aimag(emtensor(2,2)),aimag(emtensor(2,3))               
                !write(300,"(3F15.4)") aimag(emtensor(3,1)), aimag(emtensor(3,2)),aimag(emtensor(3,3))                                
                write(300,*) "##################################################################"	
	end do
	
	
	deallocate(rvec,hopmatrices,ihopmatrices,ffactor)
	deallocate(eigv)
	deallocate(autovetores)
	deallocate(kpt)
	deallocate(nbnd)
	deallocate(ovp)
	
	close(200)
	close(202)
	close(300)

end subroutine efmass
		
		
		
		
		
		
		
		
