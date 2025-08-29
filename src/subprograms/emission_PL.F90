subroutine emissionopt(nthreads,mode,outputfolder,noptics,ngrid,nc,nv,exctemp,eta,emin)

	implicit none
	integer :: erro,i,j,nthreads
	character(len=3) :: mode
	character(len=70) :: outputfolder
	character(len=5) :: cflag
	integer,dimension(3) :: ngrid
	real :: noptics
	integer :: nopt,nc,nv,dimbse
	real :: emin,eta,exctemp
	real :: aux
	real,allocatable,dimension(:,:) :: dipole
	real,allocatable,dimension(:) :: omega
	real,allocatable,dimension(:,:) :: refractive
	real,allocatable,dimension(:,:) :: plres
	
	call OMP_SET_NUM_THREADS(nthreads)
	
	nopt = int(noptics)
	dimbse = ngrid(1)*ngrid(2)*ngrid(3)*nc*nv
	
	allocate(dipole(dimbse,4))
	allocate(omega(nopt))
	allocate(refractive(nopt,3))
	allocate(plres(nopt,3))
	
	dipole = 0.0
	refractive = 0.0
	plres = 0.0
	
	select case (mode)
	
	case("IPA")
		!inputs
		OPEN(UNIT=100, FILE=trim(outputfolder)//"ipa_oscf.dat",STATUS='unknown', IOSTAT=erro)
    		if (erro/=0) stop "Error opening ipa_oscf.dat input file"
    		OPEN(UNIT=200, FILE=trim(outputfolder)//"ipa_refractive_index.dat",STATUS='unknown', IOSTAT=erro)
		if (erro/=0) stop "Error opening ipa_refractive_index output file"

		!ouputs
		OPEN(UNIT=300, FILE=trim(outputfolder)//"ipa_PL.dat",STATUS='unknown', IOSTAT=erro)
    		if (erro/=0) stop "Error opening ipa_PL.dat output file"
	
	
			read(100,*) cflag,cflag,cflag,cflag,cflag,cflag,cflag,cflag
		
		do i=1,dimbse
		
			read(100,*) dipole(i,1),dipole(i,2),dipole(i,3),dipole(i,4),aux,aux,aux
		
		end do

			read(200,*) cflag,cflag,cflag,cflag,cflag,cflag,cflag,cflag
		
		do i=1,nopt
		
			read(200,*) omega(i),refractive(i,1),refractive(i,2),refractive(i,3),aux,aux,aux
		
		end do
	
	
	case("IPP")
		!inputs	
		OPEN(UNIT=100, FILE=trim(outputfolder)//"ipa_oscf-pol.dat",STATUS='unknown', IOSTAT=erro)
    		if (erro/=0) stop "Error opening ipa_oscf-pol.dat input file"
    		OPEN(UNIT=200, FILE=trim(outputfolder)//"ipa_refractive_index-pol.dat",STATUS='unknown', IOSTAT=erro)
		if (erro/=0) stop "Error opening ipa_refractive_index-pol output file"

		!ouputs
		OPEN(UNIT=300, FILE=trim(outputfolder)//"ipa_PL-pol.dat",STATUS='unknown', IOSTAT=erro)
    		if (erro/=0) stop "Error opening ipa_PL-pol.dat output file"
    		
			read(100,*) cflag,cflag,cflag,cflag,cflag,cflag,cflag
		
		do i=1,dimbse
		
			read(100,*) dipole(i,1),aux,aux,aux,dipole(i,2),dipole(i,3)
		
		end do

			read(200,*) cflag,cflag,cflag,cflag,cflag,cflag,cflag
		
		do i=1,nopt
		
			read(200,*) omega(i),aux,aux,aux,refractive(i,1),refractive(i,2)
		
		end do	    		
	
	case("BSE")
		!inputs	
		OPEN(UNIT=100, FILE=trim(outputfolder)//"bse_oscf.dat",STATUS='unknown', IOSTAT=erro)
    		if (erro/=0) stop "Error opening bse_oscf.dat input file"    		
    		OPEN(UNIT=200, FILE=trim(outputfolder)//"bse_refractive_index.dat",STATUS='unknown', IOSTAT=erro)
		if (erro/=0) stop "Error opening bse_refractive_index output file"

		!ouputs
		OPEN(UNIT=300, FILE=trim(outputfolder)//"bse_PL.dat",STATUS='unknown', IOSTAT=erro)
    		if (erro/=0) stop "Error opening bse_PL.dat output file"
    		
 			read(100,*) cflag,cflag,cflag,cflag,cflag,cflag,cflag,cflag
		
		do i=1,dimbse
		
			read(100,*) dipole(i,1),dipole(i,2),dipole(i,3),dipole(i,4),aux,aux,aux
		
		end do

			read(200,*) cflag,cflag,cflag,cflag,cflag,cflag,cflag,cflag
		
		do i=1,nopt
		
			read(200,*) omega(i),refractive(i,1),refractive(i,2),refractive(i,3),aux,aux,aux
		
		end do   		
	
	case("BSP")
		!inputs	
		OPEN(UNIT=100, FILE=trim(outputfolder)//"bse_oscf-pol.dat",STATUS='unknown', IOSTAT=erro)
    		if (erro/=0) stop "Error opening bse_oscf-pol.dat input file"
		OPEN(UNIT=200, FILE=trim(outputfolder)//"bse_refractive_index-pol.dat",STATUS='unknown', IOSTAT=erro)
		if (erro/=0) stop "Error opening bse_refractive_index-pol output file"    		

		!ouputs
		OPEN(UNIT=300, FILE=trim(outputfolder)//"bse_PL-pol.dat",STATUS='unknown', IOSTAT=erro)
    		if (erro/=0) stop "Error opening bse_PL-pol.dat output file"
	
			read(100,*) cflag,cflag,cflag,cflag,cflag,cflag,cflag
		
		do i=1,dimbse
		
			read(100,*) dipole(i,1),aux,aux,aux,dipole(i,2),dipole(i,3)
		
		end do

			read(200,*) cflag,cflag,cflag,cflag,cflag,cflag,cflag
		
		do i=1,nopt
		
			read(200,*) omega(i),aux,aux,aux,refractive(i,1),refractive(i,2)
		
		end do	
	
	case default
	
		write(*,*) "Not available option"
		STOP
	
	end select
	

	!$omp parallel do default(shared) private(i)
	do i=1,nopt
	
		call emissionpl(dimbse,omega(i),refractive(i,:),dipole,exctemp,eta,emin,plres(i,:))
	
	end do	
	!$omp end parallel do
	
	
	select case (mode)
	
	case("IPA")
	
		write(300,*) "# energy x y z"
		
		do i=1,nopt
		
			write(300,*) omega(i),plres(i,1),plres(i,2),plres(i,3)
		
		end do
	
	case("IPP")
	
		write(300,*) "# energy sp sm"
		
		do i=1,nopt
		
			write(300,*) omega(i),plres(i,1),plres(i,2)
		
		end do		
	
	case("BSE")
	
		write(300,*) "# energy x y z"
		
		do i=1,nopt
		
			write(300,*) omega(i),plres(i,1),plres(i,2),plres(i,3)
		
		end do		
	
	case("BSP")
	
		write(300,*) "# energy sp sm"	
		
		do i=1,nopt
		
			write(300,*) omega(i),plres(i,1),plres(i,2)
		
		end do		
	
	case default
	
		write(*,*) "Not available option"
		STOP
	
	end select
	
	
	
	
	deallocate(dipole)
	deallocate(omega)
	deallocate(refractive)
	
	close(100)
	close(200)
	close(300)


end subroutine emissionopt
