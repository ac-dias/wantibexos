program main

	implicit none

	character(len=70) :: params   !parametros TB	
	character(len=70) :: outputfolder
	character(len=70) :: kpaths    !kpath


	character(len=70) :: a,b
	character(len=70) :: add
	integer :: erro
	
	
	integer :: i,j,k
	integer :: w90basis
	integer,dimension(3) :: ngrid
	
	character(len=3) :: charflag
	integer :: iflag
	real :: rflag
	
	real,dimension(3,3) :: rlat,rvec
	
	real,allocatable,dimension(:,:) :: kptmesh,kptmeshd
	real,allocatable,dimension(:,:,:) :: engw
	real,allocatable,dimension(:) :: corgwavg
	
	integer :: nks
	integer :: nkpt,nkpts
	integer :: nkpathpts
	real,allocatable,dimension(:,:) :: ks
	real,allocatable,dimension(:,:) :: kpathpts,kpathptsd
	
	real,allocatable,dimension(:,:) :: eninterp,gwcorinterp
	
	

	!default values
	
	outputfolder = "./"
	kpaths = "unknown"
	
	do 

	read(*,*,iostat=erro) a,b
	if (erro/=0) exit

	select case (a)
	
	case ("PARAMS_FILE=")

		params = b
	
	case ("OUTPUT=")
		
		outputfolder = b	
	
	case ("KPATH_FILE=")

		kpaths = b	
	
		
	case default
	 continue

	end select


	end do	

	OPEN(UNIT=300, FILE= params,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening hamiltonian input file"
	OPEN(UNIT=304, FILE=trim(outputfolder)//"gw_qp_energy_cor_avg.dat",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_qp_energy_cor_mesh_renorm input file"    		
	OPEN(UNIT=305, FILE=trim(outputfolder)//"gw_qp_energy_cor_mesh_renorm.dat",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_qp_energy_cor_mesh_renorm input file" 
	OPEN(UNIT=202, FILE= kpaths,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening kpath input file"   
    	
    	 	
	OPEN(UNIT=306, FILE=trim(outputfolder)//"gw_interpolated_bands.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_interpolated_bands output file"    
    	
	read(300,*) charflag
	read(300,*) rflag
	read(300,*) rflag


	read(300,*) rlat(1,1),rlat(1,2),rlat(1,3)
	read(300,*) rlat(2,1),rlat(2,2),rlat(2,3)
	read(300,*) rlat(3,1),rlat(3,2),rlat(3,3)    	 	
    	

	read(305,*) w90basis 
	read(305,*) nkpt
	read(305,*) ngrid(1),ngrid(2),ngrid(3) 
	
	allocate(kptmesh(nkpt,3),kptmeshd(nkpt,3))
	allocate(engw(nkpt,w90basis,2))
	
	call recvec(rlat(1,:),rlat(2,:),rlat(3,:),rvec(1,:),rvec(2,:),rvec(3,:))
	
	do i=1,nkpt
	
		read(305,*) charflag,charflag,iflag,kptmesh(i,1),kptmesh(i,2),kptmesh(i,3)
		read(305,*) charflag
		
		do j=1,w90basis
		
			read(305,*) engw(i,j,1),engw(i,j,2),rflag,rflag,rflag,rflag
			
			engw(i,j,2) = engw(i,j,2) + engw(i,j,1) 
		
		
		end do
	
	
	end do  
	
	allocate(corgwavg(w90basis))
	
	read(304,*) charflag
	
	do i=1,w90basis
	
		read(304,*) corgwavg(i),rflag,rflag,rflag,rflag
	
	end do
	
	do i=1,nkpt
	
		call kptcar2kptdir(kptmesh(i,:),rvec,kptmeshd(i,:))
	
	end do
	

	read(202,*) nks
	read(202,*) nkpts

	allocate(ks(nks,3))

	do i=1,nks

		read(202,*) ks(i,1),ks(i,2),ks(i,3)
	
	end do	
	
	nkpathpts = (nks/2)*nkpts
	
	allocate(kpathpts(nkpathpts,4),kpathptsd(nkpathpts,4))	  				

	call kpath(outputfolder,rlat(1,:),rlat(2,:),rlat(3,:),nks,ks,nkpts,kpathpts)
	
	
	
	do i=1,nkpathpts
	
		kpathptsd(i,1) = kpathpts(i,1)
	
		call kptcar2kptdir(kpathpts(i,2:4),rvec,kpathptsd(i,2:4))
		
	end do
	
	
	!interpolation of energies and G0W0 corrections

	allocate(eninterp(nkpathpts,w90basis),gwcorinterp(nkpathpts,w90basis))


	do i=1,nkpathpts
	
	 do j=1,w90basis

		call interpolate_real(kpathptsd(i,2:4),ngrid,nkpt,engw(:,j,1),kptmeshd,eninterp(i,j))	
		call interpolate_real(kpathptsd(i,2:4),ngrid,nkpt,engw(:,j,2),kptmeshd,gwcorinterp(i,j))		 
	 
	 end do
	
	end do
	

	write(306,*)"#kpt EnTB EnG0W0 EnG0W0avg"

	do i=1,w90basis
	
	 do j=1,nkpathpts
		
		write(306,"(1E18.8,4E18.4)")  kpathpts(j,1),eninterp(j,i),gwcorinterp(j,i),eninterp(j,i)+corgwavg(i)
	 
	 end do
		write(306,*)
	
	end do	


	deallocate(kptmesh,kptmeshd)
	deallocate(engw)
	deallocate(ks,kpathpts,kpathptsd)
	deallocate(eninterp,gwcorinterp)
	deallocate(corgwavg)
	
	close(300)
	close(304)
	close(305)
	close(202)
	close(306)
	
end program main



