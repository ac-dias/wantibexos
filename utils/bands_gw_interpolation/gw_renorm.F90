program main

	implicit none
	
	character(len=70) :: outputfolder
	
	character(len=70) :: a,b
	character(len=70) :: add
	integer :: erro
	
	real :: rnmgw	
	
	integer :: i,j,k
	integer :: w90basis
	integer,dimension(3) :: ngrid
	integer :: nkpt
	
	character(len=3) :: charflag
	integer :: iflag
	real :: rflag
	
	
	real,allocatable,dimension(:,:) :: kptmesh
	real,allocatable,dimension(:,:,:) :: engw
	real,allocatable,dimension(:,:) :: corgwavg	
	
	real :: eigcor,eigcoravg
	
	!default values
	
	outputfolder = "./"
	
	do 

	read(*,*,iostat=erro) a,b
	if (erro/=0) exit

	select case (a)
	
	
	case ("OUTPUT=")
		
		outputfolder = b	
	
	case ("RNM_GW=")

		read(b,*) rnmgw	
	
		
	case default
	 continue

	end select


	end do	
	
	OPEN(UNIT=304, FILE=trim(outputfolder)//"gw_qp_energy_cor_avg.dat",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_qp_energy_cor_mesh_renorm input file"		
	OPEN(UNIT=305, FILE=trim(outputfolder)//"gw_qp_energy_cor_mesh.dat",STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_qp_energy_cor_mesh input file"
    	
	OPEN(UNIT=405, FILE=trim(outputfolder)//"gw_qp_energy_cor_mesh_renorm.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening gw_qp_energy_cor_mesh_renorm output file" 
    	

	read(305,*) w90basis 
	read(305,*) nkpt
	read(305,*) ngrid(1),ngrid(2),ngrid(3) 
	
	allocate(kptmesh(nkpt,3))
	allocate(engw(nkpt,w90basis,6))
	
	
	do i=1,nkpt
	
		read(305,*) charflag,charflag,iflag,kptmesh(i,1),kptmesh(i,2),kptmesh(i,3)
		read(305,*) charflag
		
		do j=1,w90basis
		
			read(305,*) engw(i,j,1),engw(i,j,2),engw(i,j,3),engw(i,j,4),engw(i,j,5),engw(i,j,6)

		
		
		end do
	
	
	end do  
 
	allocate(corgwavg(w90basis,5))  
	
	read(304,*) charflag	
	
	do i=1,w90basis
	
		read(304,*) corgwavg(i,1),corgwavg(i,2),corgwavg(i,3),corgwavg(i,4),corgwavg(i,5)
	
	end do 



	write(405,*) w90basis
	write(405,*) nkpt	
	write(405,*) ngrid(1),ngrid(2),ngrid(3)
	write(405,*)
	
	do i=1,nkpt
	
		write(405,*) '#kpt n:',i,kptmesh(i,1),kptmesh(i,2),kptmesh(i,3)
		write(405,*) '#EnTB G0W0cor selfx Re(selfc) Im(selfc) Znk'
		write(405,*)
		
	 do j=1,w90basis
	 
	 	eigcor =  engw(i,j,2)
	 	eigcoravg = corgwavg(j,1)
	 
	 	if (abs(eigcor-eigcoravg) .gt. rnmgw ) then
	 	
	 		engw(i,j,2)= corgwavg(j,1)
	 		engw(i,j,3) = corgwavg(j,2)
	 		engw(i,j,4) = corgwavg(j,3)
	 		engw(i,j,5) = corgwavg(j,4)
	 		engw(i,j,6) = corgwavg(j,5)
	 	
	 	end if
	 
	 	write(405,'(6E18.8)') engw(i,j,1),engw(i,j,2),engw(i,j,3),engw(i,j,4),engw(i,j,5),engw(i,j,6)
	 
	 end do
	
		write(405,*)
	
	end do	

	deallocate(kptmesh)
	deallocate(engw)	
	deallocate(corgwavg)
	
	
	close(304)
	close(305)
	close(405) 	    	

end program main
