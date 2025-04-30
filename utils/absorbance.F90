program main

	implicit none
	
	integer :: i,j,k
	integer :: npt
	
	real,allocatable,dimension(:,:) :: abscoef
	real :: absorbance
	
	real :: thickness !in meters
	real :: totalabs
	
	
	!variaveis getarg

  	character*30 arg
  	integer ios
  	
       CALL getarg(1, arg)
		read(arg,*,iostat=ios) npt

      CALL getarg(2, arg)
		read(arg,*,iostat=ios) thickness 
		
		
	allocate(abscoef(npt,7))
	
	
	read(*,*) 
	
	do i=1,npt
	
		read(*,*) (abscoef(i,j), j=1,7)
		
		
		totalabs = abscoef(i,2)+abscoef(i,3)+abscoef(i,4)
		
		absorbance = 1.0 - exp(-2.0*totalabs*thickness*100.00)
		
		write(*,*) abscoef(i,1), absorbance
	
	end do
	
	deallocate(abscoef)

end program
