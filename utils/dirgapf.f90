!programa para achar o gap direto nas bandas do TB
!gfortran dirgapf.f90 -o dirgapf.x
!dirgapf.x nbands nocp nk < input > output

program main

	implicit none

	integer :: nbands,nocp,nk
	integer :: i,j,k,kflag
	double precision :: flag1,flag2
	character(len=1) :: cflag

	double precision,allocatable,dimension(:,:,:) :: bands
	double precision :: gap,gapflag
	double precision :: cbmin,vbmax

	!variaveis getarg

  	character*30 arg
  	integer ios

      CALL getarg(1, arg)
		read(arg,*,iostat=ios) nbands
      CALL getarg(2, arg)
		read(arg,*,iostat=ios) nocp
      CALL getarg(3, arg)
		read(arg,*,iostat=ios) nk

	allocate(bands(nbands,nk,2))

	do i=1,nbands

		read(*,*) cflag,kflag

		do j=1,nk
		read(*,*) bands(i,j,1),bands(i,j,2),flag1,flag2
		
		end do

		read(*,*)

	end do


	gap = 50.0
	cbmin = 50.0
	vbmax = -50.0

	do i=1,nk

		gapflag=bands(nocp+1,i,2)-bands(nocp,i,2)

		if (gapflag .lt. gap) then
			gap = gapflag

		else
			continue
		end if

		if (cbmin .gt. bands(nocp+1,i,2)) then
			cbmin = bands(nocp+1,i,2)

		else
			continue
		end if

		if (vbmax .lt. bands(nocp,i,2)) then
			vbmax = bands(nocp,i,2)

		else
			continue
		end if
	

	end do

	write(*,*) 'dirgap:',gap
	write(*,*) 'indgap:', cbmin-vbmax

	deallocate(bands)

end program main

