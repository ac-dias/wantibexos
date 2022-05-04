!programa para achar o numero de bandas de conducao e valencia para o calculo BSE
!gfortran nc_nv_finder.f90 -o nc_nv_finder.x
!nc_nv_finder.x  < input >> input

INCLUDE "../subroutines/module_input_read.f90" 

program main

	use input_variables
	use hamiltonian_input_variables
	implicit none

	integer :: nbands,nk
	integer :: i,j,k
	double precision :: flag1,flag2
	character(len=1) :: cflag

	integer :: nkpt,erro

	double precision,allocatable,dimension(:,:,:) :: bands
	double precision :: gap,gapflag,deltaen

	integer,allocatable,dimension(:) :: nocpk

	integer :: nks
	integer :: nkpts

	call input_read
	call hamiltonian_input_read(200,params)

	OPEN(UNIT=202, FILE= kpaths,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Error opening kpath input file "
	OPEN(UNIT=300, FILE=trim(outputfolder)//"bands.dat",STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Error opening bands input file"


	read(202,*) nks
	read(202,*) nkpts

	nbands = w90basis
	nk = (nks/2)*nkpts


	allocate(bands(nbands,nk,2))
	allocate(nocpk(nk))

	do i=1,nbands

		read(300,*) cflag,k

		do j=1,nk
		read(300,*) bands(i,j,1),bands(i,j,2),flag1,flag2
		
		end do

		read(300,*)

	end do


	do i=1,nk

		do j=1,nbands


		if (bands(j,i,2) .gt. 0.0) then

			nocpk(i)= j-1
			exit
		else

		continue

		end if

		end do
	end do


	gap = 50.0

	do i=1,nk

		gapflag=bands(nocpk(i)+1,i,2)-bands(nocpk(i),i,2)

		if (gapflag .lt. gap) then
			gap = gapflag
			nkpt = i

		else
			continue
		end if

	end do

	deltaen=ebsef-gap
	
	!write(*,*) nocpk(nkpt)
	!write(*,*) nkpt
	!write(*,*) deltaen

	do i=1,nocpk(nkpt)

		if (bands(i,nkpt,2) .gt. (bands(nocpk(nkpt),nkpt,2)-deltaen)) then
			nv = nocpk(nkpt)-(i-1)
			exit
		else
			nv = nocpk(nkpt)
		end if

	end do

	!write(*,*) "valence states"
	!do i=1,nocpk(nkpt)
	!	write(*,*) i,bands(i,nkpt,2)
	!end do


	do i=nocpk(nkpt)+1,nbands

		if (bands(i,nkpt,2) .gt. (bands(nocpk(nkpt)+1,nkpt,2)+deltaen)) then
			nc = i-nocpk(nkpt)-1
			exit
		else
			nc = nbands-nocpk(nkpt)

		end if

	end do
	
	if (nc .eq. 0) then
	
		nc = 1
	
	end if

	!write(*,*) "conduction states"
	!do i=nocpk(nkpt)+1,nbands
	!	write(*,*) i,bands(i,nkpt,2)
	!end do

	!write(*,*) 
	write(*,"(A9,I0)")"NBANDSC= ",nc
	write(*,"(A9,I0)")"NBANDSV= ",nv

	!write(*,*) 'dirgap:',gap

	deallocate(bands)
	deallocate(nocpk)

	close(202)
	close(303)

end program main




