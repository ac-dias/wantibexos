program main

	implicit none

	integer :: i
	logical :: test

	test= .true.
	i=1

	if (test) then

	write(*,*) i,test

	else

	write(*,*) "problem"

	end if

end program main
