subroutine interpolate_vec_complex(w90basis,qpt,ngrid,nkpt,eqijk,qptmesh,rinterp)

	implicit none
	integer :: i,j
	integer :: w90basis
	integer,dimension(3) :: ngrid
	integer :: nkpt
	
	real :: cardinal_weight
	complex,dimension(nkpt,w90basis) :: eqijk
	real,dimension(3) :: qpt
	real,dimension(nkpt,3) :: qptmesh
	complex,dimension(w90basis) :: rinterp
	
	real :: lx,ly,lz
	
	rinterp = 0.0
	
	do i=1,nkpt
	
		lx = cardinal_weight(qpt(1),ngrid(1),qptmesh(i,1))
		ly = cardinal_weight(qpt(2),ngrid(2),qptmesh(i,2))
		lz = cardinal_weight(qpt(3),ngrid(3),qptmesh(i,3))

		do j=1,w90basis
		
			rinterp(j) = rinterp(j) + eqijk(i,j)*lx*ly*lz
	
		end do
	
	end do



end subroutine interpolate_vec_complex

subroutine interpolate_vec_real(w90basis,qpt,ngrid,nkpt,eqijk,qptmesh,rinterp)

	implicit none
	integer :: i,j
	integer :: w90basis
	integer,dimension(3) :: ngrid
	integer :: nkpt
	
	real :: cardinal_weight
	real,dimension(nkpt,w90basis) :: eqijk
	real,dimension(3) :: qpt
	real,dimension(nkpt,3) :: qptmesh
	real,dimension(w90basis) :: rinterp
	
	real :: lx,ly,lz
	
	rinterp = 0.0
	
	do i=1,nkpt
	
		lx = cardinal_weight(qpt(1),ngrid(1),qptmesh(i,1))
		ly = cardinal_weight(qpt(2),ngrid(2),qptmesh(i,2))
		lz = cardinal_weight(qpt(3),ngrid(3),qptmesh(i,3))

		do j=1,w90basis
		
			rinterp(j) = rinterp(j) + eqijk(i,j)*lx*ly*lz
	
		end do
	
	end do



end subroutine interpolate_vec_real



subroutine interpolate_complex(qpt,ngrid,nkpt,eqijk,qptmesh,rinterp)

	implicit none
	integer :: i
	integer,dimension(3) :: ngrid
	integer :: nkpt
	
	real :: cardinal_weight
	complex,dimension(nkpt) :: eqijk
	real,dimension(3) :: qpt
	real,dimension(nkpt,3) :: qptmesh
	complex :: rinterp
	
	real :: lx,ly,lz
	
	rinterp = cmplx(0.0,0.0)
	
	do i=1,nkpt
	
		lx = cardinal_weight(qpt(1),ngrid(1),qptmesh(i,1))
		ly = cardinal_weight(qpt(2),ngrid(2),qptmesh(i,2))
		lz = cardinal_weight(qpt(3),ngrid(3),qptmesh(i,3))
		
		rinterp = rinterp + eqijk(i)*lx*ly*lz
	
	
	end do



end subroutine interpolate_complex





subroutine interpolate_real(qpt,ngrid,nkpt,eqijk,qptmesh,rinterp)

	implicit none
	integer :: i
	integer,dimension(3) :: ngrid
	integer :: nkpt
	
	real :: cardinal_weight
	real,dimension(nkpt) :: eqijk
	real,dimension(3) :: qpt
	real,dimension(nkpt,3) :: qptmesh
	real :: rinterp
	
	real :: lx,ly,lz
	
	rinterp = 0.0
	
	do i=1,nkpt
	
		lx = cardinal_weight(qpt(1),ngrid(1),qptmesh(i,1))
		ly = cardinal_weight(qpt(2),ngrid(2),qptmesh(i,2))
		lz = cardinal_weight(qpt(3),ngrid(3),qptmesh(i,3))
		
		rinterp = rinterp + eqijk(i)*lx*ly*lz
	
	
	end do



end subroutine interpolate_real




function cardinal_weight(qpt,ngrid,qptmesh)

	implicit none
	
	real :: qpt,qptmesh
	integer :: ngrid
	real,parameter :: pi=acos(-1.)
	
	real :: numerator
	real :: denominator
	real :: cardinal_weight
	
	if (ngrid .eq. 1) then
	
		cardinal_weight = 1.000
		
	else
	
		if (abs(qpt-qptmesh) .lt. 1.0E-13 ) then
		
			cardinal_weight = 1.000
		else
		
			numerator = sin(pi*real(ngrid)*(qpt-qptmesh))
			
			if (mod(ngrid,2) .eq. 1) then
			
				denominator = real(ngrid)*sin(pi*(qpt-qptmesh))
			
			else
			
				denominator = real(ngrid)*tan(pi*(qpt-qptmesh))
			
			end if
		
			cardinal_weight = numerator/denominator
		
		
		end if
	
	
	
	end if	



end function cardinal_weight



