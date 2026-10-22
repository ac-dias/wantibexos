program test_bse_q_optics

	use bse_q_optics
	implicit none

	type(rmn_data) :: data
	character(len=*),parameter :: filename='test_bse_q_optics_r.dat'
	character(len=*),parameter :: wsvecfile='test_bse_q_optics_wsvec.dat'
	character(len=256) :: message
	complex :: connection(3,1,1)
	complex :: phases(1,1),vc_vertex,vv_vertex
	complex :: c1(1),c2(1),v1(1),v2(1)
	complex :: phases_two(2,2),c_two(2),v_two(2),vc_reverse,vv_reverse
	real :: rlat(3,3),kpoint(3),pi
	real :: centers(3,1),kpoints(1,3)
	real :: centers_two(3,2),kpoints_two(2,3)
	integer :: unit,ios,failures
	logical :: ok

	failures=0
	pi=acos(-1.0)
	rlat=0.0
	rlat(1,1)=1.0
	rlat(2,2)=1.0
	rlat(3,3)=1.0

	open(newunit=unit,file=filename,status='replace',action='write',iostat=ios)
	if (ios /= 0) stop 'Unable to create temporary Wannier90 r-matrix fixture'
	write(unit,'(A)') 'test Wannier90 r-matrix fixture'
	write(unit,*) 1
	write(unit,*) 3
	write(unit,*) 0,0,0,1,1,1.0,0.0,0.0,0.0,3.0,0.0
	write(unit,*) 1,0,0,1,1,0.25,0.0,0.0,0.5,0.0,0.0
	write(unit,*) -1,0,0,1,1,0.25,0.0,0.0,-0.5,0.0,0.0
	close(unit)

	call rmn_read(filename,data,ok,message)
	if (.not. ok) then
		write(*,*) trim(message)
		failures=failures+1
	else
		call rmn_centers(data,centers,ok,message)
		if (.not. ok) then
			write(*,*) trim(message)
			failures=failures+1
		else
			if (maxval(abs(centers(:,1)-(/1.0,0.0,3.0/))) > 1.0e-6) failures=failures+1
			kpoints(1,:)=(/0.5*pi,0.0,0.0/)
			call center_phase_build(centers,kpoints,phases,ok,message)
			if (.not. ok) then
				write(*,*) trim(message)
				failures=failures+1
			else
				if (abs(phases(1,1)-cmplx(0.0,1.0)) > 1.0e-6) failures=failures+1
				c1=cmplx(1.0,0.0)
				c2=cmplx(1.0,0.0)
				v1=cmplx(1.0,0.0)
				v2=cmplx(1.0,0.0)
				call center_phase_direct_vertices(c1,c2,v1,v2,phases(:,1), &
					(/cmplx(1.0,0.0)/),vc_vertex,vv_vertex)
				if (abs(vc_vertex-cmplx(0.0,1.0)) > 1.0e-6) failures=failures+1
				if (abs(vv_vertex-cmplx(0.0,-1.0)) > 1.0e-6) failures=failures+1
				if (abs(vc_vertex*vv_vertex-cmplx(1.0,0.0)) > 1.0e-6) failures=failures+1
			end if
		end if

		centers_two=0.0
		centers_two(1,2)=1.0
		kpoints_two=0.0
		kpoints_two(1,1)=0.5*pi
		call center_phase_build(centers_two,kpoints_two,phases_two,ok,message)
		if (.not. ok) then
			write(*,*) trim(message)
			failures=failures+1
		else
			c_two=cmplx(1.0/sqrt(2.0),0.0)
			v_two=(/cmplx(0.0,0.0),cmplx(1.0,0.0)/)
			call center_phase_direct_vertices(c_two,c_two,v_two,v_two,phases_two(:,1), &
				phases_two(:,2),vc_vertex,vv_vertex)
			if (abs(vc_vertex-cmplx(0.5,0.5)) > 1.0e-6) failures=failures+1
			if (abs(vv_vertex-cmplx(0.0,-1.0)) > 1.0e-6) failures=failures+1
			if (abs(vc_vertex*vv_vertex-cmplx(0.5,-0.5)) > 1.0e-6) failures=failures+1
			call center_phase_direct_vertices(c_two,c_two,v_two,v_two,phases_two(:,2), &
				phases_two(:,1),vc_reverse,vv_reverse)
			if (abs(vc_reverse*vv_reverse-conjg(vc_vertex*vv_vertex)) > 1.0e-6) failures=failures+1
		end if

		kpoint=(/0.5*pi,0.0,0.0/)
		call rmn_bloch(data,kpoint,rlat,connection,ok,message)
		if (.not. ok) then
			write(*,*) trim(message)
			failures=failures+1
		else
			if (abs(connection(1,1,1)-cmplx(1.0,0.0)) > 1.0e-6) failures=failures+1
			if (abs(connection(2,1,1)-cmplx(-1.0,0.0)) > 1.0e-6) failures=failures+1
			if (abs(connection(3,1,1)-cmplx(3.0,0.0)) > 1.0e-6) failures=failures+1
		end if

		kpoint=0.0
		call rmn_bloch(data,kpoint,rlat,connection,ok,message)
		if (.not. ok) then
			write(*,*) trim(message)
			failures=failures+1
		else
			if (abs(connection(1,1,1)-cmplx(1.5,0.0)) > 1.0e-6) failures=failures+1
			if (abs(connection(2,1,1)) > 1.0e-6) failures=failures+1
			if (abs(connection(3,1,1)-cmplx(3.0,0.0)) > 1.0e-6) failures=failures+1
		end if
	end if

	call rmn_destroy(data)
	open(newunit=unit,file=wsvecfile,status='replace',action='write',iostat=ios)
	if (ios /= 0) stop 'Unable to create temporary Wannier90 wsvec fixture'
	write(unit,'(A)') 'wsvec is intentionally unsupported by this test path'
	close(unit)
	call rmn_read(filename,data,ok,message)
	if (ok) failures=failures+1
	call rmn_destroy(data)

	open(newunit=unit,file=filename,status='old',action='read',iostat=ios)
	if (ios == 0) close(unit,status='delete')
	open(newunit=unit,file=wsvecfile,status='old',action='read',iostat=ios)
	if (ios == 0) close(unit,status='delete')

	if (failures == 0) then
		write(*,*) 'bse_q_optics: PASS'
	else
		write(*,*) 'bse_q_optics: FAIL; mismatches: ',failures
		stop 1
	end if

end program test_bse_q_optics
